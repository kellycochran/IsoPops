#' Filter Truncations From Database
#'
#' Filters truncations of isoforms from a Database object. A truncated isoform
#' is defined as one which matches a longer isoform exactly, besides being
#' shortened on one or both ends (i.e. all junctions found in the shorter
#' isoform are also found in the longer isoform). Small discrepancies in
#' the 5' and 3' ends of transcripts are tolerated. When a truncation is
#' identified, its read count is added to the isoform which it is a truncation
#' of. In other words, the total sum of read counts in the database is
#' preserved even though transcript entries are removed.
#'
#' Isoform truncation is
#' determined through pairwise comparison between all isoforms of a gene.
#' Note: this filtering step is computationally intensive.
#'
#' @param database A processed \code{\link{Database}} object.
#' @return A \code{\link{Database}} object, containing only non-truncated
#' isoforms.
#' @examples
#' \dontrun{
#' rawDB <- compile_raw_db(transcript_file, abundance_file, gff_file, ORF_file)
#' DB <- process_db(rawDB, gene_ID_table)
#' DB_filtered <- filter_truncations(DB)
#' }
#' @seealso \code{\link{filter_db}}, \code{\link{filter_db_by_ID}}
#' @export
filter_truncations <- function(database) {
  # Proceed with caution.

  if ("RawDatabase" %in% class(database)) {
    # TODO: allow RawDBs
    warning("Truncations can only be filtered from processed Databases,
          not RawDatabases.")
    return(database)
  }
  truncations <- rep(NA, nrow(database$TranscriptDB))
  superset_transcripts <- rep(NA, nrow(database$TranscriptDB))
  index <- 1

  text <- "Checking truncations for "
  for (i in 1:nrow(database$GeneDB)) {
    gene <- database$GeneDB$Name[i]

    print(paste(text, gene, "...", sep = ""))
    text <- ""

    IDs <- database$TranscriptDB$PBID[database$TranscriptDB$Gene == gene]
    gene_gff <- database$GffDB[database$GffDB$TranscriptExon == "exon"
                               & database$GffDB$PBID %in% IDs, ]

    # will use these to compare start/end sites of transcripts
    # to the most extreme start/end sites for that gene
    max_start_pos_strand <- max(gene_gff$Start)
    max_start_neg_strand <- min(gene_gff$End)

    # compare each transcript pairwise with all other transcripts for gene
    for (j in seq_along(IDs)) {
      transcript_gff <- gene_gff[gene_gff$PBID == IDs[j], ]
      transcript_gff$PBID <- NULL

      other_IDs <- IDs[-j]
      matches <- sapply(other_IDs, function(other_ID) {
        check_for_match(gene_gff, transcript_gff, max_start_pos_strand,
                        max_start_neg_strand, other_ID)
      })
      if (any(matches)) {
        # transcript is indeed a truncation of another transcript
        # record both the truncated transcript and its super-transcript
        truncations[index] <- IDs[j]
        superset_transcripts[index] <- other_IDs[which(matches)[1]]
      }
      index <- index + 1
    }
  }
  truncated_IDs <- as.character(stats::na.omit(truncations))
  superset_IDs <- as.character(stats::na.omit(superset_transcripts))

  # add the abundances of truncated transcripts to their super-transcripts,
  # before filtering out truncated transcripts
  while (sum(sapply(truncated_IDs, function(ID) {
    return(database$TranscriptDB$FL_reads[database$TranscriptDB$PBID == ID])
  })) > 0) {  # end while line

    for (i in seq_along(truncated_IDs)) {
      sub_isoform <- which(database$TranscriptDB$PBID == truncated_IDs[i])
      super_isoform <- which(database$TranscriptDB$PBID == superset_IDs[i])
      new_count <- database$TranscriptDB$FL_reads[super_isoform] +
        database$TranscriptDB$FL_reads[sub_isoform]
      database$TranscriptDB$FL_reads[super_isoform] <- new_count
      database$TranscriptDB$FL_reads[sub_isoform] <- 0
    }
  }
  print("Finished finding truncations.")
  return(filter_db_by_ID(database, truncated_IDs, recalc_abundances = T))
}


# internal
check_for_match <- function(gene_gff, transcript, max_start_pos_strand,
                            max_start_neg_strand, other_ID) {
  # Proceed with a lot of caution.

  other_transcript <- gene_gff[gene_gff$PBID == other_ID, ]
  other_transcript$PBID <- NULL

  if (dim(other_transcript)[1] >= dim(transcript)[1]) {
    # Attempt to align exon starts between transcripts.
    # Assume best alignment has minimum number of unmatched *starts*
    nonmatch_starts <- min(sapply(seq_len(length(other_transcript$Start)
                                          - length(transcript$Start) + 1),
                                  function(i) {
                                    sum(!(other_transcript$Start[i:(i + length(transcript$Start) - 1)]
                                          %in% transcript$Start))
                                  }))

    # Multiple exon start positions don't line up between the two transcripts,
    # so this pair can't be a sub/super-transcript pair.
    # Exit early with a NO
    if (nonmatch_starts > 1) return(FALSE)

    # Attempt to align exon starts between transcripts.
    # Assume best alignment has minimum number of unmatched *ends*
    nonmatch_ends <- min(sapply(1:(length(other_transcript$End)
                                   - length(transcript$End) + 1),
                                function(i) {
                                  sum(!(other_transcript$End[i:(i + length(transcript$End) - 1)]
                                        %in% transcript$End))
                                }))

    # Multiple exon end positions don't line up between the two transcripts,
    # so this pair can't be a sub/super-transcript pair.
    # Exit early with a NO
    if (nonmatch_ends > 1) return(FALSE)

    # If all exon start and end positions in one transcript correspond to
    # start and end positions in the other transcript, then they are a
    # sub/super-transcript pair! Exit early with a YES
    starts_match <- (nonmatch_starts == 0)
    ends_match <- (nonmatch_ends == 0)
    if (starts_match & ends_match) return(TRUE)


    # Because PacBio/IsoSeq may misreport the outermost boundaries of
    # a transcript, we exclude only the first start / last end positions
    # and then check for a match.

    # First check if starts "almost" match
    starts <- transcript$Start[-1]
    other_starts <- other_transcript$Start[-1]
    starts_almost_match <- any(
      sapply(1:(length(other_starts) - length(starts) + 1), function(i) {
        all(other_starts[i:(i + length(starts) - 1)] %in% starts)
      }))

    # Next check ends "almost" match
    ends <- transcript$End[-length(transcript$End)]
    other_ends <- other_transcript$End[-length(other_transcript$End)]
    ends_almost_match <- any(
      sapply(1:(length(other_ends) - length(ends) + 1), function(i) {
        all(other_ends[i:(i + length(ends) - 1)] %in% ends)
      }))

    # Then see if a match can be called

    # If strand is positive:
    if (transcript$Strand[1] == "+"
        && starts_almost_match && ends_match
        # other transcript needs to be larger than main transcript:
        && nrow(other_transcript) > nrow(transcript)
        && transcript$Start[1] > other_transcript$Start[
          which(other_transcript$End %in% transcript$End)[1]])  {
      return(TRUE)
    } else {
      # If strand is negative:
      if (ends_almost_match && starts_match
          # other transcript needs to be larger than main transcript:
          && nrow(other_transcript) > nrow(transcript)
          && transcript$End[length(transcript$End)] < other_transcript$End[
            max(which(other_transcript$Start %in% transcript$Start))]) {
        return(TRUE)
      }
    }

    # if you reach here, an "almost" match could not be made, even by ignoring
    # the start of the main transcript. So now we try to ignore the end of the
    # last exon of the **gene** (not isoform), which is likely UTR, and see if
    # a match happens.
    if (transcript$Strand[1] == "+") {
      if (ends_almost_match && starts_almost_match
          && nrow(other_transcript) > nrow(transcript)
          && transcript$Start[1] > other_transcript$Start[
            which(other_transcript$End %in% transcript$End)[1]]
          && transcript$End[length(transcript$End)] > max_start_pos_strand) {
        return(TRUE)
      }
    } else {
      # If strand is negative:
      if (ends_almost_match && starts_almost_match
          && nrow(other_transcript) > nrow(transcript)
          && transcript$End[length(transcript$End)] < other_transcript$End[
            max(which(other_transcript$Start %in% transcript$Start))]
          && transcript$Start[1] < max_start_neg_strand) {
        return(TRUE)
      }
    }
  }
  # You tried everything. The transcripts are not a sub/super-transcript pair.
  return(FALSE)
}
