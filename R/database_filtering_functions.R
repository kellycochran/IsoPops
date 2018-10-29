#' Filter Database
#'
#' Filters the isoforms in a Database object using one or both of a minimum
#' exon count cutoff and a minimum full-length read count cutoff.
#'
#' @param database A processed \code{\link{Database}} object.
#' @param exon_min_count Any transcripts with fewer than \code{exon_min_count}
#' exons will be filtered out. Set to 1 to not filter by exon count.
#' @param abund_cutoff A number between 0 and 1, representing the percent
#' of isoform abundance to retain in filtering. For example, if 0.95, then only
#' the top 95% of the abundance for each gene is retained, through filtering
#' out the least abundant 5% of isoforms. Isoforms are weighted by their
#' full-length read count. In other words, the number of isoforms retained for
#' each gene is either equal to or one less than the gene's N95 value, as it
#' will not include an isoform whose abundance overlaps with the cutoff. Set to
#' 1 to not filter by abundance.
#' @param recalc_abundances Logical. If filtering is "permanent", or if the
#' resulting Database will be piped into other filtering functions in this
#' package, then it is recommended to recalculate relative abundances for each
#' gene. This only changes the values PercAbundance and CumPercAbundance (not
#' the full-length read count).
#' @return A \code{\link{Database}} object, containing only retained isoforms.
#' @examples
#' \dontrun{
#' rawDB <- compile_raw_db(transcript_file, abundance_file, gff_file, ORF_file)
#' DB <- process_db(rawDB, gene_ID_table)
#' DB_filtered <- filter_db(DB, 3, 0.90)
#' }
#' @seealso \code{\link{filter_db_by_ID}}
#' @export
filter_db <- function(database, exon_min_count = 4, abund_cutoff = 0.95,
                      recalc_abundances = T) {
  print("Filtering...")
  filters <- data.frame(PBID = database$TranscriptDB$PBID)
  # below are the implemented filtering logics:
  # columns have TRUE value if transcript should be filtered out
  filters$OffTarget <- is.na(database$TranscriptDB$Gene)
  filters$LessThan4Exons <- database$TranscriptDB$ExonCount < exon_min_count

  if (abund_cutoff < 1) {
    filters$AboveCut <- database$TranscriptDB$CumPercAbundance >= abund_cutoff
  }
  rownames(filters) <- filters$PBID
  filters$PBID <- NULL
  retain_IDs <- database$TranscriptDB$PBID[apply(filters, 1,
                                                function(x) !any(x) )]
  retain_isoforms <- database$TranscriptDB$PBID %in% retain_IDs
  filtered_transcriptDB <- database$TranscriptDB[retain_isoforms, ]

  if (recalc_abundances) {
    # TODO: move this into calc_abundances?
    filtered_transcriptDB$PercAbundance <- NULL
    filtered_transcriptDB$CumPercAbundance <- NULL
    filtered_transcriptDB <- calc_abundances(filtered_transcriptDB)
  }
  # rebuild gffDB, geneDB, and orfDB post-filtering
  gffDB <- database$GffDB[database$GffDB$PBID %in% retain_IDs,]

  geneDB <- database$GeneDB[,c("ID", "Name")]
  geneDB <- add_isoform_info(filtered_transcriptDB, geneDB)
  # remove genes with no transcripts left after filtering
  if (any(is.na(geneDB$Isoforms)) || any(geneDB$Isoforms == 0)) {
    geneDB <- geneDB[(!is.na(geneDB$Isoforms) & geneDB$Isoforms > 0), ]
  }
  geneDB <- add_transcript_lengths(filtered_transcriptDB, geneDB)
  geneDB <- add_diversity_indeces(geneDB, transcriptDB = filtered_transcriptDB)

  if ("ORF" %in% colnames(filtered_transcriptDB)) {
    orfDB <- make_ORF_db(filtered_transcriptDB)
    geneDB <- add_ORF_info(filtered_transcriptDB, geneDB, orfDB)
    geneDB <- add_diversity_indeces(geneDB, orfDB = orfDB)
  } else {
    orfDB <- NULL
  }
  print("Done.")
  return(Database(filtered_transcriptDB, gffDB, orfDB, geneDB))
}

#' Filter Database By ID
#'
#' Filters the isoforms in a Database object based on the input list of
#' isoform IDs.
#'
#' @param database A processed \code{\link{Database}} object.
#' @param IDs_to_remove List of isoform IDs to filter out.
#' @param recalc_abundances Logical. If filtering is permanent, or if the
#' resulting Database is being used with other function in this package, then
#' it is recommended to recalculate relative abundances for each gene. This
#' only changes the values PercAbundance and CumPercAbundance (not the
#' full-length read count).
#' @return A \code{\link{Database}} object, containing only retained isoforms.
#' @examples
#' \dontrun{
#' rawDB <- compile_raw_db(transcript_file, abundance_file, gff_file, ORF_file)
#' DB <- process_db(rawDB, gene_ID_table)
#' DB_filtered <- filter_db_by_ID(DB, c("PB.1.1", "PB.1.2"))
#' }
#' @seealso \code{\link{filter_db}}
#' @export
filter_db_by_ID <- function(database, IDs_to_remove, recalc_abundances = T) {
  print("Filtering...")

  filtered_transcriptDB <- database$TranscriptDB[
    !(database$TranscriptDB$PBID %in% IDs_to_remove), ]

  if (recalc_abundances) {
    filtered_transcriptDB$PercAbundance <- NULL
    filtered_transcriptDB$CumPercAbundance <- NULL
    filtered_transcriptDB <- calc_abundances(filtered_transcriptDB)
  }

  gffDB <- database$GffDB[!(database$GffDB$PBID %in% IDs_to_remove), ]
  orfDB <- make_ORF_db(filtered_transcriptDB)

  geneDB <- data.frame(ID = database$GeneDB$ID, Name = database$GeneDB$Name,
                       stringsAsFactors = F)
  geneDB <- add_isoform_info(filtered_transcriptDB, geneDB)
  # remove genes with no reads
  if (any(is.na(geneDB$Isoforms)) || any(geneDB$Isoforms == 0)) {
    geneDB <- geneDB[(!is.na(geneDB$Isoforms) & geneDB$Isoforms > 0), ]
  }
  geneDB <- add_transcript_lengths(filtered_transcriptDB, geneDB)
  geneDB <- add_ORF_info(filtered_transcriptDB, geneDB, orfDB)
  geneDB <- add_diversity_indeces(geneDB, transcriptDB = filtered_transcriptDB)

  print("Done.")
  return(Database(filtered_transcriptDB, gffDB, orfDB, geneDB))
}

