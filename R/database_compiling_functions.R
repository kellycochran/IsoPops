#' Raw Isoform Database Object
#'
#' A compiled TranscriptDB + a compiled GffDB, from the same experiment.
#' Raw Databases, unlike processed Databases, do not include any gene
#' information. Raw Databases are also different from processed Databases in
#' that they do not contain an OrfDB, or table of ORF information, although the
#' ORFs are added as part of the TranscriptDB if provided. Raw Databases are
#' unfiltered -- if your dataset includes off-target reads, they will still be
#' present in the Raw Database.
#' See \code{\link{compile_raw_db}} for more on how this object is constructed.
#'
#' @param transcriptDB A parsed transcript data frame, returned by
#' \code{get_transcript_sequences} (internal).
#' @param gffDB A parsed gff data frame, returned by
#' \code{get_gff_data} (internal).
#' @return A RawDatabase object.
#' @export
RawDatabase <- function(transcriptDB, gffDB) {
  self <- list(TranscriptDB = transcriptDB, GffDB = gffDB)
  class(self) <- append(class(self), "RawDatabase")
  return(self)
}


#' Processed Isoform Database Object
#'
#' A compiled TranscriptDB, GffDB, GeneDB, and optional OrfDB from the same
#' experiment. Data frames contain additional information compared to a Raw
#' Database object. The datasets are not automatically filtered during Database
#' creation.
#' See \code{\link{process_db}} for more on how this object is constructed.
#'
#' @param transcriptDB A parsed transcript data frame, returned by
#' \code{get_transcript_sequences} (internal).
#' @param gffDB A parsed gff data frame, returned by
#' \code{get_gff_data} (internal).
#' @param orfDB A parsed ORF data frame, returned by
#' \code{make_ORF_db} (internal).
#' @param geneDB A data frame created using gene names and IDs. See
#' \code{\link{process_db}}.
#' @return A Database object.
#' @export
Database <- function(transcriptDB, gffDB, orfDB, geneDB) {
  self <- list(TranscriptDB = transcriptDB, GffDB = gffDB, OrfDB = orfDB,
               GeneDB = geneDB)
  class(self) <- append(class(self), "Database")
  return(self)
}


#' Raw Isoform Database Creation
#'
#' Parses input files and compiles data into a TranscriptDB data frame and a
#' GffDB data frame, to be collected into a Raw Database object.
#'
#' @param transcript_filename A FASTA or tsv-formatted file containing the
#' transcript sequences for each isoform and the PB (PacBio) IDs. If a FASTA
#' file is given, a temporary tsv-formatted file will be produced as part of
#' file parsing. The preferred TSV format has the PB ID for each transcript
#' in the first column, and the corresponding sequence in the second column,
#' with columns separated by tab characters.
#' @param abundance_filename A tsv-formatted abundance file, formatted as
#' output by IsoSeq 2 or 3 -- the first column must contain isoform PB IDs
#' matching the IDs from the GFF and FASTA files, and the second column must
#' contain the number of full-length reads (or another quantification metric).
#' @param gff_filename A GFF2 format file of the isoforms, with isoform PB IDs
#' matching what is in the transcript FASTA file.
#' @param ORF_filename Optional. A FASTA or tsv-formatted file containing the
#' ORF peptide sequences for each isoform and the PB (PacBio) IDs. If a FASTA
#' file is given, a temporary tsv-formatted file will be produced as part of
#' file parsing. The preferred TSV format has the PB ID for each ORF in the
#' first column, and the corresponding AA sequence in the second column,
#' with columns separated by tab characters. IDs should match those found in
#' the FASTA/abundance/GFF files given.
#' @return A RawDatabase object.
#' @export
compile_raw_db <- function(transcript_filename, abundance_filename,
                           gff_filename, ORF_filename = NULL) {

  transcriptDB <- get_transcript_sequences(transcript_filename)
  transcriptDB <- add_abundance_data(transcriptDB, abundance_filename)
  transcriptDB <- calc_abundances(transcriptDB)

  if (!is.null(ORF_filename)) {
    transcriptDB <- add_ORF_sequences(transcriptDB, ORF_filename)
  }

  gffDB <- get_gff_data(gff_filename)
  transcriptDB <- add_gff_info(transcriptDB, gffDB)
  transcriptDB <- unique(transcriptDB)

  print("Done.")
  return(RawDatabase(transcriptDB, gffDB))
}


#' Database Gene Curation
#'
#' Processes a Raw Database object using input gene information. Assuming a
#' targeted experiment where isoforms whose IDs do not match that of a gene are
#' "off-target," this function removes off-target transcripts from the
#' TranscriptDB and GffDB and creates an OrfDB using the predicted ORFs from
#' targeted transcripts, if ORFs were provided when the Raw Database was
#' created. It also generates several statistics relevant to each gene,
#' collected in a GeneDB data frame. The processed Database can be input to any
#' plot or filtering function of this package, or its tables can be saved to
#' csv or GFF2 format for use in other analyses. Note: in order for the GeneDB
#' to calculate the Shannon Index and Gini Coefficient for isoform diversities,
#' the packages \code{ineq} and \code{vegan} are required.
#'
#' @param rawDB A \code{\link{RawDatabase}} object created by
#' \code{\link{compile_raw_db}}.
#' @param gene_ID_table A data frame consisting of a column of gene names and a
#' column of IDs (corrsponding to isoform IDs in the RawDatabase), with column
#' names "Name" and "ID".
#' @return A \code{\link{Database}} object.
#' @examples
#' \dontrun{
#' gene_ID_table <- data.frame(ID = c("PB.1", "PB.2"),
#' Name = c("Gene1", "Gene2"))
#' rawDB <- compile_raw_db(transcript_file, abundance_file, gff_file, ORF_file)
#' DB <- process_db(rawDB, gene_ID_table)
#' }
#' @export
process_db <- function(rawDB, gene_ID_table) {
  print("Processing database...")
  geneDB <- data.frame(ID = gene_ID_table$ID, Name = gene_ID_table$Name,
                       stringsAsFactors = F)
  transcriptDB <- add_gene_names_to_transcripts(rawDB$TranscriptDB, geneDB)
  # filter to only keep targeted transcripts
  # (or transcripts that we can label with a gene name)
  transcriptDB <- transcriptDB[!is.na(transcriptDB$Gene), ]

  geneDB <- add_isoform_info(transcriptDB, geneDB)
  # remove genes with no reads
  if (any(is.na(geneDB$Isoforms)) || any(geneDB$Isoforms == 0)) {
    geneDB <- geneDB[which(!is.na(geneDB$Isoforms) & geneDB$Isoforms > 0), ]
  }
  geneDB <- add_transcript_lengths(transcriptDB, geneDB)
  geneDB <- add_diversity_indeces(geneDB, transcriptDB = transcriptDB)
  gffDB <- rawDB$GffDB[rawDB$GffDB$PBID %in% transcriptDB$PBID, ]

  # add ORF info (optional)
  if ("ORF" %in% colnames(rawDB$TranscriptDB)) {
    orfDB <- make_ORF_db(transcriptDB)
    geneDB <- add_ORF_info(transcriptDB, geneDB, orfDB)
    geneDB <- add_diversity_indeces(geneDB, orfDB = orfDB)
  } else {
    orfDB <- NULL
  }
  print("Done.")
  return(Database(transcriptDB, gffDB, orfDB, geneDB))
}
