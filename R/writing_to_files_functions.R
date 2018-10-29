#' Write To GFF File
#'
#' Writes the GffDB portion of a \code{\link{RawDatabase}} or
#' \code{\link{Database}} object, or a GffDB if input by itself, to a file in
#' GFF2 format. This file can then be viewed in IGV or used with other
#' programs. See sourcecode if attempting to input your own gff-like data
#' frame to know which fields are required.
#'
#' @param database Either a \code{\link{RawDatabase}}, \code{\link{Database}},
#' or a GffDB (used by \code{\link{RawDatabase}} and \code{\link{Database}}).
#' @param dest_filename Path and name of file where the GFF file should be
#' written.
#' @examples
#' \dontrun{
#' rawDB <- compile_raw_db(transcript_file, abundance_file, gff_file)
#' write_gff_data(rawDB)
#' }
#' @export
write_gff_data <- function(database, dest_filename) {
  if ("Database" %in% class(database) || "RawDatabase" %in% class(database)) {
    # input can be a Database object...
    gff_data <- database$GffDB
  } else {
    # or just the "gffDB" data frame by itself
    gff_data <- database
  }
  # explicitly mimics original GFF2 file format
  # note: this isn't GFF3, might not be compatible with some programs :(
  df <- data.frame(Chromosome = gff_data$Chromosome)
  df$V2 <- rep("PacBio", dim(df)[1])
  df$TE <- gff_data$TranscriptExon
  df$Start <- gff_data$Start
  df$End <- gff_data$End
  df$V6 <- rep(".", dim(df)[1])
  df$Strand <- gff_data$Strand
  df$V8 <- rep(".", dim(df)[1])
  df$V9 <- sapply(gff_data$PBID, function(x) {
    return(paste(c('gene_id "', stringr::str_match(x, "PB\\.[0-9]*"),
                   '"; transcript_id "', x, '";'), collapse = ""))
  })
  utils::write.table(df, dest_filename, sep = "\t", row.names = F,
                     col.names = F, quote = F)
}


# TODO: make this function filter entire DBs
# TODO: allow this function to take in just a GffDB

#' Get Top Isoforms From Each Gene
#'
#' Filters the GffDB portion of \code{\link{Database}} object by the most
#' abundant isoforms. The result can be input to \code{\link{write_gff_data}},
#' and then viewed in IGV or used with other programs.
#'
#' @param database A processed \code{\link{Database}} object.
#' @param top_num Number of isoforms to retain for each gene (ex: top 5, 10).
#' @examples
#' \dontrun{
#' rawDB <- compile_raw_db(transcript_file, abundance_file, gff_file)
#' write_gff_data(filter_top_gff(rawDB, 5))
#' }
#' @export
filter_top_gff <- function(database, top_num) {
  # first, get IDs for the top X isoforms in each gene
  topPBIDs <- unlist(lapply(database$GeneDB$Name, function(gene) {
    one_gene <- database$TranscriptDB[database$TranscriptDB$Gene == gene, ]
    one_gene <- one_gene[rev(order(one_gene$PercAbundance)), ]
    minimum <- min(top_num, length(one_gene$PBID))
    return(one_gene$PBID[1:minimum])
  }))
  return(database$GffDB[database$GffDB$PBID %in% topPBIDs, ])
}
