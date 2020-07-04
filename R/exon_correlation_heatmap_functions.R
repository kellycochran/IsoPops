# internal
weighted_corr <- function(exon_call_df, weights) {
  # repeat each transcript's exon calls [FL_reads] times
  # to account for transcript abundance in correlation computation
  weighted_df <- data.frame(matrix(ncol = ncol(exon_call_df), byrow = T,
                              data = unlist(lapply(seq_len(nrow(exon_call_df)),
                        function(i) { rep(exon_call_df[i, ], weights[i]) }))))
  return(round(cor(weighted_df), 2))
}

#' Exon Correlation Heatmap
#'
#' Generates a 2D heatmap where each axis is the exons for a gene, and the
#' values in the heatmap correspond to the correlation between the splicing
#' events of pairs of exons. For example, the heatmap cell in row i and column
#' j contains the pearson correlation of all the observed splicing inclusions
#' and exclusions of exons i and j, according to the transcripts in the data.
#' This plot can show which exons tend to either be included or spliced out
#' together, for instance, and any exon pairs which may have mutually exclusive
#' splicing patterns. Exon presence within a transcript is determined by
#' literal string matching, so only full and completely correct matches between
#' exon sequence and transcript sequence are considered.
#'
#' @param database A compiled Database object.
#' @param exon_filename Path to a file in either FASTA or TSV format. If in
#' FASTA format, the sequences are the annotated sequences for all exons in
#' the gene, and the IDs are the exon names (will be displayed in the plot).
#' The ID line must have format ">exonname". If in TSV format, There must be
#' one column for exon names and one column for the exon sequence,
#' tab-separated.
#' @param gene The desired gene to plot. Note that the plot will be generated
#' only from exon matches to transcripts for the given gene, so no off-target
#' exon matches are possible.
#' @param weighted Logical. If TRUE, transcript abundances will be taken into
#' account when correlations are calculated (recommended).
#' @param exons_to_include Vector of exon names to subset from the input file.
#' Default is to include all exons in the inut file. This list is ordered; in
#' other words, if you would like to rearrange the order of exon names on the
#' axes of the heatmap, use this argument to do so.
#' @param weights A numeric vector specifying the weights to apply to each
#' transcript for the given gene. Default is the number of full-length reads
#' for the transcript.
#' @param plot_hist Logical. If TRUE, a histogram of all exon correlations
#' across the gene is produced.
#' @param symetric Logical. If TRUE, both sides of the symmetric heatmap are
#' shown.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' gene_ID_table <- data.frame(ID = c("PB.1", "PB.2"),
#' Name = c("Gene1", "Gene2"))
#' rawDB <- compile_raw_db(transcript_file, abundance_file, gff_file, ORF_file)
#' DB <- process_db(rawDB, gene_ID_table)
#' plot_exon_correlations(DB)
#' }
#' @seealso Database
#' @import ggplot2
#' @export
plot_exon_correlations <- function(database, exon_filename, gene, weighted = T,
                                   exons_to_include = NULL, weights = NULL,
                                   plot_hist = F, symmetric = F) {
  seqs <- database$TranscriptDB[database$TranscriptDB$Gene == gene, ]
  # default: use read counts as weights (suggested)
  if (is.null(weights)) {
    weights <- seqs[, "FL_reads"]
  }
  seqs <- seqs[, c("PBID", "Transcript")]

  # exon file should be 2 columns: name/number and nucleotide sequence
  # columns should be separated by tabs
  # if it's a fasta file, a temp. file will be made
  if (!is_tsv_format(exon_filename, PBIDs = F)) {
    tmpfile <- get_tmp_tsv_file(exon_filename, exons = T)
    exons <- read.table(tmpfile, header = F, stringsAsFactors = F)
  } else {
    exons <- read.table(exon_filename, header = F, stringsAsFactors = F)
  }

  colnames(exons) <- c("ID", "Sequence")
  exons$Sequence <- sapply(exons$Sequence, toupper)
  exons$Start <- sapply(exons$Sequence, function(exon) {
    substr(exon, 1, min(nchar(exon), 30))
  })
  exons$End <- sapply(exons$Sequence, function(exon) {
    substr(exon, max(nchar(exon) - 30 + 1, 1), nchar(exon))
  })
  tmp <- colnames(seqs)

  # exon is "in transcript" if either first or last 30bp of exon are found
  # note: you'll miss exons if you have not genome-corrected transcripts
  for (i in seq_len(nrow(exons))) {
    seqs <- cbind(seqs, as.numeric(grepl(exons$Start[i], seqs$Transcript) |
                                     grepl(exons$End[i], seqs$Transcript)))
  }
  colnames(seqs) <- c(tmp, exons$ID)
  exon.calls <- seqs[, exons$ID]
  #row.names(exon.calls) <- seqs$PBID # maybe useful later
  if (!is.null(exons_to_include)) {
    exon.calls <- exon.calls[, exons_to_include]
  }
  # generate a [# exons by # exons] correlation matrix
  if (weighted) {
    cor.matrix <- weighted_corr(exon.calls, weights)
  } else {
    cor.matrix <- round(cor(exon.calls), 2)
  }

  if (!symmetric) {
    # only want lower triangle plotted, so insert NAs in upper triangle
    for (i in seq_len(nrow(cor.matrix))) {
      for (j in seq_len(ncol(cor.matrix))) {
        if (j > i) cor.matrix[i, j] <- NA
      }
    }
  }
  if (plot_hist) {
    title <- paste("Histogram of", gene, "Exon Correlations")
    if (weighted) title <- paste(title, "(Weighted)")
    hist(cor.matrix, breaks = seq(-1, 1, 0.1), col = "cornflowerblue",
         main = title, xlab = "Pearson Correlation")
  }
  cor.melt <- reshape2::melt(cor.matrix)

  title <- paste(gene, "Exon Inclusion")
  if (weighted) title <- paste(title, "(Weighted)")

  if (symmetric) {
    plt <- ggplot(data = cor.melt, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(Var1, Var2, label = value),
                color = "black", size = 5 - ncol(exon.calls) / 10, na.rm = T) +
      scale_fill_gradient2(name = "Pearson Correlation",
                           low = "blue", high = "red", mid = "white",
                           na.value = "white", limit = c(-1, 1)) +
      theme(axis.text.x = element_text(angle = 90),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank()) +
      ggtitle(title) +
      scale_y_discrete(labels = colnames(exon.calls)) + ylab("Exons") +
      scale_x_discrete(labels = colnames(exon.calls)) + xlab("Exons")

  } else {
    plt <- ggplot(data = cor.melt, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(Var1, Var2, label = value),
                color = "black", size = 5 - ncol(exon.calls) / 10, na.rm = T) +
      scale_fill_gradient2(name = "Pearson Correlation",
                           low = "blue", high = "red", mid = "white",
                           na.value = "white", limit = c(-1, 1)) +
      theme(axis.text.x = element_text(angle = 90),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            legend.justification = c(1, 0),
            legend.position = c(0.5, 0.8),
            legend.direction = "horizontal") +
      guides(fill = guide_colorbar(barwidth = 10, barheight = 1,
                                  title.position = "top", title.hjust = 0.5)) +
      ggtitle(title) +
      scale_y_discrete(labels = colnames(exon.calls)) + ylab("Exons") +
      scale_x_discrete(labels = colnames(exon.calls)) + xlab("Exons")
  }
  return(plt)
}
