#' @export
kmer_tSNE <- function(database, kmer_counts = NULL, genes = NULL,
                      perplexity, dims = 2, iterations = 1000,
                      use_ORFs = F, verbose = T) {

  if (!(requireNamespace("Rtsne", quietly = T))) {
    stop("Error: you need the Rtsne package to perform t-SNE.")
    return(NULL)
  }
  # I require thee to input thy own choice of a perplexity value.
  print(paste("Using perplexity of ", perplexity, ".", sep = ""))

  # if kmer_counts aren't supplied (this uses all default settings!)
  # if kmer_counts aren't supplied, need to supply genes
  # otherwise will default to using entire dataset
  if (is.null(kmer_counts)) {
    if (is.null(genes)) genes <- database$GeneDB$Name
    kmer_counts <- get_kmer_counts(database, genes, use_ORFs = use_ORFs)
  }
  tsne <- Rtsne::Rtsne(kmer_counts, verbose = verbose, perplexity = perplexity,
                       max_iter = iterations, theta = 0, dims = dims)

  rownames(tsne$Y) <- rownames(kmer_counts)
  return(tsne)
}

#' @export
plot_tSNE <- function(database, tsne, use_ORFs = F, force_3D = F,
                      scale_by = NULL, insert_title = NULL) {

  if (!(ncol(tsne$Y) == 2 || ncol(tsne$Y) == 3)) {
    stop("t-SNE can only be plotted if dimensions are 2 or 3.")
    return(NULL)
  }
  if (is.null(insert_title)) {
    insert_title <- ifelse(use_ORFs, "ORF t-SNE Plot", "Isoform t-SNE Plot")
  }

  IDs <- rownames(tsne$Y)
  genes_long <- database$TranscriptDB$Gene[database$TranscriptDB$PBID %in% IDs]
  genes <- unique(genes_long)
  cols <- get_colors(length(genes) + 1)
  cols <- cols[1:length(genes)]

  if (dim(tsne$Y)[2] == 2) {
    df <- data.frame(Axis1 = tsne$Y[, 1], Axis2 = tsne$Y[, 2])
    if (force_3D) df$Axis3 <- rep(0, nrow(df))
  } else {
    df <- data.frame(Axis1 = tsne$Y[, 1], Axis2 = tsne$Y[, 2],
                                          Axis3 = tsne$Y[, 3])
  }
  threeD <- "Axis3" %in% colnames(df)
  df$Gene <- genes_long
  df$ID <- rownames(tsne$Y)

  if ("length" %in% scale_by | "abundance" %in% scale_by ) {
    if (use_ORFs) db <- database$OrfDB
    else db <- database$TranscriptDB

    isoform_subset <- sapply(db$Gene, function(gene_name) {
      return(gene_name %in% genes)
    })

    if ("length" %in% scale_by) {
      if (use_ORFs) seqs <- database$OrfDB$ORF[isoform_subset]
      else seqs <- database$TranscriptDB$Transcript[isoform_subset]
      df$Length <- sapply(seqs, function(x) nchar(x))

      if (threeD) {
        plt <- plotly::plot_ly(df, x = ~Axis1, y = ~Axis2, z = ~Axis3,
                               color = ~Gene, size = ~Length, colors = cols,   ######################
                               text = ~paste('ID:', ID)) %>%
               plotly::add_markers(opacity = 0.3)
      } else {
        print(ggplot(df, aes_string(x = "Axis1", y = "Axis2", colour = "Gene",
                                    size = "Length")) +
                geom_point(alpha = 0.3) +
                labs(title = insert_title) +
                theme_classic())
      }
    }
    if ("abundance" %in% scale_by) {
      df$Abundance <- db$FL_reads[isoform_subset]
      if (threeD) {
        plt <- plotly::plot_ly(df, x = ~Axis1, y = ~Axis2, z = ~Axis3,
                               color = ~Gene, size = ~Abundance, colors = cols,  ######################
                               text = ~paste('ID:', ID)) %>%
               plotly::add_markers(opacity = 0.3)
      } else {
        print(ggplot(df, aes_string(x = "Axis1", y = "Axis2", colour = "Gene",
                                    size = "Abundance")) +
                geom_point(alpha = 0.3) + scale_size(trans = "log2") +
                labs(title = insert_title) +
                theme_classic())
      }
    }
  } else {
    if (threeD) {
      plt <- plotly::plot_ly(df, x = ~Axis1, y = ~Axis2, z = ~Axis3, colors = cols, ######################
                             color = ~Gene, text = ~paste('ID:', ID)) %>%
             plotly::add_markers(opacity = 0.3)
    } else {
      print(ggplot(df, aes_string(x = "Axis1", y = "Axis2", colour = "Gene")) +
              geom_point(alpha = 0.3, size = 3) +
              labs(title = insert_title) +
              theme_classic())
    }
  }
  if (threeD) sink_var <- suppressWarnings(print(plt))
}
