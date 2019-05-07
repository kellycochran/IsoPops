# internal
L2_dist <- function(x) {
  m <- as.matrix(x)
  len <- nrow(m)
  msq <- matrix(rowSums(m^2), nrow = len, ncol = len)
  msq2 <- matrix(rowSums(m^2), nrow = len, ncol = len, byrow = TRUE)
  mxm <- m %*% t(m)
  result <- msq + msq2 - 2 * mxm
  diag(result) <- 0
  result[result < 0] <- 0
  return(result)
}

# internal
seq_to_string <- function(seq) {
  return(paste(as.character(strsplit(seq, " ")[[1]]), sep = "", collapse = ""))
}

# internal
string_to_seq <- function(seq) {
  return(paste(as.character(strsplit(seq, "")[[1]]), sep = "", collapse = " "))
}

#' @export
get_kmer_counts <- function(database, genes = database$GeneDB$Name,
                            kmer_size = 6, use_ORFs = F, length_normalize = T) {
  print("Creating k-mer count matrix.") # might take a while for all genes

  if (!requireNamespace("text2vec", quietly = TRUE)) {
    stop("Error: missing the text2vec package.")
    return(NULL)
  } else {
    if (use_ORFs) {
      isoform_subset <- database$OrfDB$Gene %in% genes
      # if ORF came from multiple transcripts, just use first transcript ID
      IDs <- sapply(database$OrfDB$PBID[isoform_subset], function(x) {
        x[[1]][[1]] })
      raw_sequences <- database$OrfDB$ORF[isoform_subset]
    } else {
      isoform_subset <- database$TranscriptDB$Gene %in% genes
      IDs <- database$TranscriptDB$PBID[isoform_subset]
      raw_sequences <- database$TranscriptDB$Transcript[isoform_subset]
    }
    sequences <- as.character(sapply(raw_sequences, string_to_seq))

    it <- text2vec::itoken(sequences, preprocessor = identity,
                           tokenizer = text2vec::word_tokenizer,
                           progressbar = F)
    vocab <- text2vec::create_vocabulary(it, ngram = c(kmer_size, kmer_size))
    dtm <- text2vec::create_dtm(it, text2vec::vocab_vectorizer(vocab))

    if (length_normalize) {
      dtm <- data.frame(t(apply(dtm, 1, function(row) row / sum(row))))
    }
    rownames(dtm) <- IDs
    return(dtm)
  }
}

# function for telling what the right number of clusters might be
# TODO: maybe switch to clustering stability indeces
#' @export
CH <- function(cluster, numClusters, original) {
  if (numClusters < 2) {  # not defined for a single cluster
    return(NULL)
  }
  if (!requireNamespace("clusterCrit", quietly = TRUE)) {
    stop("Error: missing the clusterCrit package.")
    return(NULL)
  }
  indeces <- rep(0, numClusters)
  for(i in 2:numClusters) {
    treecut <- stats::cutree(cluster, k = i)
    d_matrix <- data.matrix(original)
    crit <- clusterCrit::intCriteria(d_matrix, treecut, "Calinski_Harabasz")
    indeces[i] <- crit[[1]]
  }
  graphics::plot(2:length(indeces), indeces[-1],
                 main = "CH Index vs. Number of Clusters",
                 xlab = "Clusters", ylab = "Calinski-Harabasz Index",
                 col = grDevices::rgb(0.2, 0, 0.4, 0.8), pch = 19,
                 cex = exp(-numClusters) * 5 + 1)

  maximum <- max(indeces, na.rm = T)
  return(maximum)
}

#' @export
kmer_PCA <- function(database, kmer_counts = NULL, genes = NULL,
                     use_ORFs = F) {

  # if kmer_counts aren't supplied (this uses all default settings!)
  # if kmer_counts aren't supplied, need to supply genes
  # otherwise will default to using entire dataset
  if (is.null(kmer_counts)) {
    if (is.null(genes)) genes <- database$GeneDB$Name
    kmer_counts <- get_kmer_counts(database, genes, use_ORFs = use_ORFs)
  }
  IDs <- rownames(kmer_counts)
  isoform_subset <- database$TranscriptDB$PBID %in% IDs
  genes_long <- unique(database$TranscriptDB$Gene[isoform_subset])
  genes <- unique(genes_long)

  pca <- stats::prcomp(kmer_counts, scale = F)
  # in documentation, explain how to lookup points to find PBID
  rownames(pca$x) <- rownames(kmer_counts)

  eigenvalues <- (pca$sdev ^ 2)[1:max(ceiling(length(genes) / 1.5), 5)]

  graphics::plot(eigenvalues, main = "Variances Captured by PCs",
                 xlab = "Principal Component", ylab = "Variance (Eigenvalue)",
                 cex = exp(-length(genes)) * 5 + 1, pch = 19,
                 col = "darkblue")
  return(pca)
}

#' @importFrom magrittr %>%
#' @export
cluster_isoforms <- function(database, kmer_counts = NULL,
                             genes = NULL, use_ORFs = F,
                             num_clusters = -1, cut_height = -1,
                             cluster_method = "complete", return = F) {

  # if kmer_counts aren't supplied (this uses all default settings!)
  if (is.null(kmer_counts)) {
    if (is.null(genes)) genes <- database$GeneDB$Name
    kmer_counts <- get_kmer_counts(database, genes, use_ORFs = use_ORFs)
  }
  IDs <- rownames(kmer_counts)
  isoform_subset <- database$TranscriptDB$PBID %in% IDs
  genes_long <- database$TranscriptDB$Gene[isoform_subset]
  genes <- unique(genes_long)

  clust <- stats::hclust(stats::as.dist(L2_dist(kmer_counts)),
                         method = cluster_method)
  dend <- stats::as.dendrogram(clust)

  if (length(genes) > 1) {
    labs <- as.numeric(as.factor(genes_long)) + 1
    labs_ordered <- labs[stats::order.dendrogram(dend)]
    dendextend::labels_colors(dend) <- labs_ordered
  }
  IDs_ordered <- IDs[stats::order.dendrogram(dend)]

  dend <- dend %>% dendextend::set("labels_cex", 0.5) %>%
    dendextend::set("labels", IDs_ordered)

  title <- paste(ifelse(use_ORFs, "ORF", "Isoform"), "Clustering Dendrogram")

  # optional: color branches of dendrogram by # of clusters or a cut height
  # TODO: actually test cut_height parts
  if (num_clusters > 1 && cut_height > 0) {
    warning("Incompatible args num_clusters and cut_height (don't set both).")
  }
  if (num_clusters > 1 && cut_height < 0) {
    dend <- dend %>% dendextend::color_branches(k = num_clusters)
    if (length(unique(genes)) > 1) {
      clusters <- stats::cutree(clust, k = num_clusters)
      print("Gene Separation By Cluster:")
      Genes <- genes_long # just to fix header of table
      print(table(clusters, Genes))
    }
  }
  if (num_clusters < 2 && cut_height > 0) {
    dend <- dend %>% dendextend::color_branches(h = cut_height)
    if (length(unique(genes)) > 1) {
      clusters <- stats::cutree(clust, h = cut_height)
      print("Gene Separation By Cluster:")
      Genes <- genes_long # just to fix header of table
      print(table(clusters, Genes))
    }
  }
  graphics::plot(dend, main = title, ylab = "Euclidian Distance")

  if (num_clusters < 2 && cut_height < 0) {
    CH(clust, max(ceiling(length(genes) * 2), 10), kmer_counts)
  }
  if (return) return(clust)
}

#' @export
plot_PCA <- function(database, pca, use_ORFs = F, scale_by = NULL,
                     insert_title = NULL) {

  IDs <- rownames(pca$x)
  genes_long <- database$TranscriptDB$Gene[database$TranscriptDB$PBID %in% IDs]
  genes <- unique(genes_long)

  if (is.null(insert_title)) {
    insert_title <- paste(ifelse(use_ORFs, "ORF", "Isoform"), "PCA Plot")
  }

  df <- data.frame(pca$x[, 1:2])
  #colnames(df) <- c("PC1", "PC2")
  df$Gene <- genes_long

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
      # currently this only shows up if you print it...
      # ...yeah...that's not weird at all
      print(ggplot(df, aes_string(x = "PC1", y = "PC2", colour = "Gene",
                                  size = "Length")) +
              geom_point(alpha = 0.3) +
              labs(title = insert_title) + theme_classic())
    }
    if ("abundance" %in% scale_by) {
      df$Abundance <- db$FL_reads[isoform_subset]

      ggplot(df, aes_string(x = "PC1", y = "PC2", colour = "Gene",
                            size = "Abundance")) +
        geom_point(alpha = 0.3) + scale_size(trans = "log2") +
        labs(title = insert_title) + theme_classic()
    }
  } else {
    ggplot(df, aes_string(x = "PC1", y = "PC2", colour = "Gene")) +
      geom_point(alpha = 0.3, size = 3) +
      labs(title = insert_title) + theme_classic()
  }
}


#' @export
plot_3D_PCA <- function(database, pca, use_ORFs = F, scale_by = NULL,
                        insert_title = NULL) {

  if (!(requireNamespace("plotly", quietly = T))) {
    stop("Error: you need the plotly package for 3D PCA plots.")
    return(NULL)
  }
  if (is.null(insert_title)) {
    insert_title <- paste(ifelse(use_ORFs, "ORF", "Isoform"), "PCA Plot")
  }

  IDs <- rownames(pca$x)
  genes_long <- database$TranscriptDB$Gene[database$TranscriptDB$PBID %in% IDs]
  genes <- unique(genes_long)
  cols <- get_colors(length(genes) + 1)
  cols <- cols[1:length(genes)]

  df <- data.frame(pca$x[, 1:3])
  df$Gene <- genes_long
  df$ID <- IDs

  if ("length" %in% scale_by | "abundance" %in% scale_by ) {
    if (use_ORFs) db <- database$OrfDB
    else db <- database$TranscriptDB

    isoform_subset <- sapply(db$Gene, function(gene_name) {
      return(gene_name %in% genes)
    })

    if ("length" %in% scale_by) {
      if (use_ORFs) seqs <- database$OrfDB$ORF
      else seqs <- database$TranscriptDB$Transcript
      seqs <- seqs[isoform_subset]

      df$Length <- sapply(seqs, function(x) nchar(x))
      plt <- plotly::plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Gene,
                             colors = cols,  ######################
                             size = ~Length, text = ~paste('ID:', ID)) %>%
                             plotly::add_markers(opacity = 0.3)
    }
    if ("abundance" %in% scale_by) {
      df$Abundance <- db$FL_reads[isoform_subset]
      plt <- plotly::plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Gene,
                             colors = cols,  ######################
                             size = ~Abundance, text = ~paste('ID:', ID)) %>%
                             plotly::add_markers(opacity = 0.3)
    }
  } else {
    plt <- plotly::plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Gene,
                           colors = cols,  ######################
                           text = ~paste('ID:', ID)) %>%
                           plotly::add_markers(opacity = 0.3)
  }
  sink_var <- suppressWarnings(print(plt))
}
