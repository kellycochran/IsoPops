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
get_PCA_dtm <- function(raw_sequences, n, length_norm = T) {
  if (!requireNamespace("text2vec", quietly = TRUE)) {
    stop("Error: missing the text2vec package.")
    return(NULL)
  } else {
    sequences <- as.character(sapply(raw_sequences, string_to_seq))
    it <- text2vec::itoken(sequences, preprocessor = identity,
                           tokenizer = text2vec::word_tokenizer,
                           progressbar = F)
    vocab <- text2vec::create_vocabulary(it, ngram = c(n, n))
    dtm <- text2vec::create_dtm(it, text2vec::vocab_vectorizer(vocab))
    if (length_norm) {
      dtm <- data.frame(t(apply(dtm, 1, function(dtm_row) dtm_row / sum(dtm_row))))
    } else {
      # need the transformation to data frame anyways
      dtm <- data.frame(t(apply(dtm, 1, function(dtm_row) dtm_row)))
    }
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
kmer_PCA <- function(database, gene_list = database$GeneDB$Name,
                     use_ORFs = F, k = 6, length_normalize = T) {
  if (use_ORFs) {
    isoform_subset <- database$OrfDB$Gene %in% gene_list
    seqs <- database$OrfDB$ORF[isoform_subset]
  } else {
    isoform_subset <- which(database$TranscriptDB$Gene %in% gene_list)
    seqs <- database$TranscriptDB$Transcript[isoform_subset]
  }
  dtm <- get_PCA_dtm(seqs, k, length_normalize)
  pca <- stats::prcomp(dtm, scale = F)
  if (!use_ORFs) {
    rownames(pca$x) <- database$TranscriptDB$PBID[isoform_subset]
  }

  eigenvalues <- (pca$sdev ^ 2)[1:min(length(gene_list) * 3, 5)]

  graphics::plot(eigenvalues, main = "Variances Captured by PCs",
                 xlab = "Principal Component", ylab = "Variance (Eigenvalue)",
                 cex = exp(-length(gene_list)) * 5 + 1, pch = 19,
                 col = "darkblue")

  if (use_ORFs) {
    db <- database$OrfDB
  } else {
    db <- database$TranscriptDB
  }
  gene_list_long <- db$Gene[db$Gene %in% gene_list]

  pca_list <- list(pca = pca, genes = gene_list, use_ORFs = use_ORFs,
                   gene_list_long = gene_list_long)
  return(pca_list)
}

#' @importFrom magrittr %>%
#' @export
cluster_PCA <- function(pca_list, num_PCs = -1,
                        cluster_method = "complete") {
  # TODO: change from frantically debugged code to actual decent code
  gene_list <- pca_list$genes
  pca <- pca_list$pca
  use_ORFs <- pca_list$use_ORFs
  gene_list_long <- pca_list$gene_list_long

  if (num_PCs < 1) {
    num_PCs <- ceiling(log(length(gene_list), base = 2))
  }
  print(paste("Using ", num_PCs, " principal components.", sep = ""))

  pca_data <- data.frame(t(pca$x[, num_PCs]))
  colnames(pca_data) <- gene_list_long
  clust <- stats::hclust(stats::as.dist(L2_dist(t(pca_data))),
                         method = cluster_method)

  dend <- stats::as.dendrogram(clust)
  labs <- as.numeric(as.factor(gene_list_long)) + 1
  labs_ordered <- labs[stats::order.dendrogram(dend)]
  genes_ordered <- gene_list_long[stats::order.dendrogram(dend)]
  dendextend::labels_colors(dend) <- labs_ordered
  dend %>% dendextend::set("labels_cex", 0.5) %>%
           dendextend::set("labels", genes_ordered) %>%
          graphics::plot(main = "PCA Clustering Dendrogram",
                         ylab = "Euclidian Distance")

  CH(clust, max(ceiling(length(gene_list) * 2), 10), pca$x[, num_PCs])

  clust_list <- list(clustering = clust, genes = gene_list,
                     gene_list_long = gene_list_long)
  return(clust_list)
}

#' @importFrom magrittr %>%
#' @export
eval_clusters <- function(clust_list, num_clusters = length(clust_list$genes)) {
  print(paste("Number of clusters chosen: ", num_clusters, sep = ""))
  clusters <- stats::cutree(clust_list$clustering, k = num_clusters)
  print(table(clusters, clust_list$gene_list_long))

  if (requireNamespace("dendextend", quietly = TRUE)) {
    dend <- stats::as.dendrogram(clust_list$clustering)
    labs <- as.numeric(as.factor(clust_list$gene_list_long)) + 1
    labs_ordered <- labs[stats::order.dendrogram(dend)]
    genes_ordered <- clust_list$gene_list_long[stats::order.dendrogram(dend)]
    dendextend::labels_colors(dend) <- labs_ordered
    dend %>% dendextend::set("labels_cex", 0.5) %>%
             dendextend::set("labels", genes_ordered) %>%
             dendextend::color_branches(k = num_clusters) %>%
             graphics::plot(main = "PCA Clustering Dendrogram",
                            ylab = "Euclidian Distance")
  } else {
    stop("Error: missing dendextend package.")
  }
}


#' @export
plot_PCA <- function(pca_list, database, scale_by = NULL,
                     insert_title = NULL) {
  # TODO: this whole thing needs a major +/- ORFs refactor
  # TODO: less temporarily fix
  # remove: gene_list <- pca_list$genes
  pca <- pca_list$pca
  use_ORFs <- pca_list$use_ORFs
  gene_list_long <- pca_list$gene_list_long

  if (is.null(insert_title)) {
    insert_title <- ifelse(use_ORFs, "ORF PCA Plot", "Isoform PCA Plot")
  }

  df <- data.frame(pca$x[, 1:2])
  colnames(df) <- c("PC1", "PC2")
  df$Gene <- gene_list_long

  gene_list <- unique(gene_list_long)

  if ("length" %in% scale_by) {
    if (use_ORFs) {
      isoform_subset <- sapply(database$OrfDB$Gene, function(gene_name) {
        return(gene_name %in% gene_list)
      })
      seqs <- database$OrfDB$ORF[isoform_subset]
    } else {
      isoform_subset <- sapply(database$TranscriptDB$Gene, function(gene_name) {
        return(gene_name %in% gene_list)
      })
      seqs <- database$TranscriptDB$Transcript[isoform_subset]
    }
    df$Length <- sapply(seqs, function(x) nchar(x))

    ggplot(df, aes_string(x = "PC1", y = "PC2", colour = "Gene",
                          size = "Length")) + geom_point(alpha = 0.3) +
      labs(title = insert_title) + theme_classic()
  } else {
    if ("abundance" %in% scale_by) {
      if (use_ORFs) {
        db <- database$OrfDB
      } else {
        db <- database$TranscriptDB
      }
      isoform_subset <- sapply(db$Gene, function(gene_name) {
        return(gene_name %in% gene_list)
      })
      df$Abundance <- db$FL_reads[isoform_subset]

      ggplot(df, aes_string(x = "PC1", y = "PC2", colour = "Gene",
                            size = "Abundance")) +
        geom_point(alpha = 0.3) + scale_size(trans = "log2") +
        labs(title = insert_title) + theme_classic()
    } else {
      ggplot(df, aes_string(x = "PC1", y = "PC2", colour = "Gene")) +
        geom_point(alpha = 0.3, size = 3) + labs(title = insert_title) +
        theme_classic()
    }
  }
}

#' @export
plot_3D_PCA <- function(pca_list, database, scale_by = NULL,
                        insert_title = NULL) {
  if (!(requireNamespace("rgl", quietly = T) &
        requireNamespace("car", quietly = T))) {
    stop("Error: you need the rgl and car packages for 3D PCA plots.")
    return(NULL)
  }
  pca <- pca_list$pca
  use_ORFs <- pca_list$use_ORFs
  gene_list_long <- pca_list$gene_list_long
  gene_list <- unique(gene_list_long)

  if (is.null(insert_title)) {
    insert_title <- ifelse(use_ORFs, "ORF PCA Plot", "Isoform PCA Plot")
  }

  if ("length" %in% scale_by) {
    if (use_ORFs) {
      isoform_subset <- sapply(database$OrfDB$Gene, function(gene_name) {
        return(gene_name %in% unique(gene_list_long))
      })
      seqs <- database$OrfDB$ORF[isoform_subset]
    } else {
      isoform_subset <- sapply(database$TranscriptDB$Gene, function(gene_name) {
        return(gene_name %in% unique(gene_list_long))
      })
      seqs <- database$TranscriptDB$Transcript[isoform_subset]
    }
    seq_lengths <- sapply(seqs, function(x) nchar(x))

    car::scatter3d(x = pca$x[, 1], y = pca$x[, 2], z = pca$x[, 3],
                   groups = as.factor(gene_list_long),
                   radius = seq_lengths,
                   grid = FALSE, surface = FALSE)

  } else {
    if ("abundance" %in% scale_by) {
      if (use_ORFs) {
        db <- database$OrfDB
      } else {
        db <- database$TranscriptDB
      }
      isoform_subset <- sapply(db$Gene, function(gene_name) {
        return(gene_name %in% gene_list)
      })
      abundances <- db$FL_reads[isoform_subset]

      car::scatter3d(x = pca$x[, 1], y = pca$x[, 2], z = pca$x[, 3],
                     groups = as.factor(gene_list_long),
                     radius = log10(abundances),
                     grid = FALSE, surface = FALSE)
    } else {
      car::scatter3d(x = pca$x[, 1], y = pca$x[, 2], z = pca$x[, 3],
                     groups = as.factor(gene_list_long),
                     grid = FALSE, surface = FALSE)
    }
  }
}
