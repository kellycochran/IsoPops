#' @export
plot_splicing <- function(database, gene, max_isoforms = -1,
                          apply_max_to_splicing = F) {
  if(!(requireNamespace("reshape", quietly = T) |
       requireNamespace("reshape2", quietly = T))) {
    stop("Error: reshape/reshape2 package required.")
  }

  pbids <- database$TranscriptDB$PBID[database$TranscriptDB$Gene == gene]

  if (max_isoforms < 1 | max_isoforms > length(pbids)) {
    max_isoforms <- length(pbids)
  }

  if (apply_max_to_splicing) pbids <- pbids[seq_len(max_isoforms)]

  df <- database$GffDB
  df <- df[df$TranscriptExon == "exon" & df$PBID %in% pbids, ]

  real_gene_len <- max(c(df$Start, df$End)) - min(c(df$Start, df$End))
  # 250 bp / pixel just happens to be a good scaling factor, by inspection
  scale <- real_gene_len / 250

  df$Start2 <- floor(df$Start / scale)
  df$End2 <- floor(df$End / scale)

  min_coord <- min(c(df$Start2, df$End2))
  max_coord <- max(c(df$Start2, df$End2))
  gene_len <- max_coord - min_coord + 1

  df_norm <- data.frame(PBID = df$PBID, End = df$End2, Start = df$Start2)
  df_norm$Start <- df$Start2 - min_coord + 1
  df_norm$End <- df$End2 - min_coord + 1

  junc_matrix <- matrix(data = 0, nrow = gene_len + 2, ncol = gene_len + 2)

  for (i in seq_len(gene_len)) {
    # for each splice acceptor site i
    donors <- which(df_norm$End == i)
    donors <- donors[df_norm$PBID[donors] == df_norm$PBID[donors + 1]]
    donors <- stats::na.omit(donors)
    if (length(donors > 0)) {
      acceptor_coords <- df_norm$Start[donors + 1]
      acceptor_coords <- gene_len - acceptor_coords
      for (index in seq_len(length(acceptor_coords))) {
        j <- acceptor_coords[index]
        for (k in j:gene_len) {
          if (i + k <= gene_len) {
            junc_matrix[i + 1, k + 1] <- junc_matrix[i + 1, k + 1] + 1
          }
        }
        for (l in (i+1):gene_len) {
          if (j + l <= gene_len) {
            junc_matrix[l + 1, j + 1] <- junc_matrix[l + 1, j + 1] + 1
          }
        }
      }
    }
  }
  # TODO: fix this explicit dependency
  if (requireNamespace("reshape2", quietly = T)) {
    melted <- reshape2::melt(junc_matrix)
  } else {
    melted <- reshape::melt(junc_matrix)
  }
  colnames(melted) <- c("x", "y", "fill")
  melted$fill <- sqrt(melted$fill)

  if (max_isoforms == 1) {
    pbid <- pbids[1]
    df_norm <- df_norm[df_norm$PBID == pbid, ]
  } else {
    pbids <- pbids[seq_len(max_isoforms)]
    df_norm <- df_norm[df_norm$PBID %in% pbids, ]
    df_norm <- df_norm[order(match(df_norm$PBID, pbids)), ]
  }

  offset <- 5
  rect_data <- get_all_rect_data(df_norm, gene_len, offset)
  line_data <- get_all_line_data(df_norm, gene_len, offset)

  plt <- ggplot() +
    geom_tile(data = melted, aes_string(x = "x", y = "y", fill = "fill")) +
    theme_void() +
    scale_fill_gradientn(colours = c("white", "red"), na.value = "white") +
    scale_x_reverse() + coord_fixed(ratio = 1) +
    theme(legend.position = "none") +
    geom_line(data = line_data, aes_string(x = "x", y = "y", group = "groups"),
              size = 0.25) +
    geom_polygon(data = rect_data, aes_string(x = "x", y = "y",
                                              group = "groups"))

  grid::grid.newpage()
  suppressWarnings(print(plt, vp = grid::viewport(angle = 135, clip = "off")))
}



# internal
get_end_of_isoform <- function(coords) {
  if (nrow(coords) == 0 || length(unique(coords$PBID)) == 1) return(-1)
  i <- 1
  while (i < nrow(coords)) {
    if (coords$Start[i] > coords$Start[i + 1]) return(i)
    i <- i + 1
  }
  return(-1)
}

# internal
get_all_rect_data <- function(all_exons, gene_len, offset) {
  # TODO: clean this up
  coords <- data.frame(Start = all_exons$Start, End = all_exons$End)

  coords$Start[coords$End - coords$Start < 2] <- coords$Start[coords$End - coords$Start < 2] - 0.2
  coords$End[coords$End - coords$Start < 2] <- coords$End[coords$End - coords$Start < 2] + 0.2

  end_of_isoform <- get_end_of_isoform(coords)

  if (end_of_isoform == -1) {
    rect_data <- get_rect_data(coords, gene_len, offset)
  } else {
    rect_data <- get_rect_data(coords[seq_len(end_of_isoform), ], gene_len, offset)
    coords <- coords[-1 * seq_len(end_of_isoform), ]
    end_of_isoform <- get_end_of_isoform(coords)
    i <- 1
    while (end_of_isoform != -1) {
      new_rect_data <- get_rect_data(coords[seq_len(end_of_isoform), ], gene_len, offset)
      new_rect_data$x <- new_rect_data$x + offset * i
      new_rect_data$y <- new_rect_data$y + offset * i
      new_rect_data$groups <- new_rect_data$groups + i
      rect_data <- rbind(rect_data, new_rect_data)
      coords <- coords[-1 * seq_len(end_of_isoform), ]
      end_of_isoform <- get_end_of_isoform(coords)
      i <- i + 1
    }
    new_rect_data <- get_rect_data(coords, gene_len, offset)
    new_rect_data$x <- new_rect_data$x + offset * i
    new_rect_data$y <- new_rect_data$y + offset * i
    new_rect_data$groups <- new_rect_data$groups + i
    rect_data <- rbind(rect_data, new_rect_data)
  }
  rect_data$x <- rect_data$x - offset / 4
  rect_data$y <- rect_data$y - offset / 4
  return(rect_data)
}

# internal
get_rect_data <- function(coords, gene_len, offset) {
  l <- apply(coords, 1, function(exon) ex_to_rect_coords(exon, gene_len, offset))
  rect_x <- unlist(lapply(l, function(l) l[[1]]))
  rect_y <- unlist(lapply(l, function(l) l[[2]]))
  groups <- sort(rep(1, 5))
  return(data.frame(x = rect_x, y = rect_y, groups = groups))
}

# internal
ex_to_rect_coords <- function(exon_coords, gene_len, offset) {
  return(list(x = c(exon_coords[1] + offset / 2,
                    exon_coords[1] + offset,
                    exon_coords[2] + offset,
                    exon_coords[2] + offset / 2,
                    exon_coords[1] + offset / 2),
              y = c(gene_len - exon_coords[1] + offset / 2,
                    gene_len - exon_coords[1] + offset,
                    gene_len - exon_coords[2] + offset,
                    gene_len - exon_coords[2] + offset / 2,
                    gene_len - exon_coords[1] + offset / 2)))
}

# internal
get_all_line_data <- function(all_exons, gene_len, offset) {
  coords <- data.frame(Start = all_exons$Start, End = all_exons$End)

  coords$Start[coords$End - coords$Start < 2] <- coords$Start[coords$End - coords$Start < 2] - 0.2
  coords$End[coords$End - coords$Start < 2] <- coords$End[coords$End - coords$Start < 2] + 0.2

  end_of_isoform <- get_end_of_isoform(coords)

  if (end_of_isoform == -1) {
    line_data <- get_line_data(coords, gene_len)
  } else {
    line_data <- get_line_data(coords[seq_len(end_of_isoform), ], gene_len)
    coords <- coords[-1 * seq_len(end_of_isoform), ]
    end_of_isoform <- get_end_of_isoform(coords)
    i <- 1
    while (end_of_isoform != -1) {
      new_line_data <- get_line_data(coords[seq_len(end_of_isoform), ], gene_len)
      new_line_data$x <- new_line_data$x + offset * i
      new_line_data$y <- new_line_data$y + offset * i
      new_line_data$groups <- new_line_data$groups + i
      line_data <- rbind(line_data, new_line_data)
      coords <- coords[-1 * seq_len(end_of_isoform), ]
      end_of_isoform <- get_end_of_isoform(coords)
      i <- i + 1
    }
    new_line_data <- get_line_data(coords, gene_len)
    new_line_data$x <- new_line_data$x + offset * i
    new_line_data$y <- new_line_data$y + offset * i
    new_line_data$groups <- new_line_data$groups + i
    line_data <- rbind(line_data, new_line_data)
  }
  line_data$x <- line_data$x + offset / 2
  line_data$y <- line_data$y + offset / 2
  return(line_data)
}

# internal
get_line_data <- function(coords, gene_len) {
  line <- c(min(coords), max(coords))
  return(data.frame(x = line, y = gene_len - line, groups = c(1,1)))
}
