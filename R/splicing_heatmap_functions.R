#' @export
plot_splicing <- function(database, gene, max_isoforms = -1, zoom_in = c(0, 1),
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
  df_lines <- df

  if (length(zoom_in) != 2 || any(zoom_in < 0)) {
    warning("Zoom input in incorrect format, defaulting to no zoom.")
  }
  if (length(zoom_in) == 2 && !(all(zoom_in == c(0, 1)))) {
    max_coord <- max(c(df$Start, df$End))
    min_coord <- min(c(df$Start, df$End))
    if (all(zoom_in <= 1)) { # interpret zoom input as % of full gene view
      zoom_in[1] <- min_coord + round(zoom_in[1] * (max_coord - min_coord))
      zoom_in[2] <- min_coord + round(zoom_in[2] * (max_coord - min_coord))
    }
    min_limit <- min(zoom_in)
    max_limit <- max(zoom_in)

    if (all(zoom_in == as.integer(zoom_in))) {

      if (min_limit > min_coord && min_limit < max_coord) {
        df <- df[df$Start >= min_limit | df$End >= min_limit, ]
        df$Start[df$Start < min_limit] <- min_limit
        df$End[df$End < min_limit] <- min_limit
        df_lines$Start[df_lines$Start < min_limit] <- min_limit
        df_lines$End[df_lines$End < min_limit] <- min_limit
      } else {
        if (!(min_limit == min_coord)) {
          warning("Zoom coordinates out of range, defaulting to no zoom.")
        }
      }
      if (max_limit < max_coord && max_limit > min_coord) {
        df <- df[df$Start <= max_limit | df$End <= max_limit, ]
        df$Start[df$Start > max_limit] <- max_limit
        df$End[df$End > max_limit] <- max_limit
        df_lines$Start[df_lines$Start > max_limit] <- max_limit
        df_lines$End[df_lines$End > max_limit] <- max_limit
      } else {
        if (!(max_limit == max_coord)) {
          warning("Zoom coordinates out of range, defaulting to no zoom.")
        }
      }
    } else {
      warning("Zoom input in incorrect format, defaulting to no zoom.")
    }
  }

  real_gene_len <- max(c(df$Start, df$End)) - min(c(df$Start, df$End))
  # 250 bp / pixel at no zoom just happens to be a good scaling factor
  scale <- real_gene_len / 250

  df$Start2 <- floor(df$Start / scale)
  df$End2 <- floor(df$End / scale)

  min_coord_norm <- min(c(df$Start2, df$End2))
  max_coord_norm <- max(c(df$Start2, df$End2))
  gene_len <- max_coord_norm - min_coord_norm + 1

  df_norm <- data.frame(PBID = df$PBID, End = df$End2, Start = df$Start2)
  df_norm$Start <- df$Start2 - min_coord_norm + 1
  df_norm$End <- df$End2 - min_coord_norm + 1
  df_lines$Start <- floor(df_lines$Start / scale) - min_coord_norm + 1
  df_lines$End <- floor(df_lines$End / scale) - min_coord_norm + 1

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
    df_lines <- df_lines[df_lines$PBID == pbid, ]
  } else {
    pbids <- pbids[seq_len(max_isoforms)]
    df_norm <- df_norm[df_norm$PBID %in% pbids, ]
    df_norm <- df_norm[order(match(df_norm$PBID, pbids)), ]
    df_lines <- df_lines[df_lines$PBID %in% pbids, ]
    df_lines <- df_lines[order(match(df_lines$PBID, pbids)), ]
  }

  offset <- 5 # would be good not to need this
  rect_data <- get_all_rect_data(df_norm, gene_len, offset, pbids)
  line_data <- get_all_line_data(df_lines, gene_len, offset)


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

get_end_of_isoform <- function(coords) {
  if (nrow(coords) == 0 || length(unique(coords$PBID)) == 1) return(-1)
  i <- 1
  while (i < nrow(coords)) {
    if (coords$PBID[i] != coords$PBID[i + 1]) return(i)
    i <- i + 1
  }
  return(-1)
}

# internal
get_all_rect_data <- function(all_exons, gene_len, offset, pbids) {
  # TODO: clean this up
  coords <- all_exons
  # TODO: check this is ok
  coords$Start[coords$End - coords$Start < 2] <- coords$Start[coords$End - coords$Start < 2] - 0.2
  coords$End[coords$End - coords$Start < 2] <- coords$End[coords$End - coords$Start < 2] + 0.2

  rect_data <- data.frame(x = numeric(), y = numeric(), groups = numeric())
  i <- 0
  for (pbid in pbids) {
    exons <- coords[coords$PBID == pbid, c("Start", "End")]
    if (nrow(exons) > 0) {
      new_rect_data <- get_rect_data(exons, gene_len, offset)
      new_rect_data$x <- new_rect_data$x + offset * i
      new_rect_data$y <- new_rect_data$y + offset * i
      new_rect_data$groups <- new_rect_data$groups + i
      rect_data <- rbind(rect_data, new_rect_data)
    }
    i <- i + 1
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
  groups <- rep(1, 5)
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
  coords <- all_exons

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
  coords <- coords[, c("Start", "End")]
  line <- c(min(coords), max(coords))
  return(data.frame(x = line, y = gene_len - line, groups = c(1,1)))
}
