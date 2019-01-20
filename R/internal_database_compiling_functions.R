# internal
add_abundance_data <- function(transcriptDB, abundance_filename) {
  print("Loading abundances...")
  # check if file contains header
  head <- system(command = paste("head -n20 ", abundance_filename, step = ""),
                 intern = T)
  has_header <- any(grepl("pbid[[:space:]]count_fl", head))

  abundances <- utils::read.table(abundance_filename, stringsAsFactors = F,
                                  header = has_header)
  abundances <- abundances[,c(1,2)]
  names(abundances) <- c("PBID", "FL_reads")
  # abundances$RawAbundance <- abundances$FL_reads
  transcriptDB <- transcriptDB %>% dplyr::left_join(abundances, by = c("PBID"))

  # Isoseq is set to minimum 2 FL reads, so add missing 1-read counts
  # should do this better in the future
  transcriptDB$FL_reads[is.na(transcriptDB$FL_reads)] <- 1
  return(transcriptDB)
}

# internal
get_gff_data <- function(gff_filepath) {
  # GFF2-compatible only
  print("Loading GFF data...")
  gff_data <- utils::read.table(gff_filepath, header = F, stringsAsFactors = F)
  # keep only the informative columns
  gff_data <- gff_data[,c(1,3:5,7,13)]
  names(gff_data) <- c("Chromosome", "TranscriptExon", "Start", "End",
                       "Strand", "PBID")
  return(gff_data)
}



# internal
make_ORF_db <- function(transcriptDB, aa_cutoff = 250) {
  if (!("ORF" %in% colnames(transcriptDB))) return(NULL)
  # TODO: optimize speed here
  orf_db <- transcriptDB[!is.na(transcriptDB$ORF),
                         c("PBID", "Chromosome", "Gene", "FL_reads", "ORF", "ORFLength")]
  counts <- as.data.frame(table(orf_db$ORF))
  duplicate_proteins <- as.character(counts$Var1[counts$Freq > 1])

  for (i in 1:length(duplicate_proteins)) {
    duplicates <- orf_db[orf_db$ORF == duplicate_proteins[i], ]

    PBIDs <- list(duplicates$PBID)
    merged <- data.frame(ORF        = duplicates$ORF[1],
                         ORFLength  = duplicates$ORFLength[1],
                         Gene       = sort(duplicates$Gene)[1],
                         FL_reads   = sum(as.numeric(duplicates$FL_reads)),
                         Chromosome = duplicates$Chromosome[1])
    merged$PBID <- PBIDs
    orf_db <- orf_db[which(orf_db$ORF != duplicate_proteins[i]), ]
    orf_db <- rbind(orf_db, merged)
  }
  if (max(orf_db$ORFLength) < aa_cutoff) return(NULL)
  # note minimum length cutoff
  orf_db <- orf_db[which(orf_db$ORFLength >= aa_cutoff),]
  orf_db <- calc_abundances(orf_db, genes_known = T)
  return(orf_db)
}

# internal
calc_abundances <- function(df, genes_known = F) {
  # here df is either a TranscriptDB or an OrfDB
  # TODO: optimize speed here
  new_db <- data.frame(matrix(ncol = ncol(df) + 2, nrow = 0))
  colnames(new_db) <- c(names(df), "PercAbundance", "CumPercAbundance")
  if (genes_known) {
    for (gene in unique(df$Gene)) {
      one_gene <- df[df$Gene == gene, ]
      one_gene$PercAbundance <- one_gene$FL_reads / sum(one_gene$FL_reads)
      one_gene <- one_gene[rev(order(one_gene$FL_reads)), ]
      one_gene$CumPercAbundance <- cumsum(one_gene$PercAbundance)
      new_db <- rbind(new_db, one_gene)
    }
  } else {
    # if genes are not known, then PBID prefixes are used as unique gene IDs
    # TODO: apply similar logic in other functions to give RawDBs functionality
    for (pref in unique(df$Prefix)) {
      one_gene <- df[df$Prefix == pref, ]
      one_gene$PercAbundance <- one_gene$FL_reads / sum(one_gene$FL_reads)
      one_gene <- one_gene[rev(order(one_gene$FL_reads)), ]
      one_gene$CumPercAbundance <- cumsum(one_gene$PercAbundance)
      new_db <- rbind(new_db, one_gene)
    }
  }
  return(new_db)
}

# internal
add_gene_names_to_transcripts <- function(transcriptDB, geneDB) {
  # only adds a gene name to a transcript if it doesn't already have one
  # so multiple sources can be used for inputting gene names iteratively
  if (is.null(transcriptDB$Gene)) {
    transcriptDB$Gene <- rep(NA, nrow(transcriptDB))
  }
  for (i in seq_len(length(transcriptDB$Prefix))) {
    name <- geneDB$Name[which(geneDB$ID == transcriptDB$Prefix[i])]
    if (is.na(transcriptDB$Gene[i])) {
      transcriptDB$Gene[i] <- ifelse(length(name) > 0, name, NA)
    }
  }
  return(transcriptDB)
}

# internal
add_gff_info <- function(transcriptDB, gffDB) {
  counts_table <- data.frame(table(gffDB$PBID))
  # subtracting 1 to not count "transcript" lines
  counts <- data.frame(PBID = as.character(counts_table$Var1),
                       ExonCount = counts_table$Freq - 1, stringsAsFactors = F)
  transcriptDB <- transcriptDB %>% dplyr::left_join(counts, by = "PBID")
  gff_transcript_lines <- gffDB[gffDB$TranscriptExon == "transcript", ]
  chromosomes <- data.frame(PBID = gff_transcript_lines$PBID,
                            Chromosome = gff_transcript_lines$Chromosome,
                            stringsAsFactors = F)
  return(transcriptDB %>% dplyr::left_join(chromosomes, by = "PBID"))
}

# internal
add_transcript_lengths <- function(transcriptDB, geneDB) {
  geneDB$AvgTranscriptLength <- sapply(geneDB$Name, function(x) {
    mean(transcriptDB$TranscriptLength[transcriptDB$Gene == x], na.rm = T)
  })
  geneDB$MaxTranscriptLength <- sapply(geneDB$Name, function(x) {
    max(transcriptDB$TranscriptLength[transcriptDB$Gene == x], na.rm = T)
  })
  # this is a weighted average, transcript != read
  geneDB$AvgReadLength <- sapply(geneDB$Name, function(x) {
    stats::weighted.mean(transcriptDB$TranscriptLength[transcriptDB$Gene == x],
                         transcriptDB$FL_reads[transcriptDB$Gene == x], na.rm = T)
  })
  return(geneDB)
}

# internal
add_isoform_info <- function(transcriptDB, geneDB) {
  geneDB$Exons <- sapply(geneDB$Name, function(x) {
    max(transcriptDB$ExonCount[transcriptDB$Gene == x])
  })
  geneDB$N50 <- sapply(geneDB$Name, function(x) {
    sum(transcriptDB$CumPercAbundance[transcriptDB$Gene == x] < 0.5) + 1
  })
  geneDB$N75 <- sapply(geneDB$Name, function(x) {
    sum(transcriptDB$CumPercAbundance[transcriptDB$Gene == x] < 0.75) + 1
  })
  if (length(unique(transcriptDB$Prefix)) == 1) {
    counts <- data.frame(ID = transcriptDB$Prefix[1],
                         Isoforms = nrow(transcriptDB))
  } else {
    counts <- data.frame(table(transcriptDB$Prefix))
  }
  names(counts) <- c("ID", "Isoforms")
  counts$ID <- as.character(counts$ID)
  return(geneDB %>% dplyr::left_join(counts, by = "ID"))
}

# internal
add_ORF_info <- function(transcriptDB, geneDB, orfDB = NULL) {
  if ("ORFLength" %in% colnames(transcriptDB)) {
    geneDB$AvgORFLength <- sapply(geneDB$Name, function(x) {
      mean(transcriptDB$ORFLength[transcriptDB$Gene == x], na.rm = T)
    })
    geneDB$MaxORFLength <- sapply(geneDB$Name, function(x) {
      max(transcriptDB$ORFLength[transcriptDB$Gene == x], na.rm = T)
    })
  }
  if (!is.null(orfDB)) {
    if (length(unique(orfDB$Gene)) == 1) {
      counts <- data.frame(Name = c(unique(orfDB$Gene)),
                           UniqueORFs = c(length(orfDB$Gene)))
    } else {
      counts <- data.frame(table(orfDB$Gene))
      names(counts) <- c("Name", "UniqueORFs")
    }
    counts$Name <- as.character(counts$Name)
    geneDB <- geneDB %>% dplyr::left_join(counts, by = "Name")
    geneDB$ORFN50 <- sapply(geneDB$Name, function(x) {
      sum(orfDB$CumPercAbundance[orfDB$Gene == x] < 0.5) + 1
    })
    geneDB$ORFN75 <- sapply(geneDB$Name, function(x) {
      sum(orfDB$CumPercAbundance[orfDB$Gene == x] < 0.75) + 1
    })
    geneDB$ORFsPerIsoform <- geneDB$UniqueORFs / geneDB$Isoforms
  }
  return(geneDB)
}

# internal
add_diversity_indeces <- function(geneDB, indeces = c("Gini", "Shannon"),
                                  transcriptDB = NULL, orfDB = NULL) {
  if (is.null(transcriptDB) & is.null(orfDB)) {
    stop("Error: no database provided.")
    return(NULL)
  }
  if (is.null(orfDB)) {
    DB <- transcriptDB
    Shannon_str <- "Shannon"
    Gini_str <- "Gini"
  } else {
    DB <- orfDB
    Shannon_str <- "ORFShannon"
    Gini_str <- "ORFGini"
  }

  diversities <- data.frame(Name = geneDB$Name, stringsAsFactors = FALSE)

  if ("Gini" %in% indeces & requireNamespace("ineq", quietly = TRUE)) {
    diversities[Gini_str] <- sapply(diversities$Name, function(gene) {
      ineq::ineq(DB$FL_reads[DB$Gene == gene], type = "Gini")
    })
  }
  if ("Shannon" %in% indeces & requireNamespace("vegan", quietly = TRUE)) {
    diversities[Shannon_str] <- c(sapply(diversities$Name, function(gene) {
      vegan::diversity(DB$FL_reads[DB$Gene == gene])
    }))
  }
  # only merge if packages loaded ok and at least one thing was calculated
  if (dim(diversities)[2] > 1) {
    geneDB <- geneDB %>% dplyr::left_join(diversities, by = "Name")
  }
  return(geneDB)
}

# internal
is_tsv_format <- function(filename, PBIDs = T) {
  if (PBIDs) {
    head <- system(command = paste("head -n1 ", filename, step = ""),
                    intern = T)
    return(grepl("PB\\.[0-9]+\\.[0-9]+[[:space:]][acgtnACGTN]+$", head))
  } else {
    return(length(system(command = paste("awk -F '\t' 'NF < 2' ", filename,
                                         step = ""), intern = T)) == 0)
  }
}

# internal
get_orf_sequences <- function(orf_filename) {
  if (!is_tsv_format(orf_filename)) {
    tmpfile <- get_tmp_tsv_file(orf_filename)
    orfDB <- read_tsv_file(tmpfile, get_prefix = T)
  } else {
    orfDB <- read_tsv_file(orf_filename, get_prefix = T)
  }
  colnames(orfDB) <- c("PBID", "Transcript", "TranscriptLength", "Prefix")
  return(orfDB)
}

# internal
add_ORF_sequences <- function(transcriptDB, ORF_filename) {
  # the ORF file comes from SQANTI's ORF predictions
  print("Loading ORFs...")
  if (!is_tsv_format(ORF_filename)) {
    tmpfile <- get_tmp_tsv_file(ORF_filename, ORF_file = T)
    ORFs <- read_tsv_file(tmpfile, get_prefix = F)
  } else {
    ORFs <- read_tsv_file(ORF_filename, get_prefix = F)
  }
  colnames(ORFs) <- c("PBID", "ORF", "ORFLength")

  transcriptDB <- transcriptDB %>% dplyr::left_join(ORFs, by = "PBID")
  return(transcriptDB)
}

# internal
get_transcript_sequences <- function(transcript_filename) {
  print("Loading sequences...")
  if (!is_tsv_format(transcript_filename)) {
    tmpfile <- get_tmp_tsv_file(transcript_filename)
    transcriptDB <- read_tsv_file(tmpfile, get_prefix = T)
  } else {
    transcriptDB <- read_tsv_file(transcript_filename, get_prefix = T)
  }
  colnames(transcriptDB) <- c("PBID", "Transcript", "TranscriptLength", "Prefix")
  return(transcriptDB)
}

read_tsv_file <- function(tsv_filename, get_prefix = T) {
  seqDB <- utils::read.table(tsv_filename, stringsAsFactors = F)
  names(seqDB) <- c("PBID", "Sequence")
  seqDB$Sequence <- sapply(seqDB$Sequence, function(x) {
    toupper(x)
  })
  seqDB$SeqLength <- sapply(seqDB$Sequence, function(x) {
    nchar(x)
  })
  if (get_prefix) {
    seqDB$Prefix <- sapply(seqDB$PBID, function(x) {
      stringr::str_match(x, "PB\\.[0-9]*")
    })
  }
  return(seqDB)
}

get_tmp_tsv_file <- function(filename, ORFs = F, exons = F) {
  old_name <- strsplit(basename(filename), split = ".",
                       fixed = T)[[1]]
  old_name <- paste(old_name[1:length(old_name) - 1], collapse = ".")
  d_name <- dirname(filename)

  if (ORFs) ext <- "ORFs"
  else {
    if (exons) ext <- "exons"
    else ext <- "transcripts"
  }
  tmpfile <- paste(d_name, "/", old_name, ".", ext, ".tsv", sep = "")

  if (!file.exists(tmpfile)) {
    # convert fasta/fastq to tsv format for R-compatible reading in
    # example tsv file line:
    # PB.1.1   ACGTCATCATCAT
    if (!exons) {
      system(command = paste("cat ", filename,
              " | sed -E 's/^.*(PB\\.[0-9]+\\.[0-9]+).*$/DELIM1\\1DELIM2/'",
              " | tr -d '\n'  | sed 's/DELIM1/\\'$'\n/g'",
              " | sed 's/DELIM2/\'$'\t/g' | tail -n +2 > ", tmpfile, sep = ""))
    } else {
      system(command = paste("cat ", filename,
              " | sed -E 's/^>(.*)$/DELIM1\\1DELIM2/'",
              " | tr -d '\n'  | sed 's/DELIM1/\\'$'\n/g'",
              " | sed 's/DELIM2/\'$'\t/g' | tail -n +2 > ", tmpfile, sep = ""))
    }
  }
  return(tmpfile)
}
