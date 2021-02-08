################################################################################
# 20200730_m6a_project_functions.R
# Author: Mitchel Cole
# Source file for all other R/Rmd files in m6A prediction project
################################################################################

library(tidyverse)
library(plyranges)
library(pander)
library(seqTools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggpubr)

stretch_width <- 201

writeHomerInput <- function(x, y){
  #' Writes granges bed (x) to filename (y) as Homer peak file
  #' @param x Granges object
  #' @param y filename
  #' @return NA
  
  x %>%
    anchor_center() %>%
    mutate(width = stretch_width) %>%
    as_tibble() %>%
    mutate(id = sprintf("%s_%d_%d", seqnames, start, end)) %>%
    dplyr::select(id, seqnames, start, end, strand) %>%
    write_tsv(file.path("processed_data", y), col_names = FALSE)
}

createStopCodonAttr <- function(df, ...){
  #' Given biomart transcript data frame append stop codon row
  #' @param df Tibble of transcript coordinates (ie 5utr exon_chrom_start, etc)
  #' @return tibble with stop codon appended
  
  strand <- df$strand[1]
  nrows <- nrow(df)
  if (strand > 0){
    cds_end <- max(df$genomic_coding_end, na.rm = TRUE)
    df <- df %>%
      mutate(genomic_coding_end = if_else(genomic_coding_end == cds_end,
                                          cds_end - 3, as.double(genomic_coding_end)))
    df_row <- tibble_row(stop_codon_start = cds_end - 2,
                         stop_codon_end = cds_end)
  } else{
    cds_start <- min(df$genomic_coding_start, na.rm = TRUE)
    df <- df %>%
      mutate(genomic_coding_start = if_else(genomic_coding_start == cds_start,
                                            cds_start + 3, as.double(genomic_coding_start)))
    df_row <- tibble_row(stop_codon_start = cds_start,
                         stop_codon_end = cds_start + 2)
  }
  
  bind_rows(df, df_row) %>%
    return()
}
getBackground <- function(bed, granges_gene_elements, granges_drach){
  #' selects background matching gene location distribution of bed (Granges)
  #' selected background does not overlap with bed
  #' @param bed Granges object of m6A
  #' @param granges_gene_elements Granges object denoting UTRs, coding, introns
  #' @param granges_drach Granges of genome wide DRACH locations
  #' @return Granges 2.5e4 DRACH locations matched to bed distribution in
  #' granges_gene_elements
  
  chr_names <- map_chr(c(1:22,"X"), function (x) str_interp("chr${x}"))
  
  dist <- filter_by_overlaps(granges_gene_elements, bed) %>%
    filter(seqnames %in% chr_names) %>%
    as_tibble() %>%
    count(element)
  
  set.seed(1)
  
  granges_drach %>%
    filter(seqnames %in% chr_names) %>%
    filter_by_non_overlaps(bed %>% stretch(201)) %>%
    as_tibble() %>%
    left_join(dist, by="element") %>%
    distinct() %>%
    sample_n(2.5e4, weight = n) %>%
    dplyr::select(-c(n, element)) %>%
    as_granges() %>%
    return()
}

plotHomerPngs <- function(sample){
  #' Prints sample header and top3 de novo motifs from homer
  #' @param sample Path to sample
  #' @return NA
  
  cell <- str_split(sample, "\\/")[[1]][2] %>%
    str_split(., "_[dr]na_motifs") %>%
    .[[1]] %>%
    .[1] %>%
    str_to_upper()
  
  pandoc.header(cell, level=3)
  pandoc.p("\n")
  for (j in 1:3){
    img_path <- file.path("..", sample, "homerResults",
                          str_interp("motif${j}.logo.png"))
    pandoc.image(
      img_path,
      caption = str_interp("Motif ${j}"))
    pandoc.p("\n")
  }
}

getKmer <- function(seqs, k){
  #' Counts kmers for vector of seqs
  #' @param seqs vector of strings
  #' @param k length of kmer
  #' @return tibble. Rows are seqs columns are kmers
  
  seqs %>%
    map(~countDnaKmers(., k)) %>%
    purrr::reduce(cbind) %>%
    t() %>%
    as_tibble() %>%
    return()
}

plotKmer <- function(df, kmer){
  #' Plots distributions of kmer between m6A and background and significance
  #' (Mann Whitney)
  #' @param df Tibble with kmers as columns
  #' @param kmer Length of kmer
  #' @return NA
  
  if (kmer != 0){
    feature <- sprintf("%dmer", kmer)
    df2 <- df %>%
      dplyr::select(group, matches(str_interp("^[ACTG]{${kmer}}$")))
    
  } else{
    feature <- "Base Property"
    df2 <- df %>%
      dplyr::select(group, contains('content'))
  }
  plt <- df2 %>%
    tidyr::gather(key=!!sym(feature), value='value', -group) %>%
    ggplot(aes(!!sym(feature), value, fill=group)) +
    geom_boxplot() +
    theme_classic() +
    ggpubr::stat_compare_means(aes(group=group, label = ..p.signif..)) +
    labs(
      title = str_interp("${feature} Distributions"),
      subtitle = str_interp("${stretch_width} BP sequences"),
      x = feature
    ) +
    theme(axis.text.x = element_text(angle = 90, size = 16))
  
  print(plt)
}

getSurroundingSeqPlots <- function(x, name, n){
  #' Prints kmer header and plots
  #' @param x Granges object of (m6)A locations (width = 1)
  #' @param name Name of sample
  #' @param n Width of sequence to find kmers
  #' @return NA
  
  pandoc.header(name, level=2)
  pandoc.p("\n")
  
  seqs <- x %>%
    anchor_center() %>%
    mutate(width = n) %>%
    getSeq(Hsapiens, .) %>%
    as.list() %>%
    map(toString)
  
  df_temp <- tibble(seq = seqs, group=x@elementMetadata$group) %>%
    filter(str_detect(seq, str_interp("[ACTG]{${n}}")))
  
  seqs <- df_temp$seq
  
  onemer <- seqs %>%
    getKmer(1) %>%
    mutate(gc_content = (C + G) / n,
           purine_content = (A + G) / n,
           amino_content = (A + C) / n)
  twomer <- getKmer(seqs, 2)
  threemer <- getKmer(seqs, 3)
  
  df <- bind_cols(onemer, twomer, threemer) %>%
    mutate(group = df_temp$group)
  
  0:3 %>%
    walk(~plotKmer(df, .))
  
  pandoc.p('\n')
}

sample_granges <- function(granges, n=2.5e3){
  #' Sample n rows from granges object
  #' @param granges Granges object
  #' @param n Number of rows to sample. Default: 2500
  #' @return Sampled granges object
  
  set.seed(1)
  granges %>%
    as_tibble() %>%
    sample_n(n) %>%
    as_granges()
}

mergeDonors <- function(df, cell, mod){
  #' Read in files as Granges object and add cell type and histone mod metadata
  #' @param df tibble. encode metadata table
  #' @param cell label granges object
  #' @param mod chip-seq target
  df[["File accession"]] %>%
    map(function(x) file.path("raw_data", str_interp("${x}.bed.gz"))) %>%
    map(read_narrowpeaks) %>%
    purrr::reduce(bind_ranges) %>%
    reduce_ranges() %>%
    mutate(celltype = cell,
           !!mod := 1) %>%
    return()
}

makeMasterTable <- function(granges_m6a, cell, encode_beds){
  #' takes Granges objects(granges_m6a) and makes it longer
  #' each row is bp, histone mark peak is encoded as 0,1
  #' conservation score is used as is
  #' @param granges_m6a Granges object of (potential) m6A locations
  #' @param cell Name of cell
  #' @param encode_beds List of encode narrowpeak files
  #' @return tibble. Each column is annotation. Each row is bp
  
  granges_m6a <- granges_m6a %>%
    mutate(seq = sprintf("%s:%d-%d", seqnames, start, end)) %>%
    anchor_center() %>%
    mutate(width = stretch_width)
  
  granges_phylop <- read_bed_graph(
    "processed_data/hg19.phyloP100way.bedgraph",
    overlap_ranges = granges_m6a) %>%
    plyranges::select(phylop = score)
  
  granges_phast <- read_bigwig("E:/Genomes/hg19/hg19.100way.phastCons.bw",
                               overlap_ranges = granges_m6a) %>%
    plyranges::select(phastcons = score)
  
  granges_hist <- encode_beds %>%
    keep(function(x) str_detect(x@elementMetadata$celltype[1], cell)) %>%
    map(function(x) x %>% plyranges::select(-celltype))
  
  granges_m6a <-  granges_m6a %>%
    purrr::reduce(list(granges_phast, granges_phylop),
                  join_overlap_intersect_directed, .init = .)
  
  df <- granges_hist %>%
    purrr::reduce(join_overlap_left_directed, .init = granges_m6a) %>%
    as_tibble() %>%
    distinct() %>%
    group_by(seq) %>%
    arrange(start, .by_group = TRUE) %>%
    mutate(pos = as.integer(row_number() - stretch_width/2)) %>% # 1-indexed
    mutate(pos = if_else(strand == "+", pos, -pos)) %>%
    dplyr::select(-c(seqnames, start, end, starts_with("X",ignore.case = FALSE),
                     width, strand)) %>%
    replace(., is.na(.), 0)
  
  return(df)
}

histoneFisher <- function(df, hist){
  #' Performs fisher's exact test on contingency table of group and histone mark
  #' group = m6a | background
  #' @param df tibble created from masterTable function
  #' @param hist name of histone mark column in df
  #' @return tibble of fisher's test results
  
  tbl <- df %>%
    ungroup() %>%
    dplyr::select(sym(!!hist), group) %>%
    table()
  
  fisher.test(tbl) %>%
    broom::glance() %>%
    return()
}

plotHistoneOR <- function(df, hist){
  #' @param df tibble produced by masterTable function
  #' @param hist Name of annotation
  #' @return NA
  
  histone_clean <- str_extract(hist, "^[:alnum:]+")
  
  plt <- df %>%
    group_by(pos, group) %>%
    summarise(val = mean(!!sym(hist), na.rm=TRUE)) %>%
    ggplot(aes(pos, val)) +
    geom_smooth(aes(color=group), span=0.1) +
    geom_vline(xintercept = 0, linetype=2) +
    theme_classic() +
    labs(
      x = "Position Rel to m6A (BP)",
      y = histone_clean,
      title = str_interp("Position Effect of ${histone_clean} in m6A Detection")
    )
  print(plt)
  pandoc.p('\n')
  
}

plotConservation <- function(df, hist){
  #' @param df tibble produced by masterTable function
  #' @param hist Name of annotation
  #' @return NA
  
  histone_clean <- str_extract(hist, "^[:alnum:]+")
  
  plt <- df %>%
    group_by(pos, group) %>%
    summarise(val = mean(!!sym(hist), na.rm=TRUE)) %>%
    ggplot(aes(pos, val)) +
    geom_smooth(aes(color=group), span=0.1) +
    geom_vline(xintercept = 0, linetype=2) +
    theme_classic() +
    labs(
      x = "Position Rel to m6A (BP)",
      y = histone_clean,
      title = str_interp("Position Effect of ${histone_clean} in m6A Detection")
    )
  print(plt)
  pandoc.p("\n")
}

plotGroupAnnotwHeader <- function(df, cell, to_keep="human", FUN=plotHistoneOR){
  #' Prints cell header and annotation values by position group rel to m6A
  #' @param df tibble produced by masterTable function
  #' @param cell Name of cell for header
  #' @param to_keep Pattern for annotations to keep
  #' @param FUN plotting function (position group on x axis)
  #' @return NA
  
  pandoc.header(cell, level = 2)
  pandoc.p("\n")
  
  colnames(df) %>%
    keep(~str_detect(.,to_keep)) %>%
    walk(~FUN(df, .))
}

getASEnrichment <- function(granges_sample, granges_back, name){
  #' Compares AS expression between sample and background
  #' @param granges_sample Granges of m6a sample
  #' @param granges_back Granges of background
  #' @param name Name of sample
  #' @return  NA
  
  pandoc.header(name, level=2)
  pandoc.p('\n')
  granges_back2 <- granges_back %>%
    mutate(group = "background") %>%
    plyranges::select(group)
  
  granges_sample2 <- granges_sample %>%
    mutate(group = "m6A") %>%
    plyranges::select(group)
  
  tbl <- bind_ranges(granges_sample2, granges_back2) %>%
    join_overlap_left(granges_coords) %>%
    group_by(group) %>%
    reduce_ranges(Strand = abs(sum(Strand)),
                  ntxns = plyranges::n()) %>%
    mutate(AS_expr = (Strand < ntxns)) %>%
    as_tibble() %>%
    drop_na(AS_expr) %>%
    dplyr::select(AS_expr, group) 
  
  tbl %>%
    group_by(AS_expr, group) %>%
    tally() %>%
    knitr::kable(format = "markdown") %>%
    print()
  
  tbl %>%
    table() %>%
    fisher.test() %>%
    broom::glance() %>%
    knitr::kable(format = "markdown") %>%
    print()
  
  pandoc.p("\n")
}

createInputData <- function(sample, background){
  #' read in sample and background homer files
  #' extract sequences and conservation scores
  #' @param sample prefix of homer file
  #' @param background bool if file is background file
  #' @return tibble
  
  if(background){
    filename <- file.path("processed_data",
                          str_interp("${sample}_m6a.201.background.homer"))
    label <- "background"
  } else{
    filename <- file.path("processed_data",
                          str_interp("${sample}_m6a.201.homer"))
    label <- "m6a"
  }
  
  granges <-  filename %>%
    read_tsv(col_names=FALSE) %>%
    dplyr::select(2:5) %>%
    mutate(id = sprintf("%s_%d_%d", X2, X3, X4)) %>%
    as_granges(seqnames = X2, start = X3, end = X4, strand = X5) %>%
    mutate(group = label) 
  
  seqs <- getSeq(Hsapiens, granges) %>%
    as.list() %>%
    map_chr(toString)
  
  granges_phylo <- read_bed_graph("processed_data/hg19.phyloP100way.bedgraph",
                                  overlap_ranges = granges)
  
  granges <- granges %>%
    mutate(seqs = seqs) %>%
    filter(str_detect(seqs, "^[ACTG]{201}$")) # one sequence has null character
  
  df <- granges %>%
    plyranges::select(-seqs) %>%
    join_overlap_intersect(granges_phylo) %>%
    as_tibble() %>%
    distinct(id, start, .keep_all = TRUE) %>%
    group_by(id) %>%
    arrange(start, .by_group = TRUE) %>%
    mutate(pos = row_number() - 101) %>%
    mutate(pos = if_else(strand == "+", pos, -pos)) %>%
    mutate(pos = pos + 100) %>% # 0 indexed
    ungroup() %>%
    dplyr::select(-c(seqnames, start, end, width, strand)) %>%
    pivot_wider(names_from = pos, values_from = score, names_prefix= "phylo",
                names_sep = "_")
  
  df_seqs <- granges %>%
    as_tibble() %>%
    dplyr::select(id, seqs)
  
  onemer <- df_seqs$seqs %>%
    getKmer(1) %>%
    mutate(gc_content = (C + G) / 201,
           purine_content = (A + G) / 201,
           amino_content = (A + C) / 201)
  twomer <- getKmer(df_seqs$seqs, 2)
  threemer <- getKmer(df_seqs$seqs, 3)
  
  df_seqs <- bind_cols(df_seqs, onemer, twomer, threemer)
  
  bp_cols <- map_chr(0:200, function(x) str_interp("bp_${x}"))
  
  inner_join(df_seqs, df) %>%
    dplyr::select(-id) %>%
    mutate(sample = str_to_upper(sample)) %>%
    separate(seqs, bp_cols, sep = 1:200) %>%
    dplyr::select(sample, group, everything()) %>%
    return()
}