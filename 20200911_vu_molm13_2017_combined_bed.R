# Combine replicates of miclip experiments from Vu 2017 Nature Medicine
library(tidyverse)
library(plyranges)

list.files(path="raw_data/", "^GSM260*", full.names = TRUE) %>%
  map(read_bed) %>%
  reduce(join_overlap_intersect_directed) %>%
  write_bed('processed_data/vu_2017_molm13_miclip_combined.bed')
