library(tidyverse)
library(dplyr)
library(magrittr)
library(RMariaDB)
library(GenomicRanges)
source('DB_connect.R')

# load saved table if available
promoters.df = dbReadTableCached("promoters", "NishimuraLab")
promoters.df %<>% select(-GENE_NAME, -update_time, -intStrand)
head(promoters.df)  
promoters.gr = GRanges(promoters.df)
head(promoters.gr)

