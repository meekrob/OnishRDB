library(tidyverse)
library(dplyr)
library(magrittr)
library(RMariaDB)
library(GenomicRanges)
source('DB_connect.R')

# load saved table if available
evansChromatinStates = dbReadTableCached("evansChromatinStates", "NishimuraLab")
evansChromatinStates %<>% select(-id_evansChromatinStates, update_time)
head(evansChromatinStates)  
gr = GRanges(evansChromatinStates)
head(gr)
