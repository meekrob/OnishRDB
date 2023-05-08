library(tidyverse)
library(dplyr)
library(magrittr)
library(RMariaDB)
library(GenomicRanges)
source('DB_connect.R')

# load saved table if available
evansChromatinStates.df = dbReadTableCached("evansChromatinStates", "NishimuraLab")
evansChromatinStates.df %<>% select(-update_time)
head(evansChromatinStates.df)  
evansChromatinStates.gr = GRanges(evansChromatinStates.df)
head(evansChromatinStates.gr)
