library(tidyverse)
library(dplyr)
library(magrittr)
library(RMariaDB)
library(GenomicRanges)
source('DB_connect.R')

# load saved table if available
PromoterChromatinStatesOverlap.df = dbReadTableCached("PromoterChromatinStatesOverlap", "NishimuraLab")
PromoterChromatinStatesOverlap.df 
head(PromoterChromatinStatesOverlap.df)  

