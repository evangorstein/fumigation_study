# load package
library(pheatmap)
library(tidyverse)
library(dplyr)

# install required package if not already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Read in data
bac_abundance = readRDS(here("data", "clean", "bac_abundance.RDS"))

# Analyze in family level
fam_abun = bac_abundance %>% 
  mutate(Family = fct_explicit_na(Family, "Missing")) %>% 
  group_by(Family) %>% 
  summarise(across(starts_with("Samp"), sum)) %>%
  filter(Family != "Missing")

fam_abun = as.matrix(fam_abun)
fam_abun[is.na(fam_abun)] = 0
fam_abun[is.nan(fam_abun)] = 0
pheatmap(fam_abun)
