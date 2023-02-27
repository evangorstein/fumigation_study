# load package
library(pheatmap)
library(tidyverse)
library(dplyr)
library(here)

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

# change row names
fam_abun1 = fam_abun %>% 
  select(-Family)
rownames(fam_abun1) <- fam_abun$Family
fam_abun1 = as.matrix(fam_abun1)

## create simple heatmap
pheatmap(fam_abun1)


# referring to data which works on pheatmap
BiocManager::install("DESeq")
library("DESeq")
example_file <- system.file ("./../TagSeqExample.tab")
data <- read.delim("TagSeqExample.tab", header=T, row.names="gene")
data_subset <- as.matrix(data[rowSums(data)>50000,])
pheatmap(data_subset)
