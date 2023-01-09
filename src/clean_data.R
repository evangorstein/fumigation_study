## Wrangling and cleaning abundance data and sample metadata
## Saves dataframes to folder `/data/clean`

library(tidyverse)
library(here)

#### Bacteria #################

#Loads taxa and seqtab.nochim matrices
load(here("data", "16S", "Fumigation_methyl_16s.rdata"))

#Main abundance data frame
seqtab.nochim <- as_tibble(t(seqtab.nochim), rownames=NA) %>%
  #transposed so that sequences are rows, which is the format of taxa
  rownames_to_column(var="sequence") %>% 
  #Change column names so that they can be ordered by sample number
  rename_with(~str_replace(., "Fum\\-1\\-", "Samp"), !sequence) %>% 
  rename_with(~str_replace(., "\\d+", function(m) str_pad(m, 3, pad = "0")), !sequence)%>%
  #Order by sample number
  select(sequence, sort(names(.)))
  
#Metadata on sequences (i.e., phylogenetic classification)
taxa_bac <- as_tibble(taxa, rownames=NA) %>%
  rownames_to_column(var="sequence") %>%
  mutate(across(.fns = factor)) 
  
#Metadata on samples
library(lubridate)
samp_metadata <- read_csv(here("data", "Fum_metadata_16s.csv"), show_col_types = FALSE) %>%
  mutate(across(c(Time, Treatment), factor)) %>%
  mutate(Date = mdy(Date)) %>%
  arrange(samp_number)


#Prepend taxa information as additional columns to seqtab.nochim
bac_abundance <- taxa_bac %>%
  left_join(seqtab.nochim, by = "sequence") 


### Fungus #################
#Loads taxa and seqtab.nochim matrices
load(here("data", "ITS", "Fumigation_methyl_ITS.rdata"))


#Main abundance data frame
seqtab.nochim <- as_tibble(t(seqtab.nochim), rownames=NA) %>%
  #transposed so that sequences are rows, which is the format of taxa
  rownames_to_column(var="sequence") %>% 
  #Change column names so that they can be ordered by sample number
  rename_with(~str_replace(., "Fum1\\-", "Samp"), !sequence) %>% 
  rename_with(~str_replace(., "\\d+", function(m) str_pad(m, 3, pad = "0")), !sequence) %>%
  #Order by sample number
  select(sequence, sort(names(.)))

#Metadata on sequences (i.e., phylogenetic classification)
taxa_fung <- as_tibble(taxa, rownames=NA) %>%
  rownames_to_column(var="sequence") %>%
  mutate(across(.fns = factor)) 

#Prepend taxa information as additional columns to seqtab.nochim
fung_abundance <- taxa_fung %>%
  left_join(seqtab.nochim, by = "sequence") 


#Save data frames
saveRDS(bac_abundance, file = here("data", "clean", "bac_abundance.RDS"))
saveRDS(fung_abundance, file = here("data", "clean", "fung_abundance.RDS"))
saveRDS(metadata, file = here("data", "clean", "samp_metadata.RDS"))
saveRDS(taxameta, file = here("data", "clean", "taxa_metadata.RDS"))












