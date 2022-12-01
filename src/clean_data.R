library(tidyverse)
library(here)

#Loads taxa and seqtab.nochim matrices
load(here("data", "16S", "Fumigation_methyl_16s.rdata"))

#Main abundance data frame
seqtab.nochim <- as_tibble(t(seqtab.nochim), rownames=NA) %>%
  #transposed so that sequences are rows, which is the format of taxa
  rownames_to_column(var="sequence")

#Metadata on sequences (i.e., phylogenetic classification)
taxa <- as_tibble(taxa, rownames=NA) %>%
  rownames_to_column(var="sequence") %>%
  mutate(across(.fns = factor)) 
  

#Metadata on samples
library(lubridate)
metadata <- read_csv(here("data", "Fum_metadata_16s.csv"), show_col_types = FALSE) %>%
  mutate(across(c(Time, Treatment), factor)) %>%
  mutate(Date = mdy(Date)) %>%
  arrange(samp_number)

#Prepend taxa information as additional columns to seqtab.nochim
abundance <- taxa %>%
  left_join(seqtab.nochim, by = "sequence") 


#Summary information on taxa present 
table(taxa$Species, useNA = "ifany")
table(taxa$Genus, useNA = "ifany")
table(taxa$Family, useNA = "ifany")
table(taxa$Order, useNA = "ifany")
table(taxa$Class, useNA = "ifany")
table(taxa$Phylum, useNA = "ifany")
table(taxa$Kingdom, useNA = "ifany")


#Let's work at phylum and class levels

phylum_abundance = abundance %>% 
  mutate(Phylum = addNA(Phylum)) %>%
  group_by(Phylum) %>%
  summarise(across(.cols = starts_with("Fum"), .fns = sum))
  
class_abundance = abundance %>% 
  mutate(Class = addNA(Class)) %>%
  group_by(Class) %>%
  summarise(across(.cols = starts_with("Fum"), .fns = sum))


saveRDS(phylum_abundance, file = here("data", "clean", "phylum_abundance.RDS"))
saveRDS(class_abundance, file = here("data", "clean", "class_abundance.RDS"))
saveRDS(metadata, file = here("data", "clean", "class_abundance.RDS"))


