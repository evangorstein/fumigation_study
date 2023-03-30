library(here)
library(vegan)
library(goeveg) #For scree plot
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(glue)
#Read in data on the samples 
samp_metadata = readRDS(here("data", "clean", "samp_metadata.RDS"))

## Bacteria ######

#Read in cleaned data
bac_abundance = readRDS(here("data", "clean", "bac_abundance.RDS"))

# Choose to work at the family level
fam_abun = bac_abundance %>% 
  mutate(Family = fct_na_value_to_level(Family, "Unknown_Family")) %>% #Missing values included in Unknown_Family level
  group_by(Family) %>%
  summarise(across(starts_with("Samp"), sum)) 

fam = fam_abun$Family

count_data_no_missing = t(fam_abun[fam!="Unknown_Family",-1])
colnames(count_data_no_missing) = fam[-length(fam)]

#Filter to only look at families that have at least 1,000 counts in total (summing over samples) 
taxa_filt = colSums(count_data_no_missing) >= 1000 
counts_filt = count_data_no_missing[,taxa_filt]
#Left with 152 of 295 families

#Get rid of sample 173, which has much fewer counts than all the other samples and which screws up the NMDS 
counts_filt = counts_filt[-173,]


#Run nmds
set.seed(0)
nmds_fit <- metaMDS(counts_filt)

#Check nmds diagnostics
stressplot(nmds_fit) #Shepard diagram
L = dimcheckMDS(counts_filt) #Stress plot, 
#may take a minute or two to run as it fits NMDS for several different k


#Get sample scores in nice data frame with sample meta data merged
sample_scores = scores(nmds_fit)$sites %>%
  as_tibble(rownames = "samp_number") %>%
  mutate(samp_number = as.numeric(str_replace(samp_number, "Samp", ""))) %>%
  inner_join(samp_metadata) %>%
  mutate(trt_bin = factor(case_when(
    Time == "Day_0" ~ "Control",
    Treatment == "Non-fumigated chipping grass" ~ "Control",
    TRUE ~ "Fumigated"))) %>%
  mutate(`fumigation status` = factor(case_when(
    trt_bin == "Control" ~ "Never",
    trt_bin == "Fumigated" & Time == "Day_10" ~ "Recent",
    TRUE ~ "Past"
  ), levels = c("Never", "Past", "Recent")))



color_coding = c(Never = "#00BA38", Past = "#F8766D", Recent = "#619CFF")


#Plot samples
samples_plot = ggplot(sample_scores, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = `fumigation status`, shape = Treatment)) +
  #scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = color_coding) +
  ggtitle("Samples") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle("NMDS with stress 0.074")

samples_plot 

##Get NMDS scores for species
species_scores <- scores(nmds_fit)$species %>%
  as_tibble(rownames = "species") %>%
  mutate(abrev = abbreviate(species, minlength = 6))

#Alternatively, use weighted averages as given by vegan::wascores
species_scores2 = wascores(scores(nmds_fit)$sites, counts_filt, expand = TRUE) %>%
  as_tibble(rownames = "species") %>%
  mutate(abrev = abbreviate(species, minlength = 6))


#K-means cluster species scores 
set.seed(0)
kmeans_fit1 <- kmeans(select(species_scores, NMDS1, NMDS2), centers = 3)
kmeans_fit2 <- kmeans(select(species_scores2, NMDS1, NMDS2), centers = 3)
species_scores$cluster = factor(kmeans_fit1$cluster) 
species_scores2$cluster = factor(kmeans_fit2$cluster) %>%
  fct_recode(Cluster1 = "3", Cluster2 = "1", Cluster3 = "2") %>%
  fct_relevel("Cluster1", "Cluster2", "Cluster3")

color_coding_sp = c(Cluster1 = "#00BA38", Cluster3 = "#F8766D", Cluster2 = "#619CFF")
  
##Plot species NMDS results
species_plot = ggplot(species_scores2, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = cluster)) +
  scale_color_manual(values = color_coding_sp) +
  geom_text_repel(aes(label = abrev), cex = 2.5, max.overlaps = 20) +
  theme_bw() +
  ggtitle("Bacterial families") +
  theme(plot.title = element_text(hjust=0.5)) 

species_plot 

#Put two plots next to each other
grid.arrange(samples_plot, species_plot, ncol=2)


## Fungus ######

fung_abundance = readRDS(here("data", "clean", "fung_abundance.RDS"))

# Choose to work at the genus level
genus_fung_abun = fung_abundance %>% 
  mutate(Genus = fct_explicit_na(Genus, "Missing")) %>%
  group_by(Genus) %>%
  summarise(across(starts_with("Samp"), sum)) 

genuses = genus_fung_abun$Genus

# Preparing the OTU counts matrix
count_data_no_missing = t(genus_fung_abun[genus!="Missing",-1])
colnames(count_data_no_missing) = genus[-length(genus)]

#Filter to only look at families that have at least 100 counts in total (summing over samples) 
taxa_filt = colSums(count_data_no_missing) >= 100 
counts_filt = count_data_no_missing[,taxa_filt]
#Left with 150 of 408 families


#Filter to only look at samples in which at least 20 of the 408 genuses show up
samp_filt = rowSums(count_data_no_missing > 0) >= 20 
count_data_no_missing = count_data_no_missing[samp_filt,]
#Gets rid of 8 samples
#Filter to only look at samples in which at least 20 of the remaining 150 genuses show up
samp_filt = rowSums(counts_filt > 0) >= 20 
counts_filt = counts_filt[samp_filt,]
#Gets rid of 10 samples!


#Run nmds
set.seed(0)
nmds_fit <- metaMDS(counts_filt)

#Check nmds diagnostics
stressplot(nmds_fit) #Shepard diagram
L = dimcheckMDS(counts_filt) #Stress plot, 
#may take a minute or two to run as it fits NMDS for several different k


#Get sample scores in nice data frame with sample meta data merged
sample_scores = scores(nmds_fit)$sites %>%
  as_tibble(rownames = "samp_number") %>%
  mutate(samp_number = as.numeric(str_replace(samp_number, "Samp", ""))) %>%
  inner_join(samp_metadata) %>%
  mutate(trt_bin = factor(case_when(
    Time == "Day_0" ~ "Control",
    Treatment == "Non-fumigated chipping grass" ~ "Control",
    TRUE ~ "Fumigated"))) %>%
  mutate(`fumigation status` = factor(case_when(
    trt_bin == "Control" ~ "Never",
    trt_bin == "Fumigated" & Time == "Day_10" ~ "Recent",
    TRUE ~ "Past"
  ), levels = c("Never", "Past", "Recent")))


color_coding = c(Never = "#00BA38", Past = "#F8766D", Recent = "#619CFF")


#Plot samples
samples_plot = ggplot(sample_scores, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = `fumigation status`, shape = Time)) +
  #scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = color_coding) +
  ggtitle("Samples") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle(glue("NMDS with stress {round(L[2],2)}"))

samples_plot #+
  #geom_text_repel(aes(label = samp_number))




