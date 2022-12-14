library(here)
library(vegan)
library(goeveg) #For scree plot
library(tidyverse)
library(gridExtra)
library(ggrepel)

## Bacteria ######

#Read in cleaned data
bac_abundance = readRDS(here("data", "clean", "bac_abundance.RDS"))
samp_metadata = readRDS(here("data", "clean", "samp_metadata.RDS"))


# Choose to work at the family level
fam_abun = bac_abundance %>% 
  mutate(Family = fct_explicit_na(Family, "Missing")) %>%
  group_by(Family) %>%
  summarise(across(starts_with("Samp"), sum)) 

fam = fam_abun$Family

count_data_no_missing = t(fam_abun[fam!="Missing",-1])
colnames(count_data_no_missing) = fam[-length(fam)]


## Filter to only look at taxa appearing in at least 160 samples
taxa_filt = colSums(count_data_no_missing > 0) >= 160 
counts_filt = count_data_no_missing[,taxa_filt]

#Filter to only look at samples in which at least 50 of the remaining taxa show up
samp_filt = rowSums(counts_filt > 0) >= 50 
counts_filt = counts_filt[samp_filt,]
#Gets rid of one sample


#Run nmds
set.seed(0)
nmds_fit <- metaMDS(counts_filt)

#Check nmds diagnostics
stressplot(nmds_fit) #Shepard diagram
dimcheckMDS(counts_filt) #Stress plot, 
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
  mutate(trt_coarse = factor(case_when(
    trt_bin == "Control" ~ "NoFum",
    trt_bin == "Fumigated" & Time == "Day_10" ~ "RecentlyFum",
    TRUE ~ "FumMonthsAgo"
  ))) %>%
  mutate(trt_coarse = fct_relevel(trt_coarse, "NoFum", "RecentlyFum", "FumMonthsAgo"))

#Use (expanded) weighted averages as given as species scores
species_scores <- scores(nmds_fit)$species %>%
  as_tibble(rownames = "species") %>%
  mutate(abrev = abbreviate(species, minlength = 6))

#Alternatively, use weighted averages as given by wascores
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


#Create plots of NMDS results
gglayer_theme <- list(
  theme_bw(),
  scale_color_brewer(palette="Dark2")
)

samples_plot = ggplot(sample_scores, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = trt_coarse)) +
  scale_color_brewer(name = "Sample Status", 
                     palette = "Dark2", 
                     labels = c("Not\nfumigated", "Recently\nfumigated", "Fumigated\nmonths ago")) +
  ggtitle("Samples") +
  theme_bw() +
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust=0.5)) 

samples_plot +
  geom_text_repel(aes(label = samp_number))

species_plot = ggplot(species_scores2, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = cluster)) +
  geom_text_repel(aes(label = abrev), cex = 2.5, max.overlaps = 20) +
  gglayer_theme +
  ggtitle("OTUs") +
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust=0.5)) 
  

species_plot
  
grid.arrange(samples_plot, species_plot, ncol=2)


saveRDS(sample_scores, here("data", "clean", "sample_scores_bac.RDS"))
saveRDS(species_scores, here("data", "clean", "species_scores_bac.RDS"))
saveRDS(species_scores2, here("data", "clean", "species_scores2_bac.RDS"))


## Fungus ######

fung_abundance = readRDS(here("data", "clean", "fung_abundance.RDS"))

# Choose to work at the genus level
genus_fung_abun = fung_abundance %>% 
  mutate(Genus = fct_explicit_na(Genus, "Missing")) %>%
  group_by(Genus) %>%
  summarise(across(starts_with("Samp"), sum)) 

genus = genus_fung_abun$Genus

# Preparing the OTU counts matrix
count_data_no_missing = t(genus_fung_abun[genus!="Missing",-1])
colnames(count_data_no_missing) = genus[-length(genus)]

## Filter to only look at taxa appearing in at least 50 samples
taxa_filt = colSums(count_data_no_missing > 0) >= 50
counts_filt = count_data_no_missing[, taxa_filt]

#Filter to only look at samples in which at least 20 of the remaining taxa show up
samp_filt = rowSums(counts_filt > 0) >= 20 
counts_filt = counts_filt[samp_filt,]
#Gets rid of 10 samples!

#Run nmds
set.seed(0)
nmds_fit <- metaMDS(counts_filt)

#Check nmds diagnostics
stressplot(nmds_fit) #Shepard diagram
dimcheckMDS(counts_filt) #Stress plot, 
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
  mutate(trt_coarse = factor(case_when(
    trt_bin == "Control" ~ "NoFum",
    trt_bin == "Fumigated" & Time == "Day_10" ~ "RecentlyFum",
    TRUE ~ "FumMonthsAgo"
  ))) %>%
  mutate(trt_coarse = fct_relevel(trt_coarse, "NoFum", "RecentlyFum", "FumMonthsAgo"))

samples_plot = ggplot(sample_scores, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = trt_coarse, shape = Time)) +
  scale_color_brewer(name = "Sample Status", 
                     palette = "Dark2", 
                     labels = c("Not\nfumigated", "Recently\nfumigated", "Fumigated\nmonths ago")) +
  ggtitle("Samples") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        plot.title = element_text(hjust=0.5)) 

samples_plot #+
  #geom_text_repel(aes(label = samp_number))



