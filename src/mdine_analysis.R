### Performs mdine network inference ###
library(tidyverse)
library(here)

# Read in cleaned data
bac_abundance = readRDS(here("data", "clean", "bac_abundance.RDS"))
fung_abundance = readRDS(here("data", "clean", "fung_abundance.RDS"))
samp_metadata = readRDS(here("data", "clean", "samp_metadata.RDS"))


############# Fungus  #############################

###

# Choose to work at the genus level
genus_fung_abun = fung_abundance %>% 
  mutate(Genus = fct_explicit_na(Genus, "Missing")) %>%
  group_by(Genus) %>%
  summarise(across(starts_with("Samp"), sum)) 


# Preparing the OTU counts matrix

genus_fung_abun = genus_fung_abun %>%
  rowwise() %>%
  mutate(percent_0 = mean(c_across(!Genus) == 0)) %>%
  ungroup()

genus_fung_excl = genus_fung_abun %>%
  filter(percent_0 > .4 | Genus == "Missing")

genus_fung_incl = genus_fung_abun %>%
  filter(percent_0 <= .4 & Genus != "Missing")

## Reference OTU is sum of all excluded OTUs
ref_otu = colSums(genus_fung_excl[,-c(1,242)])

## Get data in format for mdine (matrix of OTU counts with last column reference OTU counts)
Y = cbind(t(genus_fung_incl[,-c(1, 242)]), ref_otu)
colnames(Y) = c(as.character(genus_fung_incl$Genus), "Reference")


# Create model matrix for modelling time, treatment, and their interaction as covariates

samp_metadata = samp_metadata %>%
  mutate(time_coarse = factor(case_when(
    Time == "Day_0" ~ "Base",
    Time == "Day_10" ~ "Seeding",
    Time %in% c("Day_43", "Day_71") ~ "Fall",
    Time %in% c("Day_255", "Day_282") ~ "Spring"))) %>%
  mutate(Treatment = relevel(Treatment, "Non-fumigated chipping grass")) %>%
  mutate(trt_coarse = factor(case_when(
    Time == "Day_0" ~ "Control",
    Treatment == "Non-fumigated chipping grass" ~ "Control",
    TRUE ~ "Fumigated")))

X = model.matrix(~ Treatment*time_coarse, data=samp_metadata)

# Binary variable, fumigated or not
fume = as.integer(samp_metadata$trt_coarse) - 1 

# Binary variable, pre-seed or post-seed
seed = c(rep(0L,80), rep(1L, 160))


library(mdine)
library(rstan)
#md.fit <- mdine(Y=Y, X=cbind(1, fume), Z=seed, mc.cores=parallel::detectCores(), 
#               iter=1000, nnet.MaxNWts=10000, chains=4)

#saveRDS(md.fit, here("models", "mdine_fung_0.4.RDS"))

md.fit = readRDS(here("models", "mdine_fung_0.4.RDS"))
stan.fit = md.fit$stan.fit

md.fit$post_mean$beta

# Estimated precision matrix for control samples (Z=0):
md.fit$post_mean$invsigma0

# Estimated precision matrix for fumigated samples (Z=1):
md.fit$post_mean$invsigma1

# Plotting the two networks
plot_networks(md.fit)

# Weighted adjacency matrices based on each precision matrix
adj <- ci2adj(md.fit, weighted = TRUE)
adj[[1]]

library(igraph)
ig0 <- adj2ig(adj$adj0)
ig1 <- adj2ig(adj$adj1)
ig1 <- set.vertex.attribute(ig1, "name", value=colnames(Y)[-23])
ig0 <- set.vertex.attribute(ig0, "label", value=colnames(Y)[-23])
lay <- layout_nicely(ig1)

par(mfrow=c(1,2))
plot.igraph(ig0, vertex_label=LETTERS[1:22])

plot.igraph(ig1, edge.size = 3, vertex.color = "white")


#Get the elemnts of the precision matrices which differ significantly


difference = sig_diff_prec(md.fit)
igdiff = adj2ig(difference)
igraph::plot.igraph(ig1)














