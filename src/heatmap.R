# load package
# reference: https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
library(pheatmap)
library(tidyverse)
library(dplyr)
library(here)
library(grid)

# install required package if not already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Read in data
bac_abundance = readRDS(here("data", "clean", "bac_abundance.RDS"))
samp_metadata = readRDS(here("data", "clean", "samp_metadata.RDS"))
# note to self- could group by day, and also treatment for visualization

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
## this does not remove any datapoints, has all 296 rows.
## plot1, before normalization, shows absolute abundance, might not be insightful.
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
plt = pheatmap(fam_abun1, show_colnames = FALSE, show_rownames = FALSE, main = "Number of OTUs in Each Sample")
setHook("grid.newpage", NULL, "replace")
grid.text("Samples 1-250", y=-0.07, gp=gpar(fontsize=16))
grid.text("Family OTUs", x=-0.07, rot=90, gp=gpar(fontsize=16))

## plot 2, normalize using library sum
library_sum <- function(x){
  x / sum(x)
}

fam_abun_norm <- apply(fam_abun1, 2, library_sum)

setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
plt2 = pheatmap(fam_abun_norm, show_colnames = FALSE, show_rownames = FALSE, main = "Normalized by Library Sum")
setHook("grid.newpage", NULL, "replace")
grid.text("Samples 1-250", y=-0.07, gp=gpar(fontsize=16))
grid.text("Family OTUs", x=-0.07, rot=90, gp=gpar(fontsize=16))


##plt3, as can be see, there are some issues within each plot, since a lot OTU's has 0 values, it is not easy to see relations.
## lots of OTUs has 0 in the dataset, the data is very sparse.
## I think it might be helpful for better visualize using OTUS that has some abundance in the soil

fam_abun_subset = as.matrix(fam_abun1[rowSums(fam_abun1)>25000,])
fam_subset_norm <- (apply(fam_abun_subset, 2, library_sum))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
plt3 = pheatmap(fam_subset_norm, show_colnames = FALSE,  main = "SUbset where Each Sample with more than 25000")
setHook("grid.newpage", NULL, "replace")
grid.text("Samples 1-250", y=-0.07, gp=gpar(fontsize=16))
grid.text("Family OTUs", x=-0.07, rot=90, gp=gpar(fontsize=16))







## Hard to see any useful information if not group by sample, also reflect staff.
## Try do some bact anal by evan
fam_abun = bac_abundance %>% 
  mutate(Family = fct_explicit_na(Family, "Missing")) %>%
  group_by(Family) %>%
  summarise(across(starts_with("Samp"), sum)) 

fam = fam_abun$Family
count_data_no_missing = t(fam_abun[fam!="Missing",-1])
colnames(count_data_no_missing) = fam[-length(fam)]

taxa_filt = colSums(count_data_no_missing > 0) >= 160
counts_filt = count_data_no_missing[,taxa_filt]

#Filter to only look at samples in which at least 50 of the remaining taxa show up
samp_filt = rowSums(counts_filt > 0) >= 50
counts_filt = counts_filt[samp_filt,]

#Get rid of single sample that doesn't meet the criteria in count_data_no_missing, as well
count_data_no_missing = t(count_data_no_missing[samp_filt,])

final_fam = as.matrix(count_data_no_missing)

## plot 4, normalize using library sum
library_sum <- function(x){
  x / sum(x)
}
fam_nomissing <- apply(final_fam, 2, library_sum)

setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
plt4 = pheatmap(fam_nomissing, show_colnames = FALSE, show_rownames = FALSE, main = "Normalized by Library Sum")
setHook("grid.newpage", NULL, "replace")
grid.text("Samples 1-240", y=-0.07, gp=gpar(fontsize=16))
grid.text("Family OTUs", x=-0.07, rot=90, gp=gpar(fontsize=16))


## other functions might be useful

# function to get dimension of heat_map, might be useful
get_plot_dims <- function(heat_map)
{
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# referring to data which works on pheatmap, does not need to be executed
BiocManager::install("DESeq")
library("DESeq")
example_file <- system.file ("./../TagSeqExample.tab")
data <- read.delim("TagSeqExample.tab", header=T, row.names="gene")
data_subset <- as.matrix(data[rowSums(data)>50000,])
pheatmap(data_subset)

##Uneeded Part for now
##plt2, normalize by row. Better Visualization
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

fam_abun_norm <- t(apply(fam_abun1, 1, cal_z_score))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
plt2 = pheatmap(fam_abun_norm, show_colnames = FALSE, show_rownames = FALSE, main = "Z-score for Each OTUs in Each Sample")
setHook("grid.newpage", NULL, "replace")
grid.text("Samples 1-250", y=-0.07, gp=gpar(fontsize=16))
grid.text("Family OTUs", x=-0.07, rot=90, gp=gpar(fontsize=16))
