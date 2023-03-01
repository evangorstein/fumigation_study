# load package
# reference: https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
library(pheatmap)
library(tidyverse)
library(dplyr)
library(here)
library(grid)

# Read in data
bac_abundance = readRDS(here("data", "clean", "bac_abundance.RDS"))
samp_metadata = readRDS(here("data", "clean", "samp_metadata.RDS"))

##only looking at family level, borrowed from Evan
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

# final data ready to use
final_fam = as.matrix(count_data_no_missing)

## plot 1, normalize using library sum, visualize over all samples
library_sum <- function(x){
  x / sum(x)
}
fam_nomissing <- apply(final_fam, 2, library_sum)

setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
plt1 = pheatmap(fam_nomissing, show_colnames = FALSE, show_rownames = FALSE, main = "Normalized Number of OTUs in each sample")
setHook("grid.newpage", NULL, "replace")
grid.text("Samples 1-240", y=-0.07, gp=gpar(fontsize=16))
grid.text("Family OTUs", x=-0.07, rot=90, gp=gpar(fontsize=16))

##plt2, as can be see, there are some issues within each plot, since a lot OTU's has 0 values, it is not easy to see relations.
## lots of OTUs has 0 in the dataset, the data is very sparse.
## I think it might be helpful for better visualize using OTUS that has some abundance in the soil
fam_abun_subset = as.matrix(final_fam[rowSums(final_fam)>25000,])
fam_subset_norm <- apply(fam_abun_subset, 2, library_sum)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
plt3 = pheatmap(fam_subset_norm, show_colnames = FALSE,  main = "SUbset where Each Sample with more than 25000")
setHook("grid.newpage", NULL, "replace")
grid.text("Samples 1-250", y=-0.07, gp=gpar(fontsize=16))
grid.text("Family OTUs", x=-0.07, rot=90, gp=gpar(fontsize=16))



## new_task, now, we should try based on different groupings.
##first- suggestion from Evan
## 3 sample groups, never fumigated, fumigated more than a month ago, recently fumigated 

## idea: - 1st - never fumigated (40 + 8*5) = 80 samples, 2nd - day10 -32 samples. Morethanmonth - any others

## find the sample number corresponding to each category
never_fumigated = samp_metadata %>% 
  filter(Time == "Day_0" | Treatment == "Non-fumigated chipping grass") %>% 
  pull(samp_number)

recently = samp_metadata %>% 
  filter(Time == "Day_10" & Treatment != "Non-fumigated chipping grass") %>% 
  pull(samp_number)

more_than_a_month = samp_metadata %>% 
  filter(Time != "Day_0" & Time != "Day_10" & Treatment != "Non-fumigated chipping grass") %>% 
  pull(samp_number)

# function that convert sample number to column names
make_sample = function(x){
  if(x>=0 & x<=9){
    return(paste0("Samp00", x))
  } else if(x>=10 & x<=99){
    return(paste0("Samp0", x))
  } else {
    return(paste0("Samp", x))
  }
}

## Create column names based on column numbers
samp_never_fumigated = c(sapply(never_fumigated, make_sample))
samp_recently = c(sapply(recently, make_sample))
samp_more = c(sapply(more_than_a_month, make_sample))

non_fumigated = rowSums(final_fam[,samp_never_fumigated])
recent_fumigated = rowSums(final_fam[,samp_recently])
more_than_a_month_fumigated = rowSums(final_fam[,more_than_a_month])

group_fam = cbind(final_fam, non_fumigated)
group_fam = cbind(group_fam, recent_fumigated)
group_fam = cbind(group_fam, more_than_a_month_fumigated)
group_fam = group_fam[,c("non_fumigated", "recent_fumigated", "more_than_a_month_fumigated")]

group_fam_normalize = group_fam[rowSums(group_fam)>0,]
group_fam_normalize <- apply(group_fam_normalize, 2, library_sum)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
plt3 = pheatmap(group_fam_normalize, show_rownames = FALSE, main = "SUbset where Each Sample with more than 25000")
setHook("grid.newpage", NULL, "replace")
grid.text("Samples 1-250", y=-0.07, gp=gpar(fontsize=16))
grid.text("Family OTUs", x=-0.07, rot=90, gp=gpar(fontsize=16))

## plot4 - look at OTUs presented the most
group_fam_normalize2 = group_fam[rowSums(group_fam)>20000,]
group_fam_normalize2 <- apply(group_fam_normalize2, 2, library_sum)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
plt4 = pheatmap(group_fam_normalize2, show_rownames = TRUE, main = "SUbset where Each Sample with more than 25000")
setHook("grid.newpage", NULL, "replace")
grid.text("Samples 1-250", y=-0.07, gp=gpar(fontsize=16))
grid.text("Family OTUs", x=-0.07, rot=90, gp=gpar(fontsize=16))

##5 - look at OTUs with big difference
group_fam_normalize2 = group_fam[rowSums(group_fam)>40000,]
group_fam_normalize2 <- apply(group_fam_normalize2, 2, library_sum)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
plt5 = pheatmap(group_fam_normalize2, show_rownames = TRUE, main = "SUbset where Each Sample with more than 25000")
setHook("grid.newpage", NULL, "replace")
grid.text("Samples 1-250", y=-0.07, gp=gpar(fontsize=16))
grid.text("Family OTUs", x=-0.07, rot=90, gp=gpar(fontsize=16))