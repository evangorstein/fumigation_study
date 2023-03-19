library(pheatmap)
library(tidyverse)
library(dplyr)
library(here)
library(grid)
Sys.setlocale("LC_ALL","English") #tianyi

#Read in cleaned data
fung_abundance = readRDS(here("data", "clean", "fung_abundance.RDS"))
samp_metadata = readRDS(here("data", "clean", "samp_metadata.RDS"))


genus_fung_abun = fung_abundance %>% 
  mutate(Genus = fct_explicit_na(Genus, "Missing")) %>%
  group_by(Genus) %>%
  summarise(across(starts_with("Samp"), sum)) 

genus = genus_fung_abun$Genus

count_data_no_missing = t(genus_fung_abun[genus!="Missing",-1])
colnames(count_data_no_missing) = genus[-length(genus)]

## Filter to only look at taxa appearing in at least 50 samples
taxa_filt = colSums(count_data_no_missing > 0) >= 150
counts_filt = count_data_no_missing[, taxa_filt]

# final data ready to use
final_data = t(as.matrix(count_data_no_missing))

# function to transform data into library sum
library_sum = function(x){
  x / sum(x)
}
# function to make a heatmap.
# argument: df(matrix): data, annotation(matrix): annotation for colname, colname, rowname(boolean): 
# show column\row name or not. ylabel, xlabel,title(String): specify labels and titles.
create_heatmap = function(df, annotation = NULL, 
                          colname = FALSE, 
                          rowname = FALSE, 
                          xlabel = NULL, 
                          ylabel = NULL, 
                          title = NULL, 
                          clustercol = TRUE){
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  heatmap = pheatmap(df, annotation_col = annotation, 
                     show_colnames = colname, 
                     show_rownames = rowname, 
                     main = title, angle_col = 0, 
                     cluster_cols = clustercol)
  setHook("grid.newpage", NULL, "replace")
  grid.text(xlabel, y=-0.07, gp=gpar(fontsize=16))
  grid.text(ylabel, x=-0.07, rot=90, gp=gpar(fontsize=16))
  return (heatmap)
}

norm_data = apply(final_data, 2, library_sum)

plt1 = create_heatmap(norm_data, 
                      title = "Normalized Number of OTUs in each sample", 
                      xlabel = "Samples 1-240", 
                      ylabel = "Genus OTUs")

data_subset = as.matrix(final_data[rowSums(final_data)>25000,])
data_subset_norm = apply(data_subset, 2, library_sum)
plt2 = create_heatmap(data_subset_norm, 
                      title = "Subset where Each Sample with more than 25000", 
                      xlabel = "Samples 1-240", 
                      ylabel = "Genus OTUs")

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

samp_never_fumigated = c(sapply(never_fumigated, make_sample))
samp_recently = c(sapply(recently, make_sample))
samp_more = c(sapply(more_than_a_month, make_sample))

all_samps = colnames(final_data)

##TODO: check to see if this is a valid way, maybe improve it.
samp_annotation = data.frame(all_samps) %>% 
  mutate(group = case_when(
    all_samps %in% samp_more ~ "more_than_a_month_fumigated",
    all_samps %in% samp_never_fumigated ~ "non_fumigated",
    all_samps %in% samp_recently ~ "recent_fumigated"
  )) %>% 
  select(-all_samps)

row.names(samp_annotation) = colnames(norm_data)
plt6 = create_heatmap(norm_data,
                      annotation = samp_annotation,
                      title = "Normalized Number of OTUs in each sample",
                      xlabel = "Samples 1-240",
                      ylabel = "Genus OTUs"
)


## filtering for graph7
plt7 = create_heatmap(data_subset_norm,
                      annotation = samp_annotation,
                      rowname = TRUE,
                      title = "Normalized Number of OTUs in each sample",
                      xlabel = "Samples 1-240",
                      ylabel = "Genus OTUs"
)

data_subset2 = as.matrix(final_data[rowSums(final_data)>10000,])
data_subset2 = apply(data_subset2, 2, library_sum)
plt8 = create_heatmap(data_subset2,
                      annotation = samp_annotation,
                      rowname= TRUE,
                      title = "Normalized Number of OTUs in each sample",
                      xlabel = "Samples 1-240",
                      ylabel = "Genus OTUs"
)