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
fam_nomissing = apply(final_fam, 2, library_sum)

# plot 1, normalize using library sum, visualize over all samples
plt1 = create_heatmap(fam_nomissing, 
                      title = "Normalized Number of OTUs in each sample", 
                      xlabel = "Samples 1-240", 
                      ylabel = "Family OTUs")

#plt2, as can be see, there are some issues within each plot, since a lot OTU's has 0 values, it is not easy to see relations.
# lots of OTUs has 0 in the dataset, the data is very sparse.
# I think it might be helpful for better visualize using OTUS that has some abundance in the soil
fam_abun_subset = as.matrix(final_fam[rowSums(final_fam)>25000,])
fam_subset_norm = apply(fam_abun_subset, 2, library_sum)
plt2 = create_heatmap(fam_subset_norm, 
                      title = "Subset where Each Sample with more than 25000", 
                      xlabel = "Samples 1-240", 
                      ylabel = "Family OTUs")


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

## construct the matrix
group_fam = cbind(final_fam, non_fumigated)
group_fam = cbind(group_fam, recent_fumigated)
group_fam = cbind(group_fam, more_than_a_month_fumigated)
group_fam = group_fam[,c("non_fumigated","recent_fumigated", "more_than_a_month_fumigated")]

#plot 3
group_fam_normalize = group_fam[rowSums(group_fam)>0,]
group_fam_normalize = apply(group_fam_normalize, 2, library_sum)
plt3 = create_heatmap(group_fam_normalize, 
                      colname = TRUE, 
                      title = "OTU Based on Fumigation Status", 
                      xlabel = "Fumigation Status", 
                      ylabel = "Family OTUs", 
                      clustercol = FALSE)


## plot4 - look at OTUs presented the most
group_fam_normalize2 = group_fam[rowSums(group_fam)>20000,]
group_fam_normalize2 <- apply(group_fam_normalize2, 2, library_sum)
plt4 = create_heatmap(group_fam_normalize2, 
                      colname = TRUE, 
                      rowname = TRUE, 
                      title = "OTU Based on Fumigation Status", 
                      xlabel = "Fumigation Status", 
                      ylabel = "Family OTUs", 
                      clustercol = FALSE)


##5 - look at OTUs that has a big weight inside the microbiome
group_fam_normalize3 = group_fam[rowSums(group_fam)>40000,]
group_fam_normalize3 <- apply(group_fam_normalize3, 2, library_sum)
plt5 = create_heatmap(group_fam_normalize3,
                      colname = TRUE,
                      rowname = TRUE,
                      title = "OTU Based on Fumigation Status",
                      xlabel = "Fumigation Status",
                      ylabel = "Family OTUs",
                      clustercol = FALSE
)

## clustering map- 6 and 7, on different sizes.

# first, create the annotation for groups
all_samps = colnames(final_fam)

##TODO: check to see if this is a valid way, maybe improve it.
samp_annotation = data.frame(all_samps) %>% 
  mutate(group = case_when(
    all_samps %in% samp_more ~ "more_than_a_month_fumigated",
    all_samps %in% samp_never_fumigated ~ "non_fumigated",
    all_samps %in% samp_recently ~ "recent_fumigated"
  )) %>% 
  select(-all_samps)
row.names(samp_annotation) = colnames(fam_nomissing)
plt6 = create_heatmap(fam_nomissing,
                      annotation = samp_annotation,
                      title = "Normalized Number of OTUs in each sample",
                      xlabel = "Samples 1-240",
                      ylabel = "Family OTUs"
)


## filtering for graph7
plt7 = create_heatmap(fam_abun_subset,
                      annotation = samp_annotation,
                      rowname = TRUE,
                      title = "Normalized Number of OTUs in each sample",
                      xlabel = "Samples 1-240",
                      ylabel = "Family OTUs"
)

fam_abun_subset2 = as.matrix(final_fam[rowSums(final_fam)>40000,])
fam_abun_subset2 = apply(fam_abun_subset2, 2, library_sum)

plt8 = create_heatmap(fam_abun_subset2,
                      annotation = samp_annotation,
                      rowname= TRUE,
                      title = "Normalized Number of OTUs in each sample",
                      xlabel = "Samples 1-240",
                      ylabel = "Family OTUs"
)




## Visualize from order OTU

order_abun = bac_abundance %>% 
  mutate(Order = fct_explicit_na(Order, "Missing")) %>%
  group_by(Order) %>%
  summarise(across(starts_with("Samp"), sum)) 
ord = order_abun$Order
count_data_no_missing2 = t(order_abun[ord!="Missing",-1])
colnames(count_data_no_missing2) = ord[-length(ord)]

taxa_filt2 = colSums(count_data_no_missing2 > 0) >= 160
counts_filt2 = count_data_no_missing2[,taxa_filt2]

#Filter to only look at samples in which at least 50 of the remaining taxa show up
samp_filt2 = rowSums(counts_filt2 > 0) >= 50
counts_filt2 = counts_filt2[samp_filt2,]

#Get rid of single sample that doesn't meet the criteria in count_data_no_missing, as well
count_data_no_missing2 = t(count_data_no_missing2[samp_filt2,])
final_ord = as.matrix(count_data_no_missing2)

ord_nomissing = apply(final_ord, 2, library_sum)
# first, create the annotation for groups
all_samps = colnames(final_ord)


samp_annotation = data.frame(all_samps) %>% 
  mutate(group = case_when(
    all_samps %in% samp_more ~ "more_than_a_month_fumigated",
    all_samps %in% samp_never_fumigated ~ "non_fumigated",
    all_samps %in% samp_recently ~ "recent_fumigated"
  )) %>% 
  select(-all_samps)
row.names(samp_annotation) = colnames(ord_nomissing)


# plot 9, normalize using library sum, visualize over all samples
plt9 = create_heatmap(ord_nomissing,
                      annotation = samp_annotation,
                      title = "Normalized Number of OTUs in each sample",
                      xlabel = "Samples 1-240",
                      ylabel = "Order OTUs"
)

# Plot 10, only look at OTUs that have over 40000 abundance across samples.
ord_abun_subset2 = as.matrix(final_ord[rowSums(final_ord)>40000,])
ord_abun_subset2 = apply(ord_abun_subset2, 2, library_sum)

plt10 = create_heatmap(ord_abun_subset2,
                      annotation = samp_annotation,
                      rowname= TRUE,
                      title = "Normalized Number of OTUs in each sample",
                      xlabel = "Samples 1-240",
                      ylabel = "Order OTUs"
)