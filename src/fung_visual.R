# Creating heatmaps of bacterial data
# reference: https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
library(pheatmap)
library(tidyverse)
library(here)

#### Load and prepare data ####

# Read in data
fung_abundance = readRDS(here("data", "clean", "fung_abundance.RDS"))
samp_metadata = readRDS(here("data", "clean", "samp_metadata.RDS"))

#Aggregate to genus level
genus.df = fung_abundance %>% 
  mutate(Genus = fct_explicit_na(Genus, "Unknown")) %>%
  group_by(Genus) %>%
  summarise(across(starts_with("Samp"), sum)) 

#Note that the "Unknown" value aggregates the counts for all the OTUs whose families are unknown
genuses = genus.df$Genus
genus.mat = as.matrix(genus.df[,-1])
rownames(genus.mat) = genuses

#Check most abundant genuses, which includes "Unknown"
sort(rowSums(genus.mat))

#Filter to only look at genuses that have at least 100 counts in total (summing over samples) 
genus_mask = rowSums(genus.mat) >= 100
genus.filt = genus.mat[genus_mask,]
#151 out of 409 genuses remain

#Filter to only look at genuses that have at least 1,000 counts in total (summing over samples) 
genus_mask_harsh = rowSums(genus.mat) >= 1000
genus.filt.harsh = genus.mat[genus_mask_harsh,]
#69 out of 409 families remain

##Filter to only look at samples in which at least 21 of the remaining 409 genuses show up
samp_mask = colSums(genus.mat > 0) >= 21
genus.mat = genus.mat[,samp_mask]
#Gets rid of 8 samples!

##Filter to only look at samples in which at least 21 of the remaining 151 genuses show up
samp_mask = colSums(genus.filt > 0) >= 21
genus.filt = genus.filt[,samp_mask]
#Gets rid of 10 samples!

##Filter to only look at samples in which at least 21 of the remaining 69 genuses show up
samp_mask = colSums(genus.filt.harsh > 0) >= 21
genus.filt.harsh = genus.filt.harsh[,samp_mask]
#Gets rid of 17 samples!


#Normalize with log ratio using Unknown_Family as reference
#In order to do this and avoid negative infinite values, have to get rid of 0's by adding pseudo counts
genus.smooth = genus.mat + 0.5
genus.filt.smooth = genus.filt + 0.5
genus.filt.harsh.smooth = genus.filt.harsh + 0.5

genus.lr = t( log(t(genus.smooth) / genus.smooth["Unknown",]) )
genus.lr.filt = t( log(t(genus.filt.smooth) / genus.filt.smooth["Unknown",]) )
genus.lr.harsh = t( log(t(genus.filt.harsh.smooth) / genus.filt.harsh.smooth["Unknown",]) )

#Now that we've normalized, safe to get rid of unknown genus
genus.lr = genus.lr[rownames(genus.lr)!="Unknown",]
genus.lr.filt = genus.lr.filt[rownames(genus.lr.filt)!="Unknown",]
genus.lr.harsh = genus.lr.harsh[rownames(genus.lr.harsh)!="Unknown",]


### Adding column annotation to heat map based on sample groups ########

#This code can be re-written in more simple/readable way but works for now
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

#Create column names based on column numbers
samp_never_fumigated = c(sapply(never_fumigated, make_sample))
samp_recently = c(sapply(recently, make_sample))
samp_more = c(sapply(more_than_a_month, make_sample))

all_samps = colnames(genus.df)[-1]

#Create dataframe for annotation
samp_annotation = data.frame(`fumigation status` = factor(case_when(
    all_samps %in% samp_never_fumigated ~ "never",
    all_samps %in% samp_recently ~ "recent",
    all_samps %in% samp_more ~ "past"
  ))) %>%
  mutate(Time = samp_metadata$Time)
row.names(samp_annotation) = all_samps

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
                     cluster_cols = clustercol, 
                     fontsize_row = 3)
  setHook("grid.newpage", NULL, "replace")
  grid.text(xlabel, y=-0.07, gp=gpar(fontsize=16))
  grid.text(ylabel, x=-0.07, rot=90, gp=gpar(fontsize=16))
  return (heatmap)
}


plt1 = create_heatmap(genus.lr,
                      annotation = samp_annotation,
                      rowname = FALSE, 
                      title = "Log ratio normalized abundances", 
                      xlabel = "232 Samples", 
                      ylabel = "408 fungal genuses"
)

plt1.filt = create_heatmap(genus.lr.filt,
                           annotation = samp_annotation,
                           title = "Log ratio normalized abundances", 
                           xlabel = "230 Samples", 
                           ylabel = "150 most abundant fungal genuses"
)


plt1.harsh = create_heatmap(genus.lr.harsh,
                            annotation = samp_annotation,
                            title = "Log ratio normalized abundances", 
                            xlabel = "223 samples", 
                            ylabel = "68 most abundant fungal family"
)





