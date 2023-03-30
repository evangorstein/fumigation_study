# Creating heatmaps of bacterial data
# reference: https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
library(pheatmap)
library(tidyverse)
library(here)

#### Load and prepare data ####

# Read in data
bac_abundance = readRDS(here("data", "clean", "bac_abundance.RDS"))
samp_metadata = readRDS(here("data", "clean", "samp_metadata.RDS"))

# Aggregate to family level
fam.df = bac_abundance %>% 
  mutate(Family = fct_na_value_to_level(Family, "Unknown_Family")) %>% #Missing values included in Unknown_Family level
  group_by(Family) %>%
  summarise(across(starts_with("Samp"), sum)) 
#Note that the "Missing" value aggregates the counts for all the OTUs whose families are unknown


#Convert to matrix with family label as rownames 
fams = fam.df$Family
fam.mat = as.matrix(fam.df[,-1])
rownames(fam.mat) = fams

#Check most abundant families, which includes "Unknown_Family"
sort(rowSums(fam.mat))

#Filter to only look at families that have at least 1,000 counts in total (summing over samples) 
fam_mask = rowSums(fam.mat) >= 1000
fam.filt = fam.mat[fam_mask,]
#153 out of 296 families remain

#Filter to only look at families which appear in at least 230 out of 240 samples
fam_mask_harsh = rowSums(fam.mat) >= 10000
fam.filt.harsh = fam.mat[fam_mask_harsh,]
#74 out of 296 families remain


#Filter to only look at samples in which at least 50 of the remaining 152 families show up
#samp_mask = colSums(fam.filt > 0) >= 50
#fam.filt = fam.filt[,samp_mask]
#Gets rid of a single sample (sample 173)

#Normalize by dividing by library sizes
#fam.ls = t( t(fam.filt) / colSums(fam.filt) )
#fam.ls.harsh = t( t(fam.filt.harsh) / colSums(fam.filt.harsh) )
#Get rid of missing family 
#fam.ls.nm = fam.ls[-nrow(fam.ls),]
#fam.ls.harsh.nm = fam.ls.harsh[-nrow(fam.ls.harsh),]


#Normalize with log ratio using Unknown_Family as reference
#In order to do this and avoid negative infinite values, have to get rid of 0's by adding pseudo counts
fam.smooth = fam.mat + 0.5
fam.filt.smooth = fam.filt + 0.5
fam.filt.harsh.smooth = fam.filt.harsh + 0.5

fam.lr = t( log(t(fam.smooth) / fam.smooth["Unknown_Family",]) )
fam.lr.filt = t( log(t(fam.filt.smooth) / fam.filt.smooth["Unknown_Family",]) )
fam.lr.harsh = t( log(t(fam.filt.harsh.smooth) / fam.filt.harsh.smooth["Unknown_Family",]) )

#Now that we've normalized, safe to get rid of unknown family
fam.lr = fam.lr[rownames(fam.lr)!="Unknown_Family",]
fam.lr.filt = fam.lr.filt[rownames(fam.lr.filt)!="Unknown_Family",]
fam.lr.harsh = fam.lr.harsh[rownames(fam.lr.harsh)!="Unknown_Family",]



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


all_samps = colnames(fam.mat)

#Create dataframe for annotation
samp_annotation = data.frame(all_samps) %>% 
  mutate(`fumigation status` = factor(case_when(
    all_samps %in% samp_never_fumigated ~ "never",
    all_samps %in% samp_recently ~ "recent",
    all_samps %in% samp_more ~ "past"
  ))) %>% 
  select(`fumigation status`)
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


plt1 = create_heatmap(fam.lr,
                      annotation = samp_annotation,
                      rowname = FALSE, 
                      title = "Log ratio normalized abundances", 
                      xlabel = "Samples 1-240", 
                      ylabel = "Bacterial family"
)

plt1.filt = create_heatmap(fam.lr.filt,
                         annotation = samp_annotation,
                         title = "Log ratio normalized abundances", 
                         xlabel = "Samples 1-240", 
                         ylabel = "Bacterial family"
)


plt1.harsh = create_heatmap(fam.lr.harsh,
                            annotation = samp_annotation,
                            title = "Log ratio normalized abundances", 
                            xlabel = "Samples 1-240", 
                            ylabel = "Bacterial family"
)




## Work at order level
order.df = bac_abundance %>% 
  mutate(Order = fct_na_value_to_level(Order, "Unknown_Order")) %>% #Missing values included in Unknown_Order level
  group_by(Order) %>%
  summarise(across(starts_with("Samp"), sum)) 
#Note that the "Missing" value aggregates the counts for all the OTUs whose orders are unknown

#Convert to matrix with family label as rownames 
ords = order.df$Order
order.mat = as.matrix(order.df[,-1])
rownames(order.mat) = ords


## Heat map of grouped samples

non_fumigated = rowMeans(fam.lr[,samp_never_fumigated])
recent_fumigated = rowMeans(fam.lr[,samp_recently])
more_than_a_month_fumigated = rowMeans(fam.lr[,more_than_a_month])

## construct the matrix
group_fam = cbind(non_fumigated, recent_fumigated, more_than_a_month_fumigated)


plt_grp = create_heatmap(group_fam, 
                      colname = TRUE, 
                      title = "Average (log-ratio normalized) abundances by sample condition", 
                      xlabel = "Sample fumigation Status", 
                      ylabel = "Bacterial family", 
                      clustercol = FALSE)


