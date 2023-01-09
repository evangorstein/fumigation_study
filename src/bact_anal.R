### Performs JRmGRN network inference ###
library(tidyverse)
library(here)
library(NetCoMi)


# Bacteria ###


#Read in cleaned data
bac_abundance = readRDS(here("data", "clean", "bac_abundance.RDS"))
samp_metadata = readRDS(here("data", "clean", "samp_metadata.RDS"))
species_nmds2 = readRDS(here("data", "clean", "species_scores2_bac.RDS")) %>%
  mutate(cluster = fct_recode(cluster, `Not fumigated` = "Cluster1", `Recently fumigated` = "Cluster2", `Fumigated months ago` = "Cluster3"))
samples_nmds = readRDS(here("data", "clean", "sample_scores_bac.RDS"))


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

#Get rid of single sample that doesn't meet the criteria in count_data_no_missing, as well
count_data_no_missing = count_data_no_missing[samp_filt,]

#NMDS results

#Create samples_plot
samples_plot = ggplot(samples_nmds, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = trt_coarse)) +
  theme_bw() +
  scale_color_brewer(name = "Sample Status", 
                     palette = "Dark2", 
                     labels = c("Not\nfumigated", "Recently\nfumigated", "Fumigated\nmonths ago"))


#Add species to samples plot

add_species <- function(p, df) {
  
  p + geom_segment(data = df, 
                   aes(x = 0, xend=NMDS1, y=0, yend=NMDS2, color = cluster), 
                   arrow = arrow(length = unit(0.25, "cm")), lwd=0.1) + #add arrows for species
    ggrepel::geom_text_repel(data = df, aes(x=NMDS1, y=NMDS2, label = abrev), cex = 2, direction = "both", segment.size = 0.25) 
  
}

add_species(samples_plot, species_nmds)
add_species(samples_plot, species_nmds2)



# Corr 1: dominant signal #############

net_pears <- netConstruct(counts_filt,
                          measure = "pearson",
                          normMethod = "clr",
                          zeroMethod = "multRepl",
                          sparsMethod = "threshold",
                          thresh = 0.70,
                          dissFunc = "signed",
                          verbose = 3)


props_pears <- netAnalyze(net_pears, 
                          clustMethod = "cluster_fast_greedy")

status = species_nmds2$cluster
names(status) = colnames(counts_filt)
library(RColorBrewer)
statcol = brewer.pal(n=3, "Dark2")
statcol_transp <- colToTransp(statcol, 60)

plot(props_pears, 
          nodeColor = "feature", 
          featVecCol = status,
          colorVec =  statcol,
          labelLength = 6,
          charToRm = "g__",
          cexLabels = 1.5,
          rmSingles = TRUE)

legend(x=-.65, y =-.4, cex=.55,
       legend=levels(status),pch=19, col = statcol_transp) 


# Spring 1 ####

net_spring <- netConstruct(counts_filt,
                           measure = "spring",
                           measurePar = list(nlambda=15, 
                                             rep.num=15, 
                                             ncores=4, 
                                             thresh=.05),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 3,
                           seed = 1)  


# net_spring <- netConstruct(data = net_spring$assoMat1,
#                                        dataType = "condDependence",
#                                        sparsMethod = "none",
#                                        dissFunc = "unsigned",
#                                        verbose = 0)

props_spring <- netAnalyze(net_spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

p <- plot(props_spring, 
          nodeColor = "feature", 
          featVecCol = status,
          colorVec =  statcol,
          #edgeFilter = "highestWeight",
          #edgeFilterPar = 150, #Include only edges with top 150 weights
          labelLength = 8,
          charToRm = "g__",
          rmSingles = TRUE,
          cexLabels = 1.5)

legend(x=.5, y = .85, cex=.55,
       legend=levels(status),pch=19, col = statcol_transp) 




## Separate Spring fits #####


### Start with no fumigant samples #####
counts_nf = count_data_no_missing[samples_nmds$trt_coarse == "NoFum",]

#Filter to only look at taxa appearing in at least two thirds of the samples (54 out of 80)
taxa_filt = colSums(counts_nf > 0) >= 54
counts_nf = counts_nf[,taxa_filt]

#Filter to only look at samples in which at least 50 of the remaining taxa show up
samp_filt = rowSums(counts_nf > 0) >= 50
counts_nf = counts_nf[samp_filt,]


net_spring_nf <- netConstruct(counts_nf,
                            measure = "spring",
                            measurePar = list(nlambda=15, 
                                              rep.num=15,
                                              ncores=4, 
                                              thresh = .05),
                            normMethod = "none", 
                            zeroMethod = "none",
                            sparsMethod = "none", 
                            dissFunc = "signed",
                            verbose = 3,
                            seed = 123456)

props_spring_nf <- netAnalyze(net_spring_nf, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)



status = species_nmds2$cluster[match(colnames(counts_nf), species_nmds2$species)] %>%
  fct_explicit_na("Missing in overall network")
names(status) = colnames(counts_nf)
statcol = brewer.pal(n=4, "Dark2")
statcol_transp <- colToTransp(statcol, 60)
plot(props_spring_nf, 
          nodeColor = "feature", 
          featVecCol = status,
          colorVec =  statcol,
          labelLength = 8,
          title1 = "Network of Bacteria Families estimated with SPRING\nfrom only the non-fumigated samples", 
          charToRm = "g__",
          rmSingles = TRUE,
          cexLabels = 1.5)
legend(x=.65, y = .88, cex=.4,
       legend=levels(status),pch=19, col = statcol_transp) 



#####Now fumigant samples, no time #########

counts_f = count_data_no_missing[samples_nmds$trt_coarse == "RecentlyFum",]

#Filter to only look at taxa appearing in at least two thirds of the samples (22 out of 32)
taxa_filt = colSums(counts_f > 0) >= 22
counts_f = counts_f[,taxa_filt]

#Filter to only look at samples in which at least 50 of the remaining taxa show up
samp_filt = rowSums(counts_f > 0) >= 50
counts_f = counts_f[samp_filt,]


net_spring_f <- netConstruct(counts_f,
                            measure = "spring",
                            measurePar = list(nlambda=15, 
                                              rep.num=15,
                                              ncores=4,
                                              thresh=.05),
                            normMethod = "none", 
                            zeroMethod = "none",
                            sparsMethod = "none", 
                            dissFunc = "signed",
                            verbose = 3,
                            seed = 123456)

props_spring_f <- netAnalyze(net_spring_f, 
                            centrLCC = TRUE,
                            clustMethod = "cluster_fast_greedy",
                            hubPar = "eigenvector",
                            weightDeg = FALSE, normDeg = FALSE)


status = species_nmds2$cluster[match(colnames(counts_f), species_nmds2$species)] %>%
  fct_explicit_na("Missing in overall network")
names(status) = colnames(counts_f)

par(xpd=TRUE)
plot(props_spring_f, 
          nodeColor = "feature", 
          featVecCol = status,
          colorVec =  statcol,
          labelLength = 8,
          edgeFilter = "highestWeight",
          edgeFilterPar = 150, #Include as many edges as in the network fit from all samples
          title1 = "Network of Bacteria Families estimated with SPRING\nfrom only the non-fumigated samples", 
          charToRm = "g__",
          rmSingles = TRUE,
          cexLabels = 1.5)
legend(x=.80, y = .88, cex=.4,
       legend=levels(status),pch=19, col = statcol_transp)  

##### Samples that have been fumigated and time has passed ####
counts_ft = count_data_no_missing[samples_nmds$trt_coarse == "FumMonthsAgo",]

#Filter to only look at taxa appearing in at least two thirds of the samples (85 out of 127)
taxa_filt = colSums(counts_ft > 0) >= 85
counts_ft = counts_ft[,taxa_filt]

#Filter to only look at samples in which at least 50 of the remaining taxa show up
samp_filt = rowSums(counts_ft > 0) >= 50
counts_ft = counts_ft[samp_filt,]


net_spring_ft <- netConstruct(counts_ft,
                             measure = "spring",
                             measurePar = list(nlambda=15, 
                                               rep.num=15,
                                               ncores=4,
                                               thresh=.05),
                             normMethod = "none", 
                             zeroMethod = "none",
                             sparsMethod = "none", 
                             dissFunc = "signed",
                             verbose = 3,
                             seed = 123456)

props_spring_ft <- netAnalyze(net_spring_ft, 
                             centrLCC = TRUE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = "eigenvector",
                             weightDeg = FALSE, normDeg = FALSE)


status = species_nmds2$cluster[match(colnames(counts_ft), species_nmds2$species)] %>%
  fct_explicit_na("Missing in overall network")
names(status) = colnames(counts_ft)
p <- plot(props_spring_ft, 
          nodeColor = "feature", 
          featVecCol = status,
          colorVec =  statcol,
          labelLength = 8,
          edgeFilter = "highestWeight",
          edgeFilterPar = 150, #Include as many edges as in the network fit from all samples
          title1 = "Network of Bacteria Families estimated with SPRING\nfrom only samples that were fumigated months prior", 
          charToRm = "g__",
          rmSingles = TRUE,
          cexLabels = 2.5)
legend(x=.80, y = .88, cex=.4,
       legend=levels(status),pch=19, col = statcol_transp) 



#Joint GLASSO
counts_filt_split = split.data.frame(counts_filt, samples_nmds$trt_coarse)

#Fit Pearson to each subset of the data to get normalized counts for JGlasso
pearson_nets = map(counts_filt_split, ~ netConstruct(.,
                                                     measure = "pearson",
                                                     normMethod = "clr",
                                                     sparsMethod = "threshold",
                                                     thresh = 0.4,
                                                     verbose = 3))
norm_data = map(pearson_nets, ~ .[["normCounts1"]])

library(JGL)
jgl.mod = JGL(norm_data, penalty = "fused", lambda1 = .25, lambda2 = .1)



net_jgl_nf = netConstruct(data = jgl.mod$theta[[1]],
                       dataType = "condDependence",
                       sparsMethod = "none",
                       verbose = 0)


props_jgl_nf = netAnalyze(net_jgl_nf, 
                         centrLCC = TRUE,
                         clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector",
                         weightDeg = FALSE, normDeg = FALSE)


status = species_nmds2$cluster
names(status) = colnames(counts_filt)
library(RColorBrewer)
statcol = brewer.pal(n=3, "Dark2")
statcol_transp <- colToTransp(statcol, 60)

p_nf = plot(props_jgl_nf, 
     nodeColor = "feature", 
     featVecCol = status,
     colorVec =  statcol,
     labelLength = 8,
     charToRm = "g__",
     cexLabels = 1.2,
     rmSingles = TRUE)

legend(x=-.90, y =-.3, cex=.55,
       legend=levels(status),pch=19, col = statcol_transp) 

net_jgl_f = netConstruct(data = jgl.mod$theta[[2]],
                          dataType = "condDependence",
                          sparsMethod = "none",
                          verbose = 0)


props_jgl_f = netAnalyze(net_jgl_f, 
                          centrLCC = TRUE,
                          clustMethod = "cluster_fast_greedy",
                          hubPar = "eigenvector",
                          weightDeg = FALSE, normDeg = FALSE)

plot(props_jgl_f, 
     nodeColor = "feature", 
     featVecCol = status,
     colorVec =  statcol,
     labelLength = 8,
     charToRm = "g__",
     cexLabels = 1.2, 
     rmSingles = TRUE)

legend(x=-.90, y =.85, cex=.55,
       legend=levels(status),pch=19, col = statcol_transp) 




net_jgl_ft = netConstruct(data = jgl.mod$theta[[3]],
                          dataType = "condDependence",
                          sparsMethod = "none",
                          verbose = 0)


props_jgl_ft = netAnalyze(net_jgl_ft, 
                          centrLCC = TRUE,
                          clustMethod = "cluster_fast_greedy",
                          hubPar = "eigenvector",
                          weightDeg = FALSE, normDeg = FALSE)

plot(props_jgl_ft, 
     nodeColor = "feature", 
     featVecCol = status,
     colorVec =  statcol,
     labelLength = 8,
     charToRm = "g__",
     cexLabels = 1.2, 
     rmSingles = TRUE)

legend(x=-1, y =.85, cex=.55,
       legend=levels(status),pch=19, col = statcol_transp) 



library(nVennR)

net_nf <- graph.data.frame(net_jgl_nf$edgelist1, directed=F)
net_f <- graph.data.frame(net_jgl_f$edgelist1, directed=F)
net_ft <- graph.data.frame(net_jgl_ft$edgelist1, directed=F)

dev.off()
myV <- plotVenn(list(`No Fumigant`=as_ids(E(net_nf)), `Recently Fumigated`=as_ids(E(net_f)), `Fumigated a while ago` = as_ids(E(net_ft))))


