### Performs JRmGRN network inference ###
library(tidyverse)
library(here)
library(NetCoMi)


# Fungus ####

#Read in cleaned data
fung_abundance = readRDS(here("data", "clean", "fung_abundance.RDS"))
samp_metadata = readRDS(here("data", "clean", "samp_metadata.RDS"))


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
taxa_filt = colSums(count_data_no_missing > 0) >= 150
counts_filt = count_data_no_missing[, taxa_filt]

# Corr 1: dominant signal #############

fume = samp_metadata$Treatment == "Non-fumigated chipping grass" | samp_metadata$Time == "Day_0"

net_pears <- netConstruct(counts_filt,
                          measure = "pearson",
                          filtSamp = "numbTaxa",
                          filtSampPar = list(numbTaxa = 20),
                          normMethod = "clr",
                          sparsMethod = "threshold",
                          thresh = 0.4,
                          verbose = 3)

props_pears <- netAnalyze(net_pears, 
                          clustMethod = "cluster_fast_greedy")

tests = apply(counts_filt, MARGIN=2, 
              FUN=function(x) wilcox.test(x~fume, alternative = "less"))
p_vals = sapply(tests, function(x) x$p.value)
abundant_in_control = head(names(sort(p_vals)), 20)
abundant_in_fume = tail(names(sort(p_vals)), 20)

status = factor(case_when(
  colnames(counts_filt) %in% abundant_in_control ~ "Abundant in Control",
  colnames(counts_filt) %in% abundant_in_fume ~ "Abundant in Fumigated",
  TRUE ~ "Neither"
))
names(status) = colnames(counts_filt)

library(RColorBrewer)
statcol = brewer.pal(n=3, "Dark2")

par(oma=c(0, 0, 0, 5))
plot(props_pears, 
     nodeColor = "feature", 
     featVecCol = status,
     colorVec =  statcol,
     title1 = "Network on Genus level based on Correlations", 
     charToRm = "g__",
     showTitle = TRUE,
     cexLabels = 1.5,
     rmSingles = TRUE)


statcol_transp <- colToTransp(statcol, 60)
legend(x=.75, y = .95, cex=.7,
       legend=levels(status),pch=19, col = statcol_transp) 



#Corr based broken into two #### 



net_pears2 <-  netConstruct(counts_filt,
                            group = fume,
                            measure = "pearson",
                            filtSamp = "numbTaxa",
                            filtSampPar = list(numbTaxa = 20),
                            normMethod = "clr",
                            sparsMethod = "threshold",
                            thresh = 0.4,
                            verbose = 3)




props_pears2 <- netAnalyze(net_pears2, 
                           clustMethod = "cluster_fast_greedy")



par(oma=c(0, 0, 0, 5))
plot(props_pears2, 
     nodeColor = "feature", 
     featVecCol = status,
     colorVec =  statcol,
     title1 = "Network on Genus level based on Correlations",
     groupNames = c("Fumigant", "Control"),
     charToRm = "g__",
     showTitle = TRUE,
     cexLabels = 1.5,
     rmSingles = "inboth",
     sameLayout = TRUE,
     layoutGroup = 2)


statcol_transp <- colToTransp(statcol, 60)
legend(x=-.2, y = 1.05, cex=.7,
       legend=levels(status),pch=19, col = statcol_transp)

# Spring 1 ####
  


net_spring <- netConstruct(counts_filt,
                           measure = "spring",
                           measurePar = list(nlambda=10, 
                                             rep.num=10),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 3,
                           seed = 123456)


props_spring <- netAnalyze(net_spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

p <- plot(props_spring, 
          nodeColor = "feature", 
          featVecCol = status,
          colorVec =  statcol,
          title1 = "Network on Genus level based on SPRING", 
          charToRm = "g__",
          showTitle = TRUE,
          cexLabels = 1.5)

legend(x=.70, y = .75, cex=.7,
       legend=levels(status),pch=19, col = statcol_transp) 

net_spring$edgelist1




#Spring broken into 2 #######


net_spring2 <- netConstruct(counts_filt,
                           group = fume,
                           filtSamp = "numbTaxa",
                           filtSampPar = list(numbTaxa = 20),
                           measure = "spring",
                           measurePar = list(nlambda=15, 
                                             rep.num=10,
                                             ncores=4),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 3,
                           seed = 123456)

props_spring2 <- netAnalyze(net_spring2, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)


par(oma=c(0, 0, 0, 5))
plot(props_spring2, 
     nodeColor = "feature", 
     featVecCol = status,
     colorVec =  statcol,
     title1 = "Network on Genus level based on Correlations",
     groupNames = c("Fumigant", "Control"),
     charToRm = "g__",
     showTitle = TRUE,
     cexLabels = 1.5,
     rmSingles = "inboth",
     sameLayout = TRUE,
     layoutGroup = 2)


statcol_transp <- colToTransp(statcol, 60)
legend(x=-.2, y = .95, cex=.7,
       legend=levels(status), pch = 19, col = statcol_transp)

summary(props_spring2)


# JGL Mod ##

library(JGL)
X1 = net_pears2$normCounts1
X2 = net_pears2$normCounts2

jgl.mod = JGL(list(X1, X2), penalty = "fused", lambda1 = .25, lambda2 = .1)
#Theta attribute of jgl.mod is list of two estimated precision matrices


net_jgl2 = netConstruct(jgl.mod$theta[[1]], jgl.mod$theta[[2]], dataType = "condDependence", 
                    sparsMethod = "none")

props_jgl2 = netAnalyze(net_jgl2, 
                            centrLCC = TRUE,
                            clustMethod = "cluster_fast_greedy",
                            hubPar = "eigenvector",
                            weightDeg = FALSE, normDeg = FALSE)


par(oma=c(0, 0, 0, 5))
plot(props_jgl2, 
     nodeColor = "feature", 
     featVecCol = status,
     colorVec =  statcol,
     title1 = "Network on Genus level based on Joint GLASSO",
     groupNames = c("Fumigant", "Control"),
     charToRm = "g__",
     showTitle = TRUE,
     cexLabels = 1.5,
     rmSingles = "inboth",
     sameLayout = TRUE,
     layoutGroup = 1)


statcol_transp <- colToTransp(statcol, 60)
legend(x=-.2, y = .95, cex=.6,
       legend=levels(status), pch = 19, col = statcol_transp)

summary(props_jgl2)


## Comparing between samples with fumigant versus without

diff_spring = netCompare(props_spring2, permTest = FALSE, 
                      verbose = TRUE,
                      seed = 123456)

diff_jgl = netCompare(props_jgl2, permTest = FALSE, 
                    verbose = TRUE,
                    seed = 123456)

summary(diff_jgl, groupNames = c("Fumigated", "Control"))

summary(diff_spring, groupNames = c("Fumigated", "Control"))



# Now work with interaction of time and fumigant treatment ######

#Fit pearson correlation network for each subset of data
#just to get normalized data for JGLASSO


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
    TRUE ~ "Fumigated"))) %>%
  mutate(f = droplevels(interaction(trt_coarse, time_coarse))) %>%
  mutate(f = fct_collapse(f, Control.Base = c("Control.Base", "Control.Seeding")))

samples_observed = data.frame(counts_filt) %>%
  mutate(f = samp_metadata$f) %>%
  group_by(f) %>%
  summarise(across(.fns = ~ sum(. > 0 ))) 

taxa_to_keep = which(colMeans(df[,-1] > 5) == 1)

counts_filt2 = counts_filt[, taxa_to_keep]
counts_filt_split = split.data.frame(counts_filt2, samp_metadata$f)


pearson_nets = map(counts_filt_split, ~ netConstruct(.,
                                      measure = "pearson",
                                      filtTax = "numbSamp",
                                      filtTaxPar = list("numbSamp"=5),
                                      filtSamp = "numbTaxa",
                                      filtSampPar = list(numbTaxa = 2),
                                      normMethod = "clr",
                                      sparsMethod = "threshold",
                                      thresh = 0.4,
                                      verbose = 3))
norm_data = map(pearson_nets, ~ .[["normCounts1"]])


jgl.mod2 = JGL(norm_data, penalty = "fused", lambda1 = .25, lambda2 = .1)



