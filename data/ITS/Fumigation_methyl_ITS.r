#Have to set the directory path to download and install packages.....Standard directory is not writable
#.libPaths(new='~/R/x86_64-pc-linux-gnu-library/3.5.0')
#dir.create('~/R/x86_64-pc-linux-gnu-library/3.5.0', showWarnings = FALSE, recursive = TRUE)
# check that the path was set correctly
#.libPaths()

#look at the list of CRAN mirrors available....you want the cloud mirror
#getCRANmirrors()

# This will set the CRAN mirror to cloud
chooseCRANmirror(ind=2)

#Have to set the directory path to download and install packages.....Standard directory is not writable
dir.create('~/R/x86_64-pc-linux-gnu-library/3.5.0', showWarnings = FALSE, recursive = TRUE)
.libPaths(new='~/R/x86_64-pc-linux-gnu-library/3.5.0')
#dir.create('~/R/x86_64-pc-linux-gnu-library/3.5.0', showWarnings = FALSE, recursive = TRUE)
# check that the path was set correctly
.libPaths()



# Install required packages
install.packages("BiocManager")
#BiocManager::install("phyloseq")
#install.packages("tidyr")
#install.packages("vegan")
#install.packages("plyr")
#install.packages("dplyr")
#BiocManager::install("DESeq2")
BiocManager::install("dada2")
BiocManager::install("DECIPHER")


#library("phyloseq")
#library("tidyr")
#library("vegan")
#library("plyr")
#library("dplyr")
#library("DESeq2")
library("dada2")
library('DECIPHER')

# After your packages have been installed you can now run your code
path = "/home/kinkelll/mmillica/Wisconsin/Fumigation/Fumigation_Methyl_Bromide/Fungal"
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,truncLen=c(250,200),maxN=0,maxEE=c(2,2),truncQ=2,rm.phix=TRUE,compress=FALSE, multithread=TRUE) 
head(out) 
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e+20, randomize = TRUE) 
errR <- learnErrors(filtRs, multithread=TRUE, nbases = 1e+20, randomize = TRUE)
plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = TRUE) 
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = TRUE) 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
save.image(file="Seqtab.nochim_16S.RData")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, file = "/home/kinkelll/mmillica/Wisconsin/Fumigation/Fumigation_Methyl_Bromide/Fungal/DADA_Fumigation_methyl_ITS_tracking.csv")
save.image(file = "/home/kinkelll/mmillica/Wisconsin/Fumigation/Fumigation_Methyl_Bromide/Fungal/DADA_PreTaxa_Fumigation_methyl_ITS.RData")

rm(fnFs) 
rm(fnRs)
rm(sample.names)
rm(filtFs)
rm(filtRs)
rm(out)
rm(errF)
rm(errR)
rm(derepFs)
rm(derepRs)
rm(dadaFs)
rm(dadaRs)
rm(mergers)
rm(seqtab)
rm(getN)
rm(track)

taxa <- assignTaxonomy(seqtab.nochim, "/home/kinkelll/mmillica/Reference_database/sh_general_release_dynamic_s_04.02.2020.fasta", multithread=TRUE, tryRC = TRUE)

save(seqtab.nochim, taxa, file= "/home/kinkelll/mmillica/Wisconsin/Fumigation/Fumigation_Methyl_Bromide/Fungal/Fumigation_methyl_ITS.RData")


################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################


