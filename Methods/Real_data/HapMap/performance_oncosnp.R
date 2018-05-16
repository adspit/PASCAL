library(gtools)
library(reshape2)
library(xlsx)
library("splitstackshape")
require(data.table)
library(plyr)
library(compareDF)
library(ggplot2)


# read Hapmap 'gold standard'
hapmap.ground.truth <- read.table("hm3_cnv_submission.txt", sep="\t", h=T, stringsAsFactors = F)
hapmap.ground.truth[1:3, 1:5]
# cnp_id chr    start      end NA06984
# 1 HM3_CNP_1   1  8105049  8112441       2
# 2 HM3_CNP_2   1 10292133 10300570       2
# 3 HM3_CNP_3   1 10466423 10467633       2
dim(hapmap.ground.truth)
# [1]  856 1188

# sample match
sample.info <- read.xlsx("HAPMAP_samples_SNP_aCGH_comp_passing_cels_sample_map.xlsx",1,h=F, stringsAsFactors=F)
colnames(sample.info) <- c("ID","Sample.id")
head(sample.info)
#         ID                                               Sample.id
# 435 NA18943 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A02_31426.CEL
# 439 NA18947 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A03_31442.CEL
# 436 NA18944 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A04_31458.CEL
# 292 NA18521 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A05_31474.CEL
# 399 NA18858 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A06_31490.CEL
# 567 NA19141 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A07_31506.CEL

dim(sample.info)
# [1] 1278    2

# select only the SHELF patients since we only have Affymetrix calls for them
sample.info = sample.info[grep("SHELF",sample.info$Sample.id),]
dim(sample.info)
# [1] 98  2
# srt samples by id
sample.info <- sample.info[order(sample.info$Sample.id),]
# remove duplicate rows
sample.info <- unique(sample.info)
dim(sample.info)
# [1] 91  2

# select only the calls and update sample names
hapmap.agilent.calls <- cbind(hapmap.ground.truth[,1:4],hapmap.ground.truth[,intersect(sample.info$ID,colnames(hapmap.ground.truth))])
dim(hapmap.agilent.calls)
# [1] 856  85

rm(hapmap.ground.truth)
# we only have Agilent calls for 85 samples - we will use these as the 'gold standard'


# preprocess the available hapmap calls
# check if we have a perfect match for the HapMap sample names
all.equal.character(colnames(hapmap.agilent.calls)[5:ncol(hapmap.agilent.calls)], sample.info$ID[match(colnames(hapmap.agilent.calls)[5:ncol(hapmap.agilent.calls)], sample.info$ID)])
# [1] TRUE

# update HapMap Agilent sample names
colnames(hapmap.agilent.calls)[5:ncol(hapmap.agilent.calls)] = sample.info$Sample.id[match(colnames(hapmap.agilent.calls)[5:ncol(hapmap.agilent.calls)], sample.info$ID)]

# create hapmap agilent regions object
hapmap.agilent.regions = data.table(hapmap.agilent.calls[,2:4])

### OncoSNP

# oncosnp result files
OncoSNP.files <- mixedsort(dir(path = "hapmap_rez/", pattern = ".1.ggs"))
length(OncoSNP.files)
# 91
# read in oncosnp results on HapMap
OncoSNP.results = mclapply(OncoSNP.files, function(x)  { aux = read.table(x, sep="\t", h=T, stringsAsFactors=F)}, mc.cores = 10)

# select only the Rank1 calls
onco.snp.output = mclapply(OncoSNP.results, function(x) {x = x[x$Rank==1,]}, mc.cores=10)

# check if the probe order matches for all elements in the list
probe_names = onco.snp.output[[1]]$ProbeID
head(probe_names)
# [1] "CN_473963" "CN_473964" "CN_473965" "CN_477984" "CN_473981" "CN_473982"
length(probe_names)
# [1] 1844399

sapply(onco.snp.output, function(x) all.equal.character(x[,1], probe_names))
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [31] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [46] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [61] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [76] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [91] TRUE

# create object that contains only the calls
names(onco.snp.output) = OncoSNP.files
oncosnp_calls = lapply(1:length(OncoSNP.files), function(x) {colnames(onco.snp.output[[x]])[7] = names(onco.snp.output)[x]; onco.snp.output[[x]][7]})
# convert list to matrix
oncosnp_c = do.call(cbind,oncosnp_calls)
rownames(oncosnp_c) = probe_names
colnames(oncosnp_c) <- gsub(".1.ggs",".CEL",colnames(oncosnp_c))

# collapse to -1, 0 and 1
oncosnp_c[oncosnp_c  < 2] <- -1
oncosnp_c[oncosnp_c == 2] <-  0
oncosnp_c[oncosnp_c  > 2] <-  1

# find probes that overlap Agilent profiled regions
oncosnp.regions = data.table(chr = onco.snp.output[[1]]$Chromosome, start = onco.snp.output[[1]]$Position-12, end = onco.snp.output[[1]]$Position+12)
oncosnp.regions = oncosnp.regions[oncosnp.regions$chr<23,]
dim(oncosnp.regions)
# [1] 1752071       3

head(oncosnp.regions)
# chr start   end
# 1:   1 51586 51610
# 2:   1 51659 51683
# 3:   1 51674 51698
# 4:   1 52003 52027
# 5:   1 52771 52795
# 6:   1 52788 52812

# save.image("oncosnp_hapmap_results_process.Rdata")


setkey(hapmap.agilent.regions)
# find the regions that were profiles by both Affymetrix SNP array and Agilent 
overlapping.regions.oncosnp.agilent = data.frame(foverlaps(oncosnp.regions, hapmap.agilent.regions, type="within", nomatch=0L))
dim(overlapping.regions.oncosnp.agilent)
# [1] 14640     5
overlapping.regions.oncosnp.agilent.indices = foverlaps(oncosnp.regions, hapmap.agilent.regions, type="within", nomatch=0L, which = T)

hapmap_oncosnp_calls = oncosnp_c[overlapping.regions.agilent.cghcall.indices$xid, intersect(colnames(i.CN), colnames(oncosnp_c))]
dim(hapmap_oncosnp_calls)
# [1] 14640    81
all.equal.character(colnames(hapmap_oncosnp_calls), colnames(hapmap.agilent.calls))
# [1] TRUE
hapmap_oncosnp_calls = as.data.frame(hapmap_oncosnp_calls)
hapmap_oncosnp_calls = data.frame(lapply(hapmap_oncosnp_calls , factor , levels = c(-1,0,1)))
dim(hapmap_oncosnp_calls)
# [1] 14640    81

confusion.matrix.oncosnp <- as.data.frame(sapply(1:ncol(i.CN), function(i) table(i.CN[,i],hapmap_oncosnp_calls[,i])))
# rows represent the initial data
# cols represent the predicted calls 
colnames(confusion.matrix.oncosnp) <- colnames(hapmap_oncosnp_calls)

confusion.matrix.oncosnp <- cbind(as.data.frame(table(i.CN[,1],hapmap_oncosnp_calls[,1]))[,1:2], confusion.matrix.oncosnp)

colnames(confusion.matrix.oncosnp)[1:2] <- c("initial", "oncoSNP")

# true states
true.i.oncosnp <- as.numeric(levels(confusion.matrix.oncosnp[,1]))[confusion.matrix.oncosnp[,1]]==as.numeric(levels(confusion.matrix.oncosnp[,2]))[confusion.matrix.oncosnp[,2]]

# accuracy method
accuracy.oncoSNP<- apply(confusion.matrix.oncosnp[,-c(1,2)], 2, function(x)sum(x[true.i.oncosnp])/sum(x))


# precision for oncosnp

precision.oncoSNP.loss <- as.numeric( confusion.matrix.oncosnp[1,-c(1,2)] / apply(confusion.matrix.oncosnp[,-c(1,2)],2,function(x) sum(x[1:3])))
precision.oncoSNP.loss[is.nan(precision.oncoSNP.loss)] <- 1

precision.oncoSNP.normal <- as.numeric(confusion.matrix.oncosnp[5,-c(1,2)] / apply(confusion.matrix.oncosnp[,-c(1,2)],2,function(x) sum(x[4:6])))
which(is.nan(precision.oncoSNP.normal))
# numeric(0)

precision.oncoSNP.gain <- as.numeric(confusion.matrix.oncosnp[9,-c(1,2)] / apply(confusion.matrix.oncosnp[,-c(1,2)],2,function(x) sum(x[7:9])))
precision.oncoSNP.gain[is.nan(precision.oncoSNP.gain)] <- 1

#recall
recall.oncoSNP.loss <- as.numeric(confusion.matrix.oncosnp[1,-c(1,2)] / apply(confusion.matrix.oncosnp[,-c(1,2)],2,function(x) sum(x[c(1,4,7)])))
recall.oncoSNP.normal <- as.numeric(confusion.matrix.oncosnp[5,-c(1,2)] / apply(confusion.matrix.oncosnp[,-c(1,2)],2,function(x) sum(x[c(2,5,8)])))
recall.oncoSNP.gain <- as.numeric(confusion.matrix.oncosnp[9,-c(1,2)] /apply(confusion.matrix.oncosnp[,-c(1,2)],2,function(x) sum(x[c(3,6,9)])))


# compute F1 score for oncoSNP
F.oncoSNP.gain <- as.numeric((2 * precision.oncoSNP.gain * recall.oncoSNP.gain) / (precision.oncoSNP.gain + recall.oncoSNP.gain))
which(is.nan(F.oncoSNP.gain))
# [1] 74
# 2 * precision.oncoSNP.gain[74] * recall.oncoSNP.gain[74] 0
# precision.oncoSNP.gain[74] + recall.oncoSNP.gain[74] 0
F.oncoSNP.gain[is.nan(F.oncoSNP.gain)] = 0

F.oncoSNP.loss <- as.numeric((2 * precision.oncoSNP.loss * recall.oncoSNP.loss) / (precision.oncoSNP.loss + recall.oncoSNP.loss))
which(is.nan(F.oncoSNP.loss))
# integer(0)

F.oncoSNP.normal <- as.numeric((2 * precision.oncoSNP.normal * recall.oncoSNP.normal) / (precision.oncoSNP.normal + recall.oncoSNP.normal))

F.score.oncoSNP = data.frame(loss = F.oncoSNP.loss, normal = F.oncoSNP.normal, gain = F.oncoSNP.gain, Method = "OncoSNP") 

# save(F.score.oncoSNP, file = 'F.score.oncosnp.Rdata')




