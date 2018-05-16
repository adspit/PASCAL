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


# CGHcall
#load CGHcalls
load("hapmap_cghcall_profiles.Rdata")

# rename file
hapmap_cghcall_profiles = result0
rm(result0)
dim(fData(hapmap_cghcall_profiles))
# [1]  1844399       3
head(fData(hapmap_cghcall_profiles))
# Chromosome Start   End
# CN_473963          1 51586 51610
# CN_473964          1 51659 51683
# CN_473965          1 51674 51698
# CN_477984          1 52003 52027
# CN_473981          1 52771 52795
# CN_473982          1 52788 52812
length(unique(rownames(fData(hapmap_cghcall_profiles))))
# [1] 1844399


# overlapping ranges with CGHCall 
# update colnames so I perform the merge
colnames(fData(hapmap_cghcall_profiles)) = colnames(hapmap.agilent.regions)
# perform merge chromosome wise
dim(hapmap_cghcall_profiles)
# Features  Samples 
# 1844399       91


cghcall.regions = data.table(fData(hapmap_cghcall_profiles))
# rmeove sex chromosomes
cghcall.regions = cghcall.regions[cghcall.regions$chr<23,]
dim(cghcall.regions)
# [1] 1752071       3

cghcall.profiles = hapmap_cghcall_profiles[fData(hapmap_cghcall_profiles)$chr<23,  ]
dim(cghcall.profiles)
# Features  Samples 
# 1752071       91        


setkey(hapmap.agilent.regions)
# find the regions that were profiled by both Affymetrix SNP array and Agilent 
overlapping.regions.cghcall.agilent = data.frame(foverlaps(cghcall.regions, hapmap.agilent.regions, type="within", nomatch=0L))
dim(overlapping.regions.cghcall.agilent)
# [1] 14640     5
overlapping.regions.agilent.cghcall.indices = foverlaps(cghcall.regions, hapmap.agilent.regions, type="within", nomatch=0L, which = T)


hapmap_cghcalls = calls(cghcall.profiles)
dim(hapmap_cghcalls) 
# [1] 1752071      91
colnames(hapmap_cghcalls) = gsub(".Log.R.Ratio", "", colnames(hapmap_cghcalls))
# match samples and overlappoing regions
hapmap_cghcalls = hapmap_cghcalls[overlapping.regions.agilent.cghcall.indices$xid,intersect(colnames(hapmap_cghcalls), colnames(hapmap.agilent.calls))]
dim(hapmap_cghcalls)
# [1] 14640    81
all.equal.character(colnames(hapmap.agilent.calls), colnames(hapmap_cghcalls))
# [1] TRUE

range(hapmap_cghcalls)
# [1] -2  2

hapmap_cghcalls[hapmap_cghcalls<0] <- -1
hapmap_cghcalls[hapmap_cghcalls>0] <- 1
range(hapmap_cghcalls)
# [1] -1  1
hapmap_cghcalls = as.data.frame(hapmap_cghcalls)
hapmap_cghcalls <- data.frame(lapply(hapmap_cghcalls , factor , levels = c(-1,0,1)))

hapmap.agilent.cghcall.matched = hapmap.agilent.calls[overlapping.regions.agilent.cghcall.indices$yid, intersect(colnames(hapmap_cghcalls), colnames(hapmap.agilent.calls))]
dim(hapmap.agilent.cghcall.matched)
# [1] 14640    81

all.equal.character(colnames(hapmap.agilent.cghcall.matched), colnames(hapmap_cghcalls))
# TRUE

i.CN = hapmap.agilent.cghcall.matched
i.CN[i.CN<2] <- -1
i.CN[i.CN>2] <- 1
i.CN[i.CN==2] = 0
i.CN = data.frame(lapply(i.CN , factor , levels = c(-1,0,1)))

confusion.matrix.cghcall <- as.data.frame(sapply(1:ncol(i.CN), function(i) table(i.CN[,i],hapmap_cghcalls[,i])))
colnames(confusion.matrix.cghcall) <- colnames(hapmap_cghcalls)

# rows represent the initial data
# cols represent the CGH calls 


confusion.matrix.cghcall <- cbind(as.data.frame(table(i.CN[,1], hapmap_cghcalls[,1]))[,1:2], confusion.matrix.cghcall)

colnames(confusion.matrix.cghcall)[1:2] <- c("initial", "CGHcall")


# compute performance measurements for each class
# for class -1: loss

# true positives
# true.pos.loss.CGHcall <- confusion.matrix[1,1]
# 
# # false positives
# false.pos.loss.CGHcall <- sum(confusion.matrix[,1])-true.pos.loss.CGHcall
# 
# # false neg
# false.neg.loss.CGHcall <- sum(confusion.matrix[1,])-confusion.matrix[1,1]

true.i <- as.numeric(levels(confusion.matrix.cghcall[,1]))[confusion.matrix.cghcall[,1]]==as.numeric(levels(confusion.matrix.cghcall[,2]))[confusion.matrix.cghcall[,2]]

# accuracy method
accuracy.CGHcall<- apply(confusion.matrix.cghcall[,-c(1,2)], 2, function(x)sum(x[true.i])/sum(x))



precision.CGHcall.loss <- as.numeric(confusion.matrix.cghcall[1,-c(1,2)] / apply(confusion.matrix.cghcall[,-c(1,2)],2,function(x) sum(x[1:3])))
which(is.nan(precision.CGHcall.loss))
# integer(0)
precision.CGHcall.normal <- as.numeric(confusion.matrix.cghcall[5,-c(1,2)] / apply(confusion.matrix.cghcall[,-c(1,2)],2,function(x) sum(x[4:6])))
which(is.nan(precision.CGHcall.normal))
# integer(0)
precision.CGHcall.gain <- as.numeric(confusion.matrix.cghcall[9,-c(1,2)] / apply(confusion.matrix.cghcall[,-c(1,2)],2,function(x) sum(x[7:9])))
which(is.nan(precision.CGHcall.gain))
# integer(0)

#recall - specificity
recall.CGHcall.loss <- as.numeric(confusion.matrix.cghcall[1,-c(1,2)] / apply(confusion.matrix.cghcall[,-c(1,2)],2,function(x) sum(x[c(1,4,7)])))
which(is.nan(recall.CGHcall.loss))
# integer(0)
recall.CGHcall.normal <- as.numeric(confusion.matrix.cghcall[5,-c(1,2)] / apply(confusion.matrix.cghcall[,-c(1,2)],2,function(x) sum(x[c(2,5,8)])))
which(is.nan(recall.CGHcall.normal))
# integer(0)
recall.CGHcall.gain <- as.numeric(confusion.matrix.cghcall[9,-c(1,2)] /apply(confusion.matrix.cghcall[,-c(1,2)],2,function(x) sum(x[c(3,6,9)])))
which(is.nan(recall.CGHcall.gain))
# integer(0)



# compute F1 score for CGHcall
F.CGHcall.gain <- as.numeric((2 * precision.CGHcall.gain * recall.CGHcall.gain) / (precision.CGHcall.gain + recall.CGHcall.gain))
which(is.nan(F.CGHcall.gain))
F.CGHcall.loss <- as.numeric((2 * precision.CGHcall.loss * recall.CGHcall.loss) / (precision.CGHcall.loss + recall.CGHcall.loss))
which(is.nan(F.CGHcall.loss))
F.CGHcall.normal <- as.numeric((2 * precision.CGHcall.normal * recall.CGHcall.normal) / (precision.CGHcall.normal + recall.CGHcall.normal))
which(is.nan(F.CGHcall.normal))

F.score.cghcall = data.frame(loss = F.CGHcall.loss, normal = F.CGHcall.normal, gain = F.CGHcall.gain, Method = "CGHcall") 
save(F.score.cghcall, file = 'F.score.cghcall.Rdata')
