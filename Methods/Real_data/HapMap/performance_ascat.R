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


# load ascat results on HapMap 
load("ascat.CN.Rdata")
# update colnames
colnames(ascat.CN) = gsub(".Log.R.Ratio","", colnames(ascat.CN))
dim(ascat.CN)
# [1] 1831313      91
ascat.CN[1:5, 1:2]
# SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A02_31426.CEL SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A03_31442.CEL
# CN_473963                                                       0                                                       0
# CN_473964                                                       0                                                       0
# CN_473965                                                       0                                                       0
# CN_477984                                                       0                                                       0
# CN_473981                                                       0                                                       0

dim(hapmap.agilent.calls)
# [1] 856  85



# we have different number of samples for ascat calls on affymetrix and agilent calls
# match samples
hapmap.agilent.calls = hapmap.agilent.calls[, intersect(colnames(ascat.CN), colnames(hapmap.agilent.calls))]
hapmap.ascat.calls = ascat.CN[, intersect(colnames(ascat.CN), colnames(hapmap.agilent.calls))]
dim(hapmap.ascat.calls)
# [1] 1831313      81
dim(hapmap.agilent.calls)
# [1] 856  81
all.equal.character(colnames(hapmap.ascat.calls), colnames(hapmap.agilent.calls))
# [1] TRUE

# samples are matched

# create ascat regions object
ascat.regions = fread("HapmapLogR.txt", select = c(1:3))
ascat.regions = as.data.frame(ascat.regions)
rownames(ascat.regions) = ascat.regions$V1
dim(ascat.regions)
# [1] 1844399       3
ascat.regions = ascat.regions[intersect(rownames(ascat.CN),  rownames(ascat.regions)),]
dim(ascat.regions)
# [1] 1831313       3
ascat.regions$start = ascat.regions$Position-12
ascat.regions$end = ascat.regions$Position+12
ascat.regions = ascat.regions[,-c(1,3)]
colnames(ascat.regions) = colnames(hapmap.agilent.regions)
unique(ascat.regions$chr)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "21" "22" "3"  "4"  "5"  "6"  "7"  "8"  "9"  "X"  "Y" 
unique(hapmap.agilent.regions$chr)
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22

# remove calls from chr X and Y
ascat.regions = ascat.regions[-which(ascat.regions$chr%in%'X' | ascat.regions$chr%in%'Y'), ]
dim(ascat.regions)
# [1] 1748029       3
ascat.regions$chr = as.numeric(ascat.regions$chr)
ascat.regions = data.table(ascat.regions)

setkey(hapmap.agilent.regions)
# find the regions that were profiles by both Affymetrix SNP array and Agilent 
overlapping.regions.ascat.agilent = data.frame(foverlaps(ascat.regions, hapmap.agilent.regions, type="within", nomatch=0L))
dim(overlapping.regions.ascat.agilent)
# [1] 14532     5
overlapping.regions.ascat.agilent.indices = foverlaps(ascat.regions, hapmap.agilent.regions, type="within", nomatch=0L, which = T)
dim(overlapping.regions.ascat.agilent.indices)
# [1] 14532     2


hapmap.matched.ascat.calls = hapmap.agilent.calls[overlapping.regions.ascat.agilent.indices$yid, ]
rownames(hapmap.matched.ascat.calls)=NULL

ascat.calls = as.data.frame(ascat.CN[overlapping.regions.ascat.agilent.indices$xid,])
ascat.calls = data.frame(lapply(ascat.calls , factor , levels = c(-1,0,1)))






# compute confusion matrix ascat and hapmap agilent
i.CN = hapmap.matched.ascat.calls
i.CN[i.CN<2] <- -1
i.CN[i.CN>2] <- 1
i.CN[i.CN==2] = 0
i.CN = data.frame(lapply(i.CN , factor , levels = c(-1,0,1)))

ascat.calls <- ascat.calls[,intersect(colnames(ascat.calls),colnames(i.CN))]
all.equal(colnames(ascat.calls),colnames(i.CN))
# [1] TRUE
dim(ascat.calls)
# [1] 14532    81

confusion.matrix.ascat <- cbind(as.data.frame(table(i.CN[,1],ascat.calls[,1]))[,1:2], confusion.matrix.ascat)
colnames(confusion.matrix.ascat)[1:2] <- c("initial", "ascat")





# true states
true.i.ascat <- as.numeric(levels(confusion.matrix.ascat[,1]))[confusion.matrix.ascat[,1]]==as.numeric(levels(confusion.matrix.ascat[,2]))[confusion.matrix.ascat[,2]]

# accuracy method
accuracy.ascat<- apply(confusion.matrix.ascat[,-c(1,2)], 2, function(x)sum(x[true.i.ascat])/sum(x))



# precision for ascat

precision.ascat.loss <- as.numeric( confusion.matrix.ascat[1,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[1:3])))
precision.ascat.loss[is.nan(precision.ascat.loss)] <- 1

precision.ascat.normal <- as.numeric(confusion.matrix.ascat[5,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[4:6])))
precision.ascat.normal[is.nan(precision.ascat.normal)] <- 1

precision.ascat.gain <- as.numeric(confusion.matrix.ascat[9,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[7:9])))
precision.ascat.gain[is.nan(precision.ascat.gain)] = 1

#recall
recall.ascat.loss <- as.numeric(confusion.matrix.ascat[1,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[c(1,4,7)])))
which(is.nan(recall.ascat.loss))
# integer(0)
recall.ascat.normal <- as.numeric(confusion.matrix.ascat[5,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[c(2,5,8)])))
which(is.nan(recall.ascat.normal))
# integer(0)
recall.ascat.gain <- as.numeric(confusion.matrix.ascat[9,-c(1,2)] /apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[c(3,6,9)])))
which(is.nan(recall.ascat.gain))
# integer(0)

# compute F1 score for ascat
F.ascat.gain <- as.numeric((2 * precision.ascat.gain * recall.ascat.gain) / (precision.ascat.gain + recall.ascat.gain))
precision.ascat.gain[is.nan(F.ascat.gain)] + recall.ascat.gain[is.nan(F.ascat.gain)]
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
F.ascat.gain[is.nan(F.ascat.gain)] = 0

F.ascat.loss <- as.numeric((2 * precision.ascat.loss * recall.ascat.loss) / (precision.ascat.loss + recall.ascat.loss))
precision.ascat.loss[is.nan(F.ascat.loss)] + recall.ascat.loss[is.nan(F.ascat.loss)]
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
F.ascat.loss[is.nan(F.ascat.loss)] = 0

F.ascat.normal <- as.numeric((2 * precision.ascat.normal * recall.ascat.normal) / (precision.ascat.normal + recall.ascat.normal))
precision.ascat.normal[is.nan(F.ascat.normal)] + recall.ascat.normal[is.nan(F.ascat.normal)]

F.score.ascat = data.frame(loss = F.ascat.loss, normal = F.ascat.normal, gain = F.ascat.gain, Method = "ASCAT") 


save(F.score.ascat, file="F.score.ascat.Rdata")




