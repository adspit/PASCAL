# performance for GenoCNA applied on HapMap
library(data.table)

# load 'gold standard'
load('hapmap.agilent.calls.Rdata')
dim(hapmap.agilent.calls)
# [1] 856  85

# load new genoCNA calls
load('genocna_calls_hapmap.Rdata')
dim(genocna.calls.hapmap)
# [1] 1844399      91

#load pfb info
load('pfb.hapmap.Rdata')
head(pfb_file)
# Name Chr Position PFB
# CN_473963 CN_473963   1    51598   2
# CN_473964 CN_473964   1    51671   2
# CN_473965 CN_473965   1    51686   2
# CN_477984 CN_477984   1    52015   2
# CN_473981 CN_473981   1    52783   2
# CN_473982 CN_473982   1    52800   2


# load sample hapmap correspondence file
load('hapmap_sample_info.Rdata')
head(sample.info)
#         ID                                               Sample.id
# 435 NA18943 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A02_31426.CEL
# 439 NA18947 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A03_31442.CEL
# 436 NA18944 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A04_31458.CEL
# 292 NA18521 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A05_31474.CEL
# 399 NA18858 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A06_31490.CEL
# 567 NA19141 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A07_31506.CEL

# check if we have a perfect match for the HapMap sample names
all.equal.character(colnames(hapmap.agilent.calls)[5:ncol(hapmap.agilent.calls)], sample.info$ID[match(colnames(hapmap.agilent.calls)[5:ncol(hapmap.agilent.calls)], sample.info$ID)])
# [1] TRUE


# update HapMap Agilent sample names
colnames(hapmap.agilent.calls)[5:ncol(hapmap.agilent.calls)] = sample.info$Sample.id[match(colnames(hapmap.agilent.calls)[5:ncol(hapmap.agilent.calls)], sample.info$ID)]

#update genocna calls sample names
colnames(genocna.calls.hapmap) = gsub('.Log.R.Ratio', '', colnames(genocna.calls.hapmap), fixed = T)


# select only samples for which we have both agilent calls and genocalls
genocna.calls.hapmap = genocna.calls.hapmap[, intersect(colnames(hapmap.agilent.calls), colnames(genocna.calls.hapmap))]
dim(genocna.calls.hapmap)
# [1] 1844399      81
hapmap.agilent.calls = cbind(hapmap.agilent.calls[,1:4], hapmap.agilent.calls[,intersect(colnames(hapmap.agilent.calls), colnames(genocna.calls.hapmap) )])


# match hapmap calls with genocna calls
hapmap.agilent.regions = data.table(hapmap.agilent.calls[,2:4])
genocna.calls.hapmap = genocna.calls.hapmap[intersect(rownames(pfb_file), rownames(genocna.calls.hapmap)), ]
all.equal.character(rownames(genocna.calls.hapmap), rownames(pfb_file))
# [1] TRUE
genocna.regions = data.table(data.frame(chr = pfb_file$Chr, start = pfb_file$Position-12, end = pfb_file$Position+12))
unique(hapmap.agilent.regions$chr)
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
# remove sex chromosomes since we don't have 'gold standard' data on them
genocna.regions = genocna.regions[genocna.regions$chr<23,]

setkey(hapmap.agilent.regions)
# find the regions that were profiles by both Affymetrix SNP array and Agilent 
overlapping.regions.agilent.genocna = data.frame(foverlaps(genocna.regions, hapmap.agilent.regions, type="within", nomatch=0L))
dim(overlapping.regions.agilent.genocna)
# [1] 14640     5
overlapping.regions.agilent.genocna.indices = foverlaps(genocna.regions, hapmap.agilent.regions, type="within", nomatch=0L, which = T)

hapmap.agilent.calls.matched.genocna = hapmap.agilent.calls[overlapping.regions.agilent.genocna.indices$yid, ]
rownames(hapmap.agilent.calls.matched.genocna)=NULL
dim(hapmap.agilent.calls.matched.genocna)
# [1] 14640    81



# compute confusion matrix ascat and hapmap agilent
i.CN = hapmap.agilent.calls.matched.genocna[,5:ncol(hapmap.agilent.calls.matched.genocna)]
i.CN[i.CN<2] <- -1
i.CN[i.CN>2] <- 1
i.CN[i.CN==2] = 0
dim(i.CN)
# [1] 14640    81

# test check
all.equal.character(colnames(i.CN), colnames(genocna.calls.hapmap))
# [1] TRUE

genocna.calls = genocna.calls.hapmap[overlapping.regions.agilent.genocna.indices$xid,]
dim(genocna.calls)
# [1] 14640    81


confusion.matrix.genocna <- as.data.frame(sapply(1:ncol(i.CN), function(i) table(i.CN[,i],genocna.calls[,i])))
colnames(confusion.matrix.genocna) <- colnames(genocna.calls)
# rows represent the initial data
# cols represent the genocna calls

confusion.matrix.genocna <- cbind(as.data.frame(table(i.CN[,1],genocna.calls[,1]))[,1:2], confusion.matrix.genocna)
colnames(confusion.matrix.genocna)[1:2] <- c("initial", "genoCNA")


# true states
true.i.genocna <- as.numeric(levels(confusion.matrix.genocna[,1]))[confusion.matrix.genocna[,1]]==as.numeric(levels(confusion.matrix.genocna[,2]))[confusion.matrix.genocna[,2]]

# accuracy genoCNA
accuracy.genocna<- apply(confusion.matrix.genocna[,-c(1,2)], 2, function(x)sum(x[true.i.genocna])/sum(x))

# precision for genocna

precision.genocna.loss <- as.numeric( confusion.matrix.genocna[1,-c(1,2)] / apply(confusion.matrix.genocna[,-c(1,2)],2,function(x) sum(x[1:3])))
which(is.nan(precision.genocna.loss))
# integer(0)

precision.genocna.normal <- as.numeric(confusion.matrix.genocna[5,-c(1,2)] / apply(confusion.matrix.genocna[,-c(1,2)],2,function(x) sum(x[4:6])))
which(is.nan(precision.genocna.normal))
# integer(0)

precision.genocna.gain <- as.numeric(confusion.matrix.genocna[9,-c(1,2)] / apply(confusion.matrix.genocna[,-c(1,2)],2,function(x) sum(x[7:9])))
which(is.nan(precision.genocna.gain))
# integer(0)

#recall
recall.genocna.loss <- as.numeric(confusion.matrix.genocna[1,-c(1,2)] / apply(confusion.matrix.genocna[,-c(1,2)],2,function(x) sum(x[c(1,4,7)])))
which(is.nan(as.numeric(recall.genocna.loss)))
# integer(0)
recall.genocna.normal <- as.numeric(confusion.matrix.genocna[5,-c(1,2)] / apply(confusion.matrix.genocna[,-c(1,2)],2,function(x) sum(x[c(2,5,8)])))
which(is.nan(as.numeric(recall.genocna.normal)))
# integer(0)
recall.genocna.gain <- as.numeric(confusion.matrix.genocna[9,-c(1,2)] /apply(confusion.matrix.genocna[,-c(1,2)],2,function(x) sum(x[c(3,6,9)])))
which(is.nan(as.numeric(recall.genocna.gain)))
# integer(0)
# compute F1 score for genocna
F.genocna.gain <- as.numeric((2 * precision.genocna.gain * recall.genocna.gain) / (precision.genocna.gain + recall.genocna.gain))

F.genocna.loss <- as.numeric((2 * precision.genocna.loss * recall.genocna.loss) / (precision.genocna.loss + recall.genocna.loss))

F.genocna.normal <- as.numeric((2 * precision.genocna.normal * recall.genocna.normal) / (precision.genocna.normal + recall.genocna.normal))


F.score.genocna = data.frame(loss = F.genocna.loss, normal = F.genocna.normal, gain = F.genocna.gain, Method = "genoCNA") 
F.score.genocna$loss[is.nan(F.score.genocna$loss)] = 0
F.score.genocna$normal[is.nan(F.score.genocna$normal)] = 0
F.score.genocna$gain[is.nan(F.score.genocna$gain)] = 0

# save(F.score.genocna, file = 'F.score.genocna.Rdata')

