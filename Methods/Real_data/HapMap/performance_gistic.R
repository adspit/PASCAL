# performance for GISTIC on HapMap
library(data.table)
library(dplyr)
library(plyr)

# load ground truth
load('hapmap.agilent.calls.Rdata')
dim(hapmap.agilent.calls)
# [1] 856  85

# read in GISTIC results
GISTIC_results = read.table('all_thresholded.by_genes.txt', h=T, sep="\t", stringsAsFactors = F)
dim(GISTIC_results)
#  1844399      94

GISTIC_results[1:3, 1:5]
#   Gene.Symbol Locus.ID Cytoband Sample.1 Sample.2
# 1       OR4F5    79501  1p36.33        0        0
# 2 OR4F16|chr1    81399  1p36.33        0        0
# 3 OR4F29|chr1   729759  1p36.33        0        0


# sample names
sample_names = read.table('hapmap_gistic_sample_names.txt', sep="\t", stringsAsFactors = F)
head(sample_names)
#                                                                  V1
# 1 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A02_31426.CEL.Log.R.Ratio
# 2 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A03_31442.CEL.Log.R.Ratio
# 3 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A04_31458.CEL.Log.R.Ratio
# 4 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A05_31474.CEL.Log.R.Ratio
# 5 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A06_31490.CEL.Log.R.Ratio
# 6 SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A07_31506.CEL.Log.R.Ratio



# update gistic sample names
colnames(GISTIC_results)[4:ncol(GISTIC_results)] = sample_names$V1

# load sample hapmap correspondence file
load('../../GenoCNA_hapmap/hapmap_sample_info.Rdata')
head(sample.info)

# check if we have a perfect match for the HapMap sample names
all.equal.character(colnames(hapmap.agilent.calls)[5:ncol(hapmap.agilent.calls)], sample.info$ID[match(colnames(hapmap.agilent.calls)[5:ncol(hapmap.agilent.calls)], sample.info$ID)])
# [1] TRUE



# update HapMap Agilent sample names
colnames(hapmap.agilent.calls)[5:ncol(hapmap.agilent.calls)] = sample.info$Sample.id[match(colnames(hapmap.agilent.calls)[5:ncol(hapmap.agilent.calls)], sample.info$ID)]

#update GISTIC calls sample names
colnames(GISTIC_results) = gsub('.Log.R.Ratio', '', colnames(GISTIC_results), fixed = T)

hapmap.agilent.calls = cbind(hapmap.agilent.calls[,1:4], hapmap.agilent.calls[,intersect(colnames(hapmap.agilent.calls), colnames(GISTIC_results) )])
hapmap.agilent.regions = data.table(hapmap.agilent.calls[,2:4])

gistic_calls = GISTIC_results[, intersect(colnames(hapmap.agilent.calls), colnames(GISTIC_results))]
gistic_regions = GISTIC_results[,1:3]


hg18_annotation = read.csv('rg.csv', h=T, sep="\t", stringsAsFactors = F) 
dim(hg18_annotation)
# [1] 26835    12
head(hg18_annotation)
# refseq                                            gene   symb locus_id   chr strand     start       end cds_start   cds_end   status chrn
# 1 NM_000014                 alpha-2-macroglobulin precursor    A2M        2 chr12      0   9111570   9159825   9111685   9159712 Reviewed   12
# 2 NM_000015                           arylamide acetylase 2   NAT2       10  chr8      1  18293034  18303003  18301793  18302666 Reviewed    8
# 3 NM_000016      acyl-Coenzyme A dehydrogenase, C-4 to C-12  ACADM       34  chr1      1  75962869  76001771  75963060  76001036 Reviewed    1
# 4 NM_000017 acyl-Coenzyme A dehydrogenase, C-2 to C-3 short  ACADS       35 chr12      1 119647953 119662194 119648071 119661634 Reviewed   12
# 5 NM_000018  acyl-Coenzyme A dehydrogenase, very long chain ACADVL       37 chr17      1   7063876   7069309   7064027   7069140 Reviewed   17
# 6 NM_000019 acetyl-Coenzyme A acetyltransferase 1 precursor  ACAT1       38 chr11      1 107497467 107523485 107497543 107523327 Reviewed   11


# check for genes with 'chr' in the symbol names
dim(gistic_regions[grep('chr', gistic_regions$Gene.Symbol), ])
# 12 3

# split the gene symbol names in two columns
gistic_regions = within(gistic_regions, Gene.Symbol <- data.frame(do.call('rbind', strsplit(as.character(Gene.Symbol), '|', fixed=TRUE))), stringsAsFactors=FALSE)
head(gistic_regions)
# Gene.Symbol.X1 Gene.Symbol.X2 Locus.ID Cytoband
# 1          OR4F5          OR4F5    79501  1p36.33
# 2         OR4F16           chr1    81399  1p36.33
# 3         OR4F29           chr1   729759  1p36.33
# 4          OR4F3           chr1    26683  1p36.33
# 5         SAMD11         SAMD11   148398  1p36.33
# 6          NOC2L          NOC2L    26155  1p36.33


# build gistic_regions df
gistic_regions = data.frame(symb = as.character(gistic_regions$Gene.Symbol$X1), chr = as.character(gistic_regions$Gene.Symbol$X2), Locus.ID = gistic_regions$Locus.ID, Cytoband = gistic_regions$Cytoband, stringsAsFactors = F)
gistic_regions$chr[-grep('chr', gistic_regions$chr)] = hg18_annotation[match(gistic_regions[-grep('chr', gistic_regions$chr),1], hg18_annotation$symb),5]
colnames(gistic_regions)[3] = 'locus_id'

dat.grp <- group_by(test, symb, chr)
head(dat.grp)
summarise(dat.grp, start = min(start), end=max(end))
dim(summarise(dat.grp, start = min(start), end=max(end))
)

resulting_annotation_2 = as.data.frame(summarise(dat.grp, start = min(start), end=max(end)))
processed_gistic_regions_2 = join(gistic_regions, resulting_annotation_2, by = c('symb', 'chr'))
# match gene order in gistic_regions and processed_gistic_regions
all.equal.character(processed_gistic_regions_2$symb, gistic_regions$symb)
# [1] TRUE

# convert reg to data table obj for performing overlap
gistic_reg_2 = data.table(processed_gistic_regions_2[, c(2,5,6)])
gistic_reg_2$chr= gsub('chr', '', gistic_reg_2$chr)
gistic_reg_2$chr= gsub('_random', '', gistic_reg_2$chr)
gistic_reg_2$chr= gsub('_cox_hap1', '', gistic_reg_2$chr)
gistic_reg_2$chr= gsub('_qbl_hap2', '', gistic_reg_2$chr)
gistic_reg_2$chr = as.numeric(gistic_reg_2$chr)


setkey(hapmap.agilent.regions)
# find the regions that were profiles by both Affymetrix SNP array and Agilent 
overlapping.regions.agilent.gistic = data.frame(foverlaps(gistic_reg_2, hapmap.agilent.regions, nomatch=0L))
dim(overlapping.regions.agilent.gistic)
# [1] 381     5
overlapping.regions.agilent.gistic.indices = foverlaps(gistic_reg_2, hapmap.agilent.regions, nomatch=0L, which = T)
dim(overlapping.regions.agilent.gistic.indices)
# [1] 381   2

hapmap.agilent.calls.matched.gistic = hapmap.agilent.calls[overlapping.regions.agilent.gistic.indices$yid, ]
rownames(hapmap.agilent.calls.matched.gistic)=NULL
dim(hapmap.agilent.calls.matched.gistic)
# [1]  381 85

# compute confusion matrix ascat and hapmap agilent
i.CN = hapmap.agilent.calls.matched.gistic[,5:ncol(hapmap.agilent.calls.matched.gistic)]
i.CN[i.CN<2] <- -1
i.CN[i.CN>2] <- 1
i.CN[i.CN==2] = 0
i.CN <- data.frame(lapply(i.CN, factor , levels = c(-1,0,1)))

dim(i.CN)
# [1] 381 81

# sample match check
all.equal.character(colnames(i.CN), colnames(gistic_calls))
# [1] TRUE

gistic_overlapping_calls = gistic_calls[overlapping.regions.agilent.gistic.indices$xid, ]
range(gistic_overlapping_calls)
gistic_overlapping_calls[gistic_overlapping_calls<2] = -1
gistic_overlapping_calls[gistic_overlapping_calls>2] = 1
gistic_overlapping_calls[gistic_overlapping_calls==2] = 0
gistic_overlapping_calls <- data.frame(lapply(gistic_overlapping_calls, factor , levels = c(-1,0,1)))


dim(gistic_overlapping_calls)
# [1] 381  81

confusion.matrix.gistic <- as.data.frame(sapply(1:ncol(i.CN), function(i) table(i.CN[,i],gistic_overlapping_calls[,i])))
colnames(confusion.matrix.gistic) <- colnames(gistic_overlapping_calls)

# rows represent the initial data
# cols represent the gistic calls

confusion.matrix.gistic <- cbind(as.data.frame(table(i.CN[,1],gistic_overlapping_calls[,1]))[,1:2], confusion.matrix.gistic)
colnames(confusion.matrix.gistic)[1:2] <- c("initial", "gistic")

true.i <- as.numeric(levels(confusion.matrix.gistic[,1]))[confusion.matrix.gistic[,1]]==as.numeric(levels(confusion.matrix.gistic[,2]))[confusion.matrix.gistic[,2]]
range(accuracy.gistic)
# [1] 0.1250000 0.2228117


# precision
precision.gistic.loss <- as.numeric(confusion.matrix.gistic[1,-c(1,2)] / apply(confusion.matrix.gistic[,-c(1,2)],2,function(x) sum(x[1:3])))
which(is.nan(precision.gistic.loss))
# integer(0)

precision.gistic.normal <- as.numeric(confusion.matrix.gistic[5,-c(1,2)] / apply(confusion.matrix.gistic[,-c(1,2)],2,function(x) sum(x[4:6])))
which(is.nan(precision.gistic.normal))
# integer(0)

precision.gistic.gain <- as.numeric(confusion.matrix.gistic[9,-c(1,2)] / apply(confusion.matrix.gistic[,-c(1,2)],2,function(x) sum(x[7:9])))
which(is.nan(precision.gistic.gain)) 
# [1] NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN

precision.gistic.gain[which(is.nan(precision.gistic.gain))] = 1

#recall
recall.gistic.loss <- as.numeric(confusion.matrix.gistic[1,-c(1,2)] / apply(confusion.matrix.gistic[,-c(1,2)],2,function(x) sum(x[c(1,4,7)])))
which(is.nan(recall.gistic.loss))
# integer(0)

recall.gistic.normal <- as.numeric(confusion.matrix.gistic[5,-c(1,2)] / apply(confusion.matrix.gistic[,-c(1,2)],2,function(x) sum(x[c(2,5,8)])))
which(is.nan(recall.gistic.normal))
# integer(0)

recall.gistic.gain <- as.numeric(confusion.matrix.gistic[9,-c(1,2)] /apply(confusion.matrix.gistic[,-c(1,2)],2,function(x) sum(x[c(3,6,9)])))
which(is.nan(recall.gistic.gain))
# integer(0)

# compute F1 score for ascat
F.gistic.gain <- as.numeric((2 * precision.gistic.gain * recall.gistic.gain) / (precision.gistic.gain + recall.gistic.gain))

F.gistic.loss <- as.numeric((2 * precision.gistic.loss * recall.gistic.loss) / (precision.gistic.loss + recall.gistic.loss))

F.gistic.normal <- as.numeric((2 * precision.gistic.normal * recall.gistic.normal) / (precision.gistic.normal + recall.gistic.normal))


F.score.gistic = data.frame(loss = F.gistic.loss, normal = F.gistic.normal, gain = F.gistic.gain, Method = "GISTIC") 
F.score.gistic$loss[is.nan(F.score.gistic$loss)] = 0
F.score.gistic$normal[is.nan(F.score.gistic$normal)] = 0
F.score.gistic$gain[is.nan(F.score.gistic$gain)] = 0

save(F.score.gistic, file = 'F.score.gistic.hapmap.Rdata')

