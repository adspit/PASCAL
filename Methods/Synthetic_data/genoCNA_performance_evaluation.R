# calculate genocna performance
library(gtools)

setwd('/storage_scratch/users/adriana.pitea/GenoCNA/GenoCNA_results/')


# calculate f-scores for t1 genocnca results
# select SNP level calls files
result_files_snp_level_genocna_t1 = dir(pattern='sample_t1_')
result_files_snp_level_genocna_t1 = mixedsort(result_files_snp_level_genocna_t1[grep('_SNP.txt', result_files_snp_level_genocna_t1)])

test_t1 =  lapply(result_files_snp_level_genocna_t1, function(x) { aux = read.table(x, sep="\t", h=T, stringsAsFactors=F, row.names = 1)})

length(test_t1)
# [1] 100
names(test_t1) = result_files_snp_level_genocna_t1
cn_calls_t1 = lapply(1:100, function(x) {colnames(test_t1[[x]])[3] = names(test_t1)[x]; test_t1[[x]][3]})
aux_t1 = Reduce(function(x, y) merge(x, y, all = TRUE), lapply(cn_calls_t1, function(y) data.table(y, keep.rownames=TRUE, key = "rn")))
predicted_states_t1 = aux_t1[,-1]
dim(predicted_states_t1)
# [1] 1844399     100
predicted_states_t1 = as.matrix(predicted_states_t1)
rownames(predicted_states_t1) = aux_t1$rn
predicted_states_t1[1:5, 1:5]
# 
#                        sample_t1_1_SNP.txt sample_t1_2_SNP.txt sample_t1_3_SNP.txt
# AFFX-SNP_10000979                  NA                  NA                  NA
# AFFX-SNP_10009702                  NA                   1                   1
# AFFX-SNP_10021569                   3                  NA                   4
# AFFX-SNP_10026879                  NA                  NA                   1
# AFFX-SNP_10034687                  NA                  NA                   1

predicted_states_t1[predicted_states_t1<2] = -1
predicted_states_t1[predicted_states_t1==2] = 0
predicted_states_t1[predicted_states_t1>2] =  1
predicted_states_t1 = as.data.frame(predicted_states_t1)
predicted_states_t1[1:5, 1:5]
#                  sample_t1_1_SNP.txt sample_t1_2_SNP.txt sample_t1_3_SNP.txt
# AFFX-SNP_10000979                  NA                  NA                  NA
# AFFX-SNP_10009702                  NA                  -1                  -1
# AFFX-SNP_10021569                   1                  NA                   1
# AFFX-SNP_10026879                  NA                  NA                  -1
# AFFX-SNP_10034687                  NA                  NA                  -1                                                             
genocna.calls.t1 <- data.frame(lapply( predicted_states_t1 , factor , levels = c(-1,0,1) ) )
genocna.calls.t1[1:5, 1:3]
#  sample_t1_1_SNP.txt sample_t1_2_SNP.txt sample_t1_3_SNP.txt
# 1                <NA>                <NA>                <NA>
# 2                <NA>                  -1                  -1
# 3                   1                <NA>                   1
# 4                <NA>                <NA>                  -1
# 5                <NA>                <NA>                  -1
rownames(genocna.calls.t1) = rownames(predicted_states_t1)
# save(genocna.calls.t1, file = 'genocna_calls_t1.Rdata')


# same for 0.7
# calculate f-scores for t1 genocnca results
# select SNP level calls files

setwd('/home/icb/adriana.pitea/GenoCNA_synthetic/')
result_files_snp_level_genocna_t07 = dir(pattern='sample_t07_')
result_files_snp_level_genocna_t07 = mixedsort(result_files_snp_level_genocna_t07[grep('_SNP.txt', result_files_snp_level_genocna_t07)])

length(result_files_snp_level_genocna_t07)
# [1] 100



test_t07 =  lapply(result_files_snp_level_genocna_t07, function(x) { aux = read.table(x, sep="\t", h=T, stringsAsFactors=F, row.names = 1)})


names(test_t07) = result_files_snp_level_genocna_t07
cn_calls_t07 = lapply(1:100, function(x) {colnames(test_t07[[x]])[3] = names(test_t07)[x]; test_t07[[x]][3]})
aux_t07 = Reduce(function(x, y) merge(x, y, all = TRUE), lapply(cn_calls_t07, function(y) data.table(y, keep.rownames=TRUE, key = "rn")))
predicted_states_t07 = aux_t07[,-1]
dim(predicted_states_t07)
# [1] 1844399     100
predicted_states_t07 = as.matrix(predicted_states_t07)
rownames(predicted_states_t07) = aux_t07$rn
predicted_states_t07[1:5, 1:2]
#                         sample_t07_101_SNP.txt sample_t07_102_SNP.txt
# AFFX-SNP_10000979                     NA                     NA
# AFFX-SNP_10009702                     NA                      2
# AFFX-SNP_10021569                      3                     NA
# AFFX-SNP_10026879                     NA                     NA
# AFFX-SNP_10034687                     NA                     NA

predicted_states_t07[predicted_states_t07<2] = -1
predicted_states_t07[predicted_states_t07==2] = 0
predicted_states_t07[predicted_states_t07>2] =  1
predicted_states_t07 = as.data.frame(predicted_states_t07)
genocna.calls.t07 <- data.frame(lapply( predicted_states_t07 , factor , levels = c(-1,0,1) ) )
rownames(genocna.calls.t07) = rownames(predicted_states_t07)

save(genocna.calls.t07, file = 'genocna_calls_t07.Rdata')
save.image('genocna.t07.Rdata')



# same for 0.5
# calculate f-scores for t1 genocnca results
# select SNP level calls files
setwd('/storage_scratch/users/adriana.pitea/GenoCNA/')
result_files_snp_level_genocna_t05 = dir(pattern='sample_t05_')
result_files_snp_level_genocna_t05 = mixedsort(result_files_snp_level_genocna_t05[grep('_SNP.txt', result_files_snp_level_genocna_t05)])
test_t05 =  lapply(result_files_snp_level_genocna_t05, function(x) { aux = read.table(x, sep="\t", h=T, stringsAsFactors=F, row.names = 1)})
length(test_t1)
# [1] 91
names(test_t05) = result_files_snp_level_genocna_t05
cn_calls_t05 = lapply(1:100, function(x) {colnames(test_t05[[x]])[3] = names(test_t05)[x]; test_t05[[x]][3]})
aux_t05 = Reduce(function(x, y) merge(x, y, all = TRUE), lapply(cn_calls_t05, function(y) data.table(y, keep.rownames=TRUE, key = "rn")))
predicted_states_t05 = aux_t05[,-1]
dim(predicted_states_t05)
# [1] 1844399     100
predicted_states_t05 = as.matrix(predicted_states_t05)
rownames(predicted_states_t05) = aux_t05$rn
predicted_states_t05[1:5, 1:2]
#                       sample_t05_201_SNP.txt sample_t05_202_SNP.txt
# AFFX-SNP_10000979                     NA                     NA
# AFFX-SNP_10009702                     NA                      2
# AFFX-SNP_10021569                      2                     NA
# AFFX-SNP_10026879                     NA                     NA
# AFFX-SNP_10034687                     NA                     NA

predicted_states_t05[predicted_states_t05<2] = -1
predicted_states_t05[predicted_states_t05==2] = 0
predicted_states_t05[predicted_states_t05>2] =  1

predicted_states_t05 = as.data.frame(predicted_states_t05)
genocna.calls.t05 <- data.frame(lapply( predicted_states_t05 , factor , levels = c(-1,0,1) ) )
rownames(genocna.calls.t05) = rownames(predicted_states_t05)
save(genocna.calls.t05, file = 'genocna_calls_t05.Rdata')
save.image('genocna.t05.Rdata')


# same for 0.3
# calculate f-scores for t1 genocnca results
# select SNP level calls files
setwd('/home/icb/adriana.pitea/GenoCNA_synthetic/')

result_files_snp_level_genocna_t03 = dir(pattern='sample_t03_')
result_files_snp_level_genocna_t03 = mixedsort(result_files_snp_level_genocna_t03[grep('_SNP.txt', result_files_snp_level_genocna_t03)])

test_t03 =  lapply(result_files_snp_level_genocna_t03, function(x) { aux = read.table(x, sep="\t", h=T, stringsAsFactors=F, row.names = 1)})


names(test_t03) = result_files_snp_level_genocna_t03
cn_calls_t03 = lapply(1:100, function(x) {colnames(test_t03[[x]])[3] = names(test_t03)[x]; test_t03[[x]][3]})
aux_t03 = Reduce(function(x, y) merge(x, y, all = TRUE), lapply(cn_calls_t03, function(y) data.table(y, keep.rownames=TRUE, key = "rn")))
predicted_states_t03 = aux_t03[,-1]
dim(predicted_states_t03)
# [1] 1844399     100
predicted_states_t03 = as.matrix(predicted_states_t03)
rownames(predicted_states_t03) = aux_t03$rn
predicted_states_t03[1:5, 1:3]
#                         sample_t03_301_SNP.txt sample_t03_302_SNP.txt
# AFFX-SNP_10000979                     NA                     NA
# AFFX-SNP_10009702                     NA                      2
# AFFX-SNP_10021569                      2                     NA
# AFFX-SNP_10026879                     NA                     NA
# AFFX-SNP_10034687                     NA                     NA

predicted_states_t03[predicted_states_t03<2] = -1
predicted_states_t03[predicted_states_t03==2] = 0
predicted_states_t03[predicted_states_t03>2] =  1
predicted_states_t03 = as.data.frame(predicted_states_t03)
genocna.calls.t03 <- data.frame(lapply( predicted_states_t03 , factor , levels = c(-1,0,1) ) )
rownames(genocna.calls.t03) = rownames(predicted_states_t03)
save(genocna.calls.t03, file = 'genocna_calls_t03.Rdata')
save.image('genocna.t03.Rdata')



genocalls = cbind(genocna.calls.t1, genocna.calls.t07, genocna.calls.t05, genocna.calls.t03)

# calculate confusion matrix
load('/home/icb/adriana.pitea/for_adriana_old/initial.regions.Rdata')

# match rownames
genocalls = genocalls[intersect(rownames(initial.regions), rownames(genocalls)), ]
all.equal.character(rownames(genocalls), rownames(initial.regions))
# [1] TRUE


confusion.matrix.genocna <- as.data.frame(sapply(1:ncol(initial.regions), function(i) table(initial.regions[,i],genocalls[,i])))
dim(confusion.matrix.genocna)
# [1]   9 400
colnames(confusion.matrix.genocna) <- colnames(initial.regions)

confusion.matrix.genocna <- cbind(as.data.frame(table(initial.regions[,1],genocalls[,1]))[,1:2], confusion.matrix.genocna)
colnames(confusion.matrix.genocna)[1:2] <- c("true", "genocna")

# true/false positions
true.i <- as.numeric(levels(confusion.matrix.genocna[,1]))[confusion.matrix.genocna[,1]]==as.numeric(levels(confusion.matrix.genocna[,2]))[confusion.matrix.genocna[,2]]

# calculate accuracy
accuracy.genocna<- apply(confusion.matrix.genocna[,-c(1,2)], 2, function(x)sum(x[true.i])/nrow(initial.regions))
#accuracy for the probes where the method didn't return NA
accuracy.genocna.nonNA<- apply(confusion.matrix.genocna[,-c(1,2)], 2, function(x)sum(x[true.i])/sum(x))


# precision for genocna

precision.genocna.loss <- as.numeric( confusion.matrix.genocna[1,-c(1,2)] / apply(confusion.matrix.genocna[,-c(1,2)],2,function(x) sum(x[1:3])))
which(is.nan(precision.genocna.loss))
# [1]  84 101 104 105 108 109 113 115 116 120 123 126 128 147 151 152 156 163 164 165 172 179 181 184 186 189 196 200
precision.genocna.loss[is.nan(precision.genocna.loss)] = 1



precision.genocna.normal <- as.numeric(confusion.matrix.genocna[5,-c(1,2)] / apply(confusion.matrix.genocna[,-c(1,2)],2,function(x) sum(x[4:6])))
which(is.nan(precision.genocna.normal))
# [1]  14  20  36  51  52  59  87  96 100
precision.genocna.normal[is.nan(precision.genocna.normal)] = 1



precision.genocna.gain <- as.numeric(confusion.matrix.genocna[9,-c(1,2)] / apply(confusion.matrix.genocna[,-c(1,2)],2,function(x) sum(x[7:9])))
which(is.nan(precision.genocna.gain))
# [1]  37 107 122 137 172 179 188 222 322 340 350 359 372 379
# confusion.matrix.genocna[7:9,which(is.nan(precision.genocna.gain))+2]
# t1sample37.txt t0.7sample7.txt t0.7sample22.txt t0.7sample37.txt
# 7              0               0                0                0
# 8              0               0                0                0
# 9              0               0                0                0
# t0.7sample72.txt t0.7sample79.txt t0.7sample88.txt t0.5sample22.txt
# 7                0                0                0                0
# 8                0                0                0                0
# 9                0                0                0                0
# t0.3sample22.txt t0.3sample40.txt t0.3sample50.txt t0.3sample59.txt
# 7                0                0                0                0
# 8                0                0                0                0
# 9                0                0                0                0
# t0.3sample72.txt t0.3sample79.txt
# 7                0                0
# 8                0                0
# 9                0                0
precision.genocna.gain[is.nan(precision.genocna.gain)] = 1



#recall
recall.genocna.loss <- as.numeric(confusion.matrix.genocna[1,-c(1,2)] / apply(confusion.matrix.genocna[,-c(1,2)],2,function(x) sum(x[c(1,4,7)])))
which(is.nan(as.numeric(recall.genocna.loss))) # a lot
recall.genocna.loss[is.nan(recall.genocna.loss)] = 1

recall.genocna.normal <- as.numeric(confusion.matrix.genocna[5,-c(1,2)] / apply(confusion.matrix.genocna[,-c(1,2)],2,function(x) sum(x[c(2,5,8)])))
which(is.nan(as.numeric(recall.genocna.normal)))
recall.genocna.normal[is.nan(recall.genocna.normal)] = 1

recall.genocna.gain <- as.numeric(confusion.matrix.genocna[9,-c(1,2)] /apply(confusion.matrix.genocna[,-c(1,2)],2,function(x) sum(x[c(3,6,9)])))
which(is.nan(as.numeric(recall.genocna.gain)))
# confusion.matrix.genocna[c(3,6,9),which(is.nan(recall.genocna.gain))+2]
# t1sample37.txt t1sample88.txt t0.7sample37.txt t0.7sample88.txt
# 3              0              0                0                0
# 6              0              0                0                0
# 9              0              0                0                0
# t0.5sample37.txt t0.5sample88.txt t0.3sample37.txt t0.3sample88.txt
# 3                0                0                0                0
# 6                0                0                0                0
# 9                0                0                0                0

recall.genocna.gain[is.nan(recall.genocna.gain)] = 1


# compute F1 score for genocna
F.genocna.gain <- as.numeric((2 * precision.genocna.gain * recall.genocna.gain) / (precision.genocna.gain + recall.genocna.gain))

F.genocna.loss <- as.numeric((2 * precision.genocna.loss * recall.genocna.loss) / (precision.genocna.loss + recall.genocna.loss))

F.genocna.normal <- as.numeric((2 * precision.genocna.normal * recall.genocna.normal) / (precision.genocna.normal + recall.genocna.normal))


F.score.genocna = data.frame(loss = F.genocna.loss, normal = F.genocna.normal, gain = F.genocna.gain, Method = "genocna", tumour.purity = c(rep(1,100),rep(.7,100),rep(.5,100),rep(.3,100)))
F.score.genocna$loss[is.nan(F.score.genocna$loss)] = 0
F.score.genocna$normal[is.nan(F.score.genocna$normal)] = 0
F.score.genocna$gain[is.nan(F.score.genocna$gain)] = 0

save.image('performance.analysis.genocna.synthetic.Rdata')
