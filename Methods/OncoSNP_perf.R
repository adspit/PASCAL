# performance OncoSNP

library(gtools)
library(reshape2)
# set path to the OncoSNP output files
# setwd("OncoSNP/oncosnp_rez/")

# save sorted OncoSNP output file names in a variable 
OncoSNP.files <- mixedsort(dir(pattern = ".1.ggs"))

# initialize an empty list for the OncoSNP calls
oncoSNP.calls <- list()

# loop that collapses the states to loss, normal and gain
for(i in seq_along(OncoSNP.files))
{
  print(i)
  onco.snp.output <- read.table(OncoSNP.files[i], sep="\t", h=T)
  # colnames(onco.snp.output) <- c("ProbeID",	"Chromosome",	"Position",	"LogRRatio",	"BAlleleFreq", "TumourState",	"CopyNumber",	"LOH",	"NumberofBAlleles",	"NormalContent",	"Amp", "Rank")
  onco.snp.output <- onco.snp.output[order(onco.snp.output$Position),]
  onco.snp.output <- onco.snp.output[onco.snp.output$Rank==1,]
  dim(onco.snp.output)
  onco.snp.output$CN <- onco.snp.output$CopyNumber
  onco.snp.output$CN[onco.snp.output$CopyNumber<2] <- -1
  onco.snp.output$CN[onco.snp.output$CopyNumber==2] <- 0
  onco.snp.output$CN[onco.snp.output$CopyNumber>2] <- 1
  oncoSNP.calls[[i]] <- onco.snp.output
  
}


#reorder initial regions so that the sample order matches the one from oncosnp

initial.regions = cbind(initial.regions[,301:400], initial.regions[,201:300], initial.regions[,101:200], initial.regions[,1:100])
names(oncoSNP.calls) = colnames(initial.regions)


# create initial list
reduced.initial = lapply(1:ncol(initial.regions), function(i) initial.regions[oncoSNP.calls[[i]]$ProbeID,i])


#compute confusion matrix for each sample
oncoSNP.c <- lapply(oncoSNP.calls, function(df) {df$CN <- factor(df$CN,levels=c(-1,0,1)); df})




confusion.matrix.oncosnp <- sapply(1:400, function(i) table(reduced.initial[[i]],oncoSNP.c[[i]]$CN))
dim(confusion.matrix.oncosnp)
# [1]   9 400

confusion.matrix.oncosnp <- cbind(as.data.frame(table(reduced.initial[[1]],oncoSNP.c[[1]]$CN))[,1:2], confusion.matrix.oncosnp)

colnames(confusion.matrix.oncosnp)[1:2] <- c("initial", "oncoSNP")

# true states
true.i.oncosnp <- as.numeric(levels(confusion.matrix.oncosnp[,1]))[confusion.matrix.oncosnp[,1]]==as.numeric(levels(confusion.matrix.oncosnp[,2]))[confusion.matrix.oncosnp[,2]]

# accuracy method
accuracy.oncoSNP<- apply(confusion.matrix.oncosnp[,-c(1,2)], 2, function(x)sum(x[true.i.oncosnp])/sum(x))

# precision

precision.oncoSNP.loss <- as.numeric( confusion.matrix.oncosnp[1,-c(1,2)] / apply(confusion.matrix.oncosnp[,-c(1,2)],2,function(x) sum(x[1:3])))


precision.oncoSNP.normal <- confusion.matrix.oncosnp[5,-c(1,2)] / apply(confusion.matrix.oncosnp[,-c(1,2)],2,function(x) sum(x[4:6]))
precision.oncoSNP.gain <- confusion.matrix.oncosnp[9,-c(1,2)] / apply(confusion.matrix.oncosnp[,-c(1,2)],2,function(x) sum(x[7:9]))

#recall
recall.oncoSNP.loss <- confusion.matrix.oncosnp[1,-c(1,2)] / apply(confusion.matrix.oncosnp[,-c(1,2)],2,function(x) sum(x[c(1,4,7)]))



recall.oncoSNP.normal <- confusion.matrix.oncosnp[5,-c(1,2)] / apply(confusion.matrix.oncosnp[,-c(1,2)],2,function(x) sum(x[c(2,5,8)]))
recall.oncoSNP.gain <- confusion.matrix.oncosnp[9,-c(1,2)] /apply(confusion.matrix.oncosnp[,-c(1,2)],2,function(x) sum(x[c(3,6,9)]))
#????


F.oncoSNP.gain <- as.numeric((2 * precision.oncoSNP.gain * recall.oncoSNP.gain) / (precision.oncoSNP.gain + recall.oncoSNP.gain))
F.oncoSNP.normal <- as.numeric((2 * precision.oncoSNP.normal * recall.oncoSNP.normal) / (precision.oncoSNP.normal + recall.oncoSNP.normal))
F.oncoSNP.loss <- as.numeric((2 * precision.oncoSNP.loss * recall.oncoSNP.loss) / (precision.oncoSNP.loss + recall.oncoSNP.loss))



F.score.oncosnp.n.10_6 = data.frame(loss = F.oncoSNP.loss, normal = F.oncoSNP.normal, gain = F.oncoSNP.gain, Method = "OncoSNP", tumour.purity = c(rep(0.3,100),rep(.5,100),rep(.7,100),rep(1,100))) 
F.score.oncosnp.n.10_6$loss[is.nan(F.score.oncosnp.n.10_6$loss)] = 0
F.score.oncosnp.n.10_6$normal[is.nan(F.score.oncosnp.n.10_6$normal)] = 0
F.score.oncosnp.n.10_6$gain[is.nan(F.score.oncosnp.n.10_6$gain)] = 0

save(F.score.oncosnp.n.10_6, file="../F.score.oncosnp.n.10_6.Rdata")
