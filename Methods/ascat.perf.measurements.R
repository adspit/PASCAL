# compute performance measurements for ASCAT
library(gtools)
library(ASCAT)
library(reshape2)

  
dim(initial.regions)
# [1] 10000   400


ascat.CN <- ascat.output$nA+ascat.output$nB

ascat.CN[ascat.CN<2] <- -1
ascat.CN[ascat.CN==2] <- 0
ascat.CN[ascat.CN>2] <- 1
ascat.CN = as.data.frame(ascat.CN)
ascat.CN = data.frame(lapply( ascat.CN , factor , levels = c(-1,0,1) ) )

#select only the samples for which ascat returned calls
i.CN = initial.regions
i.CN = i.CN[,intersect(colnames(initial.regions),colnames(ascat.CN))]

# For class x: False positive: sum of column x (without main diagonal), sum(cm(:, x))-cm(x, x).

confusion.matrix.ascat <- data.frame(sapply(1:ncol(i.CN), function(i) table(i.CN[,i],ascat.CN[,i])))
colnames(confusion.matrix.ascat) <- colnames(i.CN)

# rows represent the initial data
# cols represent the ascat calls 


confusion.matrix.ascat <- cbind(as.data.frame(table(i.CN[,1],ascat.CN[,1]))[,1:2], confusion.matrix.ascat)

colnames(confusion.matrix.ascat)[1:2] <- c("initial", "ASCAT")



#
true.i <- as.numeric(levels(confusion.matrix.ascat[,1]))[confusion.matrix.ascat[,1]]==as.numeric(levels(confusion.matrix.ascat[,2]))[confusion.matrix.ascat[,2]]




# accuracy method
accuracy.ascat<- apply(confusion.matrix.ascat[,-c(1,2)], 2, function(x)sum(x[true.i])/sum(x))



precision.ascat.loss <- as.numeric(confusion.matrix.ascat[1,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[1:3])))
#precision.ascat.loss[is.nan(precision.ascat.loss)] 


precision.ascat.normal <- as.numeric(confusion.matrix.ascat[5,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[4:6])))

precision.ascat.gain <- as.numeric(confusion.matrix.ascat[9,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[7:9])))
# precision.ascat.gain[is.nan(as.numeric(precision.ascat.gain))] <- 0


#recall - specificity
recall.ascat.loss <- as.numeric(confusion.matrix.ascat[1,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[c(1,4,7)])))
# recall.ascat.loss[is.nan(as.numeric(recall.ascat.loss))] <- 0



recall.ascat.normal <- as.numeric(confusion.matrix.ascat[5,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[c(2,5,8)])))

recall.ascat.gain <- as.numeric(confusion.matrix.ascat[9,-c(1,2)] /apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[c(3,6,9)])))



F1.ascat.gain = (2 * precision.ascat.gain * recall.ascat.gain) / (precision.ascat.gain + recall.ascat.gain)
F1.ascat.normal = (2 * precision.ascat.normal * recall.ascat.normal) / (precision.ascat.normal + recall.ascat.normal)
F1.ascat.loss = (2 * precision.ascat.loss * recall.ascat.loss) / (precision.ascat.loss + recall.ascat.loss)


length(grep("t1", colnames(ascat.CN)))
#98
length(grep("t0.7", colnames(ascat.CN)))
#99
length(grep("t0.5", colnames(ascat.CN)))
#100
length(grep("t0.3", colnames(ascat.CN)))
#100

# complex
F.score.ascat = data.frame(loss = F1.ascat.loss, normal = F1.ascat.normal, gain = F1.ascat.gain, Method = "ASCAT", tumour.purity = c(rep(1,98),rep(.7,99),rep(.5,100),rep(.3,100))) 

save(F.score.ascat,file = "F.score.ascat.Rdata")






