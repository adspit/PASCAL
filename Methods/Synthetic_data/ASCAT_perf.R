# compute performance measurements for ASCAT
library(gtools)
library(ASCAT)
library(reshape2)

  
dim(initial.regions)
# [1] 1000000   400


ascat.CN <- ascat.output$nA+ascat.output$nB

ascat.CN[ascat.CN<2] <- -1
ascat.CN[ascat.CN==2] <- 0
ascat.CN[ascat.CN>2] <- 1
ascat.CN = as.data.frame(ascat.CN)
ascat.CN = data.frame(lapply( ascat.CN , factor , levels = c(-1,0,1) ) )
rownames(ascat.CN) = rownames(ascat.output$nA)

#select only the samples and the probes for which ascat returned calls
i.CN = initial.regions[intersect(rownames(initial.regions), rownames(ascat.CN)),intersect(colnames(initial.regions),colnames(ascat.CN))]
# check if probe names and sample names match
all.equal.character(rownames(ascat.CN), rownames(i.CN))
#[1] TRUE
all.equal.character(colnames(ascat.CN), colnames(i.CN))
#[1] TRUE



# For class x: False positive: sum of column x (without main diagonal), sum(cm(:, x))-cm(x, x).

confusion.matrix.ascat <- data.frame(sapply(1:ncol(i.CN), function(i) table(i.CN[,i],ascat.CN[,i])))
colnames(confusion.matrix.ascat) <- colnames(i.CN)

# rows represent the initial data
# cols represent the ascat calls 

confusion.matrix.ascat <- cbind(as.data.frame(table(i.CN[,1],ascat.CN[,1]))[,1:2], confusion.matrix.ascat)

colnames(confusion.matrix.ascat)[1:2] <- c("initial", "ASCAT")



#
true.i <- as.numeric(levels(confusion.matrix.ascat[,1]))[confusion.matrix.ascat[,1]]==as.numeric(levels(confusion.matrix.ascat[,2]))[confusion.matrix.ascat[,2]]




# accuracy
accuracy.ascat<- apply(confusion.matrix.ascat[,-c(1,2)], 2, function(x)sum(x[true.i])/sum(x))


# precision for each class
precision.ascat.loss <- as.numeric(confusion.matrix.ascat[1,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[1:3])))
precision.ascat.loss[is.nan(precision.ascat.loss)] = 0

precision.ascat.normal <- as.numeric(confusion.matrix.ascat[5,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[4:6])))
precision.ascat.normal[is.nan(precision.ascat.normal)] = 0
                                                                               
precision.ascat.gain <- as.numeric(confusion.matrix.ascat[9,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[7:9])))
precision.ascat.gain[is.nan(as.numeric(precision.ascat.gain))] <- 0


#recall - specificity for each class
recall.ascat.loss <- as.numeric(confusion.matrix.ascat[1,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[c(1,4,7)])))
recall.ascat.loss[is.nan(as.numeric(recall.ascat.loss))] <- 0

recall.ascat.normal <- as.numeric(confusion.matrix.ascat[5,-c(1,2)] / apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[c(2,5,8)])))
recall.ascat.normal[is.nan(as.numeric(recall.ascat.normal))] <- 0
                                                                            
recall.ascat.gain <- as.numeric(confusion.matrix.ascat[9,-c(1,2)] /apply(confusion.matrix.ascat[,-c(1,2)],2,function(x) sum(x[c(3,6,9)])))
recall.ascat.gain[is.nan(as.numeric(recall.ascat.gain))] <- 0


F.ascat.gain = (2 * precision.ascat.gain * recall.ascat.gain) / (precision.ascat.gain + recall.ascat.gain)
F.ascat.normal = (2 * precision.ascat.normal * recall.ascat.normal) / (precision.ascat.normal + recall.ascat.normal)
F.ascat.loss = (2 * precision.ascat.loss * recall.ascat.loss) / (precision.ascat.loss + recall.ascat.loss)


length(grep("t1", colnames(ascat.CN)))
#97
length(grep("t0.7", colnames(ascat.CN)))
#99
length(grep("t0.5", colnames(ascat.CN)))
#97
length(grep("t0.3", colnames(ascat.CN)))
#99

# complex
F.score.ascat = data.frame(loss = F.ascat.loss, normal = F.ascat.normal, gain = F.ascat.gain, Method = "ASCAT", tumour.purity = c(rep(1,97),rep(.7,99),rep(.5,97),rep(.3,99)))
F.score.ascat$loss[is.nan(F.score.ascat$loss)] = 0
F.score.ascat$normal[is.nan(F.score.ascat$normal)] = 0
F.score.ascat$gain[is.nan(F.score.ascat$gain)] = 0
                                                                         

save(F.score.ascat,file = "F.score.ascat.Rdata")






