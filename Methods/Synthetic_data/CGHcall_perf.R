# compute CGHcall performance

# set your WD

# load CGHcall results - the Rdata object obtained after running the expandCGHcall function
load("initial.regions.Rdata")

#dim(initial.regions)
#[1] 1844399     400


# load initial data saved in initial.regions.Rdata object when running generate.synthetic.data.R
load('initial.regions.Rdata')
dim(initial.regions)
# [1] 1000000   400

load("CGH.calls.Rdata")
# colapse CGHcall output to 3 states
CGH.calls <- as.data.frame(calls(result0))
range(CGH.calls)
# [1] -2  2
CGH.calls[CGH.calls<0] <- -1
CGH.calls[CGH.calls>0] <- 1
CGH.calls <- data.frame( lapply( CGH.calls , factor , levels = c(-1,0,1) ) )
dim(CGH.calls)
#[1] 1835630     400

# reduce the initial number of probes to the number of probes resulting after running the CGHcall procedure
red.init.reg = initial.regions[intersect(rownames(initial.regions), rownames(CGH.calls)),]
dim(red.init.reg)
#[1] 1835630     400

# check if the probes match
all.equal.character(rownames(red.init.reg), rownames(CGH.calls))
#[1] TRUE


# calculate confusion matrix
confusion.matrix.cghcall <- as.data.frame(sapply(1:ncol(red.init.reg), function(i) table(red.init.reg[,i],CGH.calls[,i])))
colnames(confusion.matrix.cghcall) <- colnames(CGH.calls)
confusion.matrix.cghcall <- cbind(as.data.frame(table(red.init.reg[,1],CGH.calls[,1]))[,1:2], confusion.matrix.cghcall)
colnames(confusion.matrix.cghcall)[1:2] <- c("true", "CGHcall")

# true/false positions
true.i <- as.numeric(levels(confusion.matrix.cghcall[,1]))[confusion.matrix.cghcall[,1]]==as.numeric(levels(confusion.matrix.cghcall[,2]))[confusion.matrix.cghcall[,2]]




# accuracy method
accuracy.CGHcall<- apply(confusion.matrix.cghcall[,-c(1,2)], 2, function(x)sum(x[true.i])/sum(x))

precision.CGHcall.loss <- as.numeric(confusion.matrix.cghcall[1,-c(1,2)] / apply(confusion.matrix.cghcall[,-c(1,2)],2,function(x) sum(x[1:3])))
precision.CGHcall.loss[is.nan(precision.CGHcall.loss)] = 0


precision.CGHcall.normal <- as.numeric(confusion.matrix.cghcall[5,-c(1,2)] / apply(confusion.matrix.cghcall[,-c(1,2)],2,function(x) sum(x[4:6])))
precision.CGHcall.normal[is.nan(precision.CGHcall.normal)]= 0

precision.CGHcall.gain <- as.numeric(confusion.matrix.cghcall[9,-c(1,2)] / apply(confusion.matrix.cghcall[,-c(1,2)],2,function(x) sum(x[7:9])))
precision.CGHcall.gain[is.nan(precision.CGHcall.gain)]= 0


#recall - specificity
recall.CGHcall.loss <- as.numeric(confusion.matrix.cghcall[1,-c(1,2)] / apply(confusion.matrix.cghcall[,-c(1,2)],2,function(x) sum(x[c(1,4,7)])))
recall.CGHcall.loss[is.nan(recall.CGHcall.loss)]= 0




recall.CGHcall.normal <- as.numeric(confusion.matrix.cghcall[5,-c(1,2)] / apply(confusion.matrix.cghcall[,-c(1,2)],2,function(x) sum(x[c(2,5,8)])))
recall.CGHcall.normal[is.nan(recall.CGHcall.normal)]= 0

recall.CGHcall.gain <- as.numeric(confusion.matrix.cghcall[9,-c(1,2)] /apply(confusion.matrix.cghcall[,-c(1,2)],2,function(x) sum(x[c(3,6,9)])))
hist(as.numeric(recall.CGHcall.gain))


F.CGHcall.gain = (2 * precision.CGHcall.gain * recall.CGHcall.gain) / (precision.CGHcall.gain + recall.CGHcall.gain)
F.CGHcall.normal = (2 * precision.CGHcall.normal * recall.CGHcall.normal) / (precision.CGHcall.normal + recall.CGHcall.normal)
F.CGHcall.loss = (2 * precision.CGHcall.loss * recall.CGHcall.loss) / (precision.CGHcall.loss + recall.CGHcall.loss)

F.score.cghcall = data.frame(loss = F.CGHcall.loss, normal = F.CGHcall.normal, gain = F.CGHcall.gain, Method = "CGHcall", tumour.purity = c(rep(1,100),rep(.7,100),rep(.5,100),rep(.3,100))) 
F.score.cghcall$loss[is.nan(F.score.cghcall$loss)] = 0
F.score.cghcall$normal[is.nan(F.score.cghcall$normal)] = 0
F.score.cghcall$gain[is.nan(F.score.cghcall$gain)] = 0

save(F.score.cghcall, file = 'F.score.cghcall.Rdata')


