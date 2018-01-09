# simulate data with jointseg

library(jointseg)
#This vignette illustrates how the jointseg package may be used to generate a variety of copy-number profiles from the same biological “truth”. Such profiles have been used to compare the performance of segmentation methods


## the predefined copy number states 
regions <- c("(1,1)", "(1,2)", "(0,2)", "(0,1)") 

## microarray dataSet from which the data was generated - here we select the Affymetrix data set "GSE29172"
ds <- "GSE29172"

## load real, manually annotated copy number data
rd = loadCnRegionData(dataSet=ds)

# we observe that for tumour fraction 0.7, the total copy number has negative values for the total copy number - we limit the data to the probes that show a total copy number >= 0
regDat <- loadCnRegionData(dataSet=ds, tumorFraction=0.7)
range(regDat$c)
# [1] -0.222 13.008


## function that generates synthetic data

simFUN <- function(dataSet, tumorFraction, n) {
  regDat <- loadCnRegionData(dataSet=dataSet, tumorFraction=tumorFraction) # load data from dataSet with tumour fraction given by the tumorFraction parameter
  # we obser that for tumour purity
  regDat = regDat[regDat$c >0,]
  reg = sample(regions,sample(2:length(regions)))
  ## breakpoint positions - we samples 3 breakpoints from n probes
  bkp <- sort(sample(n, size = 3)) 
  sim <- getCopyNumberDataByResampling(n, bkp=bkp, regions=reg, regData=regDat) # simulate synthetic data
  return(sim)
}




# generate 100 samples with tumour Franction 1
a_1 = do.call(list,lapply(1:100, function(i) {var.sim = simFUN(ds, 1,n)
var.sim$profile$b = as.character(var.sim$profile$b)
write.table(var.sim$profile,paste("t", 1, "sample",i,".txt",sep=""),sep="\t", col.names = T, row.names = F, quote = F)
return(var.sim$profile)}))

# generate 100 samples with tumour Franction 0.7
a_7 = do.call(list,lapply(1:100, function(i) {var.sim = simFUN(ds, .7,n)
write.table(var.sim$profile,paste("t", .7, "sample",i,".txt",sep=""),sep="\t", col.names = T, row.names = F, quote = F)
return(var.sim$profile)}))

# generate 100 samples with tumour Franction 0.5
a_5 = do.call(list,lapply(1:100, function(i) {var.sim = simFUN(ds, .5, n)
write.table(var.sim$profile,paste("t", .5, "sample",i,".txt",sep=""),sep="\t", col.names = T, row.names = F, quote = F)
return(var.sim$profile)}))

# generate 100 samples with tumour Franction 0.3
a_3 = do.call(list,lapply(1:100, function(i) {var.sim = simFUN(ds, .3, n)
write.table(var.sim$profile,paste("t", .3, "sample",i,".txt",sep=""),sep="\t", col.names = T, row.names = F, quote = F)
return(var.sim$profile)}))


# build initial data
initial.regions = do.call(cbind, lapply(a_1, function(x) x$region))
initial.regions = cbind(initial.regions,do.call(cbind, lapply(a_7, function(x) x$region)))
initial.regions = cbind(initial.regions,do.call(cbind, lapply(a_5, function(x) x$region)))
initial.regions = cbind(initial.regions,do.call(cbind, lapply(a_3, function(x) x$region)))
colnames(initial.regions) = colnames(lrr_signals)

# collapse states to loss, normal and gain
initial.regions[initial.regions%in% "(1,1)" | initial.region%in%'(0,2)'] = 0
initial.regions[initial.regions%in% "(0,1)" ] = -1
initial.regions[initial.regions%in% "(1,2)" | initial.regions%in% "(0,3)" | initial.regions%in% "(1,3)" | initial.regions%in% "(2,3)" | initial.regions%in% "(2,2)" ] = 1
initial.regions = as.data.frame(initial.regions)
initial.regions <- data.frame( lapply( initial.regions , factor , levels = c(-1,0,1) ) )

# calculate sample-wise region lengths for each of the states 
region_lengths = do.call(cbind, lapply(initial.regions, function(x) as.data.frame(table(x))$Freq))

region_lengths[,1:5]
# t1sample1.txt t1sample2.txt t1sample3.txt t1sample4.txt t1sample5.txt
# [1,]        601308        415445        970378       1139882        314169
# [2,]       1198692        313649         99710        660118             0
# [3,]             0       1070906        729912             0       1485831
rownames(region_lengths) = as.data.frame(table(initial.regions[,3]))$Var1
save(initial.regions, file ="initial.regions.Rdata")


# simulate samples with 100000 probes

setwd("n100000/")
n = 100000
b_1 = do.call(list,lapply(1:100, function(i) {var.sim = simFUN(ds, 1,n)
write.table(var.sim$profile,paste("t", 1, "sample",i,".txt",sep=""),sep="\t", col.names = T, row.names = F, quote = F)
return(var.sim$profile)}))


b_7 = do.call(list,lapply(1:100, function(i) {var.sim = simFUN(ds, .7,n)
write.table(var.sim$profile,paste("t", .7, "sample",i,".txt",sep=""),sep="\t", col.names = T, row.names = F, quote = F)
return(var.sim$profile)}))


b_5 = do.call(list,lapply(1:100, function(i) {var.sim = simFUN(ds, .5, n)
write.table(var.sim$profile,paste("t", .5, "sample",i,".txt",sep=""),sep="\t", col.names = T, row.names = F, quote = F)
return(var.sim$profile)}))

b_3 = do.call(list,lapply(1:100, function(i) {var.sim = simFUN(ds, .3, n)
write.table(var.sim$profile,paste("t", .3, "sample",i,".txt",sep=""),sep="\t", col.names = T, row.names = F, quote = F)
return(var.sim$profile)}))


# simulate samples with 1000000 probes
setwd("/Users/adriana/Documents/SNP/JointSegSyntheticData/n1000000/")
n = 1000000
c_1 = do.call(list,lapply(1:100, function(i) {var.sim = simFUN(ds, 1,n)
write.table(var.sim$profile,paste("t", 1, "sample",i,".txt",sep=""),sep="\t", col.names = T, row.names = F, quote = F)
return(var.sim$profile)}))


c_7 = do.call(list,lapply(1:100, function(i) {var.sim = simFUN(ds, .7,n)
write.table(var.sim$profile,paste("t", .7, "sample",i,".txt",sep=""),sep="\t", col.names = T, row.names = F, quote = F)
return(var.sim$profile)}))


c_5 = do.call(list,lapply(1:100, function(i) {var.sim = simFUN(ds, .5, n)
write.table(var.sim$profile,paste("t", .5, "sample",i,".txt",sep=""),sep="\t", col.names = T, row.names = F, quote = F)
return(var.sim$profile)}))

c_3 = do.call(list,lapply(1:100, function(i) {var.sim = simFUN(ds, .3, n)
write.table(var.sim$profile,paste("t", .3, "sample",i,".txt",sep=""),sep="\t", col.names = T, row.names = F, quote = F)
return(var.sim$profile)}))



# create LRR matrix for set 10000
 a1_lrr_signals = do.call(cbind, lapply(a_1, function(x) x$c))
 a7_lrr_signals = do.call(cbind, lapply(a_7, function(x) x$c))
 a5_lrr_signals = do.call(cbind, lapply(a_5, function(x) x$c))
 a3_lrr_signals = do.call(cbind, lapply(a_3, function(x) x$c))
 lrr_signals = cbind(a1_lrr_signals, a7_lrr_signals, a5_lrr_signals, a3_lrr_signals)
 # update colnames
 colnames(lrr_signals) = c(mixedsort(dir(pattern = "t1")), mixedsort(dir(pattern = "t0.7")), mixedsort(dir(pattern = "t0.5")), mixedsort(dir(pattern = "t0.3")))
 # normalize lrr signals
 lrr_signals = log2(lrr_signals)-1
 # where log2(2) represents the normal state
 # save lrr matrix
 save(lrr_signals, file = "lrr_signals_10000_probes.Rdata")
 
 # create BAF matrix for set 10000
 a1_baf_signals = do.call(cbind, lapply(a_1, function(x) x$b))
 a7_baf_signals = do.call(cbind, lapply(a_7, function(x) x$b))
 a5_baf_signals = do.call(cbind, lapply(a_5, function(x) x$b))
 a3_baf_signals = do.call(cbind, lapply(a_3, function(x) x$b))
 # update baf matrix colnames
 baf_signals = cbind(a1_baf_signals, a7_baf_signals,a5_baf_signals,a3_baf_signals)
 colnames(baf_signals) = colnames(lrr_signals)
 # save BAF matrix for set 10000
 save(baf_signals, file = "baf_signals_10000_probes.Rdata")


# same for n=100000 and n=1000000, where n represents the number of probes



