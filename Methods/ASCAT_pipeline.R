# run ASCAT procedure
# we first load lrr and baf signals matrices

# define probe locations to the lrr signals matrix
data <- data.frame(ProbeID = as.character(seq(1,nrow(lrr_signals))),chr = rep(1,nrow(lrr_signals)), start = seq(from = 1, by = 26, length.out = nrow(lrr_signals)), end = seq(from = 26, by = 25, length.out = nrow(lrr_signals)))

data.ascat <- data.frame(chr = data$chr, pos = data$start)
# rename dataframe rows
rownames(data.ascat) <- paste("SNP",data$ProbeID,sep="")
# add probe locations to the lrr and baf signals
data.ascat.tumor.lrr <- cbind(data.ascat, lrr_signals)
data.ascat.tumor.baf <- cbind(data.ascat, baf_signals)

# save lrr and baf signals
write.table(data.ascat.tumor.lrr, "Tumor_LogR.txt", col.names = T, row.names = T, sep="\t")
write.table(data.ascat.tumor.baf, "Tumor_BAF.txt", col.names = T, row.names = T, sep="\t")

# use the ascat load function to load the data in the erquired format
ascat.bc = ascat.loadData("Tumor_LogR.txt","Tumor_BAF.txt")

# set WD for ASCAT processed files and results
setwd("ASCAT/")

# plot raw data
ascat.plotRawData(ascat.bc)

# predict germline genotypes for AffySNP6-type data
ascat.gg = ascat.predictGermlineGenotypes(ascat.bc, platform = "AffySNP6") 

# run ASPCF segmentation
ascat.bc = ascat.aspcf(ascat.bc,ascat.gg=ascat.gg)

# plot segmented data
ascat.plotSegmentedData(ascat.bc)

# run ASCAT
ascat.output = ascat.runAscat(ascat.bc)  


# save ASCAT results
save(ascat.output, file='ASCAT_output.Rdata')


