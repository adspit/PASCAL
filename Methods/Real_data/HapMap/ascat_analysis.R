# HAPMAP ASCAT analysis

library(ASCAT)


# prepare files for ASCAT

hapmap.data <- read.table("gw6.lrr_baf.txt", h=T, sep = "\t", row.names=1)

range(hapmap.data[,grep("Log.R.Ratio", colnames(hapmap.data))])
# [1] -5.6455  7.7026

range(hapmap.data[,grep("B.Allele.Freq", colnames(hapmap.data))])
# [1] 0 2


LogR <- as.matrix(hapmap.data[,grep("Log.R.Ratio",colnames(hapmap.data))])
BAF <- as.matrix(hapmap.data[,grep("B.Allele.Freq",colnames(hapmap.data))])

#replace 2's by NA
BAF[BAF==2] <- NA


# correct difference between copy number only probes and other probes
CN.probes <- grep("CN",rownames(LogR))
aux <- LogR[CN.probes,]
LogR[CN.probes,] = apply(aux,2, function(x) x-mean(x, na.rm = T))

# limit the number of digits:
LogR = round(LogR,4)

LogR <- cbind(hapmap.data[,1:2],LogR)
BAF <- cbind(hapmap.data[,1:2], BAF)

write.table(LogR,"HapmapLogR.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(BAF,"HapmapBAF.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(hapmap.data[,1:2], "SNPpos.txt", sep="\t",row.names=T, col.names = NA, quote = F)





# detect sex for each sample
gender <- read.table("birdseed.report.txt", sep="\t", skip=66, header=T)
sex <- as.vector(gender[,"computed_gender"])
sex[sex == "female"] <- "XX"
sex[sex == "male"] <- "XY"
sex[sex == "unknown"] <- "XX"


# load data
ascat.bc = ascat.loadData("HapmapLogR.txt","HapmapBAF.txt", chrs=c(1:22, "X"), gender=sex)


# correct for GC content
ascat.bc <- ascat.GCcorrect(ascat.bc, "GC_AffySNP6_102015.txt")
# plot raw data
ascat.plotRawData(ascat.bc)

ascat.gg = ascat.predictGermlineGenotypes(ascat.bc)
# run ASPCF segmentation
ascat.bc = ascat.aspcf(ascat.bc,ascat.gg=ascat.gg)

ascat.plotSegmentedData(ascat.bc)

ascat.output = ascat.runAscat(ascat.bc)
#save ASCAT results


# calc CN
ascat.CN <- ascat.output$nA+ascat.output$nB
range(ascat.CN)
# [1]  0 22

# define the three classes
ascat.CN[ascat.CN<2] <- -1
ascat.CN[ascat.CN==2] <- 0
ascat.CN[ascat.CN>2] <- 1

# ascat.CN <- data.frame( lapply( ascat.CN , factor , levels = c(-1,0,1) ) )
# save.image("HAPMAP.ASCAT.results.RData")
