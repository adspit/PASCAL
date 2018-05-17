
# genoCNA on hapmap
library(genoCN)
setwd("results/")
# apply genoCNA

hapmap.affy.signals <- read.table("/storage/users/adriana.pitea/ASCAT/gw6.lrr_baf.txt", h=T, sep="\t")
dim(hapmap.affy.signals)
# [1] 1844399     185

logRratios = hapmap.affy.signals[, grep('.CEL.Log.R.Ratio', colnames(hapmap.affy.signals))]
baf_signals = hapmap.affy.signals[, grep('.CEL.B.Allele.Freq', colnames(hapmap.affy.signals))]

pfb_file = read.table('affygw6.hg18.pfb.1', sep='\t', stringsAsFactor = FALSE, h=T)

lapply(1:91, function(x) genoCNA(pfb_file$Name, pfb_file$Chr, pfb_file$Position, logRratios[,x], baf_signals[,x], pfb_file$PFB, sampleID=paste0("sample_",x, sep=""), outputSNP=2))
