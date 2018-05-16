# prepare HapMap log2 ratios for CGHcall and CGHcall*
# read in the lrr and baf signals - the gw6.lrr_baf.txt file has been obtained after following the Affy-PennCNA pipeline

hapmap.affy.signals <- read.table("/storage/users/adriana.pitea/ASCAT/gw6.lrr_baf.txt", h=T, sep="\t")
dim(hapmap.affy.signals)
# [1] 1844399     185

# compute mean of LRR hapmap samples
ref.hapmap <- rowMeans(hapmap.affy.signals[,grep(".CEL.Log.R.Ratio",colnames(hapmap.affy.signals))])
length(ref.hapmap)
# [1] 1844399
range(ref.hapmap)
# [1] -2.063100  1.235499
 
# compute log2 ratios for CGHCall
hapmap.log2ratios <- cbind(hapmap.affy.signals[,1:3],apply(hapmap.affy.signals[,grep(".CEL.Log.R.Ratio",colnames(hapmap.affy.signals))],2, function(x) x-ref.hapmap))
dim(hapmap.log2ratios)
# [1] 1844399      94

hapmap.log2ratios[1:3, 1:4]
#           Name Chr Position SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_A02_31426.CEL.Log.R.Ratio
# 1 SNP_A-2131660   1  1145994  0.1511483516
# 2 SNP_A-1967418   1  2224111 -0.0001208791
# 3 SNP_A-1969580   1  2319424 -0.0921549451

# save(hapmap.log2ratios, file="HAPMAP.log2ratios.Rdata")
