# GenoCNA run

task.id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
setwd('/storage_scratch/users/adriana.pitea/GenoCNA/GenoCNA_results/')
load('/home/icb/adriana.pitea/for_adriana_old/input_data_cghcall.Rdata')
rownames(data) = data$ProbeID

data = data[,]

snp.info = read.table('/storage_scratch/users/adriana.pitea/GenoCNA/affygw6.hg18.pfb.1', h=T, sep="\t", stringsAsFactors=F )
snp.info = snp.info[-which(snp.info$Chr=='MT'),]
rownames(snp.info) = snp.info$Name
snp.info$Chr[snp.info$Chr=="X"]=23
snp.info$Chr[snp.info$Chr=="Y"]=24
snp.info$Chr = as.numeric(snp.info$Chr)

snp.info = snp.info[intersect(rownames(dat_tcga), rownames(snp.info)), ]



# x is going to be the task ID
genoCNA(snp.info$Name, snp.info$Chr, snp.info$Position, data[,$SGE_TASK_ID+4], baf_signals[,$SGE_TASK_ID] , snp.info$PFB, sampleID=paste0("sample_t07_",$SGE_TASK_ID, sep=""), outputSNP=2,  cnv.only=(snp.info$PFB>1))
