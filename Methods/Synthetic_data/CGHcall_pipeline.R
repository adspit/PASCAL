

# first load matrix consisting of the lrr signals saved when generating the synthetic data
# for example: load('lrr_signals.Rdata')


#change fData

# load Affymetrix coordinates
load('Affymetrix_SNP_probe_names.Rdata')
probe_names$Position = probe_names$Start+12
# analysis on the new data
snp_probes_featData = data.frame(Name = rownames(probe_names),  Chr = probe_names$Chromosome, Position = probe_names$Position)
all.equal(snp_probes_featData, snp_probes_featData[order(snp_probes_featData$Chr, snp_probes_featData$Position),])
# [1] TRUE



# add probe locations to the lrr signals matrix
data <- data.frame(ProbeID = as.character(seq(1,nrow(lrr_signals))),chr = rep(1,nrow(lrr_signals)), start = seq(from = 1, by = 26, length.out = nrow(lrr_signals)), end = seq(from = 26, by = 25, length.out = nrow(lrr_signals)))
data <- cbind(data,lrr_signals)

# create cghRaw object
cgh <- make_cghRaw(data)

# preprocess raw data - e.g. Filter out data with missing position information.
raw.data <- preprocess(cgh)

# normalize data
normalized.data = normalize(raw.data)

# apply segmentation algorithm
segmented.data <- segmentData(normalized.data, alpha=0.05)

# perform postsegmentation normalization
postsegnormalized.data <- postsegnormalize(segmented.data)

# define tumour purity levels for each sample
perc.tumor <- c(rep(1,100),rep(.7,100),rep(.5,100),rep(.3,100))

# run the CGHcall algorithm. The cellularity parameter represents the tumour purity.
CGHresult <- CGHcall(postsegnormalized.data ,cellularity=perc.tumor, ncpus = 14)

# obtain profile after running CGHcall - this return a Rdata object
CGHresult_expand <- ExpandCGHcall(CGHresult,postsegnormalized.data, memeff=T, divide = 1)

