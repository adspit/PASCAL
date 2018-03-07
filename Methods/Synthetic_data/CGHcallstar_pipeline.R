
# first load matrix consisting of the lrr signals saved when generating the synthetic data
# for example: load('lrr_signals.Rdata')

# source the updated normalization and postnormalization functions
source('norm_star.R')
source('postseg_star.R')


# add probe locations to the lrr signals matrix
data <- data.frame(ProbeID = as.character(seq(1,nrow(lrr_signals))),chr = rep(1,nrow(lrr_signals)), start = seq(from = 1, by = 26, length.out = nrow(lrr_signals)), end = seq(from = 26, by = 25, length.out = nrow(lrr_signals)))
data <- cbind(data,lrr_signals)

# create cghRaw object
cgh <- make_cghRaw(data)

# preprocess raw data - e.g. Filter out data with missing position information.
raw.data <- preprocess(cgh, nchrom = 24)

# normalize data
normalized.data = norm_star(raw.data)

# apply segmentation algorithm
segmented.data <- segmentData(normalized.data, alpha=0.05)

# perform postsegmentation normalization
postsegnormalized.data <- postseg_star(segmented.data)

# define tumour purity levels for each sample
perc.tumor <- c(rep(1,100),rep(.7,100),rep(.5,100),rep(.3,100))


# create and set the wd to CGHcallstar_results dir. 
# run the CGHcall algorithm. The cellularity parameter represents the tumour purity.
setwd('CGHcallstar_results')
CGHresult.adjusted <- CGHcall(postsegnormalized.data.adjusted,cellularity=perc.tumor, ncpus = 14)

# obtain profile after running CGHcall - this return a Rdata object
CGHresult_expand <- ExpandCGHcall(CGHresult,postsegnormalized.data, memeff=T, divide = 1)

