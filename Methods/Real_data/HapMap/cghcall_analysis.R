library(CGHcall)
require(intervals)
library(gtools)
load("HAPMAP.log2ratios.Rdata")
range(hapmap.log2ratios[,4:ncol(hapmap.log2ratios)])
# [1] -5.574670  7.627484
dim(hapmap.log2ratios)
# [1] 1844399      94

# CGHcall procedure

data <- data.frame(ProbeID = as.character(hapmap.log2ratios$Name), Chromosome = hapmap.log2ratios$chr, start = hapmap.log2ratios$Position - 12, end = hapmap.log2ratios$Position + 12)

data <- cbind(data,hapmap.log2ratios[,4:ncol(hapmap.log2ratios)])

# add probe coordinates 
hapmap.log2ratios = cbind(hapmap.log2ratios[,1:3],data[,3:ncol(data)])

# create cgh object

cgh <- make_cghRaw(data)

### This function performs the following actions on arrayCGH data:

# • Filter out data with missing position information.

# • Remove data on chromosomes larger than nchrom.

# • Remove rows with more than maxmiss percentage missing values.

# • Imputes missing values using the ‘impute.knn’ function from the impute package.

raw.data <- preprocess(cgh)



### This function normalizes arrayCGH data using the global mode or median. It can also adjust for the cellularity of your data.

normalized.data <- normalize(raw.data)
segmented.data <- segmentData(normalized.data, alpha=0.02)


# recursively searches for the interval containing the most segmented data, decreasing the interval length in each recursion.  The recursive search makes the post-segmentation
# normalization robust against local maxima. This function is particularly useful for profiles for which, after segmentation, the 0-level does not coincide with many segments.



postsegnormalized.data <- postsegnormalize(segmented.data)


# Call aberrations for array CGH data using a six state mixture
perc.tumor <- 0.8
CGHresult <- CGHcall(postsegnormalized.data,cellularity=perc.tumor, build = "GRCh36", ncpus=16)


# create data frame for loading into - see function in ads_function, try combine
# or just set divide=1
CGHresult_expand <- ExpandCGHcall(CGHresult,postsegnormalized.data, memeff=T, divide = 1)

load("Callobj_profiles1to91.Rdata")

cghcall.profiles <- result0
# match rownames
fData(cghcall.profiles)$Name = intersect(rownames(fData(cghcall.profiles)),hapmap.log2ratios$Name)

# save(cghcall.profiles, file="cghcalls.hapmap.Rdata")

