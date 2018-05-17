# process genocna hapmap results
setwd('/home/icb/adriana.pitea/for_adriana_old/HAPMAP/GenoCNA_res/')
result_files_snp_level_genocna_hapmap = mixedsort(dir(pattern="_SNP.txt"))
hapmap_genocna_results = lapply(result_files_snp_level_genocna_hapmap, function(x) { aux = read.table(x, sep="\t", h=T, stringsAsFactors=F, row.names = 1)})

length(hapmap_genocna_results)
# [1] 91
names(hapmap_genocna_results) = colnames(logRratios)
cn_calls = lapply(1:91, function(x) {colnames(hapmap_genocna_results[[x]])[3] = names(hapmap_genocna_results)[x]; hapmap_genocna_results[[x]][3]})

test = Reduce(function(x, y) merge(x, y, all = TRUE), lapply(cn_calls, function(y) data.table(y, keep.rownames=TRUE, key = "rn")))
predicted_states = test[,-1]
dim(predicted_states)
# [1] 1844399      92

predicted_states = as.matrix(predicted_states)
rownames(predicted_states) = test$rn
predicted_states[predicted_states<2] = -1
predicted_states[predicted_states==2] = 0
predicted_states[predicted_states>2] =  1
predicted_states_hapmap_genocna = as.data.frame(predicted_states)
genocna.calls.hapmap <- data.frame(lapply(predicted_states_hapmap_genocna, factor , levels = c(-1,0,1) ) )
rownames(genocna.calls.hapmap) = rownames(predicted_states_hapmap_genocna)
save(genocna.calls.hapmap, file = 'genocna_calls_hapmap.Rdata')
