library(gtools)
library(ggplot2)
library(reshape2)

load("F.score.cghcall.star.n10_6.Rdata")
load("F.score.cghcall.n10_6.Rdata")
load("F.score.ascat.n10_6.Rdata")
load("F.score.oncosnp.n10_6.Rdata")


# total F1 score

F.score = rbind(F.score.oncosnp.n.10_6, F.score.ascat.n.10_6, F.score.cghcall.n.10_6, F.score.cghcall.star.n.10_6)
F.score = melt(F.score, id.vars = c("Method","tumour.purity"))
colnames(F.score)[3:4] =c("Class","Value")

F.score$Class = factor(F.score$Class, levels = c("loss","normal","gain"))

F.score$tumour.purity = factor(F.score$tumour.purity, levels = c("0.3","0.5","0.7","1"))

# plot F1 scores

pdf("final.f.scores.class.facewrap.n10_6_new.pdf", width = 20, height = 8, useDingbats = FALSE)
ggplot(data = F.score, aes(x = tumour.purity , y = Value, fill=Method)) 
+ geom_boxplot(aes(tumour.purity, colour = Method),alpha=0.7) 
+ geom_point(position=position_jitterdodge(jitter.width = 0.1, jitter.heigh=0), aes(colour=Method)) 
+ theme_bw()
+scale_color_manual(values= c('#dc828a', "#cddc82", "#9182dc", '#82cddc'),  labels = c('OncoSNP', 'ASCAT', 'CGHcall', 'CGHcall*'))
+theme(panel.grid.major = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=16),plot.title=element_text(size=16), legend.position = "bottom", strip.text = element_text(size=16), axis.text.y=element_text(size=16), axis.text.x=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16))
+scale_fill_manual(values= c('#dc828a', "#cddc82", "#9182dc", '#82cddc'), guide=F)
+scale_y_continuous(limits = c(0, 1))
+ylab("F-score")+xlab(NULL)
+facet_wrap(~Class, ncol = 4)
+ggtitle(expression(paste("Signal length = ",10^{6},"bp")))
dev.off()


#calculate region lengths

load("simulation.results.10_6.Rdata")
region_lengths_n10_6_t1 = do.call(cbind, lapply(c_1, function(x) as.data.frame(table(x$region))$Freq))
rownames(region_lengths_n10_6_t1) = names(table(c_1[[1]]$region))
save(region_lengths_n10_6_t1, file = 'region_lengths_n10_6_t1.Rdata')

# log10 transform CNA region length
logt_region_length_106 = t(log10(region_lengths_n10_6_t1))

# add region length to F.score matrix
head(F.score.cghcall.star.n.10_6)
# loss    normal      gain                Method tumour.purity
# 1 0.9999798 0.9999943 0.9999888 CGHcall_adjusted_last             1
# 2 0.9999951 0.9999976 0.9999972 CGHcall_adjusted_last             1
# 3 0.9999982 0.9999977 1.0000000 CGHcall_adjusted_last             1
# 4 0.9976777 0.9999916 0.9999965 CGHcall_adjusted_last             1
# 5 0.9999766 1.0000000 0.9999959 CGHcall_adjusted_last             1
# 6 0.9999280 0.9999953 0.9999908 CGHcall_adjusted_last             1


sel.f.cghcallstar.106 = cbind(F.score.cghcall.star.n.10_6[1:100,], logt_region_length_106)
# save(sel.f.cghcallstar.106, file = 'cghcallstar.fscores.regions.length.106.Rdata')


# loss
mean(sel.f.cghcallstar.106$loss[sel.f.cghcallstar.106[,6]<4] )
# [1] 0.9999615
mean(sel.f.cghcallstar.106$loss[(sel.f.cghcallstar.106[,6]>=4 & sel.f.cghcallstar.106[,6]<5) | (sel.f.cghcallstar.106[,5]>=4 & sel.f.cghcallstar.106[,5]<5)] )
# [1] 0.9999567
mean(sel.f.cghcallstar.106$loss[sel.f.cghcallstar.106[,6]>5] )
# [1] 0.9999897


# gain
mean(sel.f.cghcallstar.106$gain[sel.f.cghcallstar.106[,6]<4] )
# [1] 0.9999615
mean(sel.f.cghcallstar.106$gain[(sel.f.cghcallstar.106[,6]>=4 & sel.f.cghcallstar.106[,6]<5)  | (sel.f.cghcallstar.106[,5]>=4 & sel.f.cghcallstar.106[,5]<5)] )
# [1] 0.9999742
mean(sel.f.cghcallstar.106$gain[sel.f.cghcallstar.106[,6]>5] )
# [1] 0.999946




# normal
mean(sel.f.cghcallstar.106$normal[sel.f.cghcallstar.106[,6]<4] )
# [1] 0.9999921
mean(sel.f.cghcallstar.106$normal[(sel.f.cghcallstar.106[,6]>=4 & sel.f.cghcallstar.106[,6]<5)  | (sel.f.cghcallstar.106[,5]>=4 & sel.f.cghcallstar.106[,5]<5)] )
# [1] 0.9999958
mean(sel.f.cghcallstar.106$normal[sel.f.cghcallstar.106[,6]>5] )
# [1] 0.9999897


# calculate mean f scores for different region length in samples with tumour purity 1
sel.f.mean.cghcall.star = rbind(cbind(mean(sel.f$loss[sel.f[,6]<4]), mean(sel.f$normal[sel.f[,8]<4]),mean(sel.f$gain[sel.f[,9]<4])), 
                                cbind(mean(sel.f$loss[which(sel.f[,6]>=4 & sel.f[,6]<5)]), mean(sel.f$normal[which(sel.f[,8]>=4 & sel.f[,8]<5)]),mean(sel.f$gain[which(sel.f[,9]>=4 & sel.f[,9]<5)])),
                                cbind(mean(sel.f$loss[sel.f[,6]>=5]),mean(sel.f$normal[sel.f[,8]>=5]), mean(sel.f$gain[sel.f[,9]>=5]) ))
colnames(sel.f.mean.cghcall.star) = c('loss', 'normal', 'gain')

sel.f.cghcall.106 = cbind(F.score.cghcall.n.10_6[1:100,], logt_region_length_106)
save(sel.f.cghcall.106, file = 'cghcall.fscores.regions.length.106.Rdata')

sel.f.ascat.106 = cbind(F.score.ascat.n.10_6[1:100,], logt_region_length_106)
save(sel.f.ascat.106, file = 'ascat.fscores.regions.length.106.Rdata')

sel.f.oncosnp.106 = cbind(F.score.oncosnp.n.10_6[F.score.oncosnp.n.10_6$tumour.purity==1,], logt_region_length_106)
save(sel.f.oncosnp.106, file = 'oncosnp.fscores.regions.length.106.Rdata')





# cna ratios
load("/mnt/storageGluster/users/adriana.pitea/synthetic_data_novel/n1000000/region_lengths_n10_6_t1.Rdata")
cna_ratios_106 = as.data.frame(t(region_lengths_n10_6_t1))
cna_ratios_106$ratios = (cna_ratios_106[,1]+cna_ratios_106[,4])/rowSums(cna_ratios_106)
head(cna_ratios_106)
range(cna_ratios_106$ratios)
# [1]  0.039555 0.938609
save(cna_ratios_106, file = 'cna_ratios_106.Rdata')








load('/storage/users/adriana.pitea/synthetic_data_novel/n10000/region_length_analysis/cghcall.fscores.regions.length.Rdata')
load("/mnt/storageGluster/users/adriana.pitea/synthetic_data_novel/n100000/cghcall.fscores.regions.length.105.Rdata")

# cghcall
mean.fscores.cghcall = rbind(cbind(mean(c(sel.f.cghcall.104$loss, sel.f.cghcall.105$loss[sel.f.cghcall.105[,6]<4], sel.f.cghcall.106$loss[sel.f.cghcall.106[,6]<4])), 
                                   mean(c(sel.f.cghcall.104$normal, sel.f.cghcall.105$normal[sel.f.cghcall.105[,7]<4 | sel.f.cghcall.105[,8]<4], sel.f.cghcall.106$normal[sel.f.cghcall.106[,7]<4 | sel.f.cghcall.106[,8]<4])),
                                   mean(c(sel.f.cghcall.104$gain, sel.f.cghcall.105$gain[sel.f.cghcall.105[,9]<4], sel.f.cghcall.106$gain[sel.f.cghcall.106[,9]<4]))), 
                            cbind(mean(c(sel.f.cghcall.105$loss[sel.f.cghcall.105[,6]>=4], sel.f.cghcall.106$loss[sel.f.cghcall.106[,6]<5 & sel.f.cghcall.106[,6]>=4])),
                                  mean(c(sel.f.cghcall.105$normal[(sel.f.cghcall.105[,7]>=4 & sel.f.cghcall.105[,7]<6) | (sel.f.cghcall.105[,8]>=4 & sel.f.cghcall.105[,8]<6)], sel.f.cghcall.106$normal[(sel.f.cghcall.106[,7]>=4 & sel.f.cghcall.105[,7]<6) | (sel.f.cghcall.106[,8]>=4 & sel.f.cghcall.105[,8]<6)])),
                                  mean(c(sel.f.cghcall.105$gain[sel.f.cghcall.105[,9]>=4], sel.f.cghcall.106$gain[sel.f.cghcall.106[,9]<5 & sel.f.cghcall.106[,9]>=4]))),
                            cbind(mean(sel.f.cghcall.106$loss[sel.f.cghcall.106[,6]>=5]),
                                  mean(sel.f.cghcall.106$normal[sel.f.cghcall.106[,8]>=5 | sel.f.cghcall.106[,7]>=5]), 
                                  mean(sel.f.cghcall.106$gain[sel.f.cghcall.106[,9]>=5])))

colnames(mean.fscores.cghcall) = c('loss', 'normal', 'gain')


load("/mnt/storageGluster/users/adriana.pitea/synthetic_data_novel/n10000/region_length_analysis/cna_ratios_104.Rdata")
meanf_cnaratio_cghcall = rbind(cbind(mean(c(sel.f.cghcall.104$loss[cna_ratios_104$ratio<=0.5], sel.f.cghcall.105$loss[cna_ratios_105$ratio<=0.5], sel.f.cghcall.106$loss[cna_ratios_106$ratio<=0.5])), 
                                     mean(c(sel.f.cghcall.104$normal[cna_ratios_104$ratio<=0.5], sel.f.cghcall.105$normal[cna_ratios_105$ratio<=0.5], sel.f.cghcall.106$normal[cna_ratios_106$ratio<=0.5])),
                                     mean(c(sel.f.cghcall.104$gain[cna_ratios_104$ratio<=0.5], sel.f.cghcall.105$gain[cna_ratios_105$ratio<=0.5], sel.f.cghcall.106$gain[cna_ratios_106$ratio<=0.5]))), 
                               cbind(mean(c(sel.f.cghcall.104$loss[cna_ratios_104$ratio>0.5], sel.f.cghcall.105$loss[cna_ratios_105$ratio>0.5], sel.f.cghcall.106$loss[cna_ratios_106$ratio>0.5])), 
                                     mean(c(sel.f.cghcall.104$normal[cna_ratios_104$ratio>0.5], sel.f.cghcall.105$normal[cna_ratios_105$ratio>0.5], sel.f.cghcall.106$normal[cna_ratios_106$ratio>0.5])),
                                     mean(c(sel.f.cghcall.104$gain[cna_ratios_104$ratio>0.5], sel.f.cghcall.105$gain[cna_ratios_105$ratio>0.5], sel.f.cghcall.106$gain[cna_ratios_106$ratio>0.5]))))

colnames(meanf_cnaratio_cghcall) = c('loss', 'normal', 'gain')


# oncosnp

load("/mnt/storageGluster/users/adriana.pitea/synthetic_data_novel/n100000/oncosnp.fscores.regions.length.105.Rdata")
load("/mnt/storageGluster/users/adriana.pitea/synthetic_data_novel/n10000/region_length_analysis/oncosnp.fscores.regions.length.Rdata")
load("/mnt/storageGluster/users/adriana.pitea/synthetic_data_novel/n1000000/oncosnp.fscores.regions.length.106.Rdata")
mean.fscores.oncosnp = rbind(cbind(mean(c(sel.f.oncosnp.104$loss, sel.f.oncosnp.105$loss[sel.f.oncosnp.105[,6]<4], sel.f.oncosnp.106$loss[sel.f.oncosnp.106[,6]<4])), 
                                   mean(c(sel.f.oncosnp.104$normal, sel.f.oncosnp.105$normal[sel.f.oncosnp.105[,7]<4 | sel.f.oncosnp.105[,8]<4], sel.f.oncosnp.106$normal[sel.f.oncosnp.106[,7]<4 | sel.f.oncosnp.106[,8]<4])),
                                   mean(c(sel.f.oncosnp.104$gain, sel.f.oncosnp.105$gain[sel.f.oncosnp.105[,9]<4], sel.f.oncosnp.106$gain[sel.f.oncosnp.106[,9]<4]))), 
                             cbind(mean(c(sel.f.oncosnp.105$loss[sel.f.oncosnp.105[,6]>=4], sel.f.oncosnp.106$loss[sel.f.oncosnp.106[,6]<5 & sel.f.oncosnp.106[,6]>=4])),
                                   mean(c(sel.f.oncosnp.105$normal[(sel.f.oncosnp.105[,7]>=4 & sel.f.oncosnp.105[,7]<6) | (sel.f.oncosnp.105[,8]>=4 & sel.f.oncosnp.105[,8]<6)], sel.f.oncosnp.106$normal[(sel.f.oncosnp.106[,7]>=4 & sel.f.oncosnp.105[,7]<6) | (sel.f.oncosnp.106[,8]>=4 & sel.f.oncosnp.105[,8]<6)])),
                                   mean(c(sel.f.oncosnp.105$gain[sel.f.oncosnp.105[,9]>=4], sel.f.oncosnp.106$gain[sel.f.oncosnp.106[,9]<5 & sel.f.oncosnp.106[,9]>=4]))),
                             cbind(mean(sel.f.oncosnp.106$loss[sel.f.oncosnp.106[,6]>=5]),
                                   mean(sel.f.oncosnp.106$normal[sel.f.oncosnp.106[,8]>=5 | sel.f.oncosnp.106[,7]>=5]), 
                                   mean(sel.f.oncosnp.106$gain[sel.f.oncosnp.106[,9]>=5])))

colnames(mean.fscores.oncosnp) = c('loss', 'normal', 'gain')


meanf_cnaratio_oncosnp = rbind(cbind(mean(c(sel.f.oncosnp.104$loss[cna_ratios_104$ratio<=0.5], sel.f.oncosnp.105$loss[cna_ratios_105$ratio<=0.5], sel.f.oncosnp.106$loss[cna_ratios_106$ratio<=0.5])), 
                                     mean(c(sel.f.oncosnp.104$normal[cna_ratios_104$ratio<=0.5], sel.f.oncosnp.105$normal[cna_ratios_105$ratio<=0.5], sel.f.oncosnp.106$normal[cna_ratios_106$ratio<=0.5])),
                                     mean(c(sel.f.oncosnp.104$gain[cna_ratios_104$ratio<=0.5], sel.f.oncosnp.105$gain[cna_ratios_105$ratio<=0.5], sel.f.oncosnp.106$gain[cna_ratios_106$ratio<=0.5]))), 
                               cbind(mean(c(sel.f.oncosnp.104$loss[cna_ratios_104$ratio>0.5], sel.f.oncosnp.105$loss[cna_ratios_105$ratio>0.5], sel.f.oncosnp.106$loss[cna_ratios_106$ratio>0.5])), 
                                     mean(c(sel.f.oncosnp.104$normal[cna_ratios_104$ratio>0.5], sel.f.oncosnp.105$normal[cna_ratios_105$ratio>0.5], sel.f.oncosnp.106$normal[cna_ratios_106$ratio>0.5])),
                                     mean(c(sel.f.oncosnp.104$gain[cna_ratios_104$ratio>0.5], sel.f.oncosnp.105$gain[cna_ratios_105$ratio>0.5], sel.f.oncosnp.106$gain[cna_ratios_106$ratio>0.5]))))

colnames(meanf_cnaratio_oncosnp) = c('loss', 'normal', 'gain')




# ascat

load("/mnt/storageGluster/users/adriana.pitea/synthetic_data_novel/n10000/region_length_analysis/ascat.fscores.regions.length.Rdata")
load("/mnt/storageGluster/users/adriana.pitea/synthetic_data_novel/n100000/ascat.fscores.regions.length.105.Rdata")
load("/mnt/storageGluster/users/adriana.pitea/synthetic_data_novel/n1000000/ascat.fscores.regions.length.106.Rdata")

mean.fscores.ascat = rbind(cbind(mean(c(sel.f.ascat.104$loss, sel.f.ascat.105$loss[sel.f.ascat.105[,6]<4], sel.f.ascat.106$loss[sel.f.ascat.106[,6]<4])), 
                                   mean(c(sel.f.ascat.104$normal, sel.f.ascat.105$normal[sel.f.ascat.105[,7]<4 | sel.f.ascat.105[,8]<4], sel.f.ascat.106$normal[sel.f.ascat.106[,7]<4 | sel.f.ascat.106[,8]<4])),
                                   mean(c(sel.f.ascat.104$gain, sel.f.ascat.105$gain[sel.f.ascat.105[,9]<4], sel.f.ascat.106$gain[sel.f.ascat.106[,9]<4]))), 
                             cbind(mean(c(sel.f.ascat.105$loss[sel.f.ascat.105[,6]>=4], sel.f.ascat.106$loss[sel.f.ascat.106[,6]<5 & sel.f.ascat.106[,6]>=4])),
                                   mean(c(sel.f.ascat.105$normal[(sel.f.ascat.105[,7]>=4 & sel.f.ascat.105[,7]<6) | (sel.f.ascat.105[,8]>=4 & sel.f.ascat.105[,8]<6)], sel.f.ascat.106$normal[(sel.f.ascat.106[,7]>=4 & sel.f.ascat.105[,7]<6) | (sel.f.ascat.106[,8]>=4 & sel.f.ascat.105[,8]<6)])),
                                   mean(c(sel.f.ascat.105$gain[sel.f.ascat.105[,9]>=4], sel.f.ascat.106$gain[sel.f.ascat.106[,9]<5 & sel.f.ascat.106[,9]>=4]))),
                             cbind(mean(sel.f.ascat.106$loss[sel.f.ascat.106[,6]>=5]),
                                   mean(sel.f.ascat.106$normal[sel.f.ascat.106[,8]>=5 | sel.f.ascat.106[,7]>=5]), 
                                   mean(sel.f.ascat.106$gain[sel.f.ascat.106[,9]>=5])))

colnames(mean.fscores.ascat) = c('loss', 'normal', 'gain')

meanf_cnaratio_ascat = rbind(cbind(mean(c(sel.f.ascat.104$loss[cna_ratios_104$ratio<=0.5], sel.f.ascat.105$loss[cna_ratios_105$ratio<=0.5], sel.f.ascat.106$loss[cna_ratios_106$ratio<=0.5])), 
                                     mean(c(sel.f.ascat.104$normal[cna_ratios_104$ratio<=0.5], sel.f.ascat.105$normal[cna_ratios_105$ratio<=0.5], sel.f.ascat.106$normal[cna_ratios_106$ratio<=0.5])),
                                     mean(c(sel.f.ascat.104$gain[cna_ratios_104$ratio<=0.5], sel.f.ascat.105$gain[cna_ratios_105$ratio<=0.5], sel.f.ascat.106$gain[cna_ratios_106$ratio<=0.5]))), 
                               cbind(mean(c(sel.f.ascat.104$loss[cna_ratios_104$ratio>0.5], sel.f.ascat.105$loss[cna_ratios_105$ratio>0.5], sel.f.ascat.106$loss[cna_ratios_106$ratio>0.5])), 
                                     mean(c(sel.f.ascat.104$normal[cna_ratios_104$ratio>0.5], sel.f.ascat.105$normal[cna_ratios_105$ratio>0.5], sel.f.ascat.106$normal[cna_ratios_106$ratio>0.5])),
                                     mean(c(sel.f.ascat.104$gain[cna_ratios_104$ratio>0.5], sel.f.ascat.105$gain[cna_ratios_105$ratio>0.5], sel.f.ascat.106$gain[cna_ratios_106$ratio>0.5]))))

colnames(meanf_cnaratio_ascat) = c('loss', 'normal', 'gain')



# cghcallstar

load("cghcallstar.fscores.regions.length.105.Rdata")
load("cghcallstar.fscores.regions.length.104.Rdata")
load("cghcallstar.fscores.regions.length.106.Rdata")
mean.fscores.cghcallstar = rbind(cbind(mean(c(sel.f.cghcallstar.104$loss, sel.f.cghcallstar.105$loss[sel.f.cghcallstar.105[,6]<4], sel.f.cghcallstar.106$loss[sel.f.cghcallstar.106[,6]<4])), 
                                 mean(c(sel.f.cghcallstar.104$normal, sel.f.cghcallstar.105$normal[sel.f.cghcallstar.105[,7]<4 | sel.f.cghcallstar.105[,8]<4], sel.f.cghcallstar.106$normal[sel.f.cghcallstar.106[,7]<4 | sel.f.cghcallstar.106[,8]<4])),
                                 mean(c(sel.f.cghcallstar.104$gain, sel.f.cghcallstar.105$gain[sel.f.cghcallstar.105[,9]<4], sel.f.cghcallstar.106$gain[sel.f.cghcallstar.106[,9]<4]))), 
                           cbind(mean(c(sel.f.cghcallstar.105$loss[sel.f.cghcallstar.105[,6]>=4], sel.f.cghcallstar.106$loss[sel.f.cghcallstar.106[,6]<5 & sel.f.cghcallstar.106[,6]>=4])),
                                 mean(c(sel.f.cghcallstar.105$normal[(sel.f.cghcallstar.105[,7]>=4 & sel.f.cghcallstar.105[,7]<6) | (sel.f.cghcallstar.105[,8]>=4 & sel.f.cghcallstar.105[,8]<6)], sel.f.cghcallstar.106$normal[(sel.f.cghcallstar.106[,7]>=4 & sel.f.cghcallstar.105[,7]<6) | (sel.f.cghcallstar.106[,8]>=4 & sel.f.cghcallstar.105[,8]<6)])),
                                 mean(c(sel.f.cghcallstar.105$gain[sel.f.cghcallstar.105[,9]>=4], sel.f.cghcallstar.106$gain[sel.f.cghcallstar.106[,9]<5 & sel.f.cghcallstar.106[,9]>=4]))),
                           cbind(mean(sel.f.cghcallstar.106$loss[sel.f.cghcallstar.106[,6]>=5]),
                                 mean(sel.f.cghcallstar.106$normal[sel.f.cghcallstar.106[,8]>=5 | sel.f.cghcallstar.106[,7]>=5]), 
                                 mean(sel.f.cghcallstar.106$gain[sel.f.cghcallstar.106[,9]>=5])))
colnames(mean.fscores.cghcallstar) = c('loss', 'normal', 'gain')

meanf_cnaratio_cghcallstar = rbind(cbind(mean(c(sel.f.cghcallstar.104$loss[cna_ratios_104$ratio<=0.5], sel.f.cghcallstar.105$loss[cna_ratios_105$ratio<=0.5], sel.f.cghcallstar.106$loss[cna_ratios_106$ratio<=0.5])), 
                                   mean(c(sel.f.cghcallstar.104$normal[cna_ratios_104$ratio<=0.5], sel.f.cghcallstar.105$normal[cna_ratios_105$ratio<=0.5], sel.f.cghcallstar.106$normal[cna_ratios_106$ratio<=0.5])),
                                   mean(c(sel.f.cghcallstar.104$gain[cna_ratios_104$ratio<=0.5], sel.f.cghcallstar.105$gain[cna_ratios_105$ratio<=0.5], sel.f.cghcallstar.106$gain[cna_ratios_106$ratio<=0.5]))), 
                             cbind(mean(c(sel.f.cghcallstar.104$loss[cna_ratios_104$ratio>0.5], sel.f.cghcallstar.105$loss[cna_ratios_105$ratio>0.5], sel.f.cghcallstar.106$loss[cna_ratios_106$ratio>0.5])), 
                                   mean(c(sel.f.cghcallstar.104$normal[cna_ratios_104$ratio>0.5], sel.f.cghcallstar.105$normal[cna_ratios_105$ratio>0.5], sel.f.cghcallstar.106$normal[cna_ratios_106$ratio>0.5])),
                                   mean(c(sel.f.cghcallstar.104$gain[cna_ratios_104$ratio>0.5], sel.f.cghcallstar.105$gain[cna_ratios_105$ratio>0.5], sel.f.cghcallstar.106$gain[cna_ratios_106$ratio>0.5]))))

colnames(meanf_cnaratio_cghcallstar) = c('loss', 'normal', 'gain')

# final mean fscores aggregated on region length
mean.fscores = rbind(mean.fscores.oncosnp, mean.fscores.ascat, mean.fscores.cghcall, mean.fscores.cghcallstar)
rownames(mean.fscores) = c(rep(c('short','medium', 'long'),4))

# final mean fscore aggregated on CNA burden
mean_fscores_cnaratio = rbind(meanf_cnaratio_oncosnp, meanf_cnaratio_ascat, meanf_cnaratio_cghcall, meanf_cnaratio_cghcallstar)
rownames(mean_fscores_cnaratio) = c(rep(c('CNAratio < 0.5','CNAratio > 0.5'),4))

# following code limits the lowest and highest color to 5%, and 95% of your range, respectively
quantile.range <- quantile(mean.fscores, probs = seq(0, 1, 0.01))
palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)
color.palette  <- colorRampPalette(c("white", "#fc9d8b", "#f93d18"))(10)


pdf('meanfscores_regions_length.pdf')
heatmap.2(mean.fscores, Rowv = F, Colv = F, trace = 'none', cellnote = round(mean.fscores,5), notecol="black", density.info = 'none', margins = c(9,9), dendrogram='none', RowSideColors = c(rep('#dc828a',3), rep('#cddc82',3), rep('#9182dc',3), rep('#82cddc',3)),  col=color.palette, notecex = 2,  key.title = 'Mean F-score', cexCol = 2, cexRow = 2, lhei = c(5,15), breaks = seq(0,1,0.1), srtCol= 0)
par(lend = 1)           # square line ends for the color legend
legend(.7,1.05, xpd=TRUE,     # location of the legend on the heatmap plot
       legend = c("OncoSNP", "ASCAT", "CGHcall", 'CGHcall*'), # category labels
       col = c('#dc828a', "#cddc82", "#9182dc", '#82cddc'),  # color key
       lty= 1,             # line style
       lwd = 10,  # line width
)
dev.off()


pdf('meanfscores_cna_ratio.pdf', height = 7, width = 8.5)
heatmap.2(mean_fscores_cnaratio, Rowv = F, Colv = F, trace = 'none', cellnote = round(mean_fscores_cnaratio,5), notecol="black", density.info = 'none', margins = c(9,13), dendrogram='none', RowSideColors = c(rep('#dc828a',2), rep('#cddc82',2), rep('#9182dc',2), rep('#82cddc',2)),  col=color.palette, notecex = 2,  key.title = 'Mean F-score', cexCol = 2, cexRow = 2, lhei = c(5,15), breaks = seq(0,1,0.1), srtCol= 0)
par(lend = 1)           # square line ends for the color legend
legend(.7,1.05, xpd=TRUE,     # location of the legend on the heatmap plot
       legend = c("OncoSNP", "ASCAT", "CGHcall", 'CGHcall*'), # category labels
       col = c('#dc828a', "#cddc82", "#9182dc", '#82cddc'),  # color key
       lty= 1,             # line style
       lwd = 10,  # line width
)

dev.off()



colnames(logt_region_length_106) = c('loss', '', 'normal', 'gain')
# melted_length = melt(logt_region_length[,c(1,3,4)]) 
# colnames(melted_length) = c('Sample', 'CNstate', 'CN_region_length')
# pdf('region_length_density_n106.pdf')
# ggplot(melted_length, aes(CN_region_length, fill=CNstate, colour = CNstate))+geom_density(alpha = 0.3)+theme_bw()+xlab('log10(CNA region length)')+scale_fill_manual(values = c('indianred1', 'lightgrey', '#56B4E9'))+scale_color_manual(values = c('indianred1', 'grey', '#56B4E9'))+theme(legend.position = "bottom", strip.text = element_text(size=16), axis.text.y=element_text(size=16), axis.text.x=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+ylim(0,1)+xlim(0,6)
# dev.off()


