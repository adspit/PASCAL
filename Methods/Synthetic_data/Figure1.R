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

# plot F1 scores 
pdf("f.scores.n10_6.pdf", width = 40)
ggplot(data = F.score, aes(x = Class , y = Value, fill=Method)) + geom_boxplot(aes(Class, colour = Method),alpha=0.3) + geom_point(position=position_jitterdodge(jitter.width = 0.1, jitter.heigh=0), aes(colour=Method)) + theme_bw()+scale_color_manual(values= c("mediumslateblue","yellow2", "darkturquoise","hotpink"))+theme(panel.grid.major = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=27), legend.position = "bottom")+scale_fill_manual(values= c("mediumslateblue","yellow2", "darkturquoise", "hotpink"), guide=F)+scale_y_continuous(limits = c(0, 1))+ylab("F-score")+xlab(NULL)+facet_wrap(~tumour.purity, ncol = 4)
dev.off()
