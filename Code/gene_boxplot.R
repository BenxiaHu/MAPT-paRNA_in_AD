rm(list=ls())
library(tidyr)
library(ggplot2)
library(reshape2)

expr <- read.table("AD_RNAseq_expr.bed",sep="\t",header=T)
colnames(expr)[3:ncol(expr)] <- c('Normal_brain12','Normal_brain13','Normal_brain14','Normal_brain16',
                                  'AD_brain4','AD_brain5','AD_brain7','Normal_brain10','AD_brain6',
                                  'AD_brain8','Normal_brain9','AD_brain2')
expr <- expr[,c('Geneid','Genename','Normal_brain9','Normal_brain10','Normal_brain12','Normal_brain13',
                'Normal_brain14','Normal_brain16','AD_brain2','AD_brain4','AD_brain5','AD_brain6','AD_brain7','AD_brain8')]

m6Agene <- c('YAP1')

expr <- expr[expr$Genename %in% m6Agene,]
expr$Genename <- factor(expr$Genename, levels = m6Agene)
expr <- expr[order(expr$Genename), ]


RNAseq <- melt(expr[,-1])
RNAseq <- RNAseq %>% separate(variable, c("sampleid", NA)) %>% as.data.frame()

RNAseq$sampleid <- factor(RNAseq$sampleid,levels=c("Normal","AD"))

p1 <- ggplot(RNAseq[RNAseq$Genename %in% "YAP1",], aes(x=sampleid, y=log2(value), color=sampleid)) + 
  geom_boxplot() + ggtitle("YAP1") + theme_classic() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  ylab('log2(FPKM)')

pdf("DEG_expression_boxplot.pdf")
print(p1)
dev.off()

