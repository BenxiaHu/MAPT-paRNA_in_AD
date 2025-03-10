rm(list=ls())
options(stringsAsFactors=F)
library(parallel)
library(stats4)
library(BiocGenerics)
library(Biobase)
library(IRanges)
library(S4Vectors)
library(data.table)
library(ggplot2)
library(gprofiler2)
library(forcats)

GO<-function(input,index,filename,value){
  geneExpr <- fread("MAPT_uaRNA_KD_RNAseq_FPKM.bed",sep="\t",header=T,data.table=F,nThread=10)
  geneExpr <- geneExpr[,index]
  geneExpr$ave <- apply(geneExpr[,-1],1,mean)
  ref_gene <- unique(geneExpr[geneExpr$ave>1,'gene'])
  print(head(ref_gene))
  goresult <- gost(input,organism = "hsapiens",ordered_query=F,significant=T,evcodes = TRUE,custom_bg =ref_gene,
                   domain_scope = "custom",user_threshold=0.05, correction_method="fdr",sources=c("GO"))$result
  
  goresult <- data.frame(term_id=goresult$term_id,p_value=goresult$p_value,term_name=goresult$term_name,gene=goresult$intersection,size=goresult$term_size)
  goresult <- goresult[goresult$size>15 & goresult$size<600,]
  
  GO_order<-goresult[order(goresult$p_value),][c(1:10),]
  pdf(paste0("GO_",filename,"_RNAseq_gProfiler2_1.5fold_hg19.pdf"),width=value)
  p<-ggplot(data=GO_order, aes(y=-log10(p_value), x=reorder(term_name,-log10(p_value)))) +
    geom_bar(stat="identity", width=0.5)+labs(title=paste0("GO terms for ",filename),y="-log10(padj)", x = "Term names")+ coord_flip()+
    theme_classic()+ geom_hline(yintercept =-log10(0.05),color = "red")
  print(p)
  dev.off()
}

MAPT_uaRNA_ASO1 <- read.table("MAPT_uaRNA_ASO1_MAPT_uaRNA_Ctrl_RNAseq_DESeq2_hg19.txt",sep="\t",header=T)
MAPT_uaRNA_ASO1 <- na.omit(MAPT_uaRNA_ASO1)
MAPT_uaRNA_ASO1_upgene <- MAPT_uaRNA_ASO1[MAPT_uaRNA_ASO1$padj<=0.05 & MAPT_uaRNA_ASO1$log2FoldChange > log2(1.5),]
MAPT_uaRNA_ASO1_downgene <- MAPT_uaRNA_ASO1[MAPT_uaRNA_ASO1$padj<=0.05 & MAPT_uaRNA_ASO1$log2FoldChange < (-log2(1.5)),]

GO(MAPT_uaRNA_ASO1_upgene[,2],"MAPT_uaRNA_ASO1_upregulated",7)
GO(MAPT_uaRNA_ASO1_downgene[,2],"MAPT_uaRNA_ASO1_downregulated",10)
