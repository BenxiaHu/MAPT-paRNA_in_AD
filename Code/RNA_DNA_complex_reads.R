rm(list=ls())
options(stringsAsFactors=F)
library(dplyr)
library(writexl)
library(tidyr)
library(ggplot2)
options(width = 1000)
library(data.table)
library(easyGgplot2)
library(ggpubr)

cellid <- c("Ex","In","Ast","Mic","Oli","Opc")
disid <- c("CTRL")
result <- data.frame()
for(i in 1:length(cellid)){
  for(j in 1:length(disid)){
    RNAinput <- unique(fread(paste0(cellid[i],"_",disid[j],"_RNA_DNA_interaction_RNAreads.txt"),sep="\t",header=T,nThread=10,data.table=F)) 
    RNAinput <- RNAinput %>% group_by(CB_CBMB) %>% summarise(RNAread = n()) %>% as.data.frame()
    DNAinput <- unique(fread(paste0(cellid[i],"_",disid[j],"_RNA_DNA_interaction_DNAreads.txt"),sep="\t",header=T,nThread=10,data.table=F))
    DNAinput <- DNAinput %>% group_by(CB_CBMB) %>% summarise(DNAread = n()) %>% as.data.frame()
    out <- RNAinput %>% inner_join(DNAinput, by=c("CB_CBMB")) %>% as.data.frame()
    out$sampleid <- paste0(cellid[i],"_",disid[j])
    result <- rbind(result,out)
  }
}
result$molecues <- result$RNAread + result$DNAread
result <- result %>% separate(CB_CBMB, c("CB", "CBMB"),sep="=")

result$sampleid <- factor(result$sampleid,levels=paste0(cellid,"_CTRL"))

#result2 <- result[result$RNAread<50 & result$DNAread<1000,]
result2$score <- (result2[,'RNAread'] * result2[,'DNAread'])/(result2[,'RNAread'] + result2[,'DNAread'])
cellid <- c("Ex","In","Ast","Mic","Oli","Opc")
disid <- c("CTRL")
for(i in 1:length(cellid)){
  output <- result2[result2$sampleid %in% paste0(cellid[i],"_CTRL"),]
  write.table(file=paste0(cellid[i],"_CTRL_RNA_DNA_interaction_normalization.txt"),output,sep="\t",quote=F,row.names = F)
}
