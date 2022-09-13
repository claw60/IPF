#Figure 4
rm(list = ls())
library(DESeq2)
options(stringsAsFactors = F)
countsdata=read.table("Figure 5/expr_matrix.txt",header = T,row.names = 1,sep = "\t")
head(countsdata)
colnames(countsdata)
countdata=countsdata[,c(6:12)]
head(countdata)
library(stringr)
tmp <- str_split(rownames(countdata),'\\.',simplify = T)[,1]
tmp
rownames(countdata) <- tmp
coldata <- data.frame(row.names = colnames(countdata))
sample <- c(rep("BLM",4),rep("HBO",3))
coldata$sample <- relevel(factor(sample),"BLM")
head(coldata)
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design=~sample)
nrow(dds)
rownames(dds)
##
dds <- dds[rowSums(counts(dds))>1,]
nrow(dds) 
dds <- DESeq(dds)

normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
colnames(normalized_counts)
colnames(normalized_counts)=c("BLM65_Normalized","BLM73_Normalized","BLM82_Normalized","BLM91_Normalized",
                              "HBO60_Normalized","HBO79_Normalized","HBO83_Normalized")
normalized_counts$gene_id <- rownames(normalized_counts)
resultsNames(dds)
contrast <- c("sample", "HBO", "BLM")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
library(dplyr)
library(tibble)
res <- dd1 %>% 
  data.frame() %>% 
  rownames_to_column("gene_id")
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
rs <- bitr(unique(res$gene_id), fromType = "ENSEMBL", 
           toType = c("SYMBOL","ENTREZID","GENENAME"),
           OrgDb = org.Mm.eg.db)
colnames(rs)
res_merge <- inner_join(res,rs,by=c("gene_id"="ENSEMBL"))
res_merge$change = as.factor(ifelse(res_merge$pvalue < 0.05 & res_merge$log2FoldChange > 0.57 ,"up",
                                    ifelse(res_merge$log2FoldChange < -0.57 & res_merge$pvalue<0.05 ,"down","not")))
table(res_merge$change)
countdata$gene_id <- rownames(countdata)
count_merge <- merge(res_merge,countdata,by='gene_id')

deg_merge <- merge(count_merge,normalized_counts,by='gene_id')
deg_merge %>% distinct(SYMBOL,.keep_all = T) -> deg_merge1
colnames(deg_merge1)
normalized_data <- deg_merge1[,c(8,19:25)]
#save(normalized_data,coldata,file = "Figure6_input.Rdata")
###
##Figure 4A
library(ggrepel)
library(ggplot2)
colnames(res_merge)
vol_plot <- ggplot(res_merge,aes(x=log2FoldChange,y= -log10(pvalue),color=change))+
  geom_point(data = res_merge[res_merge$pvalue<0.05&abs(res_merge$log2FoldChange)> 0.57,],size = 1.5)+ 
  geom_point(data = res_merge[res_merge$pvalue>0.05|abs(res_merge$log2FoldChange)< 0.57,],size = 1.5)+
  scale_color_manual(values=c("#4393C3","#00000033","#FC4E2A"))+
  ylab('-log10 (Pvalue)')+
  xlab('log2 (FC)')+
  geom_vline(xintercept=c(-0.57,0.57),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5) +
  theme_classic(  
    base_line_size = 1 
  )+
  theme(axis.title.x = element_text(size = 12, 
                                    color = "black"),
        axis.title.y = element_text(size = 12,
                                    color = "black", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(color="black", 
                                   size = 10, ),
        axis.text.x = element_text(size = 12,  
                                   color = "black",  
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0), 
        axis.text.y = element_text(size = 12,  
                                   color = "black",
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )+theme(legend.position = 'none')
vol_plot
ggsave("vol_lung(HBO).pdf",vol_plot,width = 4,height = 4)
###Figure 4B was done in Metascape.
###Figure 4C
library(msigdbr)
msigdbr_species() 
library(dplyr)
# mus_msigdb <- msigdbr(species="Mus musculus")
# head(mus_msigdb, 2) %>% as.data.frame
m_df = msigdbr(species = "Mus musculus") %>% dplyr::filter(gs_cat == "H")
#m_tg <- m_df %>% select(gene_symbol,human_gene_symbol)
m_t2g = m_df %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
###
geneList=deg_merge1$log2FoldChange
names(geneList)=deg_merge1$ENTREZID
geneList=sort(geneList,decreasing = T)
head(geneList)
##
Mu_GSEA <- GSEA(
  geneList,
  exponent = 1,
  minGSSize = 1,
  maxGSSize = 500,
  eps = 0,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE=m_t2g,
  verbose = F,
  seed = T
)
deg_gs<- Mu_GSEA@result
deg_gs_sig <- deg_gs[deg_gs$p.adjust < 0.25,]
#write.table(deg_gs_sig,"GO_GSEA/Hallmark_gsea.txt",sep = "\t",row.names = F,quote = F)
library(ggplot2)
deg_gs_sig <- deg_gs_sig[order(deg_gs_sig$p.adjust,decreasing = F),]
colnames(deg_gs_sig)
deg_gs_sig$Description <-  factor(deg_gs_sig$Description,levels = rev(deg_gs_sig$Description))
pf <- ggplot(deg_gs_sig,aes(NES,Description))+geom_point(aes(size=setSize,color=-log10(p.adjust)))+
  scale_color_gradient(low = "blue", high = "red")+theme_bw()+
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(color = "black",size = 12),
        legend.text = element_text(size = 12),legend.title=element_text(size=12),
        axis.title.x = element_text(size = 12))+
  theme(panel.spacing.y = unit(0,"cm"))+ylab(NULL)+xlab("NES")

pf
ggsave("gsea.pdf",pf,height = 4.5,width = 6.2)

##Figure 5
##Figure 5 A-C
rownames(normalized_data) <- normalized_data$SYMBOL
normalized_data <- normalized_data[,-1]
boxplot(normalized_data)
log_norm <- log2(normalized_data+1)
boxplot(log_norm)
hypoxia_list=list(c('Ndrg1','Pgam1','Tubb6','Tpi1','Mif','Acot7','Cdkn3',
                    'Aldoa','Ldha','Eno1','Adm','Mrps17','Slc2a1','Vegfa','P4ha1'))
library(GSVA)
gsva_Hypoxia<- gsva(as.matrix(log_norm), hypoxia_list,verbose=TRUE,mx.dif=1,kcdf = "Gaussian")
###
hypoxia <- as.data.frame(t(gsva_Hypoxia))
colnames(hypoxia) <- c("HIF_Score")
hypoxia$group <- coldata$sample
write.csv(hypoxia,file = "HIF_Score.csv",row.names = T,quote = T)
##
library(dplyr)
#install.packages("Hmisc")
library(Hmisc)
hal_gs <- lapply(readLines("gmt/h.all.v7.4.symbols.gmt"), function(x){
  y=strsplit(x,'\t')[[1]]
  y=y[3:length(y)]
  y = Hmisc::capitalize(tolower(y))
  return(y)
})
##name
library(stringr)
tmp <- unlist(lapply(readLines("gmt/h.all.v7.4.symbols.gmt"), function(x){
  y=strsplit(x,'\t')[[1]][1]
  y = unlist(strsplit(y, split="_", fixed=T))[-1]
  y = Hmisc::capitalize(tolower(y))
  y = paste(y, collapse=" ")
})
)
gsva_score <-  gsva(as.matrix(log_norm), hal_gs,verbose=TRUE,mx.dif=1)
rownames(gsva_score) <- tmp
gsva_melt <- as.data.frame(t(gsva_score))
gsva_melt$group <- coldata$sample
colnames(gsva_melt)
library(ggpubr)
EMT_Glycolysis_Score <- gsva_melt[,c(30,25,51)]
colnames(EMT_Glycolysis_Score)[1:2] <- c("EMT_Score","Glycolysis_Score")
write.csv(EMT_Glycolysis_Score,file = "gsva_melt.csv",row.names = T,quote = T)






