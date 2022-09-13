##Supplement_Figure 2 
##Affymetrix data preprocess
rm(list = ls())
library(affy)
list.files()
data <- ReadAffy() 
boxplot(data,las=2)
eset <- affy::rma(data)
class(eset)
boxplot(eset,las=2)
exp=exprs(eset)
#probe id transfer
annotation(eset)
library(GEOquery)
library(mouse430a2.db)
ls("package:mouse430a2.db")
ids=toTable(mouse4302SYMBOL)
table(rownames(exp) %in% ids$probe_id)
dim(exp)
exp = exp[rownames(exp) %in% ids$probe_id,]
dim(exp)
ids=ids[match(rownames(exp),ids$probe_id),]
ID <- rownames(exp)
eset <- exp[ID %in% ids$ID,] %>% cbind(ids) 
eset<-aggregate(x=eset[,1:i],by=list(eset$symbol),FUN=mean,na.rm=T)
write.table(eset,file = "Affy_GSE25640.txt")
###Agilent data preprocess
setwd('raw/') 
files <- dir(pattern="*\\.txt$") 
library(limma)
x <- read.maimages(files, source="agilent", green.only=TRUE, other.columns="gIsWellAboveBG") ## read ÎÄ¼þ
x <- limma::backgroundCorrect(x, method="normexp") 
dim(x)
y <- limma::normalizeBetweenArrays(x, method="quantile")
class(y)
y=as.data.frame(y)
Pos_Control <- y$ControlType==1;table(Pos_Control)
Neg_Control <- y$ControlType==-1;table(Neg_Control)
Isdup <- duplicated(y$GeneName);table(Isdup)
yfilt <- y[!Pos_Control & !Neg_Control & !Isdup, ]
dim(yfilt)
colnames(yfilt)
exp=yfilt[,c(7,11:20)];rownames(exp)=exp[,1]; exp=exp[,-1]
colnames(exp)
colnames(exp)
#probe id transfer
b <-data.table::fread("../../GPL13912_old_annotations.txt",skip ="ID")
con1 <- b$CONTROL_TYPE=='pos';table(con1)
con2 <- b$CONTROL_TYPE=='neg';table(con2)
b = b[!con1 & !con2, ]
Isdup <- duplicated(b$GENE_SYMBOL);table(Isdup)
b = b[!Isdup,]
library(dplyr)
library(tidyr)
probe_id <- b$NAME
symbol <- b$GENE_SYMBOL
ids=data.frame(probe_id=probe_id,symbol=symbol)
ids$symbol[symbol==''] =NA
ids=na.omit(ids)
table(rownames(exp) %in% ids$probe_id)
dim(exp)
exp = exp[rownames(exp) %in% ids$probe_id,]
dim(exp)
ids=ids[match(rownames(exp),ids$probe_id),]
ID <- rownames(exp)
eset <- exp[ID %in% ids$ID,] %>% cbind(ids) 
eset<-aggregate(x=eset[,1:i],by=list(eset$symbol),FUN=mean,na.rm=T)
write.table(eset,file = "Agilent_GSE34814.txt")
##Illumina data preprocess
rm(list=ls())
library(lumi)
fileName <- 'GSE37635_non-normalized.txt' 
x.lumi <- lumiR.batch(fileName) 
pData(phenoData(x.lumi))
## Do all the default preprocessing in one step
lumi.N.Q <- lumiExpresso(x.lumi)
### retrieve normalized data
dataMatrix <- exprs(lumi.N.Q)#
eSet[[1]]@annotation
##probe transfer ###
library(illuminaMousev2.db)
ls("package:illuminaMousev2.db")
ids=toTable(illuminaMousev2SYMBOL)
head(ids)
length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
table(rownames(exprSet) %in% ids$probe_id)
dim(exprSet)
exprSet = exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(exprSet)
ids=ids[match(rownames(exprSet),ids$probe_id),]
ID <- rownames(exp)
eset <- exp[ID %in% ids$ID,] %>% cbind(ids) 
eset<-aggregate(x=eset[,1:i],by=list(eset$symbol),FUN=mean,na.rm=T)
write.table(eset,file = "Agilent_GSE34814.txt")
## data merging
#Agilent 
gse112827=read.table("raw_exp/GSE112827_exp.txt",sep = "\t",header = T)
gse97825=read.table("raw_exp/GSE97825_exp.txt",sep = "\t",header = T)
gse97826=read.table("raw_exp/GSE97826_exp.txt",sep = "\t",header = T)
gse34814=read.table("raw_exp/GSE34814_exp.txt",sep = "\t",header = T)
Ag <- Reduce(intersect,list( rownames(gse112827),
                             rownames(gse97825),
                             rownames(gse97826),
                             rownames(gse34814)))
#Affymetrix 
gse40151=read.table("raw_exp/GSE40151exp.txt",sep = "\t",header = T)
gse18800=read.table("raw_exp/GSE18800exp.txt",sep = "\t",header = T)
gse16846=read.table("raw_exp/GSE16846exp.txt",sep = "\t",header = T)
gse25640=read.table("raw_exp/GSE25640.txt",sep = "\t",header = T)
gse123293=read.table("raw_exp/GSE123293_exp.txt",sep = "\t",header = T )
gse123293_young <- gse123293[,-c(10:15)]
Af <- Reduce(intersect,list( rownames(gse40151),
                             rownames(gse18800),
                             rownames(gse16846),
                             rownames(gse25640),
                             rownames(gse123293_young)))
gse37635 <- read.table("raw_exp/GSE37635.txt",sep = "\t")
merge= Reduce(intersect,list(Ag,Af,rownames(gse37635)))##
##Agilent 
m_gse112827=gse112827[merge,]
m_gse97825=gse97825[merge,]
m_gse97826=gse97826[merge,]
m_gse34814=gse34814[merge,]
##Affymetrix
m_gse40151=gse40151[merge,]
m_gse18800=gse18800[merge,]
m_gse16846=gse16846[merge,]
m_gse25640=gse25640[merge,]
##
m_gse123293=gse123293_young[merge,]
m_gse37635 = gse37635[merge,]
####raw merge
raw_exp=cbind(m_gse112827,m_gse97825,m_gse97826,m_gse34814,
              m_gse40151,m_gse18800,m_gse16846,m_gse25640,m_gse123293,m_gse37635)
##step1.Rdata is the merged data
load("step1.Rdata")
library(ggplot2)
library(umap)
ex=as.data.frame(t(raw_exp))
ex$group=Raw_group$time
#ex$group =factor(raw_group$group,levels= c('Con','BLM'))
ex$group=factor(Raw_group$time,levels = c('Con','BLM_1d','BLM_2d','BLM_7d','BLM_14d',
                                          'BLM_21d','BLM_28d','BLM_35d'))
ump <- umap(ex[,1:(ncol(ex)-1)], n_neighbors = 10, random_state = 123)
dat <- ump$layout
dat<-data.frame(dat,ex$group)
colnames(dat)=c('UMAP1','UMAP2','group')
library(ggsci)
before <- ggplot(dat,aes(UMAP1,UMAP2)) +
  geom_point(aes(colour=group,fill=group),size=1) +
  scale_color_aaas()+
  theme_classic()+
  theme(axis.line = element_line(size = 1))+
  theme(plot.title = element_text(size = 10),
        axis.text=element_text(colour = 'black',size=10),axis.title.x=element_text(colour = 'black',size=10),
        axis.title.y=element_text(colour = 'black',size=10),legend.text=element_text(color = 'black',size = 10),
  )
before
#ggsave("umap/before_time.pdf",p1,width = 4,height = 3)
table(Raw_group$time)
raw_g=Raw_group$group
group=Raw_group$time
mod = model.matrix(~as.factor(raw_g))
batch = Raw_group$batch
####remove batch effect
library(sva)
combat_norm_exp=ComBat(raw_exp,batch = batch,mod = mod)
#boxplot(combat_norm_exp,las=2)
colnames(combat_norm_exp)=paste(Raw_group$time,1:ncol(combat_norm_exp),sep = '')
d <- dist(t(combat_norm_exp),method = 'euclidean')
hc=hclust(d,method = "ward.D2")

library(factoextra)
fviz_dend(hc, k=2, cex = 0.5,horiz = F,
          k_colors = c("#2E9FDF", "#00AFBB"),
          color_labels_by_k = TRUE, rect = F)
cuts=cutree(hc,k=2)
Raw_group$cluster <- cuts
#write.csv(Raw_group,file = "combat_cluster_group.csv",row.names = F)
for (i in 1:2){
  print(paste("Sample in Cluster ",i))
  print(rownames(t(combat_norm_exp))[cuts==i])
  print(" ")
}
###remove three BLM_21d  and one Con sample
fi_exp=combat_norm_exp[,-c(107,109,175,177)]
fi_group= Raw_group[-c(107,109,175,177),]
table(fi_group$time)
###
library(ggplot2)
library(umap)
ex=as.data.frame(t(fi_exp))
ex$group=fi_group$time
#ex$group=factor(ex$group,levels = c("Con","BLM"))
ex$group =factor(fi_group$time,levels= c('Con','BLM_1d','BLM_2d','BLM_7d',
                                          'BLM_14d','BLM_21d','BLM_28d','BLM_35d'))
ump <- umap(ex[,1:(ncol(ex)-1)], n_neighbors = 30, random_state = 123)

dat <- ump$layout
dat<-data.frame(dat,ex$group)
colnames(dat)=c('UMAP1','UMAP2','group')
library(ggsci)
after <- ggplot(dat,aes(UMAP1,UMAP2)) +
  geom_point(aes(colour=group,fill=group),size=1) +
  scale_color_aaas()+
  theme_classic()+
  theme(axis.line = element_line(size = 1))+
  theme(plot.title = element_text(size = 10),
        axis.text=element_text(colour = 'black',size=10),axis.title.x=element_text(colour = 'black',size=12),
        axis.title.y=element_text(colour = 'black',size=12),legend.text=element_text(color = 'black',size = 10),
  )
after
ggsave("umap/af1.pdf",p1,width = 4,height = 3)
#save(fi_exp,fi_group,file = "Mouse_after_batch.Rdata")
###
#Supplement_Figure 3
load("GSE_matrix.Rdata")
all_gene <- Reduce(intersect,list( rownames(GSE_92592_count),
                                   rownames(GSE134692_count),
                                   rownames(GSE150910_count),
                                   rownames(GSE166036_count),
                                   rownames(GSE52463_count),
                                   rownames(GSE99621_count)))
m_gse92592=GSE_92592_count[all_gene,]
m_gse134692=GSE134692_count[all_gene,]
m_gse150910=GSE150910_count[all_gene,]
m_gse166036=GSE166036_count[all_gene,]
m_gse52463=GSE52463_count[all_gene,]
m_gse99621=GSE99621_count[all_gene,]
raw_exp=cbind(m_gse92592,m_gse134692,m_gse150910,m_gse166036,m_gse52463,m_gse99621)
##
#
raw_group <- read.csv("raw_group.csv",row.names = 1)
coldata <- data.frame(row.names = colnames(raw_exp))
sample <- raw_group$group
coldata$sample <- relevel(factor(sample),"Normal")
coldata$batch <- raw_group$batch
coldata$batch <- as.factor(coldata$batch)
head(coldata)
table(coldata$sample)
##load the merged data
load("IPF_raw.Rdata")
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = raw_exp,
                              colData = coldata,
                              design=~batch+sample)
nrow(dds)
rownames(dds)
dds <- dds[rowSums(counts(dds))>10,]
nrow(dds)     
###
vsd <- vst(dds, blind = FALSE)
library(umap)
library(ggplot2)
ex=as.data.frame(t(assay(vsd)))
ex$group=raw_group$batch1
#ex$group =factor(ex$group,levels= c('Normal','IPF'))
ump <- umap(ex[,1:(ncol(ex)-1)], n_neighbors = 45, random_state = 123)
library(ggsci)
dat <- ump$layout
dat<-data.frame(dat,ex$group)
colnames(dat)=c('UMAP1','UMAP2','group')
p1 <- ggplot(dat,aes(UMAP1,UMAP2)) +
  geom_point(aes(colour=group,fill=group),size=1) +
  scale_color_aaas()+
  theme_classic()+
  theme(axis.line = element_line(size = 1))+
  theme(plot.title = element_text(size = 10),
        axis.text=element_text(colour = 'black',size=10),axis.title.x=element_text(colour = 'black',size=12),
        axis.title.y=element_text(colour = 'black',size=12),legend.text=element_text(color = 'black',size = 10),
  )
p1
ggsave("b2.pdf",p1,width = 4,height = 3)
##
library(limma)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
mat <- assay(vsd)
##
ex1=as.data.frame(t(mat))
ex1$group=raw_group$batch1
ump <- umap(ex1[,1:(ncol(ex1)-1)], n_neighbors = 45, random_state = 123)
library(ggsci)
dat1 <- ump$layout
dat1<-data.frame(dat1,ex1$group)
colnames(dat1)=c('UMAP1','UMAP2','group')
p2 <- ggplot(dat1,aes(UMAP1,UMAP2)) +
  geom_point(aes(colour=group,fill=group),size=1) +
  scale_color_aaas()+
  theme_classic()+
  theme(axis.line = element_line(size = 1))+
  theme(plot.title = element_text(size = 10),
        axis.text=element_text(colour = 'black',size=10),axis.title.x=element_text(colour = 'black',size=12),
        axis.title.y=element_text(colour = 'black',size=12),legend.text=element_text(color = 'black',size = 10),
  )
p2
ggsave("af1.pdf",p2,width = 4,height = 3)
###
#Supplement_Figure 4
##Mice model volcano plot
rm(list = ls())
load("Sup_Figure 3/Mouse_after_batch.Rdata")
group_list=factor(fi_group$time,levels = c('Con','BLM_1d','BLM_2d','BLM_7d','BLM_14d',
                                           'BLM_21d','BLM_28d','BLM_35d'))
### do deg
library(limma)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(fi_exp)
design
contrast.matrix<-makeContrasts(BLM_1d-Con,BLM_2d-Con,BLM_7d-Con,BLM_14d-Con,
                               BLM_21d-Con,BLM_28d-Con,BLM_35d-Con,levels = design)
contrast.matrix
fit <- lmFit(fi_exp,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
##BLM_1d
deg_1d=topTable(fit2, coef = 1, n=Inf)
deg_1d$change = as.factor(ifelse(deg_1d$adj.P.Val < 0.05 & deg_1d$logFC > 0.57 ,"up",
                                 ifelse(deg_1d$logFC < -0.57 & deg_1d$adj.P.Val<0.05 ,"down","not")))
table(deg_1d$change)
##BLM_2d
deg_2d=topTable(fit2, coef = 2, n=Inf)
deg_2d$change = as.factor(ifelse(deg_2d$adj.P.Val < 0.05 & deg_2d$logFC >0.57 ,"up",
                                 ifelse(deg_2d$logFC < -0.57 & deg_2d$adj.P.Val<0.05 ,"down","not")))
table(deg_2d$change)
##BLM_7d
deg_7d=topTable(fit2, coef = 3, n=Inf)
deg_7d$change = as.factor(ifelse(deg_7d$adj.P.Val < 0.05 & deg_7d$logFC > 0.57,"up",
                                 ifelse( deg_7d$adj.P.Val<0.05 &deg_7d$logFC < -0.57 ,"down","not")))
table(deg_7d$change)
##BLM_14d
deg_14d=topTable(fit2, coef = 4, n=Inf)
deg_14d$change = as.factor(ifelse(deg_14d$adj.P.Val < 0.05 & deg_14d$logFC >0.57 ,"up",
                                  ifelse(deg_14d$logFC < -0.57 & deg_14d$adj.P.Val<0.05 ,"down","not")))
table(deg_14d$change)
##BLM_21d
deg_21d=topTable(fit2, coef = 5, n=Inf)
deg_21d$change = as.factor(ifelse(deg_21d$adj.P.Val < 0.05 & deg_21d$logFC >0.57 ,"up",
                                  ifelse(deg_21d$logFC < -0.57 & deg_21d$adj.P.Val<0.05 ,"down","not")))
table(deg_21d$change)
##BLM_28d
deg_28d=topTable(fit2, coef = 6, n=Inf)
deg_28d$change = as.factor(ifelse(deg_28d$adj.P.Val < 0.05 & deg_28d$logFC >0.57 ,"up",
                                  ifelse(deg_28d$logFC < -0.57 & deg_28d$adj.P.Val<0.05 ,"down","not")))
table(deg_28d$change)
##BLM_35d
deg_35d=topTable(fit2, coef = 7, n=Inf)
deg_35d$change = as.factor(ifelse(deg_35d$adj.P.Val < 0.05 & deg_35d$logFC >0.57 ,"up",
                                  ifelse(deg_35d$logFC < -0.57 & deg_35d$adj.P.Val<0.05 ,"down","not")))
table(deg_35d$change)
###IPF patients volcanl plot
load("Sup_Figure 3/IPF_raw.Rdata")
coldata <- data.frame(row.names = colnames(raw_exp))
sample <- raw_group$group
coldata$sample <- relevel(factor(sample),"Normal")
coldata$batch <- raw_group$batch
coldata$batch <- as.factor(coldata$batch)
head(coldata)
table(coldata$sample)
#
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = raw_exp,
                              colData = coldata,
                              design=~batch+sample)
nrow(dds)
rownames(dds)
dds <- dds[rowSums(counts(dds))>10,]
nrow(dds)                  
dds <- DESeq(dds)
resultsNames(dds)
contrast <- c("sample", "IPF", "Normal")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
#
library(ashr)
dd2 <- lfcShrink(dds,contrast=contrast, res=dd1,type="ashr")
plotMA(dd2, ylim=c(-5,5))
#
library(dplyr)
library(tibble)
IPF_deg <- dd2 %>% 
  data.frame() %>% 
  rownames_to_column("gene_id")
head(IPF_deg)
#
IPF_deg$change <- as.factor(ifelse(IPF_deg$padj < 0.05 & IPF_deg$log2FoldChange > 0.57 ,"up",
                                   ifelse(IPF_deg$log2FoldChange < -0.57 & IPF_deg$padj<0.05 ,"down","not")))
table(IPF_deg$change)
### Volcano plot
library(ggrepel)
library(ggplot2)
vol_IPF <- ggplot(res,aes(x=log2FoldChange,y= -log10(padj),color=change))+
  geom_point(data = res[res$padj<0.05&abs(res$log2FoldChange)> 0.57,],size = 1.5)+ 
  geom_point(data = res[res$padj>0.05|abs(res$log2FoldChange)< 0.57,],size = 1.5)+
  scale_color_manual(values=c("#4393C3","#00000033","#FC4E2A"))+
  ylab('-log10 (P.adjust)')+
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
vol_IPF
ggsave("vol_Lung.pdf",vol_IPF,width = 4,height = 3)
save(deg_1d,deg_2d,deg_7d,deg_14d,deg_21d,deg_28d,deg_35d,IPF_deg,file = "Mouse_IPF_GSEA_inpute.Rdata")
### 
##Figure 1
##Mouse_IPF GSEA analysis
rm(list = ls())
load("Figure 3/Mouse_IPF_GSEA_inpute.Rdata")
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
deg_1d$SYMBOL=rownames(deg_1d)
dge_1d <- bitr(unique(deg_1d$SYMBOL), fromType = "SYMBOL", 
               toType = c( "ENTREZID","GENENAME"),
               OrgDb = org.Mm.eg.db)
deg_1d <- merge(deg_1d,dge_1d,by='SYMBOL')
###
library(msigdbr)
msigdbr_species() 
mus_msigdb <- msigdbr(species="Mus musculus")
head(mus_msigdb, 2) %>% as.data.frame
m_df = msigdbr(species = "Mus musculus") %>% dplyr::filter(gs_cat == "H")
m_t2g = m_df %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
###
geneList=deg_1d$logFC
names(geneList)=deg_1d$ENTREZID
geneList=sort(geneList,decreasing = T)
head(geneList)
##
GSEA_1d <- GSEA(
  geneList,
  exponent = 1,
  minGSSize =5,
  maxGSSize = 500,
  eps = 0,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE=m_t2g,
  verbose = F,
  seed = T
)
deg_gs_1d<- GSEA_1d@result
deg_gs_1d_sig <- deg_gs_1d[deg_gs_1d$p.adjust < 0.25,]
#write.table(deg_gs_1d,"gsea_1d.txt",sep = "\t",row.names = F)
#### the data were adjusted in excel
rm(list = ls())
a=read.table("Mouse_BAL_IPF_GSEA.txt",header = T,sep = "\t")
colnames(a)
a$time
library(stringr)
#a$group <- str_split(a$time,'_',simplify = T)[,2]
a$time <- factor(a$time,levels = c('1d','2d','7d','14d','21d',
                                   '28d','35d','Lung','BAL'))#,
a$group <- factor(a$group,levels = c('Mouse','IPF'))
library(ggplot2)
pp = ggplot(a,aes(time,Pathway.Name))+ geom_point(aes(size= PP,color = NES))
pf <- pp+ theme_bw() + scale_colour_gradient(low="blue",high="red")+labs(size="-log10(p.adjust)",x="Time",y="Pathway name")+
  theme(plot.title = element_text(size = 12),
        axis.text=element_text(colour = 'black',size=12),axis.title.x=element_text(colour = 'black',size=12),
        axis.title.y=element_text(colour = 'black',size=12),legend.text=element_text(color = 'black',size = 10),
        legend.title=element_text(colour = 'black',size = 10))+
  facet_grid(cut1~group,scales = "free",space = "free")+
  theme(panel.spacing.x = unit(0,"cm"))+
  theme(panel.spacing.y = unit(0,"cm"))+ylab(NULL)+xlab(NULL)
pf
ggsave(pf, file='merge_Hypoxia_new.pdf', width=7.5, height=8)















