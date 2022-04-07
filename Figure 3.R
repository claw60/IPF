##Figure 3
rm(list = ls())
library(GEOquery)
library(survival)
library(GSVA)
#gset = getGEO('GSE70867',destdir = ".",AnnotGPL = F , getGPL = F) 
load("GSE70867_eset.Rdata")
b =gset[[2]]
phe=pData(b)
rownames(phe)
colnames(phe)
phe<-phe[23:154,c(1,43:47)]
### get expression matrix
raw_exprSet=exprs(b) 
b
### delete healthy samples
dim(raw_exprSet)
library(limma)
raw_exprSet=normalizeBetweenArrays(raw_exprSet)
boxplot(raw_exprSet)
dim(raw_exprSet)
### get expression matrix
c = gset[[3]]
c
raw_exprSet_c=exprs(c) 
library(limma)
raw_exprSet_c=normalizeBetweenArrays(raw_exprSet_c)
boxplot(raw_exprSet_c,las=2)
### get clinical information
phe_c=pData(c)
colnames(phe_c)
phe_c<-phe_c[,c(1,42:46)]
colnames(phe_c)
### merge two clinical information
phe_final<-rbind(phe,phe_c)
#write.csv(phe_final,file = "phe_final.csv")
##
d = gset[[4]]
phe_d=pData(d)
colnames(phe_d)
phe_d<-phe_d[,1:17]
raw_exprSet_d=exprs(d) 
library(limma)
raw_exprSet_d=normalizeBetweenArrays(raw_exprSet_d)
boxplot(raw_exprSet_d,las=2)
###### Probe ID transfer ####
b_anno <-data.table::fread("GPL14550-9757.txt",skip ="ID")
#
ids_b = b_anno[,c('ID','GENE_SYMBOL')]
ids_b[ids_b$GENE_SYMBOL == '',]= NA
ids_b <- na.omit(ids_b)
dim(ids_b)
library(tibble)
raw_exp_b=as.data.frame(raw_exprSet) 
colnames(ids_b) = c('ID','symbol')
table(rownames(raw_exprSet) %in% ids_b$ID)
dim(raw_exprSet)
exprSet = raw_exprSet[rownames(raw_exprSet) %in% ids_b$ID,]
dim(exprSet)
ids_b=ids_b[match(rownames(exprSet),ids_b$ID),]
table(duplicated(ids_b$symbol))
ID <- rownames(exprSet)
eset_b <- exprSet[ID %in% ids_b$ID,] %>% cbind(ids_b) 
colnames(eset_b)
eset_b<-aggregate(x=eset_b[,1:154],by=list(eset_b$symbol),FUN=mean,na.rm=T)
####
c_anno <-data.table::fread("GPL17077-17467.txt",skip ="ID")
ids_c = c_anno[,c('ID','GENE_SYMBOL')]
ids_c[ids_c$GENE_SYMBOL == '',]= NA
ids_c <- na.omit(ids_c)
exprSet_c=as.data.frame(raw_exprSet_c) 
colnames(ids_c) = c('ID','symbol')
table(rownames(exprSet_c) %in% ids_c$ID)
dim(exprSet_c)
raw_exp_c = exprSet_c[rownames(exprSet_c) %in% ids_c$ID,]
dim(raw_exp_c)
ids_cc=ids_c[match(rownames(raw_exp_c),ids_c$ID),]
table(duplicated(ids_cc$symbol))
ID <- rownames(raw_exp_c)
eset_c <- raw_exp_c[ID %in% ids_cc$ID,] %>% cbind(ids_cc) 
eset_c<-aggregate(x=eset_c[,1:64],by=list(eset_c$symbol),FUN=mean,na.rm=T)
###
d_anno <-data.table::fread("GPL570-55999.txt",skip ="ID")
colnames(d_anno)
ids_d = d_anno[,c('ID','Gene Symbol')]
ids_d[ids_d$`Gene Symbol` == '',]= NA
ids_d <- na.omit(ids_d)
library(tibble)
colnames(ids_d) = c('ID','symbol')
ids_d <- ids_d[-c(grep("///",ids_d$symbol)),]
exprSet_d = raw_exprSet_d[rownames(raw_exprSet_d) %in% ids_d$ID,]
dim(exprSet_d)
ids_d=ids_d[match(rownames(exprSet_d),ids_d$ID),]
dim(ids_d)
ID <- rownames(exprSet_d)
eset_d <- exprSet_d[ID %in% ids_d$ID,] %>% cbind(ids_d)
colnames(eset_d)
eset_d <-aggregate(x=eset_d[,1:57],by=list(eset_d$symbol),FUN=mean,na.rm=T)
##
merge_gene <- Reduce(intersect,list( eset_b$Group.1,
                                     eset_c$Group.1,
                                     eset_d$Group.1))

rownames(eset_b) <- eset_b$Group.1;eset_b <- eset_b[,-1]
colnames(eset_b)
eset_bb <- eset_b[,23:154]
rownames(eset_c) <- eset_c$Group.1;eset_c <- eset_c[,-1]
rownames(eset_d) <- eset_d$Group.1;eset_d <- eset_d[,-1]
#
mer_eset_b <- eset_bb[merge_gene,]
mer_eset_c <- eset_c[merge_gene,]
mer_eset_d <- eset_d[merge_gene,]
#
raw_merge_exp <- cbind(mer_eset_b,mer_eset_c,mer_eset_d)
colnames(raw_merge_exp)
boxplot(raw_merge_exp,las=2)
batch <- c(rep(1,82),rep(2,50),rep(3,64),rep(4,57))
### remove batch
library(sva)
norm_exp <- ComBat(raw_merge_exp,batch = batch)
boxplot(norm_exp,las=2)
save(norm_exp,phe_final,file = 'GSE70866_HDIPF_rmBatch_matrix.Rdata')
###GSVA Score caculate
library(dplyr)
#install.packages("Hmisc")
library(Hmisc)
hal_gs <- lapply(readLines("h.all.v7.4.symbols.gmt"), function(x){
  #x=readLines("gmt/h.all.v7.2.symbols.gmt")[[1]]
  y=strsplit(x,'\t')[[1]]
  y=y[3:length(y)]
  return(y)
})
##name
library(stringr)
tmp <- unlist(lapply(readLines("h.all.v7.4.symbols.gmt"), function(x){
  y=strsplit(x,'\t')[[1]][1]
  #y=str_split(y,'_',simplify = T)[2:4]
  y = unlist(strsplit(y, split="_", fixed=T))[-1]
  ### 3.¿Õ¸ñÁ¬½Ó
  y = paste(y, collapse="_")
})
)
library(GSVA)
norm_exp[1:4,1:4]
gsva_score <-  gsva(as.matrix(norm_exp), hal_gs,verbose=TRUE,mx.dif=1)
rownames(gsva_score) <- tmp
gsva_melt <- as.data.frame(t(gsva_score[,colnames(gsva_score) %in% rownames(phe_final)]))

HIF_list<-list(c('CDKN3','NDRG1','PGAM1','TUBB6','TPI1','MIF','ACOT7',
                 'ALDOA','LDHA','ENO1','ADM','MRPS17','SLC2A1','VEGFA','P4HA1'))
HIF_score <-  gsva(as.matrix(norm_exp), HIF_list,verbose=TRUE,mx.dif=1)
HIF_1 <-  as.data.frame(HIF_score[,colnames(HIF_score) %in% rownames(phe_final)])
colnames(HIF_1) <- c('HIF_score')
phe_final_IPF <- phe_final[-c(1:20),]
##
library(stringr)
tmp <- str_split(phe_final$title,'_',simplify = T)[,2:3]
tmp <- apply(tmp,1,function(x) paste(x,collapse = '_'))
tmp
row.names(phe_final) == row.names(gsva_melt)
gsva_melt$group <- tmp
rownames(HIF_1) == row.names(phe_final)
GSVA_hallmark <- cbind(gsva_melt,HIF_1)
colnames(GSVA_hallmark)
HIF_emt_Gly <- GSVA_hallmark[,c(52,30,35,51)]
colnames(HIF_emt_Gly)[2:3] <- c("EMT_Score","Glycolysis_Score")
write.csv(HIF_emt_Gly,file = "HIF_emt_Gly.csv",row.names = T)
##
Sur_Score <- HIF_emt_Gly[rownames(HIF_emt_Gly) %in% rownames(phe_final_IPF),]
colnames(phe_final_IPF)
phe_sur <- phe_final_IPF[,c(3,5,6)]
rownames(Sur_Score) == rownames(phe_sur) 
phe_Score <- cbind(Sur_Score,phe_sur)
colnames(phe_Score)
phe_Score <- phe_Score[,c(4,6:7,1:3,5)]
colnames(phe_Score)[2:3] <- c("OS","OS.time")
colnames(phe_Score)[7] <- c('GAP')
colnames(phe_Score)
### survival analysis
library(survival)
data= phe_Score
data[,2:7] <- lapply(data[,c(2:7)],as.numeric)
genes <- colnames(phe_Score)[-c(1:3)]

res <- data.frame()
for (i in 1:length(genes)) {
  print(i)
  surv =as.formula(paste('Surv(OS.time, OS)~', genes[i]))
  x = coxph(surv, data = data)
  x = summary(x)
  p.value=signif(x$wald["pvalue"], digits=2)
  HR =signif(x$coef[2], digits=2);#exp(beta)
  HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper = signif(x$conf.int[,"upper .95"],3)
  
  res[i,1] = genes[i]
  res[i,2] = HR
  res[i,3] = HR.confint.lower 
  res[i,4] = HR.confint.upper
  res[i,5] = p.value
}
names(res) <- c("ID","HR","Low.95.CI","High.95.CI","P.value")
res1 <- res[res$P.value < 0.05,]
save(res,data,file = "survival_cox_EMT.Rdata")
##
library(survival)
library(survminer)
res.cut <- surv_cutpoint(data, time = "OS.time", 
                         event = "OS", 
                         variables = names(data)[4:ncol(data)], 
                         minprop = F) 
res.cat <- surv_categorize(res.cut)    
###surv_cutpoint method choose cutoff
genes <- colnames(data)[-c(1:3)]
res2 <- data.frame()
for (i in 1:length(genes)) {
  print(i)
  surv =as.formula(paste('Surv(OS.time, OS)~', genes[i]))
  x = survdiff(surv, data = res.cat)
  pValue=1-pchisq(x$chisq,df=1)
  res2[i,1] = genes[i]
  res2[i,2] = pValue
}
names(res2) <- c("ID","pValue_log(cut_off)")

library(dplyr)
res_dat <- inner_join(res,res2,by='ID')
save(res_dat,file = "suivival_EMT.Rdata")
write.csv(res_dat,file = "sur_res_EMT.csv",row.names = F)
#
res$ID <- as.character(res$ID)
index <- res$ID
sur_data <- data[,index]
sur_data <- cbind(data[,1:3],sur_data)

library(tidyverse)
sur_data <- sur_data %>%
  select(-group)
res.cut <- surv_cutpoint(sur_data, time = "OS.time",
                         event = "OS",
                         variables = names(sur_data)[3:ncol(sur_data)],
                         minprop = F)
res.cat <- surv_categorize(res.cut)
###
genes <- colnames(res.cat)[-c(1:2)]
your.surv <- Surv(res.cat$OS.time, res.cat$OS) 
your.km.plot <- function(genes,data){
  print(genes)
  group <- res.cat[,genes] 
  survival_dat <- data.frame(group = group)
  group <- factor(group, levels = c("low", "high")) 
  fit <- survfit(your.surv ~ group)
  # log-rank test£ºpvalue
  # This function implements the G-rho family of Harrington and Fleming (1982), with weights on each death of S(t)^rho, where S is the Kaplan-Meier estimate of survival. 
  # With rho = 0 this is the log-rank or Mantel-Haenszel test, and with rho = 1 it is equivalent to the Peto & Peto modification of the Gehan-Wilcoxon test.
  sdf <- survdiff(your.surv ~ group,rho=0)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
  p.val
  photo2 <-  ggsurvplot(fit,data = survival_dat, 
                        legend.title = genes,
                        legend.labs = c("low","high"), 
                        #legend = "top",
                        #pval = T, 
                        #pval.method = TRUE,
                        #conf.int = TRUE,
                        risk.table = T, 
                        #risk.table.col = "strata", 
                        risk.table.y.text = F,
                        #linetype = "strata", 
                        #surv.median.line = "hv", 
                        xlab = "Time to death(days)", 
                        ylab="Survival",
                        #xlim = c(0,max(res.cat$OS.time)), 
                        #break.time.by = 10, 
                        size = 1.5, 
                        break.x.by = 500,
                        #ggtheme = theme_bw(), 
                        palette = c("red", "blue"),
                        font.main = 18,     
                        font.x = c("bold",16),        
                        font.y = c("bold",16),        
                        font.tickslab = 14,legend=c(0.15,0.2)
                        
  ) 
  photo2  
} 
###Figure 3 A/B
your.km.plot("Glycolysis_Score",data = res.cat)
####multivariate proportional hazards models 
##Figure 3C 
library(survminer)
single_line<-Surv(time = data$OS.time,event = data$OS)
data$single_line<- with(data,single_line)
##
Gcox<-coxph(single_line~Glycolysis_Score+EMT_Score+GAP,
            data = data)
x <- summary(Gcox)
x
mul_cox <- ggforest(Gcox,data = data,main = 'Hazard ratio',
                    fontsize =1.2,refLabel = 'reference',noDigits = 2)
mul_cox
ggsave("multi_Gly_GAP.pdf",mul_cox,height = 3.5,width = 10.5)
###Figure 3 DE was done in prism


