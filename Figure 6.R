###Figure 6A¡¢B¡¢C
rm(list = ls())
load("Figure 6/Figure6_input.Rdata")
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
