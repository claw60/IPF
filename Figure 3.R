###Figure 3 A
rm(list = ls())
load("IPF_HIF_EMT_Gly.Rdata")
library(ggplot2)
library(ggforce)
library(ggsci)
colnames(IPF_HIF_EMT_Gly)
pf <-  ggplot(IPF_HIF_EMT_Gly,aes(x=group,y=HIF_Score))+
  geom_violin(aes(fill=group,color=group),trim = F,show.legend = F)+
  geom_boxplot(aes(x=group,y=HIF_Score),outlier.shape = NA,fill = "#a6a7ac",color = "#a6a7ac",width = 0.08)+
  scale_fill_manual(values=c("#d1d2d2","#E64B35B2"))+
  scale_color_manual(values=c("#d1d2d2","#E64B35B2"))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line())+
  theme(plot.title = element_text(size = 10),
        axis.text=element_text(colour = 'black',size=10),axis.title.x=element_text(colour = 'black',size=12),
        axis.title.y=element_text(colour = 'black',size=12),legend.text=element_text(color = 'black',size = 10),
        legend.title=element_text(colour = 'black',size = 10))+
  scale_y_continuous(limits = c(-1.2,1.2))+
  
  labs(x=NULL,
       y="HIF Score")
pf
ggsave("HIF_Score_Lung_tissue.pdf",pf,height = 2.3,width = 2.5)
##Figure 3 F/I
library(ggpubr)
colnames(IPF_HIF_EMT_Gly)
p <- ggplot(data = IPF_HIF_EMT_Gly, aes(x = HIF_Score, y = Glycolysis_Score)) + 
  geom_point() + geom_smooth(method = lm,color = "black",se=F) +
  xlab("HIF Score")+
  ylab("Glycolysis Score")+
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson')+
  theme(axis.line = element_line(size = 1))+
  theme(plot.title = element_text(size = 12),
        axis.text=element_text(colour = 'black',size=12),axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12))
p
cor.test(as.numeric(IPF_HIF_EMT_Gly$EMT_Score),as.numeric(IPF_HIF_EMT_Gly$HIF_Score))

######Figure 3 B
load("BAL_HIF_EMT_Gly.Rdata")
cor_HIF_EMT_Gly$group1 <- c(rep("Healthy",20),rep("IPF",176))
library(ggplot2)
library(ggforce)
library(ggsci)
pf <-  ggplot(cor_HIF_EMT_Gly,aes(x=group1,y=HIF_score))+
  geom_violin(aes(fill=group1,color=group1),trim = F,show.legend = F)+
  geom_boxplot(aes(x=group1,y=HIF_score),outlier.shape = NA,fill = "#a6a7ac",color = "#a6a7ac",width = 0.08)+
  scale_fill_manual(values=c("#d1d2d2","#E64B35B2"))+
  scale_color_manual(values=c("#d1d2d2","#E64B35B2"))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line())+
  theme(plot.title = element_text(size = 10),
        axis.text=element_text(colour = 'black',size=10),axis.title.x=element_text(colour = 'black',size=12),
        axis.title.y=element_text(colour = 'black',size=12),legend.text=element_text(color = 'black',size = 10),
        legend.title=element_text(colour = 'black',size = 10))+
  scale_y_continuous(limits = c(-1.2,1.2))+
  
  labs(x=NULL,
       y="HIF Score")
pf+ stat_compare_means()
ggsave("HIF_Score_BAL.pdf",pf,height = 2.8,width = 3)
##Figure 3 G/J
library(ggpubr)
colnames(cor_HIF_EMT_Gly)
p <- ggplot(data = cor_HIF_EMT_Gly, aes(x = HIF_score, y = Glycolysis_Score)) + 
  geom_point() + geom_smooth(method = lm,color = "black",se=F) +
  xlab("HIF Score")+
  ylab("Glycolysis Score")+
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson')+
  theme(axis.line = element_line(size = 1))+
  theme(plot.title = element_text(size = 12),
        axis.text=element_text(colour = 'black',size=12),axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12))
p
cor.test(cor_HIF_EMT_Gly$Glycolysis_Score,cor_HIF_EMT_Gly$HIF_score)

###Figure 3 C
load("Mouse_HIF_EMT_Gly.Rdata")
library(ggplot2)
library(ggforce)
library(ggsci)
pf <-  ggplot(Mouse_HIF_EMT_Gly,aes(x=group,y=HIF_Score))+
  geom_violin(aes(fill=group,color=group),trim = F,show.legend = F)+
  geom_boxplot(aes(x=group,y=HIF_Score),outlier.shape = NA,fill = "#a6a7ac",color = "#a6a7ac",width = 0.08)+
  scale_fill_manual(values=c("#d1d2d2","#E64B35B2","#fbd3b9",
                             "goldenrod1",
                             "#a1c9e5",
                             "#4DBBD5B2",
                             "#00A087B2",
                             "#417bb9"
  ))+
  scale_color_manual(values=c("#d1d2d2","#E64B35B2","#fbd3b9",
                              "goldenrod1",
                              "#a1c9e5",
                              "#4DBBD5B2",
                              "#00A087B2",
                              "#417bb9"))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line())+
  theme(plot.title = element_text(size = 10),
        axis.text=element_text(colour = 'black',size=10),axis.title.x=element_text(colour = 'black',size=12),
        axis.title.y=element_text(colour = 'black',size=12),legend.text=element_text(color = 'black',size = 10),
        legend.title=element_text(colour = 'black',size = 10))+
  scale_y_continuous(limits = c(-1.2,1.2))+
  
  labs(x=NULL,
       y="HIF Score")
pf
ggsave("HIF_Score_Mouse.pdf",pf,height = 2.8,width = 4)
##Figure 3 H/K
colnames(Mouse_HIF_EMT_Gly)
p <- ggplot(data = Mouse_HIF_EMT_Gly, aes(x = HIF_Score, y = Glycolysis_Score)) + 
  geom_point() + geom_smooth(method = lm,color = "black",se=F) +
  xlab("HIF Score")+
  ylab("Glycolysis Score")+
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson')+
  theme(axis.line = element_line(size = 1))+
  theme(plot.title = element_text(size = 12),
        axis.text=element_text(colour = 'black',size=12),axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12))
p
p1 <- p+scale_y_continuous(limits = c(-0.4,0.5))+scale_x_continuous(limits = c(-0.73,0.65))
p1
ggsave(p1,file="Mouse_HIF_Gly_1.pdf",width = 4.5,height = 3.5)
