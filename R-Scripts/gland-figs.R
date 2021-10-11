rm(list=ls())
library(lattice)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(naniar)
library(ggpubr)
library(GGally)
library(scatterplot3d) 


#########################################################

# Figures to complete data

#######################################################



gland<-read.table("gland.data", sep=",", header=F)
gland<-as.data.frame(gland)

gland[,1]<-ifelse(gland[,1]==1, "normal", gland[,1])
gland[,1]<-ifelse(gland[,1]==2, "hyper", gland[,1])
gland[,1]<-ifelse(gland[,1]==3, "hypo", gland[,1])

colnames(gland)<-c("Class", "T3", "TTS", "TST", "TSH", "DTSH")
head(gland)
gland$Class<-as.factor(gland$Class)
ggpairs(gland, columns = 2:6, ggplot2::aes(colour=Classes))


#Use box plot as marginal plots
pdf("T3TST.pdf",width=7, height=7)
fig1 <- ggscatterhist(
  gland, x = "T3", y = "TST", xlab = "T3", ylab = "TST",
  color = "Class", size = 5, alpha = 0.7,
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  margin.plot = "boxplot",
  ggtheme = theme_classic(base_size = 17))

ggsave(fig1, filename = "T3TST.pdf",  bg = "transparent", width=7, height=7)


fig2 <- ggscatterhist(
  gland, x = "TST", y = "TTS", xlab = "TST", ylab = "TTS",
  color = "Class", size = 5, alpha = 0.7,
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  margin.plot = "boxplot",
  ggtheme = theme_classic(base_size = 17))

#dev.off()
ggsave(fig2, filename = "TSTTTS.pdf",  bg = "transparent", width=7, height=7)

#<>--------------------------------------
#    Missings mapping 
#<>--------------------------------------


dados5<-read.table("glandNA5.data",  header=F)
dados10<-read.table("glandNA10.data", header = F)
dados15<-read.table("glandNA15.data", header = F)
dados20<-read.table("glandNA20.data", header = F)
dados5[dados5==9999]<-NA
dados10[dados10==9999]<-NA
dados15[dados15==9999]<-NA
dados20[dados20==9999]<-NA

colnames(dados5)<-c("Class", "T3", "TTS", "TST", "TSH", "DTSH")
colnames(dados10)<-c("Class", "T3", "TTS", "TST", "TSH", "DTSH")
colnames(dados15)<-c("Class", "T3", "TTS", "TST", "TSH", "DTSH")
colnames(dados20)<-c("Class", "T3", "TTS", "TST", "TSH", "DTSH")
head(dados5)
head(dados10)
head(dados15)
head(dados20)

ggplot_missing <- function(x){
  
  x %>% 
    is.na %>%
    melt %>%
    ggplot(data = .,
           aes(x = Var2,
               y = Var1)) +
    geom_raster(aes(fill = value)) +
    scale_fill_grey(name = "" ,start=0.85, end=0.2, labels = c("Presente","Missing")) +
    theme_classic(base_size = 15) + 
    theme(axis.text.x  = element_text(angle=59, vjust=0.5, colour="black"),axis.text.y = element_text(angle=0, colour="black")) + 
    labs(x = "",
         y = "Observações")
}


a<-ggplot_missing(dados5[,-1])+ ylim(0,216)
b<-gg_miss_var(dados5[,-1])+labs(y = "Number of missing values",x='Variables')+theme_classic(base_size = 15)+theme(axis.text.x = element_text(angle=0, colour="black"),axis.text.y = element_text(angle=0, colour="black"))
grid.arrange(a,b, nrow=2, ncol=1, heights=c(2,1))
dev.off()



postscript("/gmap10.eps", horizontal=FALSE, paper="special",width=7, height=8)
c<-ggplot_missing(dados10[,-1])+ ylim(0,216)
d<-gg_miss_var(dados10[,-1])+labs(y = "Number of missings values",x='Variables')+theme_classic(base_size = 15)+theme(axis.text.x = element_text(angle=0, colour="black"),axis.text.y = element_text(angle=0, colour="black"))
grid.arrange(c,d, nrow=2, ncol=1, heights=c(2,1))
dev.off()  


postscript("gmap15.eps", horizontal=FALSE, paper="special",width=7, height=8)
e<-ggplot_missing(dados15[,-1])+ ylim(0,215)
f<-gg_miss_var(dados15[,-1])+ labs(y = "Number of missings values",x='Variables')+theme_classic(base_size = 15)+theme(axis.text.x = element_text(angle=0, colour="black"),axis.text.y = element_text(angle=0, colour="black"))
grid.arrange(e,f, nrow=2, ncol=1, top='', heights=c(2,1))
dev.off()  


postscript("gmap20.eps", horizontal=FALSE, paper="special",width=7, height=8)
g<-ggplot_missing(dados20[,-1])+ ylim(0,215)
h<-gg_miss_var(dados20[,-1])+ labs(y = "Number of missings values",x='Variables')+theme_classic(base_size = 15)+theme(axis.text.x = element_text(angle=0, colour="black"),axis.text.y = element_text(angle=0, colour="black"))
grid.arrange(g,h, nrow=2, ncol=1, top='', heights=c(2,1))
dev.off()  

#<>--------------------------------------
# Performance Charts
#<>--------------------------------------

dados<-read.table("Tgmissclassication.txt", header = T)
names(dados)
attach(dados)
dados<-as.data.frame(dados)
dados$met.<-as.factor(dados$met.)
dados$miss<-as.factor(dados$miss)

dados[,3:5]<-dados[,3:5]*100
postscript("gdesempenho.eps", horizontal=FALSE, paper="special",width=12, height=6)
ggplot(dados, aes(x=miss, y=mean, group=met.), col='black') + theme_bw(base_size = 17)+
  theme(axis.text.x = element_text(angle=0, colour="black"),axis.text.y = element_text(angle=0, colour="black"))+
  geom_errorbar(aes(ymin=dados$mean-dados$sd, ymax=dados$mean+dados$sd), width=0.6) +
  geom_point(data=dados) + geom_line(data=dados,linetype="dashed")+
labs(x=" (%) de Missings", y=" Taxa de Erro (%)")+
  facet_grid(. ~ met.)
dev.off()

#<>--------------------------------------

#Comparison between Methods - Accuracy

#<>--------------------------------------

dados<-read.table("gcompara.txt", header = T)
names(dados)
attach(dados)
dados<-as.data.frame(dados)
dados$met.<-as.factor(dados$met.)
dados$miss<-as.factor(dados$miss)

dados$erro<-(1-erro)*100

postscript("/gacuracia.eps", horizontal=FALSE, paper="special",width=11, height=6)
ggplot(dados, aes(x=dados$miss, y=dados$erro, colour=factor(met.), group=met.)) +
  geom_point() + geom_line(size=2) +
  scale_colour_manual(values = c("Média" = "#004586", "Mediana" = "#ff420e", "OCS" = "#ffd320"), name="Métodos", labels=c("Média", "Mediana", "ECO")) +
  scale_y_continuous(breaks=seq(75, 90, 5), name="Acurácia (%)") +
  #expand_limits(y=0) +
  scale_x_discrete(name="(%) Missings") +
  theme_classic(base_size = 17)+
  theme(axis.text.x = element_text(angle=0, hjust=.5, colour="black", size=rel(1.2)), axis.text.y = element_text(hjust=1, colour="black", size=rel(1.2)),
        legend.title = element_text(colour="black", size=rel(1.2), face="bold"), legend.text = element_text(colour="black", size = rel(1.2)), legend.position="top" )#+
  #ggtitle("Assertividade") + theme(plot.title = element_text(lineheight=.8, face="bold"))
#Fechar arquivo
dev.off()


#<>------------------------------------------------------------------------------------------------------

##                            Principal component analysis

#<>-------------------------------------------------------------------------------------------------------

gland<-read.table("gland.data", sep=",", header = F)
colnames(gland)<-c("Classes", "T3", "TTS", "TST", "TSH", "DTSH")
attach(gland)

gland[,1]<-ifelse(gland[,1]==1, "normal", gland[,1])
gland[,1]<-ifelse(gland[,1]==2, "hyper", gland[,1])
gland[,1]<-ifelse(gland[,1]==3, "hypo", gland[,1])

dados.pca <- gland[-1] 
pca.obj <- prcomp(dados.pca)

postscript("/glandPCA.eps", horizontal=FALSE, paper="special",width=11, height=6)
dtp <- data.frame('Classes' = gland$Classes, pca.obj$x[,1:2]) # the first two componets are selected (NB: you can also select 3 for 3D plottings or 3+)
ggplot(data = dtp) + 
  geom_point(size=5,aes(x = PC1, y = PC2, col =Classes)) + 
  theme_classic(base_size = 16)+
theme(axis.text.x = element_text(angle=0, hjust=.5, colour="black", size=rel(1.2)), axis.text.y = element_text(hjust=1, colour="black", size=rel(1.2)),
      legend.title = element_text(colour="black", size=rel(1.2), face="bold"), legend.text = element_text(colour="black", size = rel(1.2)), legend.position="top" )



dev.off()


