
rm(list=ls())
library(lattice)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(naniar)
library(ggpubr)
library(scatterplot3d) 

#<>--------------------------------------

# Figures - Complete Data

#<>--------------------------------------


head(iris)

pdf("3dIris.pdf",width=7, height=7)
shapes = c(16, 17, 18) 
shapes <- shapes[as.numeric(iris$Species)]
colors <- c("red", "blue", "green")
colors <- colors[as.numeric(iris$Species)]
scatterplot3d(iris[,1:3], pch = shapes, color=colors)
legend( 'right' ,legend = levels(iris$Species),
        col = c("red", "blue", "green") ,box.lty = 0, pch = c(16,17,18), inset = 0.08, cex=0.7)
dev.off()

pdf("PSlength.pdf",width=7, height=7)
ggplot(iris) + geom_point(aes(x=Sepal.Length, y=Petal.Length)) + geom_point(size=4, aes(x=Sepal.Length, y=Petal.Length, col=Species, shape=Species))  +
  theme_classic()
dev.off()

names(iris)
pdf("PSwidth.pdf",width=7, height=7)
ggplot(iris) + geom_point(aes(x=Sepal.Width, y=Petal.Width)) + geom_point(size=4, aes(x=Sepal.Width, y=Petal.Width, col=Species, shape=Species))  +
  theme_classic()
dev.off()


box <- ggplot(data=iris, aes(x=Species, group=Species))
box + geom_boxplot(aes(fill=Species)) + 
  ylab("Sepal Length") + ggtitle("Iris Boxplot") +
  stat_summary(fun.y=mean, geom="point", shape=5, size=4) 


colnames(iris) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width", "Species")

fig2 <- ggscatterhist(
  iris, x = "Sepal.Length", y = "Petal.Width", xlab = "SL", ylab = "PL",
  color = "Species", size = 5, alpha = 0.7,
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  margin.plot = "boxplot",
  ggtheme = theme_classic(base_size = 17))
#dev.off()

ggsave(fig2, filename = "PSlength1.pdf",  bg = "transparent", width=7, height=7)

#pdf("boxPSwidth.pdf",width=7, height=7)

fig1 <- ggscatterhist(
  iris, x = "Sepal.Width", y = "Petal.Width", xlab = "SW", ylab = "PW",
  color = "Species", size = 5, alpha = 0.7,
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  margin.plot = "boxplot",
  ggtheme = theme_classic(base_size = 17))
#dev.off()
ggsave(fig1, filename = "PSwidth1.pdf",  bg = "transparent", width=7, height=7)


#<>--------------------------------------
#    Missings mapping 
#<>--------------------------------------

dados5<-read.table("irisNA5.data",  header=F)
dados10<-read.table("irisNA10.data", header = F)
dados15<-read.table("irisNA15.data", header = F)
dados20<-read.table("irisNA20.data", header = F)
dados5[dados5==9999]<-NA
dados10[dados10==9999]<-NA
dados15[dados15==9999]<-NA
dados20[dados20==9999]<-NA

colnames(dados5)<-c('Species','SL', 'SW', 'PL', 'PW')
colnames(dados10)<-c('Species','SL', 'SW', 'PL', 'PW')
colnames(dados15)<-c('Species','SL', 'SW', 'PL', 'PW')
colnames(dados20)<-c('Species','SL', 'SW', 'PL', 'PW')
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
    scale_fill_grey(name = "" ,start=0.85, end=0.2, labels = c("Observed","Missing")) +
    theme_classic(base_size = 15) + 
    theme(axis.text.x  = element_text(angle=59, vjust=0.5, colour="black"),axis.text.y = element_text(angle=0, colour="black")) + 
    labs(x = "",
         y = "Observations")
}


#pdf("mapping5.pdf",width=6, height=8)
postscript("mapping5.eps", horizontal=FALSE, paper="special",width=7, height=8)
a<-ggplot_missing(dados5[,-1])+ ylim(0,151)
b<-gg_miss_var(dados5[,-1]) +labs(y = "Número de Missings",x='Variáveis')+
theme_classic(base_size = 15)+theme(axis.text.x = element_text(angle=0, colour="black"),axis.text.y = element_text(angle=0, colour="black"))
grid.arrange(a,b, nrow=2, ncol=1, heights=c(2,1))
dev.off()

postscript("mapping10.eps", horizontal=FALSE, paper="special",width=7, height=8)
c<-ggplot_missing(dados10[,-1]) + 
  ylim(0,151)
d<-gg_miss_var(dados10[,-1])+labs(y = "Number of missing values",x='Variables') + 
  theme_classic(base_size = 15) + 
  theme(axis.text.x = element_text(angle=0, colour="black"),axis.text.y = element_text(angle=0, colour="black"))
grid.arrange(c,d, nrow=2, ncol=1, heights=c(2,1))

dev.off()  


postscript("mapping15.eps", horizontal=FALSE, paper="special",width=7, height=8)
e<-ggplot_missing(dados15[,-1])+ ylim(0,151)
f<-gg_miss_var(dados15[,-1])+ labs(y = "Number of missing values",x='Variables')+theme_classic(base_size = 15)+theme(axis.text.x = element_text(angle=0, colour="black"),axis.text.y = element_text(angle=0, colour="black"))
grid.arrange(e,f, nrow=2, ncol=1, top='', heights=c(2,1))
dev.off()  


postscript("mapping20.eps", horizontal=FALSE, paper="special",width=7, height=8)
g<-ggplot_missing(dados20[,-1])+ ylim(0,151)
h<-gg_miss_var(dados20[,-1])+ labs(y = "Number of missing values",x='Variables')+theme_classic(base_size = 15)+theme(axis.text.x = element_text(angle=0, colour="black"),axis.text.y = element_text(angle=0, colour="black"))
grid.arrange(g,h, nrow=2, ncol=1, top='', heights=c(2,1))
dev.off()  



#<>--------------------------------------
# Performance Charts
#<>--------------------------------------

dados<-read.table("missclassication.txt", header = T)
names(dados)
attach(dados)
dados<-as.data.frame(dados)
dados$met.<-as.factor(dados$met.)
dados$miss<-as.factor(dados$miss)

dados[,3:5]<-dados[,3:5]*100
#pdf("desempenho.pdf",width=10, height=6)
postscript("desempenho.eps", horizontal=FALSE, paper="special",width=12, height=6)
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


dados<-read.table("Iriscompara.txt", header = T)
names(dados)
attach(dados)
dados<-as.data.frame(dados)
dados$met.<-as.factor(dados$met.)
dados$miss<-as.factor(dados$miss)

dados$erro<-(1-erro)*100

postscript("Irisacuracy.eps", horizontal=FALSE, paper="special",width=11, height=6)
ggplot(dados, aes(x=dados$miss, y=dados$erro, colour=factor(met.), group=met.)) +
  geom_point() + geom_line(size=2) +
  scale_colour_manual(values = c("Média" = "#004586", "Mediana" = "#ff420e", "OCS" = "#ffd320"), name="Métodos", labels=c("Média", "Mediana", "ECO")) +
  scale_y_continuous(breaks=seq(75, 90, 5), name="Acurácia (%)") +
  #expand_limits(y=0) +
  scale_x_discrete(name="(%) Missings") +
  theme_classic(base_size = 17)+
  theme(axis.text.x = element_text(angle=0, hjust=.5, colour="black", size=rel(1.2)), axis.text.y = element_text(hjust=1, colour="black", size=rel(1.2)),
        legend.title = element_text(colour="black", size=rel(1.2), face="bold"), legend.text = element_text(colour="black", size = rel(1.2)), legend.position="top") #+
  #ggtitle("Assertividade") + theme(plot.title = element_text(lineheight=.8, face="bold"))
#Fechar arquivo
dev.off()

#<>------------------------------------------------------------------------------------------------------

##                            Principal component analysis

#<>-------------------------------------------------------------------------------------------------------

iris.pca <- iris[c(1, 2, 3, 4)] 
pca.obj <- prcomp(iris.pca)

postscript("IriPCA.eps", horizontal=FALSE, paper="special",width=11, height=6)
dtp <- data.frame('Species' = iris$Species, pca.obj$x[,1:2]) # the first two componets are selected (NB: you can also select 3 for 3D plottings or 3+)
ggplot(data = dtp, groups=Species) + 
  geom_point(size=5, aes(x = PC1, y = PC2, col =Species)) + 
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=0, hjust=.5, colour="black", size=rel(1.2)), axis.text.y = element_text(hjust=1, colour="black", size=rel(1.2)),
        legend.title = element_text(colour="black", size=rel(1.2), face="bold"), legend.text = element_text(colour="black", size = rel(1.2)), legend.position="top" )
dev.off()




