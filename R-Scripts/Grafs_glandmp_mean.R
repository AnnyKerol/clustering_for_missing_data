library(ggpubr)
library(VIM)
library(tibble)

# Grouped Scatter plot with marginal density plots
ggscatterhist(
  iris, x = "Sepal.Length", y = "Sepal.Width",
  color = "Species", size = 3, alpha = 0.6,
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  margin.params = list(fill = "Species", color = "black", size = 0.2)
)

# Use box plot as marginal plots
ggscatterhist(
  iris, x = "Sepal.Length", y = "Petal.Width", xlab = "SL", ylab = "PL",
  color = "Species", size = 3, alpha = 0.6,
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  margin.plot = "boxplot",
  ggtheme = theme_classic()
)

ggscatterhist(
  iris, x = "Sepal.Width", y = "Petal.Width", xlab = "SL", ylab = "PL",
  color = "Species", size = 3, alpha = 0.6,
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  margin.plot = "boxplot",
  ggtheme = theme_classic()
)


#######################################################################################

                           ### Imputação da Média -Thyroid Gland

#####################################################################################



#install.packages('e1071')
rm(list = ls())
library(e1071)
library(ggplot2)
dat<-read.table("Thyroid gland/Gmedia/glandNA15.data",  header=F)
dat[dat==9999]<-NA
 
colnames(dat)<-c("Classes", "T3", "TTS", "TST", "TSH", "DTSH")
head(dat) 
# mean imputation
imp<- as.data.frame(impute(dat, what = "mean"))
d<-write.table(imp, "Gmedia/glandMNA15.data",row.names = FALSE, col.names=FALSE, sep= "\t")
colnames(imp)<-c("Classes", "T3", "TTS", "TST", "TSH", "DTSH")
head(imp)

library(dplyr)

dt1 = dat %>% 
  select(TST, T3) %>% 
  rename(T3_imp = T3) %>% 
  mutate(
    T3_imp = as.logical(ifelse(is.na(T3_imp), "TRUE", "FALSE"))
  ) %>% 
  rename(TST_imp = TST) %>% 
  mutate(
    TST_imp = as.logical(ifelse(is.na(TST_imp), "TRUE", "FALSE"))
  )%>%
  rownames_to_column()

dt2 = imp %>% 
  select(T3, TST) %>% 
  rownames_to_column()

dt = left_join(dt2,dt1)
head(dt)




pdf("impM15_T3TST.pdf",width=7, height=7)
vars <- c("T3", "TST", "TST_imp", "T3_imp")
marginplot(dt[,vars], delimiter="imp",alpha=0.9, pch=c(19), cex = 1.5, cex.numbers=0.9, col= c("#00AFBB", "#E7B800", "#FC4E07"))
#ggsave("a.eps", width = 4, height = 4)
dev.off()


dt1 = dat %>% 
  select(TST, TTS) %>% 
  rename(TST_imp = TST) %>% 
  mutate(
    TST_imp = as.logical(ifelse(is.na(TST_imp), "TRUE", "FALSE"))
  ) %>% 
  rename(TTS_imp = TTS) %>% 
  mutate(
    TTS_imp = as.logical(ifelse(is.na(TTS_imp), "TRUE", "FALSE"))
  )%>%
  rownames_to_column()


dt2 = imp %>% 
  select(TTS, TST) %>% 
  rownames_to_column()

dt = left_join(dt1, dt2)
head(dt)

pdf("Thyroid gland/Gmedia/impM20_TTSTST.pdf",width=7, height=7)
vars <- c("TST_imp", "TTS_imp",  "TTS", "TST")
marginplot(dt[,vars], delimiter="imp",alpha=0.7, pch=c(19), cex = .9, cex.numbers=0.6, col= c("#00AFBB", "#E7B800", "#FC4E07"))
dev.off()

