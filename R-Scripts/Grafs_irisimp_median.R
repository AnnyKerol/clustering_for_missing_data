library(ggpubr)
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

                           ### Median Imputation

#####################################################################################


rm(list = ls())
library(e1071)
library(dplyr)
library(tibble)

dat<-read.table("irisNA20.data",  header=F)
dat[dat==9999]<-NA
 
colnames(dat)<-c('Species','SL', 'SW', 'PL', 'PW')
head(dat)
#median imputation
imp<- as.data.frame(impute(dat, what = "median"))
d<-write.table(imp, "grafs_impMediana/irisimp20.data",row.names = FALSE, col.names=FALSE, sep= "\t")
colnames(imp)<-c('Species','SL', 'SW', 'PL', 'PW')
head(imp)



dt1 = dat %>% 
  select(PL, SL) %>% 
  rename(PL_imp = PL) %>% 
  mutate(
    PL_imp = as.logical(ifelse(is.na(PL_imp), "TRUE", "FALSE"))
  ) %>% 
  rename(SL_imp = SL) %>% 
  mutate(
    SL_imp = as.logical(ifelse(is.na(SL_imp), "TRUE", "FALSE"))
  )%>%
  rownames_to_column()


dt2 = imp %>% 
  select(PL, SL) %>% 
  rownames_to_column()

dt = left_join(dt1, dt2)
head(dt)

pdf("grafs_impMediana/imp20_PSL.pdf",width=7, height=7)
vars <- c("SL","PL", "SL_imp","PL_imp")
marginplot(dt[,vars], delimiter="imp",alpha=0.7, pch=c(19), cex = .9, cex.numbers=0.6, col= c("#00AFBB", "#E7B800", "#FC4E07"))
dev.off()


dt1 = dat %>% 
  select(PW, SW) %>% 
  rename(PW_imp = PW) %>% 
  mutate(
    PW_imp = as.logical(ifelse(is.na(PW_imp), "TRUE", "FALSE"))
  ) %>% 
  rename(SW_imp = SW) %>% 
  mutate(
    SW_imp = as.logical(ifelse(is.na(SW_imp), "TRUE", "FALSE"))
  )%>%
  rownames_to_column()


dt2 = imp %>% 
  select(PW, SW) %>% 
  rownames_to_column()

dt = left_join(dt1, dt2)
head(dt)

pdf("grafs_impMediana/imp20_PSWidth.pdf",width=7, height=7)
vars <- c("SW","PW", "SW_imp","PW_imp")
marginplot(dt[,vars], delimiter="imp",alpha=0.7, pch=c(19), cex = .9, cex.numbers=0.6, col= c("#00AFBB", "#E7B800", "#FC4E07"))
dev.off()

