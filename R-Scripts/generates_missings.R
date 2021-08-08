
##################################################################

#  Generating missings from a Bernoulli distribution (0, p)

#################################################################

rm(list=ls())

#m= sample numbers
#p=probability of success 

bernoulli.1 <- function(m,p){
  set.seed(13) # iris 
  x <- rep(NA,m)
  for(i in 1:m){
    k <- 0
    if(runif(1)<=p) k <- 1
    x[i] <- k
  }
  return(x)
}

# Reading dataset iris
dados<-read.table("iris.data")
dados1<-dados[-1] 
dados1<-as.matrix(dados1)
dados1<-c(dados1) # converting dataframe in vector
length(dados1)
x<-length(dados1)

# generating the positions that will be missing values
posicao<-bernoulli.1(x,0.05)

# When posicao == 1 the values is replacing with NA

for(i in 1:length(dados1)){
  for(j in 1:length(posicao)){
    if(i==j && posicao4[j]==1)
      dados1[i]=9999
  }
}

dados1<-matrix(dados1, nrow=150,ncol=4)
dados1<-data.frame(dados1)
dados1<-cbind(dados[1],dados1)
dados1<-as.data.frame(dados1)
write.table(dados1, "irisNA5.data",row.names = FALSE, col.names=FALSE, sep= "\t")

# Function to estimate sigma
sigma2est = function(x, frac = .5) {
  x = as.matrix(x)
  n = nrow(x)
  m = floor(n*frac)
  idx1 = sample(1:n, m, replace = T)
  idx2 = sample(1:n, m, replace = T)
  tmp = (x[idx1,, drop = F] - x[idx2,, drop = F])^2
  dist = apply(tmp, 1, sum)
  mean(quantile(dist[dist != 0], probs = c(.9, .1)))
}


#####################################################

# Generate NA for dataset Thyroid Gland

####################################################

bernoulli.1 <- function(m,p){
  set.seed(347) # gland
  x <- rep(NA,m)
  for(i in 1:m){
    k <- 0
    if(runif(1)<=p) k <- 1
    x[i] <- k
  }
  return(x)
}



dados<-read.table("gland.data", sep = ',')
attach(dados)
dim(dados)
dados1<-dados[-1] 
dados1

# estimating sigma 
set.seed(347)
sigma2est(dados1[,1])
sigma2est(dados1[,2])
sigma2est(dados1[,3])
sigma2est(dados1[,4])
sigma2est(dados1[,5])


dados1<-as.matrix(dados1)
dados1<-c(dados1) 
length(dados1)
x<-length(dados1)

posicao4<-bernoulli.1(x,0.20) 

for(i in 1:length(dados1)){
  for(j in 1:length(posicao4)){
    if(i==j && posicao4[j]==1)
      dados1[i]=9999
  }
}

dados1<-matrix(dados1, nrow=215,ncol=5)

dados1<-data.frame(dados1)
dados1<-cbind(dados[1],dados1)
dados1<-as.data.frame(dados1)

write.table(dados1, "glandNA20.data",row.names = FALSE, col.names=FALSE, sep= "\t")

