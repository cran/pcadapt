## ----echo=TRUE,include=FALSE---------------------------------------------
library(pcadapt)

## ------------------------------------------------------------------------
data <- read4pcadapt(x="geno3pops",option="example")
print(dim(data))

## ------------------------------------------------------------------------
x <- pcadapt(data,K=20)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="screeplot")

## ------------------------------------------------------------------------
poplist <- c(rep(1,50),rep(2,50),rep(3,50))
print(poplist)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="scores",pop=poplist)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="scores",i=3,j=4,pop=poplist)

## ------------------------------------------------------------------------
x <- pcadapt(data,K=2)

## ------------------------------------------------------------------------
summary(x)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="manhattan")

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="qqplot",threshold=0.1)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
hist(x$pvalues,xlab="p-values",main=NULL)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x,option="stat.distribution")

## ----echo=FALSE----------------------------------------------------------
library(qvalue)

## ------------------------------------------------------------------------
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval<alpha)
print(outliers)

## ------------------------------------------------------------------------
pooldata <- read4pcadapt("pool3pops",option="example",transpose=FALSE)

## ------------------------------------------------------------------------
xpool <- pcadapt(pooldata,data.type="pool")

## ------------------------------------------------------------------------
x_com <- pcadapt(data,K=2,method="communality")

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x_com,option="stat.distribution")

## ------------------------------------------------------------------------
x_eucl <- pcadapt(data,K=2,method="euclidean")

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x_eucl,option="stat.distribution")

## ------------------------------------------------------------------------
x_cw <- pcadapt(data,K=2,method="componentwise")
summary(x_cw$pvalues)

## ----fig.width=7,fig.height=5,fig.align='center'-------------------------
plot(x_cw,option="stat.distribution",K=2)

