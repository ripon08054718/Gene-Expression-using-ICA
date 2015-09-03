#### Gene Expression Data Analysis using ICA....

## Directory setting..
setwd("H:\\BreastICA")


## Library Function...
library(MineICA)
library(fastICA)
library(moments)
library(ggplot2)
library(gplots)
library(plotrix)
library(blockmodeling)
library(heatmap3)
library(genefilter)
library(ALL)
library(dplyr)
library(reshape)
library(scales)
library(reshape2)
library(plyr)
library(made4)
library(d3heatmap)
library(scales)
library(mygene)
library(hgu133a.db)

## Read txt data
rawData<-read.table("Breast_Cancer.txt",header=TRUE)

## Selecting 12 samples (3 normal...)..
subData12<-rawData[,c("ID_REF","GSM940873", "GSM940874", "GSM940875", "GSM940882", "GSM940883", "GSM940884", 
"GSM940879", "GSM940880", "GSM940881", "GSM940888", "GSM940889", "GSM940890")]



# Gene Filtering using IQR (Inter Quartile Range)...
#subData<- selectFeatures_IQR(subData12,4000)
#dim(subData)

# Gene Filtering using CV (Coefficient of Variation and ...
ffun <- filterfun(pOverA(0.05,500), cv(0.15,10))
t.fil <- genefilter(subData12[,-1],ffun)
subData <- log2(subData12[,-1][t.fil,])
rownames(subData)<-subData12[as.numeric(rownames(subData)),1]

head(subData)
dim(subData)

## Heatmap...
jpeg("rawHeatmap.jpeg")
heatmap(as.matrix(subData),col=redgreen(100),margins = c(10,20),
keep.dendro = FALSE,scale="none",labRow = NULL, labCol = NULL)
dev.off()


jpeg("rawHeatmap.jpeg")
heatplot(as.matrix(subData),scaleKey=F, dend="none",
lowcol = "white", highcol = "magenta",cols.default =T)
dev.off()


rc <- rainbow(nrow(subData), start = 0, end = .7)
cc <- rainbow(ncol(subData), start = 0, end = .7)

hv <- heatmap(as.matrix(subData), col = cm.colors(100), scale = "column",
              RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
              xlab = "Samples", ylab =  "Genes",
              main = "")
utils::str(hv) # the two re-ordering index vectors


#geneNames<-subData12[as.numeric(rownames(subData)),1]
#subData<-cbind(geneNames,subData)

#subData$geneNames <- with(subData, reorder(geneNames, GSM940873))
#str(subData)
#nba.m <- melt(subData)
#head(nba.m)
#str(nba.m)

#nba.m <- ddply(nba.m, .(variable), transform,
#               rescale = scales::rescale(value))


#(p <- ggplot(nba.m, aes(variable, geneNames)) + geom_tile(aes(fill = rescale),
#      colour = "white") + scale_fill_gradient(low = "white",high = "steelblue"))



############ -------------  Fast ICA -------------################
### Working with fastICA algorithm...
res <- clusterFastICARuns(X=as.matrix(subData), nbComp=12, alg.type="deflation",
                          nbIt=10, funClus="hclust", method="ward")


## Hinton Diagram for identifying gene symbol...
plot.mat(res$A)
border.col<-color.scale(res$A,extremes=2:3)
color2D.matplot(res$A*0.01,extremes=c(2,3),main="Hinton diagram (green +, red -)",
                Hinton=TRUE,border=border.col)



heatmap(t(res$A),col=redgreen(100))
heatmap(as.matrix(res$S),col=redgreen(100),keep.dendro = F)
kurtosis(res$A)
kurtosis(res$S)

heatplot(res$S)


#####----- Significant Gene Identification using ICA ------########

M<-res$S; dim(M);
M1<-M[,rank(kurtosis(res$A))==1]
M2<-M[,rank(kurtosis(res$A))==2]


mean(M2)
2*sd(M2)
median(M2)

mean(M1)
2*sd(M1)
median(M1)

IQR(M1)
quantile(M1,0.025)




library(fitdistrplus)
library(logspline)
descdist(M1, discrete=F,obs.col = "darkblue", boot = 1000,
obs.pch = 16, boot.col="red")


Thresold<-2.0

par(mfrow=c(2,1))
barplot(t(M1))
abline(h = Thresold, v=NULL, col = "red")
abline(h=-Thresold, v=NULL, col="red")

barplot(t(M2))
abline(h = Thresold, v=NULL, col = "red")
abline(h=-Thresold, v=NULL, col="red")



cond<-abs(M1)<Thresold
MM1<-M1[abs(M1)>=Thresold]
write.csv(MM1,"gene.csv")

cond<-abs(M2)<Thresold
M2[cond]<-0.0
MM2<-M2
write.csv(MM2,"gene2.csv")


par(mfrow=c(2,1))
barplot(t(MM1))
abline(h = Thresold, v=NULL, col="red")
abline(h=-Thresold, v=NULL, col="red")

barplot(t(MM2))
abline(h = Thresold, v=NULL, col = "red")
abline(h=-Thresold, v=NULL, col="red")


genes<-MM1[which(MM1>0)]

sort(genes, T)
count(MM1>0)
count(MM1<0)

count(MM2>0)
count(MM2<0)


################### ------------- PCA ---------------- ###################.

source("H:/BreastICA/Necessary Functions/genePCA.R")
gene_expr_pca("")

source("http://bioconductor.org/biocLite.R")
biocLite("hgu133a.db")


library(mygene)
library(hgu133a.db)
geneList<-names(MM1)
res <- queryMany(geneList, scopes='symbol', fields=c('entrezgene', 'go'), 
species='human')
res[1, 'go.CC'][[1]]

