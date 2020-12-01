library(edgeR)

### read counts

fDir <-  file.choose()
### fName <- "WCFS1_cnts.txt"

### "WCFS1_cnts.txt"
### "WCFS1_RNA-Seq.txt"

cnts <- read.delim(paste0(fDir), comment.char = "#") ### ,fName), comment.char="#")

### used for topTags identification
row.names(cnts) <- cnts[,"ID"]

### Create DGEList

exp <- c("WCFS1.glc","WCFS1.glc","WCFS1.rib","WCFS1.rib", "NC8.glc.1", "NC8.glc.2", "NC8.rib.1", "NC8.rib.2")

group <- factor(exp)
y <- DGEList(counts=cnts[,2:9],group=group)

### Normalise counts
### Trimmed mean of M values : remove lowest and highest values
### (percentile) and calculate mean

y <- calcNormFactors(y, method="TMM" )

### Check statistics

print("Count statistics")
print(summary(y$counts))
print(y$samples)

### Create design matrix

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
print(design)

### Estimate Dispersion

#y <- estimateDisp(y, design)

y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design, method="power")
y <- estimateGLMTagwiseDisp(y,design)

### Plot results

pdf(paste0(fDir,"test.pdf"))
plotMDS(y)
plotBCV(y)
dev.off()

### Fit data

fit <- glmFit(y,design)

### Determine fold changes

mc <- makeContrasts(exp.r=WCFS1.glc-WCFS1.rib, levels=design)

fit <- glmLRT(fit, contrast=mc)

### Print top tags

res<-topTags(fit)
print(res)

### Correlatie van de data

cor.test(y$counts[,1], y$counts[,2])
cor.test(y$counts[,3], y$counts[,4])

### plot van diviance en dispersion

plot(deviance(fit))
plot(fit$dispersion)

### Alle p-values

PV <- fit[["table"]][["PValue"]]
head(PV, 100)

### Alle fold changes in log

logFC <- fit[["table"]][["logFC"]]
head(logFC, 100)

### Alle fold changes normaal

FC <- logFC^2
head(FC, 33)

### Cluster

cntsCluster <- kmeans(cnts[, 2:5], 4, nstart = 20)
cntsCluster
cntsCluster$cluster <- as.factor(cntsCluster$cluster)
plot(cntsCluster$cluster)

bing =kmeans(cnts[, 2:9], 4)
bing
plot(cnts[c("WCFS1.glc.1", "NC8.glc.1")], col = bing$cluster)#"WCFS1.glc.2", "WCFS1.rib.1", "WCFS1.rib.2", "NC8.glc.1", "NC8.glc.2", "NC8.rib.1", "NC8.rib.2")], col = bing$cluster)
bing$size

install.packages("dplyr")
install.packages("ggfortify")
library(stats)
library(dplyr)
library(ggplot2)
library(ggfortify)
view(cnts)

wssplot <- function(data, nc=15, seed=1234)
{
  wss <- (nrow(data)-1)*sum(apply(data, 2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type='b', xlab = "number of clusters", 
       ylab = "within groups sum of squares")
}

mydata = select(cnts, c(2,3,4,5,6,7,8,9))
wssplot(mydata)
KM = kmeans(mydata, 4)
autoplot(KM, mydata, frame=TRUE)
KM$centers
###

summary(de <- decideTestsDGE(fit))

detags <- rownames(y)[as.logical(de)]
plotSmear(fit, de.tags = detags)
abline(h=c(-1,1), col="blue")

?kmeans

res <- aggregate(.~group,df,mean) 

res["fc",]=c("A.vs.B",as.numeric(res[1,-1])/as.numeric(res[2,-1])) 

bing <- "bong"
bing