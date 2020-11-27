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

exp <- c("WCFS1.glc","WCFS1.glc","WCFS1.rib","WCFS1.rib")

group <- factor(exp)
y <- DGEList(counts=cnts[,2:5],group=group)

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
head(FC, 100)

###

summary(de <- decideTestsDGE(fit))

detags <- rownames(y)[as.logical(de)]
plotSmear(fit, de.tags = detags)
abline(h=c(-1,1), col="blue")

