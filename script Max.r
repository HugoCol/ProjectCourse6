# Huidige werkmap instellen naar map waar script in staat.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Benodigde bibiotheken laden.
require(edgeR)

# Bestand RNA-seq data inlezen.
inhoud = read.table("./WCFS1_cnts.txt", header=TRUE)
row.names(inhoud) = inhoud[,"ID"]
# Annotatiedata toevoegen.
annotatiebestand = "./WCFS1_anno.txt"
annotatie = read.table(annotatiebestand, header=TRUE, sep="\t"
                  , quote="")
names(annotatie)[names(annotatie) == "name"] = "ID"
data = merge(inhoud, annotatie, by="ID")

# Data omzetten naar DGEList.
groep = factor(c("WCFS1.glc.1",	"WCFS1.glc.2",	"WCFS1.rib.1"
                 ,	"WCFS1.rib.2"))
specialeLijst = DGEList(counts=data[,2:5], group=groep)

# Met edgeR de data normaliseren (edgeRUsersGuide()
# voor gebruiksaanwijzing).
x = calcNormFactors(specialeLijst, method="TMM")
keep = filterByExpr(specialeLijst)
x = specialeLijst[keep, , keep.lib.sizes=FALSE]
x$samples$lib.size = colSums(x$counts)

# Dispersie berekenen
#matrix = model.matrix(~0+group, data=x$samples)
#colnames(matrix) = levels(x$samples$group)
#dispersie = estimateDisp(x, matrix)
#matrix

summary(x$counts)
plotMDS(x)
#plotBCV(dispersie)


