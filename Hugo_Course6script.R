# Huidige werkmap instellen naar map waar script in staat. 
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) 

# Benodigde bibiotheken laden. 
require(edgeR) 
require(limma)

# Bestand RNA-seq data inlezen. 
inhoud = read.table("./WCFS1_cnts.txt", header=TRUE) 
row.names(inhoud) = inhoud[,"ID"] 

# Annotatiedata toevoegen. 
annotatiebestand = "./WCFS1_anno.txt" 
annotatie = read.table(annotatiebestand, header=TRUE, sep="\t" , quote="") 
names(annotatie)[names(annotatie) == "name"] = "ID" 
data = merge(inhoud, annotatie, by="ID") 

# Data omzetten naar DGEList. 
groep = factor(c("WCFS1.glc.1",	"WCFS1.glc.2",	"WCFS1.rib.1" ,	"WCFS1.rib.2")) 
specialeLijst = DGEList(counts=data[,2:5], group=groep) 

# Met edgeR de data normaliseren (edgeRUsersGuide() 
# voor gebruiksaanwijzing). 

x = calcNormFactors(specialeLijst, method="TMM") 
keep = filterByExpr(specialeLijst) 
x = specialeLijst[keep, , keep.lib.sizes=FALSE] 
x$samples$lib.size = colSums(x$counts) 
x

# Dispersie berekenen 

matrix = model.matrix(~0+group, data=x$samples) 
colnames(matrix) = levels(x$samples$group) 
dispersie = estimateDisp(x, matrix) 
matrix 

# fold change berekenen

summary(x$counts) 
plotMDS(x)


# Hier maak je een dataframe met de informatie die je wilt vergelijken 
df=data.frame(group=rep(c("A","B"),5), value1=1:10,value2=21:30,value3=41:50,stringsAsFactors = F)
df
# Dit is de code om een Fold Change te berekening uit de limma package 
res=aggregate(.~group,df,mean) 
res["fc",]=c("A.vs.B",as.numeric(res[1,-1])/as.numeric(res[2,-1])) 
res
# welke genen zijn interessant om te kiezen adhv hun Fold Change? 
# Een histogram maken endan kijken hoe het utpakt 

# Welke P-waarde gebruik je om je gen expressie te controleren? 
# Dit heb je niet echt nodig. Genen die niet veel tot expressie gekomen zijn worden door filteren en normaliseren al uit je data gehaald. Je zou eventueel naar de FDR (FoldChange Discovery Rate) kunnen kijken als je een cut-off wilt voor de frequentie waarin genen voorkomen.  

# Wat is K-mean? 

# Is het berekenen van afstand. Je begint met twee samples en je kijkt dan welke samples het dichtste bij vallen. Je maakt zo groepjes van bij elkaar vallende samples.
# Je geeft aan de functie het aantal centers mee. Voor ieder center word data gegroepeerd die het meest met elkaar overeenkomt. 

# Hoe cluster je genen? 

# Doe je met de functie hclust 
# Groepeer niet op de genen maar de samples. Dan is het duidelijk 

# wanneer is dit hierarchisch?
# Je kan dit hierarchisch maken door de samples die het meest op elkaar lijken te koppelen en dan een soort fylogenetisch boom te maken.  