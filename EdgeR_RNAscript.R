
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR", force= TRUE)
BiocManager::install("limma", force = TRUE)
BiocManager::install("AnnotationDbi", force = TRUE)
BiocManager::install("org.Hs.eg.db", force = TRUE) 
BiocManager::install("RColorBrewer", force= TRUE)
BiocManager::install("vidger", force= TRUE)
install.packages("gplots")

if (!requireNamespace("GO.db", quietly = TRUE)) {
  BiocManager::install("GO.db")
}


library(limma)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(RColorBrewer) 
library(gplots)
library(vidger)
library(GO.db)


setwd("/home/sajid/Downloads/NIB_projects/RNA_seq/edger")  #set your working directory 


seqdata <- read.delim("final_htcount.txt", stringsAsFactors=FALSE)
sampleinfo <- read.delim("sampleinfo.txt", stringsAsFactors=FALSE)

head(seqdata)
countdata0 <- seqdata[,-(1:2)]
head(countdata0)
rownames(countdata0) <- seqdata[,1]
colnames(countdata0) <- sampleinfo$sampleid
countdata <- countdata0[rowSums(countdata0[])>0,]
head(countdata)

group = factor(sampleinfo$condition)
y <- DGEList(countdata, group=group)

columns(org.Hs.eg.db)
ENTREZID <- mapIds(org.Hs.eg.db,rownames(y),
                   keytype="SYMBOL",column="ENTREZID")   #'select()' returned 1:many mapping between keys and columns

rownames(y$counts) <- ENTREZID
ann<-select(org.Hs.eg.db,keys=rownames(y$counts),
            columns=c("ENTREZID","SYMBOL","GENENAME"))  #'select()' returned 1:many mapping between keys and columns
head(ann)
y$genes <- ann
i <- is.na(y$genes$ENTREZID)
y <- y[!i, ]
condition <- factor(sampleinfo$condition)
design <- model.matrix(~ 0 + condition)
design

keep <- filterByExpr(y, design)
class(countdata)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes=FALSE]
yNorm <- calcNormFactors(y)
yNorm <- estimateDisp(yNorm, design)
names(yNorm)
yNorm$common.dispersion
head(yNorm$trended.dispersion)
head(yNorm$tagwise.dispersion)
jpeg('myBCVPlot.jpg')
plotBCV(yNorm, pch=16, cex=1.2)
dev.off()



png(file="libsizeplot.png")
x <-barplot(yNorm$samples$lib.size/1000000,
           names=colnames(yNorm),
           las=2, ann=FALSE,
           cex.names=0.75, 
           col="lightskyblue",
           space = .5)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")
dev.off()
x




jpeg('mdplots.jpg')
par(mfrow=c(2,2))
for (i in 1:4) {
  plotMD(yNorm, column=i,
         xlab="Average log CPM (all samples)",
         ylab="log-ratio (this vs others)")
  abline(h=0, col="red", lty=2, lwd=2)
}
dev.off()




png(file="logcpmboxplot.png")   
logcounts <- cpm(y, log=TRUE) 
# Boxplot for Unnormalized Data
boxplot(logcounts, xlab="Samples", ylab="Log2 counts per million", las=2, main="Unnormalized logCPMs") 
abline(h=median(logcounts), col="blue") 
# Boxplot for Normalized Data
logcountsNorm <- cpm(yNorm, log=TRUE)
boxplot(logcountsNorm, xlab="Samples", ylab="Log2 counts per million", las=2, main="Normalized logCPMs") 
abline(h=median(logcountsNorm), col="blue") 
dev.off()




# If I want to keep the uunnormalized data
png(file="logcpmboxplot.png")
par(mfrow=c(1,2))
logcounts <- cpm(y,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Unnormalized logCPMs")
logcountsNorm <- cpm(yNorm,log=TRUE)
boxplot(logcountsNorm, xlab="", ylab="Log2 counts per
        million",las=2)
abline(h=median(logcountsNorm),col="blue")
title("Normalized logCPMs")
dev.off()




plotMDS(y, top=500, gene.selection = "pairwise", method = "logFC")
plotMDS(yNorm$counts, top=500, gene.selection = "pairwise", method = "logFC")
str(x)
plotMDS(yNorm$counts, top=500, gene.selection = "pairwise", method = "logFC") 
y
x
yNorm

# Generate an MDS Plot
png(file="MDSPlot.png")
pseudoCounts <- log2(yNorm$counts + 1)
# Color palette based on conditions
colConditions <- brewer.pal(3, "Set2") 
colConditions <- colConditions[match(sampleinfo$condition, 
                                     levels(factor(sampleinfo$condition)))]
# Patient symbols
patients <- c(8, 15, 16)[match(sampleinfo$patient, 
                               levels(factor(sampleinfo$patient)))]
plotMDS(pseudoCounts, 
        pch = patients, 
        col = colConditions, 
        xlim = c(-2, 2)) 
# Legends
legend("topright", 
       lwd = 2, 
       col = brewer.pal(3, "Set2")[1:2], 
       legend = levels(factor(sampleinfo$condition)))
legend("bottomright", 
       pch = c(8, 15, 16), 
       legend = levels(factor(sampleinfo$patient)))
dev.off() 





#png(file="MDSPlot.png")
#pseudoCounts <- log2(yNorm$counts + 1)
#colConditions <- brewer.pal(3, "Set2")
#colConditions <- colConditions[match(sampleinfo$condition,
#                                     levels(factor(sampleinfo$condition)))]
#patients <- c(8, 15, 16)[match(sampleinfo$patient,
#                               levels(factor(sampleinfo$patient)))]
#plotMDS(pseudoCounts, pch = patients, col = colConditions, xlim =
#          c(-2,2))
#legend("topright", lwd = 2, col = brewer.pal(3, "Set2")[1:2],
#        legend = levels(factor(sampleinfo$condition)))
#legend("bottomright", pch = c(8, 15, 16),
#        legend = levels(factor(sampleinfo$patient)))
#dev.off()


#Modified one
png(file = "MDSPlot.png")
# Assuming yNorm$counts contains your data
pseudoCounts <- log2(yNorm$counts + 1)
# Assuming you have four samples, specifying colors manually
colConditions <- c("red", "blue", "green", "orange")
# Assuming you have four patients
patients <- c(8, 15, 16, 4)
# Plotting MDS
plotMDS(pseudoCounts, pch = patients, col = colConditions, xlim = c(-2, 2))
# Creating legends
legend("topright", lwd = 2, col = colConditions[1:2],
       legend = c("Condition 1", "Condition 2"))
legend("bottomright", pch = patients,
       legend = c("Patient 8", "Patient 15", "Patient 16", "Patient 4"))
dev.off()




#png(file="heatmap1.png")
#logcountsNorm <- cpm(yNorm, log = TRUE)
#var_genes <- apply(logcountsNorm, 1, var)
#select_var <- names(sort(var_genes, decreasing = TRUE))[1:10]
#highly_variable_lcpm <- logcountsNorm[select_var,]
#mypalette <- brewer.pal(11, "RdYlBu")
#morecols <- colorRampPalette(mypalette)
#col.con <- c(rep("purple", 2), rep("orange", 4))[factor(sampleinfo$condition)]
#heatmap.2(highly_variable_lcpm,
#          col = rev(morecols(50)), trace = "none",
#          main = "Top 10 most variable genes",
#          ColSideColors = col.con, scale = "row",
#          margins = c(12, 8), srtCol = 45)
# Retrieve and print the names of the highly variable genes
#select_var_names <- rownames(highly_variable_lcpm)
#cat("Highly variable genes selected:\n")
#print(select_var_names)
#dev.off()

#Modifiied one
# Begin PNG device for plot
png(file = "heatmap1.png")
# Assuming yNorm contains your gene expression data
# Compute log-transformed counts per million (CPM)
logcountsNorm <- cpm(yNorm, log = TRUE)
# Calculate variance for each gene across samples
var_genes <- apply(logcountsNorm, 1, var)
# Select the top 10 genes with the highest variance
select_var <- names(sort(var_genes, decreasing = TRUE))[1:10]
# Extract highly variable genes from the log-transformed counts
highly_variable_lcpm <- logcountsNorm[select_var,]
# Define a color palette for the heatmap
mypalette <- brewer.pal(11, "RdYlBu")
morecols <- colorRampPalette(mypalette)

# Assuming sampleinfo contains information about the conditions
# Adjust the col.con based on the number of conditions and samples
col.con <- c(rep("purple", 2), rep("orange", 4))[factor(sampleinfo$condition)]

# Generate heatmap
heatmap.2(highly_variable_lcpm,
          col = rev(morecols(50)), trace = "none",
          main = "Top 10 most variable genes",
          ColSideColors = col.con, scale = "row",
          margins = c(12, 8), srtCol = 45)

# Retrieve and print the names of the highly variable genes
select_var_names <- rownames(highly_variable_lcpm)
cat("Highly variable genes selected:\n")
print(select_var_names)
# End PNG device
dev.off()

fitq <- glmQLFit(yNorm, design)
names(fitq)

jpeg("qlDispplots.jpg")
fitq <- glmQLFit(yNorm, design)
plotQLDisp(fitq, pch=16, cex=1.2)
dev.off()




my.contrasts <- makeContrasts(conditioncontrol-conditiontreated,levels=design)  
fitq <- glmQLFit(yNorm, design)
qlfq<- glmQLFTest(fitq,contrast=my.contrasts)  
topTags(qlfq, n=10, adjust.method="BH", sort.by="PValue",
        p.value=0.05)  


my.contrasts <- makeContrasts(conditioncontrol-conditiontreated,levels=design)
fitg<- glmFit(yNorm, design)
qlfg<- glmLRT(fitg,contrast=my.contrasts)
topTags(qlfg, n=10, adjust.method="BH", sort.by="PValue",
        p.value=0.05)

my.contrasts <- makeContrasts(conditioncontrol-conditiontreated,levels=design)
fitq <- glmQLFit(yNorm, design)
qlfq<- glmTreat(fitq, contrast=my.contrasts, lfc=1)
topTags(qlfq, n=10, adjust.method="BH", sort.by="PValue",
        p.value=0.05)


my.contrasts <- makeContrasts(conditioncontrol-conditiontreated,levels=design)
fitq<-glmQLFit(yNorm, design)
qlfq<-glmQLFTest(fitq, contrast=my.contrasts)
#DEGenes<-decideTestsDGE(qlfq, adjust.method="BH", p.value=0.05, lfc=2)
DEGenes<-decideTests(qlfq, adjust.method="BH", p.value=0.05, lfc=2)
summary(DEGenes)



jpeg('mdPlotfitted.jpg')
my.contrasts <- makeContrasts(conditioncontrol-conditiontreated,levels=design)
fitq<-glmQLFit(yNorm, design)
qlfq<-glmQLFTest(fitq, contrast=my.contrasts)
DEGenes<-decideTestsDGE(qlfq,
                        adjust.method="BH", p.value=0.05, lfc=2)
plotMD(qlfq, status=DEGenes, values=c(1,0,-1), col=c("red","black","blue"), legend="topright")
dev.off()

jpeg('volcano.jpg')
fitq <- glmQLFit(yNorm, design)
qlfq <- glmTreat(fitq, contrast=my.contrasts, lfc=2)
resFilt <- topTags(qlfq, n=100, adjust.method="BH", sort.by="PValue", p.value=1)
volcanoData <- cbind(resFilt$table$logFC, -log2(resFilt$table$PValue))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19)
dev.off()

jpeg('heatmap2.jpg')

# Define contrasts
my.contrasts <- makeContrasts(conditioncontrol-conditiontreated, levels=design)

# Perform differential expression analysis
fitq <- glmQLFit(yNorm, design)
qlfq <- glmQLFTest(fitq, contrast=my.contrasts)
DEGenes <- decideTestsDGE(qlfq, adjust.method="BH", p.value=0.05, lfc=2)

# Prepare data for heatmap
logCPM <- cpm(yNorm, prior.count=2, log=TRUE)
rownames(logCPM) <- yNorm$genes$SYMBOL
colnames(logCPM) <- paste(yNorm$samples$group, 1:3, sep="-")
o <- order(qlfq$table$PValue)
logCPM <- logCPM[o[1:20],]
logCPM <- t(scale(t(logCPM)))

# Define colors for heatmap
col.pan <- colorpanel(100, "blue", "white", "red")

# Generate heatmap
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none", trace="none",
          dendrogram="both", cexRow=1, cexCol=1.4, margin=c(10,9),
          lhei=c(2,10), lwid=c(2,6))

dev.off()




# Define contrasts
my.contrasts <- makeContrasts(conditioncontrol-conditiontreated, levels=design)

# Perform differential expression analysis
fitq <- glmQLFit(yNorm, design)
qlfq <- glmTreat(fitq, contrast=my.contrasts, lfc=1)

# Perform GO enrichment analysis
go <- goana(qlfq, species="Hs")

# Get top 20 enriched GO terms
topGO20 <- topGO(go, sort="down", n=20)

# Write results to a CSV file
write.csv(topGO20, file="topGO20.csv")

my.contrasts <- makeContrasts(conditioncontrol-conditiontreated,levels=design)
fitq<-glmQLFit(yNorm, design)

# Perform differential expression analysis
fitq <- glmQLFit(yNorm, design)
qlfq <- glmTreat(fitq, contrast=my.contrasts, lfc=2)

# Perform KEGG pathway enrichment analysis
keg <- kegga(qlfq, species="Hs")

# Get top 20 enriched KEGG pathways
keg20 <- topKEGG(keg, sort="up", n=20)

# Write results to a CSV file
write.csv(keg20, file="keg20.csv")

#Visualization of the RNA-seq data

#Preparing data
countdata0 <- seqdata[,-(1:2)]
rownames(countdata0) <- seqdata[,1]
colnames(countdata0) <- sampleinfo$sampleid
countdata <- countdata0[rowSums(countdata0[])>0,]
group = factor(sampleinfo$condition)
#Creating DGEList object
y <- DGEList(countdata, group=group)
#Adding annotation
ENTREZID <- mapIds(org.Hs.eg.db,rownames(y),
                   keytype="SYMBOL",column="ENTREZID")
rownames(y$counts) <- ENTREZID
ann<-select(org.Hs.eg.db,keys=rownames(y$counts),
            columns=c("ENTREZID","SYMBOL","GENENAME"))
y$genes <- ann
#Removing rows without Entrez ids
i <- is.na(y$genes$ENTREZID)
y <- y[!i, ]
#Creating the design matrix
condition <- factor(sampleinfo$condition)
design <- model.matrix(~ 0 + condition)
#Filtering genes with low abundance
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes=FALSE]
# Normalizing count data
yNorm <- calcNormFactors(y)
#Estimating dispersions:
yNorm <- estimateDisp(yNorm, design)


jpeg('boxplot.jpg')

# Creating box plot
vsBoxPlot(
  data = yNorm,
  d.factor = NULL,
  type = 'edger',
  title = "Box plots of normal and tumor count data",
  legend = TRUE,
  grid = TRUE
)

dev.off()

jpeg('CPMscatterplot.jpg')

# Creating scatter plot
vsScatterPlot(
  x = 'control',
  y = 'treated',
  data = yNorm,
  type = 'edger',
  d.factor = NULL,
  title = TRUE,
  grid = TRUE
)

dev.off()


jpeg('MAPlot.jpg')

# Creating MA plot
vsMAPlot(
  x = 'control',
  y = 'treated',
  data = yNorm,
  d.factor = NULL,
  type = 'edger',
  padj = 0.05,
  y.lim = NULL,
  lfc = 2,
  title = TRUE,
  legend = TRUE,
  grid = TRUE
)

dev.off()


jpeg('volcanoPlot2.jpg')

# Creating volcano plot
jpeg('volcanoPlot2.jpg')
vsVolcano(
  x = 'control',
  y = 'treated',
  data = yNorm,
  d.factor = NULL,
  type = 'edger',
  padj = 0.05,
  x.lim = NULL,
  lfc = 2,
  title = TRUE,
  legend = TRUE,
  grid = TRUE,
  data.return = FALSE
)

dev.off()





