#Input SSc dataset
df <- read.table("SSc dataset.txt", header = TRUE)

#Change column names
names(df)[names(df) == "FDG1_TPM"] <- "Healthy_Fibroblast_1"
names(df)[names(df) == "FDG2_TPM"] <- "Healthy_Fibroblast_2"
names(df)[names(df) == "FDG3_TPM"] <- "Healthy_Fibroblast_3"
names(df)[names(df) == "FDG4_TPM"] <- "SSc_Clone_1"
names(df)[names(df) == "SSC_Clone1_TPM"] <- "SSc_Clone_2"
names(df)[names(df) == "SSC_Clone2_TPM"] <- "SSc_Clone_3"
names(df)[names(df) == "SSC_Clone3_TPM"] <- "SSc_Clone_4"
names(df)[names(df) == "SSC_Clone4_TPM"] <- "SSc_Clone_5"
names(df)[names(df) == "SSC_Clone5_TPM"] <- "SSc_Clone_6"
names(df)[names(df) == "SSC_Clone6_TPM"] <- "SSc_Clone_7"

#Make gene names unique and make them the rowname
rownames(df) <- make.names(df[,1], unique = TRUE)
row.names ( df ) <- df $ Name

#Delete the two gene name/engs columns
df <- df [ , -c (1:2) ]

#Delete ssc clones 2 and 3
df <- df[ , -c (5:6) ]

#Add in metadata and change the rownames to the sample-id
metadata <- read.delim(file.choose (),
                       stringsAsFactor = FALSE) 

#Remove ssc clones 2 and 3 from metadata
metadata <- metadata[-c(5:6), ]
row.names ( metadata ) <- metadata $ Sample_ID 
metadata <- metadata [,-1]

#Make the type column a factor
metadata$Type <- as.factor(metadata$Type)

#Make the Batch column a character
metadata$Batch <- as.character(metadata$Batch)

#Check that the rownames of metadata match main data
rownames(metadata) == colnames(df)
#All TRUE

#Make the type column a factor
metadata$Type <- as.factor(metadata$Type)

#Make the Batch column a character
metadata$Batch <- as.character(metadata$Batch)


#log the tpm value
logTPMs <- df
logTPMs[,1:8] <- log2(logTPMs[,1:8]+1)

#Make sce 
library(SingleCellExperiment)

sce <- SingleCellExperiment(
  assays = list(tpm = df), 
  colData = metadata
)

#remove rows that are zeros
keep_feature <- rowSums(tpm(sce) > 0) > 0
sce <- sce[keep_feature,]

#Add logTPM to sce
assay(sce, "logTPM") <- log2(tpm(sce)+1)
exprs_values="logTPM"

#make sca
library(MAST)

sca = SceToSingleCellAssay(sce)
SceToSingleCellAssay(sce, class = "SingleCellAssay", check_sanity = TRUE)

#Do DEG analysis
cond<-factor(colData(sca)$Type)
cond<-relevel(cond,"WT")
colData(sca)$Type<-cond
zlmCond <- zlm(~ Type, sca)

# summary with likelihood test ratio
summary_lrt <- summary(zlmCond, doLRT= 'TypeSSc') 
summary_lrt #TO view it

#getsummarytable
fit <- summary_lrt$datatable

#pvalue df
pvalue <- fit[contrast == 'TypeSSc' & component == 'H', .(primerid, `Pr(>Chisq)`)]

#logFC df
logFC <- fit[contrast == 'TypeSSc' & component == 'logFC', .(primerid, coef)]

#pvalue and logFC
fit <- merge(pvalue, logFC, by = 'primerid')

#print head of fit
head(fit)

#adjusted pvalues
fit[, padjusted:=p.adjust(`Pr(>Chisq)`, 'fdr')]

#result table
res <- data.frame(gene = fit$primerid,
                  pvalue = fit[,'Pr(>Chisq)'],
                  padjusted = fit$padjusted,
                  logFC = fit$coef)

#export results
library(openxlsx)
write.xlsx(res, "DEG genes without clone 2 and 3.xlsx")

#run PCA
library(scater)
PCAsce <- sce
names(assays(PCAsce))[names(assays(PCAsce)) == "logTPM"] <- "logcounts"
PCAsce <- runPCA(PCAsce)
str(reducedDim(PCAsce, "PCA"))
plotPCA(PCAsce,ncomponents=2,colour_by="Type",shape_by="Batch")
plotPCA(PCAsce,ncomponents=4,colour_by="Type",shape_by="Batch")

#Heatmap of all genes
library(pheatmap)
norm <- assay(sce, "logTPM")
mat <- as.matrix(norm)
anno <- as.data.frame(colData(sce)[,c("Type","Batch")])
pheatmap(mat, cluster_rows = F, show_rownames=F,
         annotation= anno, border_color=NA, fontsize = 10, scale="none",
         fontsize_row = 10, height=20)

#hetamap of report 2 genes
sce


#logTPM spread of Top 20 significant genes
res <- res[order(res$padjusted), ]
mostDE <- res$gene[1:20]
mostDE10 <- res$gene[1:10]
top_20 <- mostDE[1:20]

#normalized counts for top 20 significant genes
top_20_norm <- data.frame(logTPMs[top_20, ])

# Create a column with the gene names (from row names)
top_20_norm$gene <- rownames(top_20_norm)

#need to gather counts for all samples into one column
library(reshape)
melted_top_20 <- melt(top_20_norm)

# check the column header in the "melted" data frame
View(melted_top_20)

## add column names that make sense
colnames(melted_top_20) <- c("gene", "CellType", "logTPM")

# add metadata to "melted" dataframe
metadata$sample_ID <- rownames(metadata)
melted_top_20 <- merge(melted_top_20, metadata)

#plot using ggplot2
ggplot(melted_top_20) +
  geom_point(aes(x = gene, y = logTPM, color = CellType)) +
  scale_y_log10() +
  xlab("Gene_ID") +
  ylab("logTPM") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))

#Heatmap of Top 20 significant genes
library(pheatmap)
norm20 <- assay(sce[mostDE, ], "logTPM")
mat20 <- as.matrix(norm20)
anno20 <- as.data.frame(colData(sce)[,c("Type","Batch")])
pheatmap(mat20, cluster_rows = T, show_rownames=T,
         annotation= anno20, border_color=NA, fontsize = 10, scale="none",
         fontsize_row = 10, height=20)

#Violin plots pf top 10
scevoilin1 <- sce
sceviolin1 <- sce[mostDE10,]
plotExpression(sceviolin1, rownames(sceviolin1) [1:10],
               x = "Type", exprs_values = "logTPM", 
               colour = "Batch",log=TRUE)


plotExpression(sceviolin, rownames(sceviolin) [11:20],
               x = "Type", exprs_values = "logTPM", 
               colour = "Batch",log=TRUE)

#Too squished if you plot all twenty genes at once

#See significnat genes wrt thresholds
#set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

threshold <- res$padjusted < padj.cutoff & abs(res$logFC) > lfc.cutoff
#only got two genes when we apply threshold;ENSG00000171246 and  ENSG00000153885

#add this vector to our results table
res$threshold <- threshold 
significantgenes <- data.frame(subset(res, threshold==TRUE))

#Volcano plots
library(ggrepel)



p <- ggplot(data=res, aes(x=logFC, y=-log10(padjusted))) + geom_point()+ theme_minimal()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.1), col="red")
res$diffexpressed <- "NO"
res$diffexpressed[res$logFC > 0.6 & res$padjusted < 0.1] <- "UP"
res$diffexpressed[res$logFC < -0.6 & res$padjusted < 0.1] <- "DOWN"
p <- ggplot(data=res, aes(x=logFC, y=-log10(padjusted), col=diffexpressed)) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.1), col="red")
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3

library(ggrepel)




res$Key <- "Not differentially expressed"
res$Key[res$logFC > 0.6 & res$padjusted < 0.05] <- "Upregulated"
res$Key[res$logFC < -0.6 & res$padjusted < 0.05] <- "Downregulated"
p <- ggplot(data=res, aes(x=logFC, y=-log10(padjusted), col=Key)) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("Downregulated", "Upregulated", "Not differentially expressed")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3
#cant label bc there is too many genes

#TABLE OF SIG GENES

FCTHRESHOLD <- log2(1.5)
summaryCond <- summary(zlmCond, doLRT='TypeSSc') 

print(summaryCond, n=4)
print(summaryCond, n=4, by='D')
print(summaryCond, n=4, by='C')

summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='TypeSSc' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='TypeSSc' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]

library(data.table)
fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
head(fcHurdleSig)

entrez_to_plot <- fcHurdleSig[1:50,primerid]
symbols_to_plot <- fcHurdleSig[1:50,primerid]

flat <- as(sca[entrez_to_plot,], 'data.table')
head(flat)

library(openxlsx)
write.xlsx(fcHurdleSig, "fchurdleSigDEGgeneswithoutclone2and3.xlsx")



