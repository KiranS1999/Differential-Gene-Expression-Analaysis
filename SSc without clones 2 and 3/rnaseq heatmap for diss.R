#Input SSc tpm for qpcr
dfq <- read.table("rna seq heatmap for report 2.txt", header = TRUE)
row.names ( dfq ) <- dfq $ Name
dfq <- dfq [ , -c (1:2) ]
logTPMsq <- dfq
logTPMsq[,1:10] <- log2(logTPMsq[,1:10]+1)
metadataq <- read.delim(file.choose (),
                       stringsAsFactor = FALSE) 
row.names ( metadataq ) <- metadataq $ Sample_ID 
metadataq <- metadataq [,-1]
#Make the type column a factor
metadataq$Type <- as.factor(metadataq$Type)

#Make the Batch column a character
metadataq$Batch <- as.character(metadataq$Batch)

#Check that the rownames of metadata match main data
rownames(metadataq) == colnames(dfq)
#All TRUE


library(SingleCellExperiment)

sceq <- SingleCellExperiment(
  assays = list(tpm = dfq), 
  colData = metadataq
)

#Add logTPM to sce
assay(sceq, "logTPMq") <- log2(tpm(sceq)+1)
exprs_values="logTPMq"

#Heatmap of the genes
library(pheatmap)
normq <- assay(sceq, "logTPMq")
matq <- as.matrix(normq)
annoq <- as.data.frame(colData(sceq)[,c("Type","Batch")])
pheatmap(matq, cluster_rows = T, show_rownames= T,
         annotation= annoq, border_color=NA, fontsize = 10, scale="none",
         fontsize_row = 10, height=20)

normq1 <- assay(sceq, "tpm")
matq1 <- as.matrix(normq1)
annoq1 <- as.data.frame(colData(sceq)[,c("Type","Batch")])
pheatmap(matq1, cluster_rows = T, show_rownames= T,
         annotation= annoq, border_color=NA, fontsize = 10, scale="none",
         fontsize_row = 10, height=20)




