# Differential-Gene-Expression-Analysis
Differential gene expression analysis in R to identify over/under expressed genes in SSc fibroblasts.


# Introduction
Systemic sclerosis (SSc) is a rare but severe chronic autoimmune disease characterised by a 
triad of symptoms: vasculopathy, autoimmunity and visceral/cutaneous fibrosis. The extent of 
skin sclerosis has been used to define two subtypes of the disease: limited cutaneous SSc 
(lcSSc) and diffuse cutaneous SSc (dcSSc). In lcSSc skin thickening is distal to elbows and 
knees, with or without face involvement whereas in dcSSc it is proximal i.e., above the elbows 
and knees.

# Single cell RNA-seq
Single cell RNA-seq (scRNA-Seq), a next generation sequencing method, maintains several 
advantages over microarray technology in terms of resolution, accuracy and sensitivity, 
however it has had limited applications with respect to SSc. Therefore, this study aimed
to perform scRNA-seq on dermal fibroblasts cultured from dcSSc patients; skin is commonly 
used as a marker for organ involvement of SSc therefore this cell type was chosen for this 
analysis.

# Methodology
Differentially expressed genes (DEGs) were identified using a differential expression hurdle 
model with the MAST package.The 
p-values for the differential expression were found using the likelihood ratio test and multiple 
hypothesis testing correction was through setting the false discovery rate as <0.05 via 
p.adjust. 


# Results (data_visualisation)
dcSSc clones clustered into two subsets: ‘dcSSc’ and ‘normal-like’. When the subset (clones 2/3)
which resembled the gene expression of healthy fibroblasts was removed from analysis, 5123
DEGs were identified, of which 3665 were upregulated with a fold change greater than 1.5
(P<0.05). 

# PCA 
To gain an insight into the heterogeneity of the SSc and WT cells, PCA analysis was 
performed. 

The WT samples (orange) and the SSc 
samples (blue) in Batch 1 clustered closely together in the positive direction of the PC1 axis 
whereas the SSc samples from Batch 2 clustered in the negative direction.

# Heatmap/ volcano 
To identify if clones should be removed from further analysis,  a heatmap of the top twenty significant genes with respect to the adjusted p value, 
was performed. Samples did cluster by cell type however the SSc clones 2 and 3 appeared to 
branch off early from the other SSc clones and clustered closely to the healthy fibroblasts.

A volcano plot was used to visualise genes that are significantly upregulated/downregulated in SSc.


#Future Study
This analysis provides a set of gene sets and the pathways implicit for future analysis. n addition to this, the 
study corroborated previous evidence of a ‘normal-like’ subtype of dcSSc which has gene 
expression signatures of both healthy and dcSSc cells. 

Overall, the ability of scRNA-seq in elucidating the heterogeneity within cell populations has been shown to be 
crucial in understanding the complex pathogenesis of SSc and is a vital tool for further 
research

