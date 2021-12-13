#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("pasilla")
#BiocManager::install("apeglm")
#BiocManager::install("BiocParallel")
#BiocManager::install("pheatmap")
#BiocManager::install("DEGreport")
#BiocManager::install("pathview")
#BiocManager::install("clusterProfiler")
library(DESeq2)
library(DEGreport)
library("pasilla")
library(apeglm)
library("BiocParallel")
register(MulticoreParam(4))
library("ggplot2")
library(ggrepel)
library("pheatmap")
library("RColorBrewer")
library(dplyr)
library(tibble)
library(pathview)
library(clusterProfiler)
require(biomaRt)
library(ggnewscale)
library(ggpubr)

###Count matrix input
cts <- as.matrix(read.csv("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/counts.txt",sep="\t",row.names="GeneID"))
coldata<- read.csv("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/sample_annotation.csv")
colnames(cts)[match(coldata[,"Index"], colnames(cts))] = coldata[,"Sample"]

##format sample annotation file
rownames(coldata)<-coldata$Sample
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts)) ###not in the same order
coldata$Transfers<-factor(coldata$Transfers)

cts <- cts[, rownames(coldata)] ###order the counts data by sample annotation file order
all(rownames(coldata) == colnames(cts)) ###in same order

###create DESeq data set
##compare condition effect accross isolates
###m<-model.matrix(~ Isolate+Isolate:condition,coldata) ###uses one of the isolates 8115 as the reference
m<-model.matrix(~ condition+Isolate:condition,coldata) #Control is the reference level (Intercept)
colnames(m) 

### Complex model design:
####condition:Isolate studies the effect of the Isolate on the condition effect 
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~Isolate+condition+Isolate:condition)

###minimal filtering on genes with few reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

###indicate the comparison
dds$condition <- factor(dds$condition, levels = c("Control","Irradiated"))
dds$Isolate <- factor(dds$Isolate, 
                      levels = c("WT0","WT.C1","PKS.C1",
                                 "WT15.1.2","WT15.2.1","WT15.2.2","WT15.4.1","WT15.4.2",
                                 "PKS15.1.1","PKS15.2.1","PKS15.2.2","PKS15.3.2","PKS15.4.2",
                                 "8115","R52"))

##############################################
###Data Quality Assessment and Exploration####
##############################################
####transformation on variance
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
####heatmap of count matrix
df <- as.data.frame(colData(dds)[,c("condition","Isolate")])
pheatmap(assay(ntd), cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

####Heatmap of sample clustering
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$Isolate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pcaData<-plotPCA(vsd, intgroup=c("condition", "Isolate"),returnData=TRUE)
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition), shape = factor(Isolate))) + 
  geom_point(size =3, aes(fill=factor(Isolate), alpha=as.character(Isolate))) + 
  scale_shape_manual(values=1:nlevels(pcaData$Isolate))

pcaData<-plotPCA(vsd, intgroup=c("condition", "genotype"),returnData=TRUE)
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition), shape = factor(genotype))) + 
  geom_point(size =3, aes(fill=factor(genotype), alpha=as.character(genotype)))  

pcaData<-plotPCA(vsd, intgroup=c("condition", "Sample"),returnData=TRUE)
pcaData$name<- sub("_.*", "", pcaData$name)
pcaData$name<- sub("PKS C.*", "PKS C1", pcaData$name)
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition), shape = factor(name))) + 
  geom_point(size =3, aes(fill=factor(name), alpha=as.character(name))) + 
  geom_text(aes(label=Sample),hjust=0, vjust=0)

##Two WT15.1.2 libraries (one irradiated WT15.1.2_3 and one not irradiated WT15.1.2_6) do not cluster with other WT but with R52, toss out
idx<-which(colnames(vsd)==c("WT15.1.2_3","WT15.1.2_6"))
vsd2<-vsd[,-idx]

##create DESeqDataSet for only WT and PKS (evolution dataset)
idx_evol<-which(grepl("WT.*|PKS.*", colnames(vsd2)))
vsd_evol<-vsd2[,idx_evol]
pheatmap(assay(vsd_evol), cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)

sampleDists <- dist(t(assay(vsd_evol)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_evol$condition, vsd_evol$Isolate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pcaData<-plotPCA(vsd_evol, intgroup=c("condition", "Sample"),returnData=TRUE)
pcaData$name<- sub("_.*", "", pcaData$name)
pcaData$name<- sub("PKS C.*", "PKS C1", pcaData$name)
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition), shape = factor(name))) + 
  geom_point(size =3, aes(fill=factor(Sample), alpha=as.character(Sample))) + 
  geom_text(aes(label=Sample),hjust=0, vjust=0)+
  theme(legend.title = element_blank()) 

##create DESeqDataSet for only 8115 and R52
idx_genes<-which(grepl("WT0.*|8115.*|R52.*", colnames(vsd)))
vsd_genes<-vsd[,idx_genes]

sampleDists <- dist(t(assay(vsd_genes)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_genes$condition, vsd_genes$Isolate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pcaData<-plotPCA(vsd_genes, intgroup=c("condition", "Sample"),returnData=TRUE)
pcaData$name<- sub("_.*", "", pcaData$name)
pcaData$name<- sub("PKS C.*", "PKS C1", pcaData$name)
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition), shape = factor(name))) + 
  geom_point(size =3, aes(fill=factor(Sample), alpha=as.character(Sample))) + 
  geom_text(aes(label=Sample),hjust=0, vjust=0)


###############################################################
#####Dispersion during differential expression analysis########
###############################################################
####estimates of variation for each gene are often unreliable, 
####To address this problem, DESeq2 shares information across genes to generate more accurate 
####estimates of variation based on the mean expression level of the gene using a method called ‘shrinkage’. 
####DESeq2 assumes that genes with similar expression levels should have similar dispersion.
idx<-which(colnames(cts)==c("WT15.1.2_3","WT15.1.2_6"))
cts2<-cts[,-idx]
idx_evol<-which(grepl("WT.*|PKS.*", colnames(cts2)))
cts_evol<-cts2[,idx_evol]

idx_genes<-which(grepl("WT0.*|8115.*|R52.*", colnames(cts2)))
cts_genes<-cts2[,idx_genes]

dds_evol <- DESeqDataSetFromMatrix(countData = cts_evol,
                              colData = coldata[c(idx_evol),],
                              design = ~Isolate+condition+Isolate:condition)
dds_genes <- DESeqDataSetFromMatrix(countData = cts_genes,
                                   colData = coldata[c(idx_genes),],
                                   design = ~Isolate+condition)

dds <- DESeq(dds)
dds_evol <- DESeq(dds_evol)
dds_genes <- DESeq(dds_genes)

plotDispEsts(dds)
plotDispEsts(dds_evol)
plotDispEsts(dds_genes)
#Red Curve: expected dispersion value for genes of a given expression strength
#Black dot: a gene with an associated mean expression level and maximum likelihood estimation (MLE) of the dispersion (Step 1).
#BlueCircle: genes with extremely high dispersion values are not shrunk towards the cuve; the gene does not follow the modeling assumptions and has higher variability than others for biological or technical reasons

#Genes with low dispersion estimates are shrunken towards the curve, and the more accurate, higher shrunken values are output for fitting of the model and differential expression testing
#Dispersion estimates that are slightly above the curve are also shrunk toward the curve for better dispersion estimation

#######################################################################
###Example 4: Two conditions, multiple genotpes with interaction terms#
#######################################################################

###The effect of treatment in wild-type/WT0 (the main effect)
###This is for WT, treated compared with untreated. 
###Note that WT is not mentioned, because it is the reference level.
res_WT_irr_vs_ctl = results(dds, contrast=c("condition","Irradiated","Control"),alpha = 0.05)
###The effect of treatment in mutant (8115 and R52)
###The main effect plus the interaction term (the extra condition effect in genotype Mutant compared to genotype WT).
res_R52_vs_WT_irr_vs_ctl <- results(dds, list( c("condition_Irradiated_vs_Control","IsolateR52.conditionIrradiated") ),alpha = 0.05)
res_8115_vs_WT_irr_vs_ctl <- results(dds, list( c("condition_Irradiated_vs_Control","Isolate8115.conditionIrradiated") ),alpha = 0.05)

###The effect of R52 vs 8115, under condition Control.
###This is the effect of “genotype_R52_vs_WT0” minus “genotype_8115_vs_WT0”
#resultsNames(dds) ### look at the order of the coefficient names
#res_R52_vs_8115_ctl <- results(dds, contrast=c(0,0,0,0,0,0,0,0,0,0,
#                                               0,0,0,-1,1,0,0,0,0,0,
#                                               0,0,0,0,0,0,0,0,0,0),alpha = 0.05) #### same as list("Isolate_R52_vs_WT0","Isolate_8115_vs_WT0")
res_R52_vs_8115_ctl <- results(dds, contrast=list("Isolate_R52_vs_WT0","Isolate_8115_vs_WT0"),alpha = 0.05)
###The effect of R52 vs 8115, under condition Irradiation.
res_R52_vs_8115_irr <- results(dds, contrast=list("IsolateR52.conditionIrradiated","Isolate8115.conditionIrradiated"),alpha = 0.05)

###The different response in genotypes (interaction term)
###Is the effect of treatment different between mutant and WT? This is the interaction term.
res_R52_vs_WT_irr = results(dds, name="IsolateR52.conditionIrradiated",alpha = 0.05)
res_8115_vs_WT_irr = results(dds, name="Isolate8115.conditionIrradiated",alpha = 0.05)

###What is the difference between mutant and wild-type without treatment?
###As Ctrl is the reference level, we can just retrieve the “genotype_MU_vs_WT”
res_R52_vs_WT_ctl = results(dds,name="Isolate_R52_vs_WT0",alpha = 0.05) ####same as res_R52_vs_WT_ctl = results(dds, contrast=c("Isolate","R52","WT0"),alpha = 0.05) 
res_8115_vs_WT_ctl = results(dds,name="Isolate_8115_vs_WT0",alpha = 0.05)

Plot_results<-function(res,name){
  ix = which.min(res$padj) # most significant
  res <- res[order(res$padj),] # sort
  res$id<-as.character(gsub("jgi.p.*Exode1.","",row.names(res)))
  pdf(paste0("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/",name,".pdf"))
  par(mfrow=c(2,1))
  barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]) #Here we show the most significant gene.
  plotMA(res)
  graphics.off()
  ###Volcano plots
  pdf(paste0("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/",name,".volcano.pdf"))
  p<-DEGreport::degVolcano(
    data.frame(res[,c("log2FoldChange","padj")]), # table - 2 columns
    plot_text = data.frame(res[1:10,c("log2FoldChange","padj","id")])) # table to add names
  print(p)
  graphics.off()
  #Create a tibble of results (all the genes with transcripts)
  res_tableOE_tb <- res%>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  res_tableOE_tb <- res_tableOE_tb %>% 
    mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)
  # Compare to numbers we had from Wald test
  sigOE <- res_tableOE_tb %>%
    filter(padj < 0.05)
    ## Volcano plot
  ###output in table format significant genes 
  write.csv(as.data.frame(sigOE), 
            file=paste0("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/",name,".csv"))
  }

Plot_results(res_R52_vs_WT_irr,"res_R52_vs_WT_irr")
Plot_results(res_8115_vs_WT_irr,"res_8115_vs_WT_irr") 
Plot_results(res_WT_irr_vs_ctl,"res_WT_irr_vs_ctl")
Plot_results(res_R52_vs_WT_irr_vs_ctl,"res_R52_vs_WT_irr_vs_ctl")
Plot_results(res_8115_vs_WT_irr_vs_ctl,"res_8115_vs_WT_irr_vs_ctl")
Plot_results(res_R52_vs_8115_ctl,"res_R52_vs_8115_ctl")
Plot_results(res_R52_vs_8115_irr,"res_R52_vs_8115_irr") 
Plot_results(res_R52_vs_WT_ctl,"res_R52_vs_WT_ctl")
Plot_results(res_8115_vs_WT_ctl,"res_8115_vs_WT_ctl")



######################################################################
### GO Enrichment Analysis ###########################################
######################################################################
ED_gff_GO<-read.table("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/Exode1_GeneCatalog_proteins_20160902_GO.genename.tab",
                      sep='\t')
colnames(ED_gff_GO)<-c("gene","gotermId","goName", "keytypes","goAcc")

ED_gff_GO$gene<-as.character(gsub("jgi.p.*Exode1.","",ED_gff_GO$gene))

ED_gff_GO$keytypes[ED_gff_GO$keytypes=="biological_process"]<-"BP"
ED_gff_GO$keytypes[ED_gff_GO$keytypes=="cellular_component"]<-"CC"
ED_gff_GO$keytypes[ED_gff_GO$keytypes=="molecular_function"]<-"MF"

GO_results<-function(res,name){
  #Create a tibble of results (all the genes with transcripts)
  row.names(res)<-as.character(gsub("jgi.p.*Exode1.","",row.names(res)))
  res_tableOE_tb <- res%>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  res_tableOE_tb <- res_tableOE_tb %>% 
    mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)
  ## Merge the AnnotationHub dataframe with the results 
  res_ids <- left_join(res_tableOE_tb, ED_gff_GO, by=c("gene"="gene"))  
  ## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
  allOE_genes <- as.character(res_ids$gene)
  ## Extract significant results
  sigOE <- res_tableOE_tb %>%
    filter(padj < 0.05)
  ## Volcano plot
  sigOE_genes <- as.character(sigOE$gene)
  
  term2gene<-ED_gff_GO[,c("goName","gene", "keytypes")]
  #term2gene_BP<-term2gene[which(term2gene$keytypes=="BP"),]
  #term2gene_CC<-term2gene[which(term2gene$keytypes=="CC"),]
  #term2gene_MF<-term2gene[which(term2gene$keytypes=="MF"),]
  
  ########### ALL Processes #########################
  ego <- enricher(gene = sigOE_genes, 
                     universe = allOE_genes,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05,
                     TERM2GENE = term2gene[,c("goName","gene")])
  ## Dotplot 
  a<-dotplot(ego,orderBy="GeneRatio")
  ## Output results from GO analysis to a table
  cluster_summary <- data.frame(ego)
  ## Add similarity matrix to the termsim slot of enrichment result
  ego <- enrichplot::pairwise_termsim(ego)
  ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
  b<-emapplot(ego, showCategory = 50)
  ## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
  OE_foldchanges <- sigOE$log2FoldChange
  
  names(OE_foldchanges) <- sigOE$gene
  
  ## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
  c<-cnetplot(ego, 
           categorySize="pvalue", 
           showCategory = 5, 
           foldChange=OE_foldchanges, 
           vertex.label.font=6)
  pdf(paste0("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/",name,".GO.pdf"),width=14,height=18)
  print(ggarrange(a,b,c,nrow = 3,ncol=1))
  graphics.off()
  
  # ########### Biological Process #########################
  # ego_BP <- enricher(gene = sigOE_genes, 
  #                 universe = allOE_genes,
  #                 pAdjustMethod = "BH", 
  #                 qvalueCutoff = 0.05,
  #                 TERM2GENE = term2gene_BP[,c("goName","gene")])
  # ## Dotplot 
  # dotplot(ego_BP,orderBy="GeneRatio")
  # ## Output results from GO analysis to a table
  # cluster_summary <- data.frame(ego_BP)
  # ## Add similarity matrix to the termsim slot of enrichment result
  # ego_BP <- enrichplot::pairwise_termsim(ego_BP)
  # ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
  # emapplot(ego_BP, showCategory = 50)
  # ## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
  # OE_foldchanges <- sigOE$log2FoldChange
  # 
  # names(OE_foldchanges) <- sigOE$gene
  # 
  # ## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
  # cnetplot(ego_BP, 
  #          categorySize="pvalue", 
  #          showCategory = 5, 
  #          foldChange=OE_foldchanges, 
  #          vertex.label.font=6)
  # 
  # ########### Cellular Component #########################
  # ego_CC <- enricher(gene = sigOE_genes, 
  #                    universe = allOE_genes,
  #                    pAdjustMethod = "BH", 
  #                    qvalueCutoff = 0.05,
  #                    TERM2GENE = term2gene_CC[,c("goName","gene")])
  # ## Dotplot 
  # dotplot(ego_CC,orderBy="GeneRatio")
  # ## Output results from GO analysis to a table
  # cluster_summary <- data.frame(ego_CC)
  # ## Add similarity matrix to the termsim slot of enrichment result
  # ego_CC <- enrichplot::pairwise_termsim(ego_CC)
  # ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
  # emapplot(ego_CC, showCategory = 50)
  # ## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
  # OE_foldchanges <- sigOE$log2FoldChange
  # 
  # names(OE_foldchanges) <- sigOE$gene
  # 
  # ## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
  # cnetplot(ego_CC, 
  #          categorySize="pvalue", 
  #          showCategory = 5, 
  #          foldChange=OE_foldchanges, 
  #          vertex.label.font=6)
  # 
  # ########### Molecular Function #########################
  # ego_MF <- enricher(gene = sigOE_genes, 
  #                    universe = allOE_genes,
  #                    pAdjustMethod = "BH", 
  #                    qvalueCutoff = 0.05,
  #                    TERM2GENE = term2gene_MF[,c("goName","gene")])
  # ## Dotplot 
  # dotplot(ego_MF,orderBy="GeneRatio")
  # ## Output results from GO analysis to a table
  # cluster_summary <- data.frame(ego_MF)
  # ## Add similarity matrix to the termsim slot of enrichment result
  # ego_MF <- enrichplot::pairwise_termsim(ego_MF)
  # ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
  # emapplot(ego_MF, showCategory = 50)
  # ## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
  # OE_foldchanges <- sigOE$log2FoldChange
  # 
  # names(OE_foldchanges) <- sigOE$gene
  # 
  # ## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
  # cnetplot(ego_MF, 
  #          categorySize="pvalue", 
  #          showCategory = 5, 
  #          foldChange=OE_foldchanges, 
  #          vertex.label.font=6)
}

GO_results(res_R52_vs_WT_irr_vs_ctl,"res_R52_vs_WT_irr_vs_ctl")
GO_results(res_8115_vs_WT_irr_vs_ctl,"res_8115_vs_WT_irr_vs_ctl")
GO_results(res_R52_vs_8115_irr,"res_R52_vs_8115_irr") 
######################################################################
###### GSEA:coordinated differential expression over gene sets #######
###### is tested instead of changes of individual genes ##############
######################################################################


# GSEA using gene sets associated with BP Gene Ontology terms
gseaGO <- gseGO(geneList = foldchanges, 
                OrgDb = test,
                ont = 'BP', 
                nPerm = 1000, 
                minGSSize = 20, 
                pvalueCutoff = 0.05,
                verbose = FALSE) 

gseaGO_results <- gseaGO@result



# # Apply fold change shrinkage (apeglm less bias than normal)
# res_tableOE <- lfcShrink(dds, coef="IsolateR52.conditionIrradiated", type="apeglm")
# plotMA(res_tableOE, ylim=c(-2,2))
# ## Summarize results
# summary(res_tableOE, alpha = 0.05)
# ##Volcano plot (DEGreport)
# DEGreport::degPlot(dds = dds, res = res, n = 20, xs = "type", group = "condition") # dds object is output from DESeq2
# res$id<-as.character(gsub("jgi.p.*Exode1.","",row.names(res)))
# DEGreport::degVolcano(
#   data.frame(res[,c("log2FoldChange","padj")]), # table - 2 columns
#   plot_text = data.frame(res[1:10,c("log2FoldChange","padj","id")])) # table to add names
# # Available in the newer version for R 3.4
# DEGreport::degPlotWide(dds, genes = row.names(res)[1:5], group = "condition")

######################################################################
### Likelihood ratio test: the ‘reduced’ model (removing factor(s))###
######################################################################
# dds_lrt <- DESeq(dds, test="LRT", reduced = ~ Isolate+condition) 
# #dds_lrt <- DESeq(dds, test="LRT", reduced = ~ Isolate+Isolate:condition)
# #dds_lrt <- DESeq(dds, test="LRT", reduced = ~ condition+Isolate:condition)
# 
# dds_lrt <- DESeq(dds, test="LRT", reduced = ~ Isolate)
# dds_lrt <- DESeq(dds, test="LRT", reduced = ~ condition)
# dds_lrt <- DESeq(dds, test="LRT", reduced = ~ Isolate:condition)
# dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)
# # Extract results for LRT
# res_LRT <- results(dds_lrt)
# 
# res_LRT_tb <- res_LRT %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>% 
#   as_tibble()
# 
# # Subset to return genes with padj < 0.05
# sigLRT_genes <- res_LRT_tb %>% 
#   filter(padj < padj.cutoff)
# 
# # Get number of significant genes
# nrow(sigLRT_genes)










