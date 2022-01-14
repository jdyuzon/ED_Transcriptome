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
#BiocManager::install("WGCNA")
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
library(VennDiagram)
library(WGCNA)
library(enrichplot)

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
row.names(cts)<-as.character(gsub("jgi.p.*Exode1.","",row.names(cts)))
### Read in Gene Names and format
genenames<-read.table("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/transcript_genenames.txt",
                      header=TRUE)
colnames(genenames)<-c("gene","gene_name")
genenames$gene<-as.character(genenames$gene)

genenames_order<-genenames[,"gene_name"][match(rownames(cts),genenames[,"gene"])] 
###check missing data
rownames(cts)[is.na(genenames_order)]
### rename RAD52, 8115, PKS transcripts and genes, and rename 1000 transcript
genenames$gene[genenames$gene=="2864"]<-"2864_rad52" 
genenames$gene[genenames$gene=="3244"]<-"3244_pks1"  
genenames$gene[genenames$gene=="8372"]<-"8372_8115"

genenames$gene_name[genenames$gene=="2864_rad52"]<-"rad52" 
genenames$gene_name[genenames$gene=="3244_pks1"]<-"pks1"  
genenames$gene_name[genenames$gene=="8372_8115"]<-"8115"

genenames_order<-genenames[,"gene_name"][match(rownames(cts),genenames[,"gene"])] 
genenames_order[is.na(genenames_order)]<-"transcript_1000"

### Replace transcript IDs with gene names
rownames(cts)<-genenames_order

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
###This would say, e.g. filter out genes where there are less than 3 samples with normalized counts greater than or equal to 5.
dds <- estimateSizeFactors(dds) #estimates the size factors using the "median ratio method" 
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]

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
  geom_point(size =3, aes(fill=factor(Isolate))) + 
  scale_shape_manual(values=1:nlevels(pcaData$Isolate))

pcaData<-plotPCA(vsd, intgroup=c("condition", "genotype"),returnData=TRUE)
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition), shape = factor(genotype))) + 
  geom_point(size =3, aes(fill=factor(genotype), alpha=as.character(genotype)))  

pcaData<-plotPCA(vsd, intgroup=c("condition", "Sample"),returnData=TRUE)
pcaData$name<- sub("_.*", "", pcaData$name)
pcaData$name<- sub("PKS C.*", "PKS C1", pcaData$name)
pdf("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/ED_transcriptome_WT_PKS_R52_8115.pca.pdf", width=14,height=14)
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition))) + 
  geom_text(aes(label=Sample, alpha=0.5),hjust=0, vjust=0)
graphics.off()

##Two WT15.1.2 libraries (one irradiated WT15.1.2_3 and one not irradiated WT15.1.2_6) do not cluster with other WT but with R52, toss out
##WTC1_4 and WT15.4.1_4 Control conditions cluster with Irradiated
idx<-which(colnames(vsd)=="WT C1_4" |
             colnames(vsd)=="WT15.1.2_3" |
             colnames(vsd)=="WT15.1.2_6" |
             colnames(vsd)=="WT15.4.1_4")

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
pdf("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/ED_transcriptome_WT_PKS.pca.pdf", width=14,height=14)
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition))) + 
  geom_text(aes(label=Sample,alpha=0.5),hjust=0, vjust=0,size = 3)+
  theme(legend.title = element_blank()) 
graphics.off()

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
pdf("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/ED_transcriptome_R52_8115.pca.pdf", width=14,height=14)
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition))) + 
  geom_text(aes(label=Sample),hjust=0, vjust=0)
graphics.off()


###############################################################
#####Dispersion during differential expression analysis########
###############################################################
####estimates of variation for each gene are often unreliable, 
####To address this problem, DESeq2 shares information across genes to generate more accurate 
####estimates of variation based on the mean expression level of the gene using a method called ‘shrinkage’. 
####DESeq2 assumes that genes with similar expression levels should have similar dispersion.
idx<-which(colnames(vsd)=="WT C1_4" |
             colnames(vsd)=="WT15.1.2_3" |
             colnames(vsd)=="WT15.1.2_6" |
             colnames(vsd)=="WT15.4.1_4")
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
res_8115_vs_R52_irr <- results(dds, contrast=list("Isolate8115.conditionIrradiated","IsolateR52.conditionIrradiated"),alpha = 0.05)  ####same results

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

### Read in GO annotation and format
ED_gff_GO<-read.csv("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/Exode1_GeneCatalog_proteins_20160902_GO.genename.csv",
                     header=TRUE)
colnames(ED_gff_GO)<-c("gene","gotermId","goName", "keytypes","goAcc")
ED_gff_GO$gene<-as.character(gsub(".*Exode1.","",ED_gff_GO$gene))
ED_gff_GO<-left_join(ED_gff_GO, genenames, by=c("gene"="gene")) 

ED_gff_GO$keytypes[ED_gff_GO$keytypes=="biological_process"]<-"BP"
ED_gff_GO$keytypes[ED_gff_GO$keytypes=="cellular_component"]<-"CC"
ED_gff_GO$keytypes[ED_gff_GO$keytypes=="molecular_function"]<-"MF"

ED_gff_GO$gene_name[ED_gff_GO$gene=="2864"]<-"rad52" 
ED_gff_GO$gene_name[ED_gff_GO$gene=="3244"]<-"pks1"  
ED_gff_GO$gene_name[ED_gff_GO$gene=="8372"]<-"8115"

### Read in Interproscan annotation
ED_gff_IPR<-read.csv("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/Exode1_GeneCatalog_proteins_20160902_IPR.csv",
                      header=TRUE)
colnames(ED_gff_IPR)[1]<-"gene"
ED_gff_IPR$gene<-as.character(ED_gff_IPR$gene)
ED_gff_IPR<-left_join(ED_gff_IPR, genenames, by=c("gene"="gene")) 

ED_gff_IPR$gene_name[ED_gff_IPR$gene=="2864"]<-"rad52" 
ED_gff_IPR$gene_name[ED_gff_IPR$gene=="3244"]<-"pks1"  
ED_gff_IPR$gene_name[ED_gff_IPR$gene=="8372"]<-"8115"
### Read in KEGG annotation
ED_gff_KEGG<-read.csv("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/Exode1_GeneCatalog_proteins_20160902_KEGG.csv",
                     header=TRUE)
colnames(ED_gff_KEGG)[1]<-"gene"
ED_gff_KEGG$gene<-as.character(ED_gff_KEGG$gene)
ED_gff_KEGG<-left_join(ED_gff_KEGG, genenames, by=c("gene"="gene")) 

ED_gff_KEGG$gene_name[ED_gff_KEGG$gene=="2864"]<-"rad52" 
ED_gff_KEGG$gene_name[ED_gff_KEGG$gene=="3244"]<-"pks1"  
ED_gff_KEGG$gene_name[ED_gff_KEGG$gene=="8372"]<-"8115"
### Read in KOG annoation
ED_gff_KOG<-read.csv("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/Exode1_GeneCatalog_proteins_20160902_KOG.csv",
                     header=TRUE)
colnames(ED_gff_KOG)[1]<-"gene"
ED_gff_KOG$gene<-as.character(ED_gff_KOG$gene)
ED_gff_KOG<-left_join(ED_gff_KOG, genenames, by=c("gene"="gene")) 

ED_gff_KOG$gene_name[ED_gff_KOG$gene=="2864"]<-"rad52" 
ED_gff_KOG$gene_name[ED_gff_KOG$gene=="3244"]<-"pks1"  
ED_gff_KOG$gene_name[ED_gff_KOG$gene=="8372"]<-"8115"

GO_results<-function(res,name){
  ## assign up/down regulation
  res$regulation<-NULL
  res$regulation[res$log2FoldChange<0]<-"down"
  res$regulation[res$log2FoldChange>0]<-"up"
  #Create a tibble of results (all the genes with transcripts)
  row.names(res)<-as.character(gsub("jgi.p.*Exode1.","",row.names(res)))
  res_tableOE_tb <- res%>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  res_tableOE_tb <- res_tableOE_tb %>% 
    mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)
  ## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
  allOE_genes <- as.character(res_ids$gene)
  ## Extract significant results
  sigOE <- res_tableOE_tb %>%
    filter(padj < 0.05)
  ## Volcano plot
  sigOE_genes <- as.character(sigOE$gene)
  
  term2gene<-ED_gff_GO[,c("goName","gene_name", "keytypes")]
  #term2gene_BP<-term2gene[which(term2gene$keytypes=="BP"),]
  #term2gene_CC<-term2gene[which(term2gene$keytypes=="CC"),]
  #term2gene_MF<-term2gene[which(term2gene$keytypes=="MF"),]
  
  ########### GO Enrichment #########################
  ego <- enricher(gene = sigOE_genes, 
                     universe = allOE_genes,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05,
                     TERM2GENE = term2gene[,c("goName","gene_name")])
  ## Dotplot 
  ego_up <- enricher(gene = as.character(sigOE$gene[sigOE$regulation=="up"]), 
                       universe = allOE_genes,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       TERM2GENE = term2gene[,c("goName","gene_name")])
  a1<-dotplot(ego_up,orderBy="GeneRatio")+ggtitle("activated")
  ego_down <- enricher(gene = as.character(sigOE$gene[sigOE$regulation=="down"]), 
                         universe = allOE_genes,
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         TERM2GENE = term2gene[,c("goName","gene_name")])
  a2<-dotplot(ego_down,orderBy="GeneRatio")+ggtitle("suppressed")
  a3<-ggarrange(a1,a2,ncol=2,common.legend = TRUE)
  #a<-dotplot(ego,orderBy="GeneRatio")
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
           node_label_size=2)
  pdf(paste0("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/",name,".GO.pdf"),width=14,height=18)
  print(ggarrange(a3,b,c,nrow = 3,ncol=1))
  graphics.off()

  ########### KEGG pathway analysis ######################### 
  term2gene<-ED_gff_KEGG[,c("gene_name", "pathway", "pathway_class")]
  ekegg <- enricher(gene = sigOE_genes, 
                  universe = allOE_genes,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  TERM2GENE = term2gene[,c("pathway","gene_name")])

  ekegg_up <- enricher(gene = as.character(sigOE$gene[sigOE$regulation=="up"]), 
                     universe = allOE_genes,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05,
                     TERM2GENE = term2gene[,c("pathway","gene_name")])
  a<-dotplot(ekegg_up,orderBy="GeneRatio")
  ekegg_down <- enricher(gene = as.character(sigOE$gene[sigOE$regulation=="down"]), 
                       universe = allOE_genes,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       TERM2GENE = term2gene[,c("pathway","gene_name")])
  b<-dotplot(ekegg_down,orderBy="GeneRatio")
  c<-cnetplot(ekegg, 
              categorySize="pvalue", 
              showCategory = 5, 
              foldChange=OE_foldchanges, 
              vertex.label.font=6)
  pdf(paste0("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/",name,".KEGG.pdf"),width=14,height=18)
  print(ggarrange(a,b,c)) 
  graphics.off()
  ########### KOG Enrichment #########################
  term2gene<-ED_gff_KOG[,c("gene_name", "kogdefline", "kogClass")]
  ekog <- enricher(gene = sigOE_genes, 
                    universe = allOE_genes,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    TERM2GENE = term2gene[,c("kogdefline","gene_name")])
  
  ekog_up <- enricher(gene = as.character(sigOE$gene[sigOE$regulation=="up"]), 
                       universe = allOE_genes,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       TERM2GENE = term2gene[,c("kogdefline","gene_name")])
  a<-dotplot(ekog_up,orderBy="GeneRatio")
  ekog_down <- enricher(gene = as.character(sigOE$gene[sigOE$regulation=="down"]), 
                         universe = allOE_genes,
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         TERM2GENE = term2gene[,c("kogdefline","gene_name")])
  b<-dotplot(ekog_down,orderBy="GeneRatio")
  c<-cnetplot(ekog, 
              categorySize="pvalue", 
              showCategory = 5, 
              foldChange=OE_foldchanges, 
              vertex.label.font=6)
  pdf(paste0("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/",name,".KOG.pdf"),width=14,height=18)
  print(ggarrange(a,b,c)) 
  graphics.off()
  
  #### Table of all the data
  ## Merge the AnnotationHub dataframe with the results 
  ed_gff_go<-mutate(ED_gff_GO[c("gene_name","goName")])%>%
    group_by(gene_name)%>% summarise(goName = paste0(goName, collapse = ";")) %>%
    ungroup()
  ed_gff_kegg<-mutate(ED_gff_KEGG[c("gene_name","pathway")])%>%
    group_by(gene_name)%>% summarise(pathway = paste0(pathway, collapse = ";")) %>%
    ungroup()
  ed_gff_kog<-mutate(ED_gff_KOG[c("gene_name","kogdefline")])%>%
    group_by(gene_name)%>% summarise(kogdefline = paste0(kogdefline, collapse = ";")) %>%
    ungroup()
  ed_gff_ipr<-mutate(ED_gff_IPR[c("gene_name","iprDesc")],)%>%
    group_by(gene_name)%>% summarise(iprDesc = paste0(iprDesc, collapse = ";")) %>%
    ungroup()
  #### Annotation of all GO,KEGG,KOG,IPR
  res_ids <- left_join(res_tableOE_tb, ed_gff_go, by=c("gene"="gene_name"))  
  res_ids <- left_join(res_ids, ed_gff_kegg, by=c("gene"="gene_name")) 
  res_ids <- left_join(res_ids, ed_gff_kog, by=c("gene"="gene_name")) 
  res_ids <- left_join(res_ids, ed_gff_ipr, by=c("gene"="gene_name")) 
  res_sigOEids <- res_ids %>%
    filter(padj < 0.05)
  write.csv(as.data.frame(res_sigOEids), 
            file=paste0("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/",name,".GO_KEGG_KOG_IPR_annotation.csv"))
  ####
}

GO_results(res_R52_vs_WT_irr_vs_ctl,"res_R52_vs_WT_irr_vs_ctl") ##### no KEGG enriched results
GO_results(res_8115_vs_WT_irr_vs_ctl,"res_8115_vs_WT_irr_vs_ctl")
GO_results(res_R52_vs_8115_irr,"res_R52_vs_8115_irr") 
GO_results(res_R52_vs_WT_ctl,"res_R52_vs_WT_ctl")
GO_results(res_8115_vs_WT_ctl,"res_8115_vs_WT_ctl")

######################################################################
###### WGCNA:weighted correlation network analysis ###################
###### weighted gene co-expression network analysis ##################
######################################################################

###From: https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#45_Format_normalized_data_for_WGCNA
options(stringsAsFactors = FALSE)
###WGCNA requires at least 15 samples to produce a meaningful result
# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)

# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) %>%
  t() # Transpose this data
##power parameter will affect the number of modules identified, and the WGCNA modules provides the pickSoftThreshold() function to help identify good choices for this parameter
sft <- pickSoftThreshold(normalized_counts,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

###removing features whose counts are consistently low (for example, removing all features that have a count of less than say 10 in more than 90% of the samples) because such low-expressed features tend to reflect noise 










######################################################################
####### Venn Diagram: Shared and Unique Genes ########################
######################################################################
# Calculate the intersection of the two sets

test<-res_R52_vs_WT_irr_vs_ctl[!is.na(res_R52_vs_WT_irr_vs_ctl$padj),]
dim(test[test$padj<0.05,]) #1852

test<-res_8115_vs_WT_irr_vs_ctl[!is.na(res_8115_vs_WT_irr_vs_ctl$padj),]
dim(test[test$padj<0.05,]) #2152

test<-res_R52_vs_8115_irr[!is.na(res_R52_vs_8115_irr$padj),]
test<-test[test$padj<0.05,] #245

!res_sigOE_R52_vs_WT_irr_vs_ctl$gene %in% row.names(test)
!row.names(test) %in% res_sigOE_R52_vs_WT_irr_vs_ctl$gene


res_sigOE_R52_vs_WT_irr_vs_ctl <- as.data.frame(res_R52_vs_WT_irr_vs_ctl[!is.na(res_R52_vs_WT_irr_vs_ctl$padj),])
res_sigOE_R52_vs_WT_irr_vs_ctl <- res_sigOE_R52_vs_WT_irr_vs_ctl[res_sigOE_R52_vs_WT_irr_vs_ctl$padj<0.05,] #1852

res_sigOE_8115_vs_WT_irr_vs_ctl <- as.data.frame(res_8115_vs_WT_irr_vs_ctl[!is.na(res_8115_vs_WT_irr_vs_ctl$padj),])
res_sigOE_8115_vs_WT_irr_vs_ctl <- res_sigOE_8115_vs_WT_irr_vs_ctl[res_sigOE_8115_vs_WT_irr_vs_ctl$padj<0.05,] #2152

res_sigOE_R52_vs_WT_irr_vs_ctl$gene<-row.names(res_sigOE_R52_vs_WT_irr_vs_ctl)
res_sigOE_8115_vs_WT_irr_vs_ctl$gene<-row.names(res_sigOE_8115_vs_WT_irr_vs_ctl)

deg.intersect = length(intersect(res_sigOE_R52_vs_WT_irr_vs_ctl$gene, res_sigOE_8115_vs_WT_irr_vs_ctl$gene)) 
deg.venn <- list('intersect' = deg.intersect,
                 'irr_R52' = length(row.names(res_sigOE_R52_vs_WT_irr_vs_ctl)),
                 'irr_8115' = length(row.names(res_sigOE_8115_vs_WT_irr_vs_ctl)))

# Arguments for a pairwise (two-sets) venn-diagram are sizes for set1, set2 and overlap (intersect)
# Many more functions are available for triple, quad and quantuple diagrams (starting with 'draw.***')
venn.plot <- draw.pairwise.venn(deg.venn$irr_R52, deg.venn$irr_8115, deg.venn$intersect,
                                category = c("irr_R52", "irr_8115"), scaled = F,
                                fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                                cat.pos = c(0, 0))

# Actually plot the plot
grid.draw(venn.plot)

#### Create a table of DES common between expression profiles

res_common_R52_8115_irr_vs_ctl<-merge(as.data.frame(res_sigOE_R52_vs_WT_irr_vs_ctl), as.data.frame(res_sigOE_8115_vs_WT_irr_vs_ctl), by=0) 
colnames(res_common_R52_8115_irr_vs_ctl)[1]<-"gene"
res_common_R52_8115_irr_vs_ctl$gene<-as.character(gsub(".*Exode1.","",res_common_R52_8115_irr_vs_ctl$gene))
res_common_R52_8115_irr_vs_ctl<-inner_join(res_common_R52_8115_irr_vs_ctl, ED_gff_GO, by=c("gene"="gene")) ###length of rows increases because multiple GO terms
write.csv(as.data.frame(res_common_R52_8115_irr_vs_ctl), 
          file=paste0("~/Desktop/NRL_fungalprojects/JGI/JGI_transcriptome/res_common_R52_8115_irr_vs_ctl.GOterms.csv"))

module_19_df <- module_eigengenes %>%
  tibble::rownames_to_column("accession_code") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metadata %>%
                      dplyr::select(Sample, condition),
                    by = c("accession_code" = "Sample")
  )
ggplot(
  module_19_df,
  aes(
    x = condition,
    y = ME19,
    color = condition
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()













