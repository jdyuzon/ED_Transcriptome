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
#BiocManager::install("biomaRt")
#install.packages("ggplot2")
#install.packages("ggrepel")
#install.packages("ggnewscale")
#install.packages("ggpubr")
#install.packages("VennDiagram")
#install.packages("gplots")

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
library(tidyr)
library('igraph')
library(tidyverse)
library( "genefilter" )
library("gplots")
library("ggplotify")

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/NRL_Postdoc/ED_Transcriptome_evol/")


###Count matrix input
cts <- as.matrix(read.csv("counts.txt",sep="\t",row.names="GeneID"))
coldata<- read.csv("sample_annotation_evol.csv")
colnames(cts)[match(coldata[,"Index"], colnames(cts))] = coldata[,"Sample"]

##format sample annotation file
rownames(coldata)<-coldata$Sample
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts)) ###not in the same order
coldata$Transfers<-factor(coldata$Transfers)
coldata$Group<-gsub("15..*","15",coldata$genotype)

cts <- cts[, rownames(coldata)] ###order the counts data by sample annotation file order
all(rownames(coldata) == colnames(cts)) ###in same order
row.names(cts)<-as.character(gsub("jgi.p.*Exode1.","",row.names(cts)))
### Read in Gene Names and format
genenames<-read.table("transcript_genenames.txt",
                      header=FALSE)
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

### Read in GO annotation and format
ED_gff_GO<-read.csv("../ED_Transcriptome/Exode1_GeneCatalog_proteins_20160902_GO.genename.csv",
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
ED_gff_IPR<-read.csv("../ED_Transcriptome/Exode1_GeneCatalog_proteins_20160902_IPR.csv",
                     header=TRUE)
colnames(ED_gff_IPR)[1]<-"gene"
ED_gff_IPR$gene<-as.character(ED_gff_IPR$gene)
ED_gff_IPR<-left_join(ED_gff_IPR, genenames, by=c("gene"="gene")) 

ED_gff_IPR$gene_name[ED_gff_IPR$gene=="2864"]<-"rad52" 
ED_gff_IPR$gene_name[ED_gff_IPR$gene=="3244"]<-"pks1"  
ED_gff_IPR$gene_name[ED_gff_IPR$gene=="8372"]<-"8115"
### Read in KEGG annotation
ED_gff_KEGG<-read.csv("../ED_Transcriptome/Exode1_GeneCatalog_proteins_20160902_KEGG.csv",
                      header=TRUE)
colnames(ED_gff_KEGG)[1]<-"gene"
ED_gff_KEGG$gene<-as.character(ED_gff_KEGG$gene)
ED_gff_KEGG<-left_join(ED_gff_KEGG, genenames, by=c("gene"="gene")) 

ED_gff_KEGG$gene_name[ED_gff_KEGG$gene=="2864"]<-"rad52" 
ED_gff_KEGG$gene_name[ED_gff_KEGG$gene=="3244"]<-"pks1"  
ED_gff_KEGG$gene_name[ED_gff_KEGG$gene=="8372"]<-"8115"
### Read in KOG annoation
ED_gff_KOG<-read.csv("../ED_Transcriptome/Exode1_GeneCatalog_proteins_20160902_KOG.csv",
                     header=TRUE)
colnames(ED_gff_KOG)[1]<-"gene"
ED_gff_KOG$gene<-as.character(ED_gff_KOG$gene)
ED_gff_KOG<-left_join(ED_gff_KOG, genenames, by=c("gene"="gene")) 

ED_gff_KOG$gene_name[ED_gff_KOG$gene=="2864"]<-"rad52" 
ED_gff_KOG$gene_name[ED_gff_KOG$gene=="3244"]<-"pks1"  
ED_gff_KOG$gene_name[ED_gff_KOG$gene=="8372"]<-"8115"

term2gene_go<-ED_gff_GO[,c("goName","gene_name", "keytypes")]
term2gene_kegg<-ED_gff_KEGG[,c("gene_name", "pathway", "pathway_class")]
term2gene_kog<-ED_gff_KOG[,c("gene_name", "kogdefline", "kogClass")]

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
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
###This would say, e.g. filter out genes where there are less than 3 samples with normalized counts greater than or equal to 5.
dds <- estimateSizeFactors(dds) #estimates the size factors using the "median ratio method" 
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]

###indicate the comparison (WT0 is for 8115 and R52, WTC is for evolution experiment)
dds$condition <- factor(dds$condition, levels = c("Control","Irradiated"))
dds$Isolate <- factor(dds$Isolate, 
                      levels = c("WT.C1","PKS.C1",
                                 "WT15.1.2","WT15.2.1","WT15.2.2","WT15.4.1","WT15.4.2",
                                 "PKS15.1.1","PKS15.2.1","PKS15.2.2","PKS15.3.2","PKS15.4.2","WT0",
                                 "8115","R52"))
dds$genotype<- factor(dds$genotype, levels = c("WT","PKS",
                                               "WT15.1.2","WT15.2.1","WT15.2.2","WT15.4.1","WT15.4.2",
                                               "PKS15.1.1","PKS15.2.1","PKS15.2.2","PKS15.3.2","PKS15.4.2",
                                               "8115","R52"))


##############################################
###Data Quality Assessment and Exploration####
##############################################
####transformation on variance
ntd <- normTransform(dds)
vsd <- vst(dds)
rld <- rlog(dds)
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

pcaData<-plotPCA(vsd, intgroup=c("condition", "Isolate"),returnData=TRUE)
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition), shape = factor(Isolate))) + 
  geom_point(size =3, aes(fill=factor(Isolate), alpha=as.character(Isolate)))  

pcaData<-plotPCA(vsd, intgroup=c("condition", "Sample"),returnData=TRUE)
pcaData$name<- sub("_.*", "", pcaData$name)
pcaData$name<- sub("PKS C.*", "PKS C1", pcaData$name)

##Two WT15.1.2 libraries (one irradiated WT15.1.2_3 and one not irradiated WT15.1.2_6) do not cluster with other WT but with R52, toss out
##WTC1_4 and WT15.4.1_4 Control conditions cluster with Irradiated
idx<-which(colnames(vsd)=="WTC1_4" |
             colnames(vsd)=="WT15.1.2_3" |
             colnames(vsd)=="WT15.1.2_6" |
             colnames(vsd)=="WT15.4.1_4")

vsd2<-vsd[,-idx]

##create DESeqDataSet for only WT and PKS (evolution dataset)
idx_evol<-which(grepl("WT.*|PKS.*", colnames(vsd2)))
vsd_evol<-vsd2[,idx_evol]
vsd_evol$condition <- factor(vsd_evol$condition, levels = c("Irradiated", "Control"))

pheatmap(assay(vsd_evol), cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)

sampleDists <- dist(t(assay(vsd_evol)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_evol$condition, vsd_evol$Isolate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
evol_heatmap<-pheatmap(sampleDistMatrix,
                       clustering_distance_rows=sampleDists,
                       clustering_distance_cols=sampleDists,
                       col=colors)

pcaData<-plotPCA(vsd_evol, intgroup=c("condition", "Sample"),returnData=TRUE)
pcaData$name<- sub("_.*", "", pcaData$name)
pcaData$name<- sub("PKS C.*", "PKS C1", pcaData$name)
pdf("output/ED_transcriptome_WT_PKS.pca.pdf", width=14,height=14)
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition))) + 
  geom_text(aes(label=Sample,alpha=0.5),hjust=0, vjust=0,size = 3)+
  theme(legend.title = element_blank()) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
graphics.off()


pcaData$name_color<- sub("PKS.*", "PKS", pcaData$name)
pcaData$name_color<- sub("WT.*", "WT", pcaData$name_color)
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(name_color))) + 
  geom_text(aes(label=Sample,alpha=0.5),hjust=0, vjust=0,size = 3)+
  theme(legend.title = element_blank()) +
  coord_flip()

percentVar <- round(100 * attr(pcaData, "percentVar"))
evol_pca<-ggplot(pcaData, aes(x = PC1, y = PC2, label =Sample)) + 
  geom_text_repel(aes(color=factor(condition)),segment.color = 'transparent',
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 50))+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
###############################################################
#####Dispersion during differential expression analysis########
###############################################################
####estimates of variation for each gene are often unreliable, 
####To address this problem, DESeq2 shares information across genes to generate more accurate 
####estimates of variation based on the mean expression level of the gene using a method called ‘shrinkage’. 
####DESeq2 assumes that genes with similar expression levels should have similar dispersion.
idx<-which(colnames(vsd)=="WTC1_4" |
             colnames(vsd)=="WT15.1.2_3" |
             colnames(vsd)=="WT15.1.2_6" |
             colnames(vsd)=="WT15.4.1_4")
cts2<-cts[,-idx]
cts_evol<-cts2[,grepl("WT.*|PKS.*", colnames(cts2))]

dds_evol <- DESeqDataSetFromMatrix(countData = cts_evol,
                                   colData = coldata[colnames(cts_evol),],
                                   design = ~Group+condition+Group:condition)

dds_evol_individual<- DESeqDataSetFromMatrix(countData = cts_evol,
                                             colData = coldata[colnames(cts_evol),],
                                             design =  ~Isolate+condition+Isolate:condition)

###indicate the comparison (WT0 is for 8115 and R52, WTC is for evolution experiment)
dds_evol$condition <- factor(dds_evol$condition, levels = c("Control","Irradiated"))
dds_evol$Group <- factor(dds_evol$Group, levels = c("WT","PKS","WT15","PKS15"))
dds_evol$Isolate <- factor(dds_evol$Isolate, 
                           levels = c("WT.C1","PKS.C1",
                                      "WT15.1.2","WT15.2.1","WT15.2.2","WT15.4.1","WT15.4.2","WT0",
                                      "PKS15.1.1","PKS15.2.1","PKS15.2.2","PKS15.3.2","PKS15.4.2"))
dds_evol_individual$condition <- factor(dds_evol_individual$condition, levels = c("Control","Irradiated"))
dds_evol_individual$Group <- factor(dds_evol_individual$Group, levels = c("WT","PKS","WT15","PKS15"))
dds_evol_individual$Isolate <- factor(dds_evol_individual$Isolate, 
                           levels = c("WT.C1","PKS.C1",
                                      "WT15.1.2","WT15.2.1","WT15.2.2","WT15.4.1","WT15.4.2","WT0",
                                      "PKS15.1.1","PKS15.2.1","PKS15.2.2","PKS15.3.2","PKS15.4.2"))


dds <- DESeq(dds)
dds_evol <- DESeq(dds_evol)
dds_evol_individual <- DESeq(dds_evol_individual)

plotDispEsts(dds)
plotDispEsts(dds_evol)
plotDispEsts(dds_evol_individual)
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
resultsNames(dds_evol)
res_WT_irr_vs_ctl = results(dds_evol, contrast=c("condition","Irradiated","Control"),alpha = 0.05)
res_PKS_irr_vs_ctl <- results(dds_evol, list( c("condition_Irradiated_vs_Control","GroupPKS.conditionIrradiated") ),alpha = 0.05)
res_WT15common_irr_vs_ctl<- results(dds_evol, list( c("condition_Irradiated_vs_Control","GroupWT15.conditionIrradiated") ),alpha = 0.05)
res_PKS15common_irr_vs_ctl<- results(dds_evol, list( c("condition_Irradiated_vs_Control","GroupPKS15.conditionIrradiated") ),alpha = 0.05)

res_PKS15common_vs_PKS_irr_vs_ctl<-results(dds_evol, contrast= c(0,-1,0,1,
                                                                 0,-1,0,1),alpha = 0.05)
res_PKS15common_vs_WT15common_irr_vs_ctl<-results(dds_evol, contrast= c(0,0,-1,1,
                                                                        0,0,-1,1),alpha = 0.05)

res_WT15common_vs_WT_irr <-results(dds_evol, list( c("Group_WT15_vs_WT","GroupWT15.conditionIrradiated") ),alpha = 0.05)
res_WT15common_vs_WT_ctl<-results(dds_evol, contrast=c("Group","WT15","WT"),alpha = 0.05) 
res_PKS15common_vs_PKS_irr<-results(dds_evol, contrast=list("GroupPKS15.conditionIrradiated","GroupPKS.conditionIrradiated"),alpha = 0.05)
res_PKS15common_vs_PKS_ctl<-results(dds_evol, contrast=c("Group","PKS15","PKS"),alpha = 0.05) 
res_PKS15_vs_WT15_irr<-results(dds_evol, contrast=list("GroupPKS15.conditionIrradiated","GroupWT15.conditionIrradiated"),alpha = 0.05)
res_PKS15_vs_WT15_ctl<-results(dds_evol, contrast=c("Group","PKS15","WT15"),alpha = 0.05) 

resultsNames(dds_evol_individual)
res_WT15.1.2_vs_WT_irr <- results(dds_evol_individual, list( c("Isolate_WT15.1.2_vs_WT.C1","IsolateWT15.1.2.conditionIrradiated")),alpha = 0.05)  
res_WT15.2.1_vs_WT_irr <- results(dds_evol_individual, list( c("Isolate_WT15.2.1_vs_WT.C1","IsolateWT15.2.1.conditionIrradiated")),alpha = 0.05)  
res_WT15.2.2_vs_WT_irr <- results(dds_evol_individual, list( c("Isolate_WT15.2.2_vs_WT.C1","IsolateWT15.2.2.conditionIrradiated")),alpha = 0.05) 
res_WT15.4.1_vs_WT_irr <- results(dds_evol_individual, list( c("Isolate_WT15.4.1_vs_WT.C1","IsolateWT15.4.1.conditionIrradiated")),alpha = 0.05)  
res_WT15.4.2_vs_WT_irr <- results(dds_evol_individual, list( c("Isolate_WT15.4.2_vs_WT.C1","IsolateWT15.4.2.conditionIrradiated")),alpha = 0.05)
res_PKS15.1.1_vs_PKS_irr <- results(dds_evol_individual, contrast=list("Isolate_PKS15.1.1_vs_WT.C1","IsolatePKS15.1.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.1_vs_PKS_irr <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.1_vs_WT.C1","IsolatePKS15.2.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.2_vs_PKS_irr <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.2_vs_WT.C1","IsolatePKS15.2.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.3.2_vs_PKS_irr <- results(dds_evol_individual, contrast=list("Isolate_PKS15.3.2_vs_WT.C1","IsolatePKS15.3.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.4.2_vs_PKS_irr <- results(dds_evol_individual, contrast=list("Isolate_PKS15.4.2_vs_WT.C1","IsolatePKS15.4.2.conditionIrradiated"),alpha = 0.05)
res_WT15.1.2_vs_WT_ctl <- results(dds_evol_individual, contrast=c("Isolate","WT15.1.2","WT.C1"),alpha = 0.05)  
res_WT15.2.1_vs_WT_ctl <- results(dds_evol_individual, contrast=c("Isolate","WT15.2.1","WT.C1"),alpha = 0.05)  
res_WT15.2.2_vs_WT_ctl <- results(dds_evol_individual, contrast=c("Isolate","WT15.2.2","WT.C1"),alpha = 0.05) 
res_WT15.4.1_vs_WT_ctl <- results(dds_evol_individual, contrast=c("Isolate","WT15.4.1","WT.C1"),alpha = 0.05)  
res_WT15.4.2_vs_WT_ctl <- results(dds_evol_individual, contrast=c("Isolate","WT15.4.2","WT.C1"),alpha = 0.05)
res_PKS15.1.1_vs_PKS_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.1.1_vs_WT.C1","Isolate_PKS.C1_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.1_vs_PKS_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.1_vs_WT.C1","Isolate_PKS.C1_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.2_vs_PKS_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.2_vs_WT.C1","Isolate_PKS.C1_vs_WT.C1"),alpha = 0.05)
res_PKS15.3.2_vs_PKS_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.3.2_vs_WT.C1","Isolate_PKS.C1_vs_WT.C1"),alpha = 0.05)
res_PKS15.4.2_vs_PKS_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.4.2_vs_WT.C1","Isolate_PKS.C1_vs_WT.C1"),alpha = 0.05)

###Within differences of the evolved WT15 lines
res_WT15.1.2_vs_WT15.2.1_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.1.2.conditionIrradiated","IsolateWT15.2.1.conditionIrradiated"),alpha = 0.05)
res_WT15.1.2_vs_WT15.2.2_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.1.2.conditionIrradiated","IsolateWT15.2.2.conditionIrradiated"),alpha = 0.05)
res_WT15.1.2_vs_WT15.4.1_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.1.2.conditionIrradiated","IsolateWT15.4.1.conditionIrradiated"),alpha = 0.05)
res_WT15.1.2_vs_WT15.4.2_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.1.2.conditionIrradiated","IsolateWT15.4.2.conditionIrradiated"),alpha = 0.05)

res_WT15.2.1_vs_WT15.1.2_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.2.1.conditionIrradiated","IsolateWT15.1.2.conditionIrradiated"),alpha = 0.05)
res_WT15.2.1_vs_WT15.2.2_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.2.1.conditionIrradiated","IsolateWT15.2.2.conditionIrradiated"),alpha = 0.05)
res_WT15.2.1_vs_WT15.4.1_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.2.1.conditionIrradiated","IsolateWT15.4.1.conditionIrradiated"),alpha = 0.05)
res_WT15.2.1_vs_WT15.4.2_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.2.1.conditionIrradiated","IsolateWT15.4.2.conditionIrradiated"),alpha = 0.05)

res_WT15.2.2_vs_WT15.1.2_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.2.2.conditionIrradiated","IsolateWT15.1.2.conditionIrradiated"),alpha = 0.05)
res_WT15.2.2_vs_WT15.2.1_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.2.2.conditionIrradiated","IsolateWT15.2.1.conditionIrradiated"),alpha = 0.05)
res_WT15.2.2_vs_WT15.4.1_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.2.2.conditionIrradiated","IsolateWT15.4.1.conditionIrradiated"),alpha = 0.05)
res_WT15.2.2_vs_WT15.4.2_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.2.2.conditionIrradiated","IsolateWT15.4.2.conditionIrradiated"),alpha = 0.05)

res_WT15.4.1_vs_WT15.1.2_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.4.1.conditionIrradiated","IsolateWT15.1.2.conditionIrradiated"),alpha = 0.05)
res_WT15.4.1_vs_WT15.2.1_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.4.1.conditionIrradiated","IsolateWT15.2.1.conditionIrradiated"),alpha = 0.05)
res_WT15.4.1_vs_WT15.2.2_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.4.1.conditionIrradiated","IsolateWT15.2.2.conditionIrradiated"),alpha = 0.05)
res_WT15.4.1_vs_WT15.4.2_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.4.1.conditionIrradiated","IsolateWT15.4.2.conditionIrradiated"),alpha = 0.05)

res_WT15.4.2_vs_WT15.1.2_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.4.2.conditionIrradiated","IsolateWT15.1.2.conditionIrradiated"),alpha = 0.05)
res_WT15.4.2_vs_WT15.2.1_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.4.2.conditionIrradiated","IsolateWT15.2.1.conditionIrradiated"),alpha = 0.05)
res_WT15.4.2_vs_WT15.2.2_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.4.2.conditionIrradiated","IsolateWT15.2.2.conditionIrradiated"),alpha = 0.05)
res_WT15.4.2_vs_WT15.4.1_irr <- results(dds_evol_individual, contrast=list("IsolateWT15.4.2.conditionIrradiated","IsolateWT15.4.1.conditionIrradiated"),alpha = 0.05)

###Within differences of the evolved PKS15 lines
res_PKS15.1.1_vs_PKS15.2.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.1.1.conditionIrradiated","IsolatePKS15.2.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.1.1_vs_PKS15.2.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.1.1.conditionIrradiated","IsolatePKS15.2.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.1.1_vs_PKS15.3.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.1.1.conditionIrradiated","IsolatePKS15.3.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.1.1_vs_PKS15.4.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.1.1.conditionIrradiated","IsolatePKS15.4.2.conditionIrradiated"),alpha = 0.05)

res_PKS15.2.1_vs_PKS15.1.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.1.conditionIrradiated","IsolatePKS15.1.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.1_vs_PKS15.2.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.1.conditionIrradiated","IsolatePKS15.2.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.1_vs_PKS15.3.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.1.conditionIrradiated","IsolatePKS15.3.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.1_vs_PKS15.4.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.1.conditionIrradiated","IsolatePKS15.4.2.conditionIrradiated"),alpha = 0.05)

res_PKS15.2.2_vs_PKS15.1.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.2.conditionIrradiated","IsolatePKS15.1.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.2_vs_PKS15.2.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.2.conditionIrradiated","IsolatePKS15.2.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.2_vs_PKS15.3.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.2.conditionIrradiated","IsolatePKS15.3.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.2_vs_PKS15.4.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.2.conditionIrradiated","IsolatePKS15.4.2.conditionIrradiated"),alpha = 0.05)

res_PKS15.3.2_vs_PKS15.1.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.3.2.conditionIrradiated","IsolatePKS15.1.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.3.2_vs_PKS15.2.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.3.2.conditionIrradiated","IsolatePKS15.2.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.3.2_vs_PKS15.2.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.3.2.conditionIrradiated","IsolatePKS15.2.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.3.2_vs_PKS15.4.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.3.2.conditionIrradiated","IsolatePKS15.4.2.conditionIrradiated"),alpha = 0.05)

res_PKS15.4.2_vs_PKS15.1.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.4.2.conditionIrradiated","IsolatePKS15.1.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.4.2_vs_PKS15.2.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.4.2.conditionIrradiated","IsolatePKS15.2.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.4.2_vs_PKS15.2.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.4.2.conditionIrradiated","IsolatePKS15.2.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.4.2_vs_PKS15.3.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.4.2.conditionIrradiated","IsolatePKS15.3.2.conditionIrradiated"),alpha = 0.05)

###Between differences of the evolved WT15 and PKS15 lines
res_PKS15.1.1_vs_WT15.1.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.1.1.conditionIrradiated","IsolateWT15.1.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.1.1_vs_WT15.2.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.1.1.conditionIrradiated","IsolateWT15.2.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.1.1_vs_WT15.2.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.1.1.conditionIrradiated","IsolateWT15.2.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.1.1_vs_WT15.4.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.1.1.conditionIrradiated","IsolateWT15.4.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.1.1_vs_WT15.4.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.1.1.conditionIrradiated","IsolateWT15.4.2.conditionIrradiated"),alpha = 0.05)

res_PKS15.2.1_vs_WT15.1.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.1.conditionIrradiated","IsolateWT15.1.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.1_vs_WT15.2.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.1.conditionIrradiated","IsolateWT15.2.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.1_vs_WT15.2.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.1.conditionIrradiated","IsolateWT15.2.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.1_vs_WT15.4.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.1.conditionIrradiated","IsolateWT15.4.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.1_vs_WT15.4.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.1.conditionIrradiated","IsolateWT15.4.2.conditionIrradiated"),alpha = 0.05)

res_PKS15.2.2_vs_WT15.1.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.2.conditionIrradiated","IsolateWT15.1.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.2_vs_WT15.2.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.2.conditionIrradiated","IsolateWT15.2.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.2_vs_WT15.2.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.2.conditionIrradiated","IsolateWT15.2.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.2_vs_WT15.4.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.2.conditionIrradiated","IsolateWT15.4.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.2.2_vs_WT15.4.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.2.2.conditionIrradiated","IsolateWT15.4.2.conditionIrradiated"),alpha = 0.05)

res_PKS15.3.2_vs_WT15.1.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.3.2.conditionIrradiated","IsolateWT15.1.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.3.2_vs_WT15.2.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.3.2.conditionIrradiated","IsolateWT15.2.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.3.2_vs_WT15.2.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.3.2.conditionIrradiated","IsolateWT15.2.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.3.2_vs_WT15.4.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.3.2.conditionIrradiated","IsolateWT15.4.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.3.2_vs_WT15.4.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.3.2.conditionIrradiated","IsolateWT15.4.2.conditionIrradiated"),alpha = 0.05)

res_PKS15.4.2_vs_WT15.1.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.4.2.conditionIrradiated","IsolateWT15.1.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.4.2_vs_WT15.2.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.4.2.conditionIrradiated","IsolateWT15.2.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.4.2_vs_WT15.2.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.4.2.conditionIrradiated","IsolateWT15.2.2.conditionIrradiated"),alpha = 0.05)
res_PKS15.4.2_vs_WT15.4.1_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.4.2.conditionIrradiated","IsolateWT15.4.1.conditionIrradiated"),alpha = 0.05)
res_PKS15.4.2_vs_WT15.4.2_irr <- results(dds_evol_individual, contrast=list("IsolatePKS15.4.2.conditionIrradiated","IsolateWT15.4.2.conditionIrradiated"),alpha = 0.05)

###Within differences of the evolved WT15 lines
res_WT15.1.2_vs_WT15.2.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.1.2_vs_WT.C1","Isolate_WT15.2.1_vs_WT.C1"),alpha = 0.05)
res_WT15.1.2_vs_WT15.2.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.1.2_vs_WT.C1","Isolate_WT15.2.2_vs_WT.C1"),alpha = 0.05)
res_WT15.1.2_vs_WT15.4.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.1.2_vs_WT.C1","Isolate_WT15.4.1_vs_WT.C1"),alpha = 0.05)
res_WT15.1.2_vs_WT15.4.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.1.2_vs_WT.C1","Isolate_WT15.4.2_vs_WT.C1"),alpha = 0.05)

res_WT15.2.1_vs_WT15.1.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.2.1_vs_WT.C1","Isolate_WT15.1.2_vs_WT.C1"),alpha = 0.05)
res_WT15.2.1_vs_WT15.2.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.2.1_vs_WT.C1","Isolate_WT15.2.2_vs_WT.C1"),alpha = 0.05)
res_WT15.2.1_vs_WT15.4.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.2.1_vs_WT.C1","Isolate_WT15.4.1_vs_WT.C1"),alpha = 0.05)
res_WT15.2.1_vs_WT15.4.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.2.1_vs_WT.C1","Isolate_WT15.4.2_vs_WT.C1"),alpha = 0.05)

res_WT15.2.2_vs_WT15.1.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.2.2_vs_WT.C1","Isolate_WT15.1.2_vs_WT.C1"),alpha = 0.05)
res_WT15.2.2_vs_WT15.2.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.2.2_vs_WT.C1","Isolate_WT15.2.1_vs_WT.C1"),alpha = 0.05)
res_WT15.2.2_vs_WT15.4.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.2.2_vs_WT.C1","Isolate_WT15.4.1_vs_WT.C1"),alpha = 0.05)
res_WT15.2.2_vs_WT15.4.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.2.2_vs_WT.C1","Isolate_WT15.4.2_vs_WT.C1"),alpha = 0.05)

res_WT15.4.1_vs_WT15.1.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.4.1_vs_WT.C1","Isolate_WT15.1.2_vs_WT.C1"),alpha = 0.05)
res_WT15.4.1_vs_WT15.2.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.4.1_vs_WT.C1","Isolate_WT15.2.1_vs_WT.C1"),alpha = 0.05)
res_WT15.4.1_vs_WT15.2.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.4.1_vs_WT.C1","Isolate_WT15.2.2_vs_WT.C1"),alpha = 0.05)
res_WT15.4.1_vs_WT15.4.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.4.1_vs_WT.C1","Isolate_WT15.4.2_vs_WT.C1"),alpha = 0.05)

res_WT15.4.2_vs_WT15.1.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.4.2_vs_WT.C1","Isolate_WT15.1.2_vs_WT.C1"),alpha = 0.05)
res_WT15.4.2_vs_WT15.2.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.4.2_vs_WT.C1","Isolate_WT15.2.1_vs_WT.C1"),alpha = 0.05)
res_WT15.4.2_vs_WT15.2.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.4.2_vs_WT.C1","Isolate_WT15.2.2_vs_WT.C1"),alpha = 0.05)
res_WT15.4.2_vs_WT15.4.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_WT15.4.2_vs_WT.C1","Isolate_WT15.4.1_vs_WT.C1"),alpha = 0.05)

###Within differences of the evolved PKS15 lines
res_PKS15.1.1_vs_PKS15.2.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.1.1_vs_WT.C1","Isolate_PKS15.2.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.1.1_vs_PKS15.2.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.1.1_vs_WT.C1","Isolate_PKS15.2.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.1.1_vs_PKS15.3.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.1.1_vs_WT.C1","Isolate_PKS15.3.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.1.1_vs_PKS15.4.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.1.1_vs_WT.C1","Isolate_PKS15.4.2_vs_WT.C1"),alpha = 0.05)

res_PKS15.2.1_vs_PKS15.1.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.1_vs_WT.C1","Isolate_PKS15.1.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.1_vs_PKS15.2.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.1_vs_WT.C1","Isolate_PKS15.2.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.1_vs_PKS15.3.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.1_vs_WT.C1","Isolate_PKS15.3.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.1_vs_PKS15.4.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.1_vs_WT.C1","Isolate_PKS15.4.2_vs_WT.C1"),alpha = 0.05)

res_PKS15.2.2_vs_PKS15.1.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.2_vs_WT.C1","Isolate_PKS15.1.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.2_vs_PKS15.2.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.2_vs_WT.C1","Isolate_PKS15.2.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.2_vs_PKS15.3.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.2_vs_WT.C1","Isolate_PKS15.3.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.2_vs_PKS15.4.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.2_vs_WT.C1","Isolate_PKS15.4.2_vs_WT.C1"),alpha = 0.05)

res_PKS15.3.2_vs_PKS15.1.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.3.2_vs_WT.C1","Isolate_PKS15.1.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.3.2_vs_PKS15.2.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.3.2_vs_WT.C1","Isolate_PKS15.2.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.3.2_vs_PKS15.2.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.3.2_vs_WT.C1","Isolate_PKS15.2.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.3.2_vs_PKS15.4.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.3.2_vs_WT.C1","Isolate_PKS15.4.2_vs_WT.C1"),alpha = 0.05)

res_PKS15.4.2_vs_PKS15.1.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.4.2_vs_WT.C1","Isolate_PKS15.1.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.4.2_vs_PKS15.2.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.4.2_vs_WT.C1","Isolate_PKS15.2.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.4.2_vs_PKS15.2.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.4.2_vs_WT.C1","Isolate_PKS15.2.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.4.2_vs_PKS15.3.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.4.2_vs_WT.C1","Isolate_PKS15.3.2_vs_WT.C1"),alpha = 0.05)

###Between differences of the evolved WT15 and PKS15 lines
res_PKS15.1.1_vs_WT15.1.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.1.1_vs_WT.C1","Isolate_WT15.1.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.1.1_vs_WT15.2.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.1.1_vs_WT.C1","Isolate_WT15.2.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.1.1_vs_WT15.2.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.1.1_vs_WT.C1","Isolate_WT15.2.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.1.1_vs_WT15.4.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.1.1_vs_WT.C1","Isolate_WT15.4.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.1.1_vs_WT15.4.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.1.1_vs_WT.C1","Isolate_WT15.4.2_vs_WT.C1"),alpha = 0.05)

res_PKS15.2.1_vs_WT15.1.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.1_vs_WT.C1","Isolate_WT15.1.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.1_vs_WT15.2.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.1_vs_WT.C1","Isolate_WT15.2.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.1_vs_WT15.2.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.1_vs_WT.C1","Isolate_WT15.2.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.1_vs_WT15.4.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.1_vs_WT.C1","Isolate_WT15.4.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.1_vs_WT15.4.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.1_vs_WT.C1","Isolate_WT15.4.2_vs_WT.C1"),alpha = 0.05)

res_PKS15.2.2_vs_WT15.1.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.2_vs_WT.C1","Isolate_WT15.1.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.2_vs_WT15.2.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.2_vs_WT.C1","Isolate_WT15.2.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.2_vs_WT15.2.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.2_vs_WT.C1","Isolate_WT15.2.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.2_vs_WT15.4.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.2_vs_WT.C1","Isolate_WT15.4.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.2.2_vs_WT15.4.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.2.2_vs_WT.C1","Isolate_WT15.4.2_vs_WT.C1"),alpha = 0.05)

res_PKS15.3.2_vs_WT15.1.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.3.2_vs_WT.C1","Isolate_WT15.1.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.3.2_vs_WT15.2.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.3.2_vs_WT.C1","Isolate_WT15.2.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.3.2_vs_WT15.2.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.3.2_vs_WT.C1","Isolate_WT15.2.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.3.2_vs_WT15.4.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.3.2_vs_WT.C1","Isolate_WT15.4.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.3.2_vs_WT15.4.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.3.2_vs_WT.C1","Isolate_WT15.4.2_vs_WT.C1"),alpha = 0.05)

res_PKS15.4.2_vs_WT15.1.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.4.2_vs_WT.C1","Isolate_WT15.1.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.4.2_vs_WT15.2.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.4.2_vs_WT.C1","Isolate_WT15.2.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.4.2_vs_WT15.2.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.4.2_vs_WT.C1","Isolate_WT15.2.2_vs_WT.C1"),alpha = 0.05)
res_PKS15.4.2_vs_WT15.4.1_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.4.2_vs_WT.C1","Isolate_WT15.4.1_vs_WT.C1"),alpha = 0.05)
res_PKS15.4.2_vs_WT15.4.2_ctl <- results(dds_evol_individual, contrast=list("Isolate_PKS15.4.2_vs_WT.C1","Isolate_WT15.4.2_vs_WT.C1"),alpha = 0.05)


Plot_results<-function(res){
  name<-deparse(substitute(res))
  ix = which.min(res$padj) # most significant
  res <- res[order(res$padj),] # sort
  res$id<-as.character(gsub("jgi.p.*Exode1.","",row.names(res)))
  pdf(paste0("output_group/",name,".pdf"))
  par(mfrow=c(2,1))
  barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]) #Here we show the most significant gene.
  plotMA(res)
  graphics.off()
  ###Volcano plots
  pdf(paste0("output_group/",name,".volcano.pdf"))
  p<-DEGreport::degVolcano(
    data.frame(res[,c("log2FoldChange","padj")]), # table - 2 columns
    plot_text = data.frame(res[1:5,c("log2FoldChange","padj","id")])) # table to add names
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
            file=paste0("output_group/",name,".csv"))
}

Plot_results(res_WT_irr_vs_ctl)
Plot_results(res_PKS_vs_WT_irr_vs_ctl)
Plot_results(res_WT15common_vs_WT_irr_vs_ctl)
Plot_results(res_PKS15common_vs_PKS_irr_vs_ctl)
Plot_results(res_PKS15common_vs_WT15common_irr_vs_ctl)

Plot_results(res_WT15common_vs_WT_irr) 
Plot_results(res_WT15common_vs_WT_ctl)
Plot_results(res_PKS15common_vs_PKS_irr)
Plot_results(res_PKS15common_vs_PKS_ctl)
Plot_results(res_PKS15_vs_WT15_irr)
Plot_results(res_PKS15_vs_WT15_ctl)

Plot_results(res_WT15.1.2_vs_WT_irr) 
Plot_results(res_WT15.2.1_vs_WT_irr) 
Plot_results(res_WT15.2.2_vs_WT_irr) 
Plot_results(res_WT15.4.1_vs_WT_irr) 
Plot_results(res_WT15.4.2_vs_WT_irr) 
Plot_results(res_PKS15.1.1_vs_PKS_irr) 
Plot_results(res_PKS15.2.1_vs_PKS_irr) 
Plot_results(res_PKS15.2.2_vs_PKS_irr) 
Plot_results(res_PKS15.3.2_vs_PKS_irr) 
Plot_results(res_PKS15.4.2_vs_PKS_irr) 
Plot_results(res_WT15.1.2_vs_WT_ctl) 
Plot_results(res_WT15.2.1_vs_WT_ctl) 
Plot_results(res_WT15.2.2_vs_WT_ctl) 
Plot_results(res_WT15.4.1_vs_WT_ctl) 
Plot_results(res_WT15.4.2_vs_WT_ctl) 
Plot_results(res_PKS15.1.1_vs_PKS_ctl) 
Plot_results(res_PKS15.2.1_vs_PKS_ctl) 
Plot_results(res_PKS15.2.2_vs_PKS_ctl) 
Plot_results(res_PKS15.3.2_vs_PKS_ctl) 
Plot_results(res_PKS15.4.2_vs_PKS_ctl) 

######################################################################
### GO Enrichment Analysis ###########################################
######################################################################

GO_results<-function(res){
  ## assign up/down regulation
  name<-deparse(substitute(res))
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
  #######allOE_genes <- as.character(res_tableOE_tb$gene)
  allOE_genes <- row.names(res)
  ## Extract significant results
  sigOE <- res_tableOE_tb %>%
    filter(padj < 0.05 & threshold_OE==TRUE)
  ## Volcano plot
  sigOE_genes <- as.character(sigOE$gene)
  #term2gene_BP<-term2gene[which(term2gene$keytypes=="BP"),]
  #term2gene_CC<-term2gene[which(term2gene$keytypes=="CC"),]
  #term2gene_MF<-term2gene[which(term2gene$keytypes=="MF"),]
  counts<-ggplot(data=sigOE, aes(x=regulation,fill=regulation)) +
    geom_bar(stat="count")+
    scale_x_discrete(limits=rev)+xlab("Direction of regulation")+
    scale_y_continuous(expand = c(0,0))+ ylab("No. of DEGs")+
    scale_fill_brewer(palette="Paired")+
    geom_text(stat='count', aes(label=..count..), vjust=1.7,size=2)+
    theme_classic()
  gene.counts[[name]]<<-counts
  gene.counts.numbers[[name]]<<-table(sigOE$regulation)
  ########### GO Enrichment #########################
  ego <- enricher(gene = sigOE_genes, 
                  universe = allOE_genes,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  TERM2GENE = term2gene_go[,c("goName","gene_name")])
  ## Dotplot 
  ego_up <- enricher(gene = as.character(sigOE$gene[sigOE$regulation=="up"]), 
                     universe = allOE_genes,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05,
                     TERM2GENE = term2gene_go[,c("goName","gene_name")])
  Title=gsub("WT ","WTC ",gsub("PKS ","PKSC ",gsub("_"," ",gsub("res"," ",name))))
  a1<-dotplot(ego_up,orderBy="GeneRatio",showCategory = 5,title=Title)
  a1<-tryCatch({print(a1)}, error = function(e){
    a1<-NA
    return(a1)})
  GO.plot.up[[name]]<<-a1
  ego_down <- enricher(gene = as.character(sigOE$gene[sigOE$regulation=="down"]), 
                       universe = allOE_genes,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       TERM2GENE = term2gene_go[,c("goName","gene_name")])
  a2<-dotplot(ego_down,orderBy="GeneRatio",showCategory = 5,title=Title)+scale_y_discrete(position = "right")
  a2<-tryCatch({print(a2)}, error = function(e){
    a2<-NA
    return(a2)})
  GO.plot.down[[name]]<<-a2
  a3<-ggarrange(a1,a2,ncol=2,common.legend = TRUE, legend = "right")
  # #a<-dotplot(ego,orderBy="GeneRatio")
  # ## Output results from GO analysis to a table
  # cluster_summary <- data.frame(ego)
  # ## Add similarity matrix to the termsim slot of enrichment result
  # ego <- enrichplot::pairwise_termsim(ego)
  # ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
  # b<-emapplot(ego, showCategory = 50)
  # ## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
  # OE_foldchanges <- sigOE$log2FoldChange
  # names(OE_foldchanges) <- sigOE$gene
  # ## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
  # c<-cnetplot(ego, 
  #          categorySize="pvalue", 
  #          showCategory = 5, 
  #          foldChange=OE_foldchanges,
  #          node_label_size=2)
  pdf(paste0("output_group/",name,".GO.pdf"),width=10,height=8)
  print(a3)
  graphics.off()
  write.table(summary(ego_up)[,-6][,-8], paste0("output_group/Supplementary_",name,"_go_up.txt"), sep="\t", row.names=F, quote=F)
  write.table(summary(ego_down)[,-6][,-8], paste0("output_group/Supplementary_",name,"_go_down.txt"), sep="\t", row.names=F, quote=F)
  
  # ########### KEGG pathway analysis ######################### 
  # ekegg <- enricher(gene = sigOE_genes, 
  #                 universe = allOE_genes,
  #                 pAdjustMethod = "BH",
  #                 qvalueCutoff = 0.05,
  #                 TERM2GENE = term2gene_kegg[,c("pathway","gene_name")])
  # 
  # ekegg_up <- enricher(gene = as.character(sigOE$gene[sigOE$regulation=="up"]), 
  #                    universe = allOE_genes,
  #                    pAdjustMethod = "BH",
  #                    qvalueCutoff = 0.05,
  #                    TERM2GENE = term2gene_kegg[,c("pathway","gene_name")])
  # a<-dotplot(ekegg_up,orderBy="GeneRatio")
  # ekegg_down <- enricher(gene = as.character(sigOE$gene[sigOE$regulation=="down"]), 
  #                      universe = allOE_genes,
  #                      pAdjustMethod = "BH",
  #                      qvalueCutoff = 0.05,
  #                      TERM2GENE = term2gene_kegg[,c("pathway","gene_name")])
  # b<-dotplot(ekegg_down,orderBy="GeneRatio")
  # # c<-cnetplot(ekegg, 
  # #             categorySize="pvalue", 
  # #             showCategory = 5, 
  # #             foldChange=OE_foldchanges, 
  # #             vertex.label.font=6)
  # pdf(paste0("output_group/",name,".KEGG.pdf"),width=14,height=18)
  # print(ggarrange(a,b)) 
  # graphics.off()
  # ########### KOG Enrichment #########################
  # ekog <- enricher(gene = sigOE_genes, 
  #                   universe = allOE_genes,
  #                   pAdjustMethod = "BH",
  #                   qvalueCutoff = 0.05,
  #                   TERM2GENE = term2gene_kog[,c("kogdefline","gene_name")])
  # 
  # ekog_up <- enricher(gene = as.character(sigOE$gene[sigOE$regulation=="up"]), 
  #                      universe = allOE_genes,
  #                      pAdjustMethod = "BH",
  #                      qvalueCutoff = 0.05,
  #                      TERM2GENE = term2gene_kog[,c("kogdefline","gene_name")])
  # a<-dotplot(ekog_up,orderBy="GeneRatio")
  # ekog_down <- enricher(gene = as.character(sigOE$gene[sigOE$regulation=="down"]), 
  #                        universe = allOE_genes,
  #                        pAdjustMethod = "BH",
  #                        qvalueCutoff = 0.05,
  #                        TERM2GENE = term2gene_kog[,c("kogdefline","gene_name")])
  # b<-dotplot(ekog_down,orderBy="GeneRatio")
  # # c<-cnetplot(ekog, 
  # #             categorySize="pvalue", 
  # #             showCategory = 5, 
  # #             foldChange=OE_foldchanges, 
  # #             vertex.label.font=6)
  # pdf(paste0("output_group/",name,".KOG.pdf"),width=14,height=18)
  # print(ggarrange(a,b)) 
  # graphics.off()
  
  #### Table of all the data
  #### Annotation of all GO,KEGG,KOG,IPR
  res_ids <- left_join(res_tableOE_tb, ed_gff_go, by=c("gene"="gene_name"))  
  res_ids <- left_join(res_ids, ed_gff_kegg, by=c("gene"="gene_name")) 
  res_ids <- left_join(res_ids, ed_gff_kog, by=c("gene"="gene_name")) 
  res_ids <- left_join(res_ids, ed_gff_ipr, by=c("gene"="gene_name")) 
  res_sigOEids <- res_ids %>%
    filter(padj < 0.05 & threshold_OE==TRUE)
  write.csv(as.data.frame(res_sigOEids), 
            file=paste0("output_group/",name,".GO_KEGG_KOG_IPR_annotation.csv"))
  ####
}

GO.plot.up=list()
GO.plot.down=list()
gene.counts=list()
gene.counts.numbers=list()

GO_results(res_WT15common_vs_WT_irr) #up=NULL
GO_results(res_WT15common_vs_WT_ctl)
GO_results(res_PKS15common_vs_PKS_irr)
GO_results(res_PKS15common_vs_PKS_ctl)
GO_results(res_PKS15_vs_WT15_irr)
GO_results(res_PKS15_vs_WT15_ctl)

GO_results(res_WT_irr_vs_ctl)
GO_results(res_PKS_irr_vs_ctl)
GO_results(res_WT15common_irr_vs_ctl)
GO_results(res_PKS15common_irr_vs_ctl)
GO_results(res_PKS15common_vs_WT15common_irr_vs_ctl)

GO_results(res_WT15.1.2_vs_WT_irr) 
GO_results(res_WT15.2.1_vs_WT_irr) 
GO_results(res_WT15.2.2_vs_WT_irr) 
GO_results(res_WT15.4.1_vs_WT_irr) 
GO_results(res_WT15.4.2_vs_WT_irr) 
GO_results(res_PKS15.1.1_vs_PKS_irr) 
GO_results(res_PKS15.2.1_vs_PKS_irr) 
GO_results(res_PKS15.2.2_vs_PKS_irr) 
GO_results(res_PKS15.3.2_vs_PKS_irr) 
GO_results(res_PKS15.4.2_vs_PKS_irr) 
GO_results(res_WT15.1.2_vs_WT_ctl) 
GO_results(res_WT15.2.1_vs_WT_ctl) 
GO_results(res_WT15.2.2_vs_WT_ctl) 
GO_results(res_WT15.4.1_vs_WT_ctl) 
GO_results(res_WT15.4.2_vs_WT_ctl) 
GO_results(res_PKS15.1.1_vs_PKS_ctl) 
GO_results(res_PKS15.2.1_vs_PKS_ctl) 
GO_results(res_PKS15.2.2_vs_PKS_ctl) 
GO_results(res_PKS15.3.2_vs_PKS_ctl) 
GO_results(res_PKS15.4.2_vs_PKS_ctl) 


###Within differences of the evolved WT15 lines
GO_results(res_WT15.1.2_vs_WT15.2.1_irr)
GO_results(res_WT15.1.2_vs_WT15.2.2_irr)
GO_results(res_WT15.1.2_vs_WT15.4.1_irr)
GO_results(res_WT15.1.2_vs_WT15.4.2_irr)

GO_results(res_WT15.2.1_vs_WT15.1.2_irr)
GO_results(res_WT15.2.1_vs_WT15.2.2_irr)
GO_results(res_WT15.2.1_vs_WT15.4.1_irr)
GO_results(res_WT15.2.1_vs_WT15.4.2_irr)

GO_results(res_WT15.2.2_vs_WT15.1.2_irr)
GO_results(res_WT15.2.2_vs_WT15.2.1_irr)
GO_results(res_WT15.2.2_vs_WT15.4.1_irr)
GO_results(res_WT15.2.2_vs_WT15.4.2_irr)

GO_results(res_WT15.4.1_vs_WT15.1.2_irr)
GO_results(res_WT15.4.1_vs_WT15.2.1_irr)
GO_results(res_WT15.4.1_vs_WT15.2.2_irr)
GO_results(res_WT15.4.1_vs_WT15.4.2_irr)

GO_results(res_WT15.4.2_vs_WT15.1.2_irr)
GO_results(res_WT15.4.2_vs_WT15.2.1_irr)
GO_results(res_WT15.4.2_vs_WT15.2.2_irr)
GO_results(res_WT15.4.2_vs_WT15.4.1_irr)

###Within differences of the evolved PKS15 lines
GO_results(res_PKS15.1.1_vs_PKS15.2.1_irr)
GO_results(res_PKS15.1.1_vs_PKS15.2.2_irr)
GO_results(res_PKS15.1.1_vs_PKS15.3.2_irr)
GO_results(res_PKS15.1.1_vs_PKS15.4.2_irr)

GO_results(res_PKS15.2.1_vs_PKS15.1.1_irr)
GO_results(res_PKS15.2.1_vs_PKS15.2.2_irr)
GO_results(res_PKS15.2.1_vs_PKS15.3.2_irr)
GO_results(res_PKS15.2.1_vs_PKS15.4.2_irr)

GO_results(res_PKS15.2.2_vs_PKS15.1.1_irr)
GO_results(res_PKS15.2.2_vs_PKS15.2.1_irr)
GO_results(res_PKS15.2.2_vs_PKS15.3.2_irr)
GO_results(res_PKS15.2.2_vs_PKS15.4.2_irr)

GO_results(res_PKS15.3.2_vs_PKS15.1.1_irr)
GO_results(res_PKS15.3.2_vs_PKS15.2.1_irr)
GO_results(res_PKS15.3.2_vs_PKS15.2.2_irr)
GO_results(res_PKS15.3.2_vs_PKS15.4.2_irr)

GO_results(res_PKS15.4.2_vs_PKS15.1.1_irr)
GO_results(res_PKS15.4.2_vs_PKS15.2.1_irr)
GO_results(res_PKS15.4.2_vs_PKS15.2.2_irr)
GO_results(res_PKS15.4.2_vs_PKS15.3.2_irr)

###Between differences of the evolved WT15 and PKS15 lines
GO_results(res_PKS15.1.1_vs_WT15.1.2_irr)
GO_results(res_PKS15.1.1_vs_WT15.2.1_irr)
GO_results(res_PKS15.1.1_vs_WT15.2.2_irr)
GO_results(res_PKS15.1.1_vs_WT15.4.1_irr)
GO_results(res_PKS15.1.1_vs_WT15.4.2_irr)

GO_results(res_PKS15.2.1_vs_WT15.1.2_irr)
GO_results(res_PKS15.2.1_vs_WT15.2.1_irr)
GO_results(res_PKS15.2.1_vs_WT15.2.2_irr)
GO_results(res_PKS15.2.1_vs_WT15.4.1_irr)
GO_results(res_PKS15.2.1_vs_WT15.4.2_irr)

GO_results(res_PKS15.2.2_vs_WT15.1.2_irr)
GO_results(res_PKS15.2.2_vs_WT15.2.1_irr)
GO_results(res_PKS15.2.2_vs_WT15.2.2_irr)
GO_results(res_PKS15.2.2_vs_WT15.4.1_irr)
GO_results(res_PKS15.2.2_vs_WT15.4.2_irr)

GO_results(res_PKS15.3.2_vs_WT15.1.2_irr)
GO_results(res_PKS15.3.2_vs_WT15.2.1_irr)
GO_results(res_PKS15.3.2_vs_WT15.2.2_irr)
GO_results(res_PKS15.3.2_vs_WT15.4.1_irr)
GO_results(res_PKS15.3.2_vs_WT15.4.2_irr)

GO_results(res_PKS15.4.2_vs_WT15.1.2_irr)
GO_results(res_PKS15.4.2_vs_WT15.2.1_irr)
GO_results(res_PKS15.4.2_vs_WT15.2.2_irr)
GO_results(res_PKS15.4.2_vs_WT15.4.1_irr)
GO_results(res_PKS15.4.2_vs_WT15.4.2_irr)

###Within differences of the evolved WT15 lines
GO_results(res_WT15.1.2_vs_WT15.2.1_ctl)
GO_results(res_WT15.1.2_vs_WT15.2.2_ctl)
GO_results(res_WT15.1.2_vs_WT15.4.1_ctl)
GO_results(res_WT15.1.2_vs_WT15.4.2_ctl)

GO_results(res_WT15.2.1_vs_WT15.1.2_ctl)
GO_results(res_WT15.2.1_vs_WT15.2.2_ctl)
GO_results(res_WT15.2.1_vs_WT15.4.1_ctl)
GO_results(res_WT15.2.1_vs_WT15.4.2_ctl)

GO_results(res_WT15.2.2_vs_WT15.1.2_ctl)
GO_results(res_WT15.2.2_vs_WT15.2.1_ctl)
GO_results(res_WT15.2.2_vs_WT15.4.1_ctl)
GO_results(res_WT15.2.2_vs_WT15.4.2_ctl)

GO_results(res_WT15.4.1_vs_WT15.1.2_ctl)
GO_results(res_WT15.4.1_vs_WT15.2.1_ctl)
GO_results(res_WT15.4.1_vs_WT15.2.2_ctl)
GO_results(res_WT15.4.1_vs_WT15.4.2_ctl)

GO_results(res_WT15.4.2_vs_WT15.1.2_ctl)
GO_results(res_WT15.4.2_vs_WT15.2.1_ctl)
GO_results(res_WT15.4.2_vs_WT15.2.2_ctl)
GO_results(res_WT15.4.2_vs_WT15.4.1_ctl)

###Within differences of the evolved PKS15 lines
GO_results(res_PKS15.1.1_vs_PKS15.2.1_ctl)
GO_results(res_PKS15.1.1_vs_PKS15.2.2_ctl)
GO_results(res_PKS15.1.1_vs_PKS15.3.2_ctl)
GO_results(res_PKS15.1.1_vs_PKS15.4.2_ctl)

GO_results(res_PKS15.2.1_vs_PKS15.1.1_ctl)
GO_results(res_PKS15.2.1_vs_PKS15.2.2_ctl)
GO_results(res_PKS15.2.1_vs_PKS15.3.2_ctl)
GO_results(res_PKS15.2.1_vs_PKS15.4.2_ctl)

GO_results(res_PKS15.2.2_vs_PKS15.1.1_ctl)
GO_results(res_PKS15.2.2_vs_PKS15.2.1_ctl)
GO_results(res_PKS15.2.2_vs_PKS15.3.2_ctl)
GO_results(res_PKS15.2.2_vs_PKS15.4.2_ctl)

GO_results(res_PKS15.3.2_vs_PKS15.1.1_ctl)
GO_results(res_PKS15.3.2_vs_PKS15.2.1_ctl)
GO_results(res_PKS15.3.2_vs_PKS15.2.2_ctl)
GO_results(res_PKS15.3.2_vs_PKS15.4.2_ctl)

GO_results(res_PKS15.4.2_vs_PKS15.1.1_ctl)
GO_results(res_PKS15.4.2_vs_PKS15.2.1_ctl)
GO_results(res_PKS15.4.2_vs_PKS15.2.2_ctl)
GO_results(res_PKS15.4.2_vs_PKS15.3.2_ctl)

###Between differences of the evolved WT15 and PKS15 lines
GO_results(res_PKS15.1.1_vs_WT15.1.2_ctl)
GO_results(res_PKS15.1.1_vs_WT15.2.1_ctl)
GO_results(res_PKS15.1.1_vs_WT15.2.2_ctl)
GO_results(res_PKS15.1.1_vs_WT15.4.1_ctl)
GO_results(res_PKS15.1.1_vs_WT15.4.2_ctl)

GO_results(res_PKS15.2.1_vs_WT15.1.2_ctl)
GO_results(res_PKS15.2.1_vs_WT15.2.1_ctl)
GO_results(res_PKS15.2.1_vs_WT15.2.2_ctl)
GO_results(res_PKS15.2.1_vs_WT15.4.1_ctl)
GO_results(res_PKS15.2.1_vs_WT15.4.2_ctl)

GO_results(res_PKS15.2.2_vs_WT15.1.2_ctl)
GO_results(res_PKS15.2.2_vs_WT15.2.1_ctl)
GO_results(res_PKS15.2.2_vs_WT15.2.2_ctl)
GO_results(res_PKS15.2.2_vs_WT15.4.1_ctl)
GO_results(res_PKS15.2.2_vs_WT15.4.2_ctl)

GO_results(res_PKS15.3.2_vs_WT15.1.2_ctl)
GO_results(res_PKS15.3.2_vs_WT15.2.1_ctl)
GO_results(res_PKS15.3.2_vs_WT15.2.2_ctl)
GO_results(res_PKS15.3.2_vs_WT15.4.1_ctl)
GO_results(res_PKS15.3.2_vs_WT15.4.2_ctl)

GO_results(res_PKS15.4.2_vs_WT15.1.2_ctl)
GO_results(res_PKS15.4.2_vs_WT15.2.1_ctl)
GO_results(res_PKS15.4.2_vs_WT15.2.2_ctl)
GO_results(res_PKS15.4.2_vs_WT15.4.1_ctl)
GO_results(res_PKS15.4.2_vs_WT15.4.2_ctl)
######################################################################
###### WGCNA:weighted correlation network analysis ###################
###### weighted gene co-expression network analysis ##################
######################################################################

###From: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.pdf
options(stringsAsFactors = FALSE)
###WGCNA requires at least 15 samples to produce a meaningful result
# Choose a set of soft-thresholding powers
powers = c(c(1:20))

# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)

# Retrieve the normalized data from the `DESeqDataSet`
datExpr<- assay(dds_norm) %>%
  t() # Transpose this data
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,
                        corFnc = cor,corOptions = list(use = 'p'))
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.7;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.93,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

###We now calculate the adjacencies, using the soft thresholding power: the lowest power for which the scale-free topology fit index reaches 0.85 or 0.90
softPower = 17;
direction<-"signed"
adjacency = adjacency(datExpr, power = softPower, type=direction);

###To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType=direction);
dissTOM = 1-TOM

### Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
### Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

### We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
### Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

### Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
### Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

### Calculate eigengenes; excludeGrey=True because has NA values (missing data)
MEList = moduleEigengenes(datExpr, colors = dynamicColors, excludeGrey = TRUE)
MEs = MEList$eigengenes
### Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
### Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
### Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
###We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
MEDissThres = 0.25
### Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
### Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
### The merged module colors
mergedColors = merge$colors;
### Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
pdf(file = paste0("output_softpower",softPower,"/Figure_geneDendro-3_",direction,".pdf"), wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
graphics.off()

table(mergedColors)
### Rename to moduleColors
moduleColors = mergedColors
### Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
### Save module colors and labels for use in subsequent parts
#save(MEs, moduleLabels, moduleColors, geneTree, file = paste0(paste0("output_softpower",softPower,"/WGCNA_8115_Rad52-02-networkConstruction-stepByStep","_",direction,".RData")))


### Relating modules to external information and identifying important genes ###
### Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

### Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
###get "trait" or "radiation" data
datTraits<-coldata[rownames(MEs),c("genotype","Transfers","condition","radiation_level")]
### Condition as numeric values
datTraits$condition[datTraits$condition=="Control"]<-0
datTraits$condition[datTraits$condition=="Irradiated"]<-1
datTraits$condition<-as.numeric(datTraits$condition)
### Genotype as numeric values: assign 0 to deletion mutants because missing gene
datTraits$genotype_PKS<-gsub("WT.*|8115|R52","0",datTraits$genotype)
datTraits$genotype_PKS<-gsub("PKS.*","1",datTraits$genotype_PKS) 
datTraits$genotype_WT<-gsub("PKS.*|8115|R52",0,datTraits$genotype)
datTraits$genotype_WT<-gsub("WT.*","1",datTraits$genotype_WT) 
datTraits$genotype<-NULL

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
### Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
### Display the correlation values within a heatmap plot: correlation is color, parentheses are p-values
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
View(t(rbind(colnames(datExpr),moduleColors)))
###Gene relationship to trait and important modules: Gene Significance and Module Membership
### Define variable condition containing the condition column of datTrait
#radiation_level = as.data.frame(datTraits$Transfers);
#names(radiation_level) = "Transfers"
radiation_level = as.data.frame(datTraits$radiation_level);
names(radiation_level) = "radiation_level"
### names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, radiation_level, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(radiation_level), sep="");
names(GSPvalue) = paste("p.GS.", names(radiation_level), sep="");

##### Test different modules and their correlation with Irradiation
scatter_plot<-vector('list',length(unique(moduleColors)[unique(moduleColors)!="grey"]))

scatter_plot<-function(module){
  print(module)
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  selected<-c("8115","rad52","pks1")
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for Irradiation",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  text(abs(geneModuleMembership[moduleGenes,][selected,3]),
       abs(geneTraitSignificance[moduleGenes,,drop=FALSE][selected, 1]),
       labels=rownames(geneModuleMembership))
  
  p<-recordPlot()
  
}

scatter_plots <- lapply(unique(moduleColors)[unique(moduleColors)!="grey"], scatter_plot)
###will return all gene IDs included in the analysis
probes=colnames(datExpr)
###will return all gene IDs in the color module
colnames(datExpr)[moduleColors=="darkorange2"]
colnames(datExpr)[moduleColors=="greenyellow"]

annot<-Reduce(function(...) merge(..., by='gene_name', all=TRUE), 
              list(geneTraitSignificance%>%
                     mutate(gene_name = rownames(geneTraitSignificance)),
                   ed_gff_go,ed_gff_kegg,ed_gff_kog,ed_gff_ipr))


module_colors<-cbind(moduleColors,probes)
colnames(module_colors)<-c("moduleColors","gene_name")

geneInfo0 <- Reduce(
  function(x, y, ...) merge(x, y, all = TRUE, ...),
  list(annot, module_colors,
       tibble::rownames_to_column(geneTraitSignificance, "gene_name"),
       tibble::rownames_to_column(GSPvalue, "gene_name"))
)

colnames(geneInfo0)<-c("genes","GS.condition","GO","KEGG","KOG","IPR","moduleColors","p.GS.condition")


### Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, radiation_level, use = "p")));
### Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
### Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.condition));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = paste0("output_softpower",softPower,"/ModuleMembership_",names(radiation_level),"_geneInfo_annotation","_",direction,".csv"))

######## Module connectivity/correlation between modules and Trait ###############

MET = (orderMEs(cbind(subset(MEs, select=-(MEgrey)), radiation_level)))
MET = (orderMEs(cbind(MEs, radiation_level)))

sizeGrWindow(5,10);
par(cex = 0.9)
pdf(paste0("output_softpower",softPower,"/Figure_betweenModuleCorrelation_trait_signed.pdf"))
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(5,5,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
graphics.off()

#################### Gene Network Visualization ###############
### get the negative/positive correlation between genes
### identify hub genes for all modules
hubs = chooseTopHubInEachModule(datExpr, moduleColors,power = softPower,type = direction)
# Select modules
top5_sub_hubs<-data.frame(matrix(ncol=2,nrow=0))
colnames(top5_sub_hubs)<-c("module","hub_score")
cytoscape_tables<-function(modules){
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule];
  modGenes = annot$gene_name[match(modProbes, annot$gene_name)];
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  ### top hub genes (most connections) and 8115, rad52, pks1
  nTop = 50;
  IMConn = softConnectivity(datExpr[, modProbes]);
  top = (rank(-IMConn) <= nTop | colnames(modTOM) %in%c("8115","rad52","pks1"))
  adj <- modTOM[top, top]
  network <- graph.adjacency(adj,weighted=TRUE)
  network <- simplify(network)  # removes self-loops
  V(network)$color <- modules
  par(mar=c(0,0,0,0))
  # remove unconnected nodes
  network <- delete.vertices(network, degree(network)==0)
  ### identify Hub genes
  sub_hubs=sort(hub_score(network)$vector,decreasing=TRUE)
  top5_sub_hubs_module<-data.frame(module=modules,
                                   hub_score=sub_hubs[1:5])
  top5_sub_hubs<<-rbind(top5_sub_hubs,top5_sub_hubs_module)
  pdf(paste0("output_softpower",softPower,"/Figure-network-", paste(modules, collapse="-"),"_",direction,".pdf", sep=""),paper="us")
  plot(network,directed=FALSE,  
       edge.color=adjustcolor("grey", alpha.f = 0.3),
       edge.arrow.mode=0,
       vertex.size= ifelse(names(V(network)) %in% c(hubs), 10, 4), #4,
       vertex.color = ifelse(names(V(network)) %in% c("8115","rad52","pks1",hubs,names(sub_hubs)[1:5]), adjustcolor(modules, alpha.f = 1), adjustcolor(modules, alpha.f = 0.4)),#adjustcolor(modules, alpha.f = 0.5),
       vertex.frame.color	= adjustcolor(modules, alpha.f = 0),
       vertex.label.color = "black",
       vertex.label.font = 2,
       vertex.label.dist=1.5,
       vertex.label.cex = 0.4, 
       vertex.label.degree = -pi/2,
       vertex.label = ifelse(names(V(network)) %in% c("8115","rad52","pks1",hubs,names(sub_hubs)[1:5]), names(V(network)), NA)
  )
  graphics.off() 
  
  # Export the network into edge and node list files Cytoscape can read
  #cyt = exportNetworkToCytoscape(modTOM,
  #                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"),"_",direction, ".txt", sep=""),
  #                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"),"_",direction, ".txt", sep=""),
  #                               weighted = TRUE,
  #                               threshold = 0.2,
  #                               nodeNames = modProbes,
  #                               altNodeNames = modGenes,
  #                               nodeAttr = moduleColors[inModule]);
  
}
lapply(unique(moduleColors)[unique(moduleColors)!="grey"], cytoscape_tables)
#cytoscape_tables(modules = c("greenyellow"))

top5_sub_hubs$gene_name<-rownames(top5_sub_hubs)
top5_sub_hubs<-left_join(top5_sub_hubs,annot, by="gene_name")
write.table(top5_sub_hubs, paste0("output_softpower",softPower,"/Table_top5_hubgenes_module","_",direction, ".txt"), sep="\t", row.names=T, quote=F)

############ WGCNA GO enrichment per module ##############
########### GO Enrichment #########################
GO_wgcna<-function(modules){
  inModule = is.finite(match(moduleColors, modules))
  modProbes = probes[inModule]
  modGenes = annot$gene_name[match(modProbes, annot$gene_name)]
  #############GO Enrichment######################
  ego <- enricher(gene = modGenes, 
                  universe = term2gene_go$gene_name,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  TERM2GENE = term2gene_go[,c("goName","gene_name")])
  ## Dotplot 
  a1<-dotplot(ego,orderBy="GeneRatio",showCategory = 5)+ggtitle(modules)+
    theme(plot.title = element_text(size = 7, face = "bold"),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7))+
    font("xy.text", size = 6)
  a1<-tryCatch({print(a1)}, error = function(e){
    a1<-NA
    return(a1)})
  GO.plot.wgcna[[modules]]<<-a1
  #pdf(paste0("output_softpower",softPower,"/",modules,"_",direction,".GO.pdf"),width=14,height=18)
  #print(a1)
  #graphics.off()
  write.table(summary(ego)[,-6][,-8], paste0("output_softpower",softPower,"/","Supplementary_",modules,"_",direction,"_go.txt"), sep="\t", row.names=F, quote=F)
  #############KEGG Enrichment######################
  #ekegg <- enricher(gene = modGenes, 
  #                  universe = term2gene_kegg$gene_name,
  #                  pAdjustMethod = "BH",
  #                  qvalueCutoff = 0.05,
  #                  TERM2GENE = term2gene_kegg[,c("pathway","gene_name")])
  #a2<-dotplot(ekegg,orderBy="GeneRatio")
  #pdf(paste0("output_softpower",softPower,"/",modules,"_",direction,".KEGG.pdf"),width=14,height=18)
  #print(a2)
  #graphics.off()
  #write.table(summary(ekegg)[,-6][,-8], paste0("output_softpower",softPower,"/","Supplementary_",modules,"_",direction,"_kegg.txt"), sep="\t", row.names=F, quote=F)
  #############KOG Enrichment######################
  #ekog <- enricher(gene = modGenes,
  #                 universe = term2gene_kog$gene_name,
  #                 pAdjustMethod = "BH",
  #                 qvalueCutoff = 0.05,
  #                 TERM2GENE = term2gene_kog[,c("kogdefline","gene_name")])
  #a3<-dotplot(ekog,orderBy="GeneRatio")
  #pdf(paste0("output_softpower",softPower,"/",modules,".KOG.pdf"),width=14,height=18)
  #print(a3)
  #graphics.off()
  #write.table(summary(ekog)[,-6][,-8], paste0("output_softpower",softPower,"/","Supplementary_",modules,"_",direction,"_kog.txt"), sep="\t", row.names=F, quote=F)
}

#### Identify modules with significant correlation with radiation
significant_modules<-moduleTraitPvalue[,"radiation_level"][moduleTraitPvalue[,"radiation_level"]<0.05]
significant_modules<-names(significant_modules)
####Positive Correlation with Irradiation
GO_wgcna(modules="cyan") #no KOG
GO_wgcna(modules="orangered4") #no GO,KEGG,KOG
GO_wgcna(modules="plum1") #no KEGG
GO_wgcna(modules="darkorange") 
GO_wgcna(modules="darkgreen") 
GO_wgcna(modules="paleturquoise") #no GO,KEGG,KOG
####Negative Correlation with Irradiation
GO_wgcna(modules="lightyellow")
GO_wgcna(modules="violet")
GO_wgcna(modules="brown")
GO_wgcna(modules="yellow")

GO_wgcna(modules="black")
GO_wgcna(modules="lightcyan1")
GO_wgcna(modules="green")
GO_wgcna(modules="grey")

GO.plot.wgcna<-list()
GO_wgcna(modules="green")
GO_wgcna(modules="turquoise")
GO_wgcna(modules="brown")
GO_wgcna(modules="salmon")
GO_wgcna(modules="lightcyan1")
GO_wgcna(modules="blue")
GO_wgcna(modules="violet")    
GO_wgcna(modules="lightyellow")
GO_wgcna(modules="black")
GO_wgcna(modules="pink")
GO_wgcna(modules="darkgreen")
GO_wgcna(modules="orange")
GO_wgcna(modules="paleturquoise")
GO_wgcna(modules="darkorange")   
GO_wgcna(modules="cyan")
GO_wgcna(modules="skyblue3")
GO_wgcna(modules="yellow")
GO_wgcna(modules="orangered4")
GO_wgcna(modules="plum1")
GO_wgcna(modules="sienna3")

############## Candidate Genes: DEG and WGCNA analyses ###############
###### PKS15 vs PKS control conditions enriched for oxidoreductase activity
###### Top positive correlated Transfer module (black) enriched for oxidoreductase activity
Transfers_df<-read.csv(paste0("output_softpower",softPower,"/","ModuleMembership_Transfers_geneInfo_annotation_signed.csv"))[2:9]
Transfers_df<-subset(Transfers_df,Transfers_df$p.GS.condition<=0.0001)
PKS15common_vs_PKS_ctl_df<-read.csv("output/res_PKS15common_vs_PKS_ctl.GO_KEGG_KOG_IPR_annotation.csv")[2:10]
PKS15common_vs_PKS_ctl_df<-subset(PKS15common_vs_PKS_ctl_df,PKS15common_vs_PKS_ctl_df$padj<=0.0001 & PKS15common_vs_PKS_ctl_df$regulation=="up")

DEG_WGCNA_Transfers<-inner_join(Transfers_df,PKS15common_vs_PKS_ctl_df,by=c("genes"="gene"))
###### PKS15 vs PKS control irradiated enriched for ribosome and translation activity
###### Top negative correlated Transfer module (yellow) enriched for translation and ribosome activity

################################################################################
####################### WT15 and PKS15 Manuscript Figures ######################
################################################################################

############################### Figure S3 ########################################
FigS3A<-ggarrange(gene.counts[["res_WT_irr_vs_ctl"]]+ rremove("xlab")+font("xy.text", size = 7)+font("ylab", size = 7),
                 gene.counts[["res_PKS_irr_vs_ctl"]]+ rremove("xlab")+font("xy.text", size = 7)+font("ylab", size = 7), 
                 gene.counts[["res_WT15common_irr_vs_ctl"]]+ rremove("xlab")+font("xy.text", size = 7)+font("ylab", size = 7), 
                 gene.counts[["res_PKS15common_irr_vs_ctl"]]+font("xlab", size = 7)+font("xy.text", size = 7)+font("ylab", size = 7), 
                 nrow =4, ncol = 1, common.legend = TRUE, legend = "none",labels ="A")

FigS3B<-ggarrange(GO.plot.up[["res_WT_irr_vs_ctl"]]+ rremove("xlab")+font("xy.text", size = 7),
                 GO.plot.up[["res_PKS_irr_vs_ctl"]]+ rremove("xlab")+font("xy.text", size = 7),
                 GO.plot.up[["res_WT15common_irr_vs_ctl"]]+ rremove("xlab")+font("xy.text", size = 7),
                 GO.plot.up[["res_PKS15common_irr_vs_ctl"]]+font("xlab", size = 7)+font("xy.text", size = 7),
                 nrow = 5, ncol = 1, labels = c("B"),align = "hv",common.legend = TRUE,legend="none")

FigS3C<-ggarrange(GO.plot.down[["res_WT_irr_vs_ctl"]]+ rremove("xlab")+font("xy.text", size = 7)+
                   theme(legend.key.size = unit(0.5, 'cm'),
                         legend.text = element_text(size = 7),
                         legend.title = element_text(size = 7)),
                 GO.plot.down[["res_PKS_irr_vs_ctl"]]+ rremove("xlab")+font("xy.text", size = 7)+
                   theme(legend.key.size = unit(0.5, 'cm'),
                         legend.text = element_text(size = 7),
                         legend.title = element_text(size = 7)),
                 GO.plot.down[["res_WT15common_irr_vs_ctl"]]+ rremove("xlab")+font("xy.text", size = 7)+
                   theme(legend.key.size = unit(0.5, 'cm'),
                         legend.text = element_text(size = 7),
                         legend.title = element_text(size = 7)),
                 GO.plot.down[["res_PKS15common_irr_vs_ctl"]]+font("xlab", size = 7)+font("xy.text", size = 7)+
                   theme(legend.key.size = unit(0.5, 'cm'),
                         legend.text = element_text(size = 7),
                         legend.title = element_text(size = 7)),
                 nrow = 5, ncol = 1, labels = c("C"),align = "hv",common.legend = TRUE,legend="right")


FigS3<-ggarrange(FigS3A,FigS3B,FigS3C,ncol=3, widths=c(1,3,3.5),align = "hv")
pdf("output_group/evolFigS3_WT15_PKS15_WTC_PKSC_irr_vs_ctl_GOenrich.pdf",width = 11,height=8.5)
FigS3
graphics.off()

############################### Figure 2 ########################################
Fig3A<-ggarrange(gene.counts[["res_WT15common_vs_WT_ctl"]]+ rremove("xlab")+font("xy.text", size = 6)+font("ylab", size = 6),
                 gene.counts[["res_PKS15common_vs_PKS_ctl"]]+ rremove("xlab")+font("xy.text", size = 6)+font("ylab", size = 6),
                 gene.counts[["res_PKS15_vs_WT15_ctl"]]+ rremove("xlab")+font("xy.text", size = 6)+font("ylab", size = 6), 
                 gene.counts[["res_WT15common_vs_WT_irr"]]+ rremove("xlab")+font("xy.text", size = 6)+font("ylab", size = 6), 
                 gene.counts[["res_PKS15common_vs_PKS_irr"]]+ rremove("xlab")+font("xy.text", size = 6)+font("ylab", size = 6), 
                 gene.counts[["res_PKS15_vs_WT15_irr"]]+font("xlab", size = 6)+font("xy.text", size = 6)+font("ylab", size = 6),
                 nrow =6, ncol = 1, common.legend = TRUE, legend = "none",labels ="A")

Fig3B<-ggarrange(GO.plot.up[["res_WT15common_vs_WT_ctl"]]+ rremove("xlab")+font("xy.text", size = 5),
                 GO.plot.up[["res_PKS15common_vs_PKS_ctl"]]+ rremove("xlab")+font("xy.text", size = 5),
                 GO.plot.up[["res_PKS15_vs_WT15_ctl"]]+ rremove("xlab")+font("xy.text", size = 5),
                 GO.plot.up[["res_WT15common_vs_WT_irr"]]+ rremove("xlab")+font("xy.text", size = 5),
                 GO.plot.up[["res_PKS15common_vs_PKS_irr"]]+ rremove("xlab")+font("xy.text", size = 5),
                 GO.plot.up[["res_PKS15_vs_WT15_irr"]]+font("xlab", size = 6)+font("xy.text", size = 5),
                 nrow = 6, ncol = 1, labels = c("B"),align = "hv",common.legend = TRUE,legend="none")

Fig3C<-ggarrange(GO.plot.down[["res_WT15common_vs_WT_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+
                   theme(legend.key.size = unit(0.5, 'cm'),
                         legend.text = element_text(size = 6),
                         legend.title = element_text(size = 6)),
                 GO.plot.down[["res_PKS15common_vs_PKS_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+
                   theme(legend.key.size = unit(0.5, 'cm'),
                         legend.text = element_text(size = 6),
                         legend.title = element_text(size = 6)),
                 GO.plot.down[["res_PKS15_vs_WT15_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+
                   theme(legend.key.size = unit(0.5, 'cm'),
                         legend.text = element_text(size = 6),
                         legend.title = element_text(size = 6)),
                 GO.plot.down[["res_WT15common_vs_WT_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+
                   theme(legend.key.size = unit(0.5, 'cm'),
                         legend.text = element_text(size = 6),
                         legend.title = element_text(size = 6)),
                 GO.plot.down[["res_PKS15common_vs_PKS_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+
                   theme(legend.key.size = unit(0.5, 'cm'),
                         legend.text = element_text(size = 6),
                         legend.title = element_text(size = 6)),
                 GO.plot.down[["res_PKS15_vs_WT15_irr"]]+font("xlab", size = 6)+font("xy.text", size = 5)+
                   theme(legend.key.size = unit(0.05, 'cm'),
                         legend.text = element_text(size = 6),
                         legend.title = element_text(size = 6)),
                 nrow = 6, ncol = 1, labels = c("C"),align = "hv",common.legend = TRUE,legend="none")


Fig3<-ggarrange(Fig3A,Fig3B,Fig3C,ncol=3, widths=c(1.5,4,6),align = "hv")
pdf("output_group/evolFig3_WT15_PKS15_irr_ctl_GOenrich.pdf",width = 8.5,height=11)
Fig3
graphics.off()

############################### Figure 3 ########vjust =#########################
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
net_cluster <- recordPlot()  
plot.new() 

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

pdf("evolFig3_WGCNA_corrModuleTrait.pdf", width=8.5,height=11)
par(mar = c(6,12, 1, 1.5))
labeledHeatmap(Matrix = moduleTraitCor[,c(1,2,4,5)],
               xLabels = c("generations","irradiation","PKS","WT"),
               yLabels = names(MEs),
               ySymbols = c("regulation of transcription 
                            DNA dependent",
                            "ATP synthesis coupled 
                 proton transport",
                            "binding",
                            "transport",
                            "signal transduction",
                            "DNA binding/
                 replication/repair",
                            "protein binding",
                            "zinc finger",
                            "unfolded protein binding",
                            "oxidoreductase activity/
                 pyridoxal phosphate 
                 binding",
                            "membrane",
                            "oxidoreductase activity/
                 FAD binding",
                            "tricarboxylic acid cycle",
                            "intracellular 
                 protein transport",
                            "oxidoreductase activity/
                 carbohydrate metabolism",
                            "Starch and 
                 sucrose metabolism",
                            "WD40 repeat protein",
                            "nucleic acid binding",
                            "cytoplasm",
                            "translation",
                            "integral to membrane"),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix[,c(1,2,4,5)],
               setStdMargins = FALSE,
               cex.text = 0.5,
               font.lab.x=0.01,
               font.lab.y=0.01,
               zlim = c(-1,1))
graphics.off()
#"MEsalmon" = regulation of transcription, DNA−dependent
#"MEsienna3" = ATP synthesis coupled proton transport
#"MEviolet" = binding
#"MEdarkgreen" = transport
#"MEpink" = signal transduction
#"MEcyan" = DNA binding/replication/repair
#"MEdarkorange" = protein binding
#"MEorangered4" =*hub zinc finger
#"MEplum1" = unfolded protein binding
#"MEorange" = oxidoreductase activity/pyridoxal phosphate binding
#"MEgreen" = membrane
#"MEblack" = oxidoreductase activity/FAD binding
#"MElightcyan1" =tricarboxylic acid cycle
#"MEblue" = intracellular protein transport
#"MEturquoise" = oxidoreductase activity/carbohydrate metabolic process
#"MEpaleturquoise" =*hub Starch and sucrose metabolism
#"MEskyblue3" =*hub WD40 repeat protein
#"MEbrown" = nucleic acid binding
#"MElightyellow" = cytoplasm
#"MEyellow" = translation
#"MEgrey" = integral to membrane
############################### Figure S1 #######################################

pdf("evolFigS1_WT15_PKS15.heatmap.pca.pdf",width = 11,height=8.5)
ggarrange(as.grob(evol_heatmap),
          evol_pca+
            theme(legend.position="none"),
          ncol=2,labels=c("A","B"))
graphics.off()
############################### Figure S2 #######################################
df<-as.data.frame(unlist(gene.counts.numbers))
df$regulation<-gsub(".*irr.|.*ctl.","",rownames(df))
df$condition<-gsub(".*_|.up|.down","",rownames(df))
df$condition<-gsub("irr","irradiation",df$condition)
df$condition<-gsub("ctl","non-irradiated",df$condition)
df$Isolate<-gsub("res_|_vs.*","",rownames(df))
colnames(df)[1]<-"counts"
pdf("output_group/evolFigS2_WT15_PKS15_irr_ctl_genecounts.pdf",width = 8.5,height=11)
ggplot(data=df, aes(x=Isolate,y=counts, fill=regulation)) +
  facet_wrap(~condition)+
  geom_bar(stat="identity", position=position_dodge())+
  scale_x_discrete(limits=rev)+xlab("Direction of regulation")+
  scale_y_continuous(expand = c(0,0))+ ylab("No. of DEGs")+
  scale_fill_brewer(palette="Paired")+
  theme_classic()+
  geom_text(stat="identity", aes(label=counts), vjust=1.7,size=2,position = position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
graphics.off()


############################### Figure S3 #######################################
Fig3SA<-ggarrange(
  GO.plot.up[["res_WT15.1.2_vs_WT_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_PKS_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.1_vs_WT_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_PKS_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.2_vs_WT_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_PKS_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7), 
  GO.plot.up[["res_WT15.4.1_vs_WT_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_PKS_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.2_vs_WT_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_PKS_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")
Fig3SB<-ggarrange( 
  GO.plot.up[["res_WT15.1.2_vs_WT_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_PKS_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7), 
  GO.plot.up[["res_WT15.2.1_vs_WT_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_PKS_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.2_vs_WT_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_PKS_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.1_vs_WT_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_PKS_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7), 
  GO.plot.up[["res_WT15.4.2_vs_WT_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_PKS_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7), 
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")
pdf("output_group/evolFigS3_WT15PKS15ind_irr_ctl_GOenrich.up.pdf",width = 8.5,height=11)
Fig3SA
Fig3SB
graphics.off()
############################### Figure S4 #######################################
Fig4SA<-ggarrange(
  GO.plot.down[["res_WT15.1.2_vs_WT_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_PKS15.1.1_vs_PKS_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_WT15.2.1_vs_WT_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_PKS15.2.1_vs_PKS_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_WT15.2.2_vs_WT_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_PKS15.2.2_vs_PKS_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7), 
  GO.plot.down[["res_WT15.4.1_vs_WT_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_PKS15.3.2_vs_PKS_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_WT15.4.2_vs_WT_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_PKS15.4.2_vs_PKS_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")
Fig4SB<-ggarrange( 
  GO.plot.down[["res_WT15.1.2_vs_WT_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_PKS15.1.1_vs_PKS_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7), 
  GO.plot.down[["res_WT15.2.1_vs_WT_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_PKS15.2.1_vs_PKS_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_WT15.2.2_vs_WT_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_PKS15.2.2_vs_PKS_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_WT15.4.1_vs_WT_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_PKS15.3.2_vs_PKS_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7), 
  GO.plot.down[["res_WT15.4.2_vs_WT_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.down[["res_PKS15.4.2_vs_PKS_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7), 
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")
pdf("output_group/evolFigS4_WT15PKS15ind_irr_ctl_GOenrich.down.pdf",width = 11,height=8.5)
Fig4SA
Fig4SB
graphics.off()
############################### Figure S5 #######################################
FigS3<-ggarrange(plotlist = GO.plot.wgcna,
                 ncol = 2,nrow=5,align = "hv",common.legend = TRUE,legend="bottom")
pdf("evolFigS5_GOenrich_wgcna.pdf",width=8.5,height=11)
FigS3
graphics.off()
############################### Figure S Heatmap #######################################
topVarGenes <- head( order( rowVars( assay(vsd_evol) ), decreasing=TRUE ), 35 )
pdf("evolFigS2_WT15_PKS15.heatmap.genes.pdf",width = 11,height=8.5)
heatmap.2( assay(vsd_evol)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="both", 
           col = colorRampPalette(c("green", "black", "red"))(100),
           ColSideColors = c( Control="blue", Irradiated="red")[
             colData(vsd_evol)$condition ],margins=c(5,10) )
graphics.off()

############################### Figure S Individual #######################################
pdf("evolFigS_WT15_PKS15.individual_pairwise.GOenrich.pdf",width = 11,height=8.5)
ggarrange(	
  GO.plot.up[["res_WT15.1.2_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.1.2_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.1.2_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.1.2_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.1.2_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.1.2_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.1.2_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.1.2_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_WT15.2.1_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.1_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.1_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.1_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.1_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.1_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.1_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.1_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_WT15.2.2_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.2_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.2_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.2_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.2_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.2_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.2_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.2_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_WT15.4.1_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.1_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.1_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.1_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.1_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.1_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.1_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.1_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_WT15.4.2_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.2_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.2_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.2_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.2_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.2_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.2_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.2_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.1.1_vs_PKS15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_PKS15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_PKS15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_PKS15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_PKS15.3.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_PKS15.3.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_PKS15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_PKS15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.2.1_vs_PKS15.1.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_PKS15.1.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_PKS15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_PKS15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_PKS15.3.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_PKS15.3.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_PKS15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_PKS15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.2.2_vs_PKS15.1.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_PKS15.1.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_PKS15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_PKS15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_PKS15.3.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_PKS15.3.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_PKS15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_PKS15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.3.2_vs_PKS15.1.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_PKS15.1.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_PKS15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_PKS15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_PKS15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_PKS15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_PKS15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_PKS15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.4.2_vs_PKS15.1.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_PKS15.1.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_PKS15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_PKS15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_PKS15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_PKS15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_PKS15.3.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_PKS15.3.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.1.1_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.2.1_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.2.2_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.3.2_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.4.2_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_WT15.1.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_WT15.2.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_WT15.2.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_WT15.4.1_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_WT15.4.2_irr"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_WT15.1.2_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.1.2_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.1.2_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.1.2_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.1.2_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.1.2_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.1.2_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.1.2_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_WT15.2.1_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.1_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.1_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.1_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.1_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.1_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.1_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.1_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_WT15.2.2_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.2_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.2_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.2_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.2_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.2_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.2.2_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.2.2_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_WT15.4.1_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.1_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.1_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.1_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.1_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.1_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.1_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.1_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_WT15.4.2_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.2_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.2_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.2_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.2_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.2_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_WT15.4.2_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_WT15.4.2_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.1.1_vs_PKS15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_PKS15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_PKS15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_PKS15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_PKS15.3.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_PKS15.3.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_PKS15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_PKS15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.2.1_vs_PKS15.1.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_PKS15.1.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_PKS15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_PKS15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_PKS15.3.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_PKS15.3.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_PKS15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_PKS15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.2.2_vs_PKS15.1.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_PKS15.1.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_PKS15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_PKS15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_PKS15.3.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_PKS15.3.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_PKS15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_PKS15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.3.2_vs_PKS15.1.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_PKS15.1.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_PKS15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_PKS15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_PKS15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_PKS15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_PKS15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_PKS15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.4.2_vs_PKS15.1.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_PKS15.1.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_PKS15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_PKS15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_PKS15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_PKS15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_PKS15.3.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_PKS15.3.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.1.1_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.1.1_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.1.1_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.2.1_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.1_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.1_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.2.2_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.2.2_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.2.2_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.3.2_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.3.2_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.3.2_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")	
ggarrange(	
  GO.plot.up[["res_PKS15.4.2_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_WT15.1.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_WT15.2.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_WT15.2.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_WT15.4.1_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  GO.plot.up[["res_PKS15.4.2_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),	GO.plot.down[["res_PKS15.4.2_vs_WT15.4.2_ctl"]]+ rremove("xlab")+font("xy.text", size = 5)+font("title", size = 7),
  nrow=5,ncol = 2,align = "hv",common.legend = TRUE,legend="right")		
graphics.off()