# singlecell-analysis-G

[Step 1: Demultiplexing of Samples](#Identification-of-singlecells-based-on-the-genotype)

[Step 2: QC and Filtering the data](#QC-and-filtering-data)

[Step 3: Normalization of the cell data](#Data-normalization-of-data-by-samples)

[Step 4: Cell type identification](#Cell-type-identification)

[Step 5: Cell Composition](#Cell-Composition)

[Step 6: Aggregation of single-cell perl celltype to pseudobulk data](#Pseudo-bulk)

[Step 7: scRNA eQTL](#sceQTL)

[Software Requirments](#Software-requirments)


# Identification of singlecells based on the genotype

**Step (a): Demuxlet** 

Deconvoluation and Identification of sample identity of multiplexed singlecell run is done using program Demuxlet. The program uses the genotype information in form of VCF along with the bam output from cellranger to predict the likelihood to assign a barcode to a specific sample.
```
demuxlet --sam pool001/outs/possorted_genome_bam.bam \
--vcf Filtered_2.Sorted.vcf.gz  --field GT \
--group-list pool001/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
--min-callrate 0.30 \
--out Pool001-GT

```

The output from the above command produces a file called \[prefix\].best. The 6th column 'BEST' contains the barcodes likelihood to a given sample. Only those with a prefix tag 'SNG' is utilized for the downstream analysis.


**Step (b): Subsetting** singlecells per sample and creating pool wise subset

We created individual files containing the barcodes assigned for each of the samples. The list of the barcodes for each of the samples are provided along with the *GEOID*. 

```{r}
library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)

pool.data <-NULL;
pool.data <- Read10X(data.dir ="pool01/outs/filtered_gene_bc_matrices/GRCh38/")
pool<-NULL;
pool<- CreateSeuratObject(counts = pool.data, project = "SC", min.cells = 5)
SAMBCPAG008_V1<-read.table("PAG00_V1_BC.txt",header=FALSE)
pool.vectorPAG00_V1<-SAMBCPAG008_V1[,"V1"]
SAM_poolPAG00_V1<-subset(pool,cell=pool.vectorPAG00_V1)
SAM_poolPAG00_V1$Visit <- "Noninfected"
SAM_poolPAG00_V1$SAMP <- "PAG00"
SAM_poolPAG00_V1$PAGID <- "PAG00_V1"
SAM_poolPAG00_V1$Ethnic <- "Gouin"
SAM_poolPAG00_V1$Batch <- "GB1"
saveRDS (SAM_poolPAG008_V1, file="PAG008_V1.Rds")


CombinedPoolPool1 <- merge(SAM_poolPAG00_V1,y=c(SAM_poolSAM2_V1,SAM_poolSAM3_V1,SAM_poolPSAM4_V1),add.cell.ids=c("GB1-Pool1-SAM1_V1","GB1-Pool1-SAM2_V1" ,"GB1-Pool1-PAGSAM3_V1" ,"GB1-Pool1-PAGSAM4_V1"))

```

**Step (c) : Combine all pools and run seurat**
```{r}
pag.combined<- merge(Pool1.1, y=c(Pool2.1,Pool3.1,Pool4.1,Pool5.1,Pool6.1,Pool7.1,Pool8.1,Pool9.1,Pool10.1,Pool11.1,Pool12.1,Pool13.1,Pool14.1,Pool15.1,Pool16.1), project = "AllSamplesCombined")

```
# QC and filtering data 
```{r}
pag.combined<- subset(pag.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 5200 & percent.mt < 10)
````

# Data normalization by samples

```{r}
pag.combined_batches <- SplitObject(pag.combined, split.by = "PAGID")
pag.combined_batches <- lapply(X = pag.combined_batches, FUN = SCTransform, verbose = FALSE)

```

# Cell type identification

**Step (a) : Multimodal Reference mapping PBMC**
```{r}
reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")
anchors <- list()
for (i in 1:length(pag.combined_batches)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = reference,
    query = pag.combined_batches[[i]],
   k.filter = NA,
  reference.assay = "SCT",
  query.assay = "SCT",
  reference.reduction = "spca",
  normalization.method = "SCT",
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
  )
}

for (i in 1:length(pag.combined_batches)) {
  pag.combined_batches[[i]] <- MapQuery(
  anchorset = anchors[[i]],
  query = pag.combined_batches[[i]],
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)
}
```

**Step (b) : Merging the Dataset and filtering out the cells with low prediction score**

```{r}
Merged_Combined_Batch<- merge(pag.combined_batches[[1]], pag.combined_batches[2:length(pag.combined_batches)], merge.dr = "ref.umap")
Data.filtered <- subset(Merged_Combined_Batch, subset = predicted.celltype.l1.score > 0.7)
```

# Cell Composition

**(i)Removing Samples** 

```
Idents(Merged_Combined_Batch) <- 'PAGID'
MyData_Subset<-subset(Merged_Combined_Batch, idents = c("PAG145_V3","PAG157_V1","PAG158_V1","PAG159_V1","PAG028_V1"), invert = TRUE)

```

**(ii)Calculation number of cells per Sample overall**
```
CountTable<-table(MyData_Subset@meta.data$PAGID)
write.table(CountTable,file="Filtered_CountTable_Per_Sample.txt",sep="\t")
```

**(iii)Calculation number of cells for each celltype per Sample**
```
Idents(MyData_Subset) <- 'predicted.celltype.l1'

## Predicted LeveL: 1
CountTable<-table(MyData_Subset@meta.data$predicted.celltype.l1, MyData_Subset@meta.data$PAGID)
write.table(CountTable,file="Filtered_CountTable_Per_Sample_CellType_Predicted_Level1.txt",sep="\t")
```
**Ploting**
```
CountTable<-t(CountTable)
P1melt<-melt(CountTable)
colnames(P1melt)<-c("PAGID","CellType","Count")
P1melt$Individual<-as.factor(substr(P1melt$PAGID,1,6))
P1melt$Visit<-as.factor(substr(P1melt$PAGID,8,9))
ggplot(P1melt,aes(fill=Visit,y=Count,x=CellType))+geom_bar(position='dodge',stat='identity')+theme_bw()+scale_fill_manual(values=c("#009e73","#E69F00"))+theme(axis.text.x = element_text(angle = 90))   

```
**(iv)Calculation cell proportion for each celltype per Sample**
```
PropCountTable<-prop.table(table(MyData_Subset@meta.data$predicted.celltype.l1, MyData_Subset@meta.data$PAGID), margin = 2)
write.table(PropCountTable,file="Filtered_proportion_Per_Sample_CellType_Predicted_Level1.txt",sep="\t")
```
**PCA plot with grouping based on Infection Status**

```
PropCountTable<-prop.table(table(MyData_Subset@meta.data$predicted.celltype.l1, MyData_Subset@meta.data$PAGID), margin = 2)
PropCountTable1<-as.data.frame.matrix(t(PropCountTable))
PropCountTable1$Visit<-as.factor(substr(rownames(PropCountTable1),8,9))
PropCountTable1melt.pca<-prcomp(PropCountTable1[,c(1:8)], center = TRUE,scale. = TRUE)
ggbiplot(PropCountTable1melt.pca, groups = PropCountTable1$Visit,ellipse=TRUE,varname.adjust=1,repel=TRUE)+xlim(-2.5,2.5)+scale_color_manual(values=c("#009e73","#E69F00"))+theme_bw()
dev.off()
```

# Pseudo bulk

**Step (a) : Converting to SingleCellExperiment**
```
Idents(Merged_Combined_Batch) <- 'predicted.celltype.l1'
GData.SE<-as.SingleCellExperiment(Merged_Combined_Batch)
GData.SE$cluster_id <- GData.SE$predicted.celltype.l2
GData.SE$sample_id <- GData.SE$PAGID
GData.SE$group_id <- GData.SE$Visit
GData.SE$id <-  GData.SE$PAGID
(GData.SE <- prepSCE(GData.SE, kid = "cluster_id",gid = "group_id",  sid = "sample_id",drop = TRUE)) 

nk <- length(kids <- levels(GData.SE$cluster_id))
ns <- length(sids <- levels(GData.SE$sample_id))
names(kids) <- kids; names(sids) <- sids
assayNames(GData.SE)
```

**Step (b) : Per celltype mean aggregate**
```
pbCMean <-NULL
pbCMean <- aggregateData(GData.SE, assay = "counts", fun = "mean",  by = c("cluster_id", "sample_id"))
CellTypes<-assayNames(pbCMean)
for(i in 1:length(CellTypes)) {# Head of for-loop
  Fname<-NULL
  ROutput<-NULL
  Fname<-paste(CellTypes[i],"_aggregate_Count_mean.txt",sep="")
  ROutput<-t((assay(pbCMean[],i)))
  write.table(ROutput,file=Fname,sep = "\t",row.names = F)
}

```



# sceQTL

**Step: Linear Dominant Model**

```
plink --all-pheno --linear dominant interaction  --bfile AnalysisSamples  --no-sex  --pheno $i --covar Covariate_Filtered.txt --out INT  --parameters 1-3 --vif 9999
```


# Software requirments
1: R version 4.0.3

2: Seurat : 4.0.0

3: sctransform : 0.3.2

4: SeuratObject : 4.0.0

5: SeuratDisk : 0.0.0.9019

6: plink : beta5.3

7: [Demuxlet](https://github.com/statgen/demuxlet)

8: VCFtools : 0.1.14

9: cellranger : 2.2.0

10: Muscat :1.2.1 
