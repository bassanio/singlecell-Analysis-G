# singlecell-analysis-G
Singlecell Analysis

[Step 1: Demultiplexing of Samples](#Identification-of-singlecells-based-on-the-genotype)

[Step 2: QC and Filtering the data](#QC-and-filtering-data)

[Step 3: Normalization of the cell data](#Data-normalization-of-data-by-samples)

[Step 4: Cell type identification](#Cell-type-identification)

[Step 5: scRNA eQTL](#sceQTL)

[Software Requirments](#Software-requirments)


# Identification of singlecells based on the genotype

**Step (a): Demuxlet** 
```
demuxlet --sam pool001/outs/possorted_genome_bam.bam \
--vcf Filtered_2.Sorted.vcf.gz  --field GT \
--group-list pool001/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
--min-callrate 0.30 \
--out Pool2-GT

```

**Step (b): Subsetting** singlecells per sample and creating pool wise subset

```{r}
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

# Data normalization of data by samples

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
