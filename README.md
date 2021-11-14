# singlecell-analysis-G

[Step 1: Demultiplexing of Samples](#Identification-of-singlecells-based-on-the-genotype)

[Step 2: QC and Filtering the data](#QC-and-filtering-data)

[Step 3: Normalization of the cell data](#Data-normalization-of-data-by-samples)

[Step 4: Cell type identification](#Cell-type-identification)

[Step 5: Cell Composition](#Cell-Composition)

[Step 6: Aggregation of single-cell perl celltype to pseudobulk data](#Pseudo-bulk)

[Step 7: scRNA eQTL](#sceQTL)

[Software Requirments](#Software-requirments)

[Session Info](#Session-Info)



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
Downstream analysis was performed only on the cells identified with predictedscore of 0.7 and above 

**(i)Removing Unmatched Samples** 

```
Idents(Data.filtered) <- 'PAGID'
MyData_Subset<-subset(Data.filtered, idents = c("PAG145_V3","PAG157_V1","PAG158_V1","PAG159_V1","PAG028_V1"), invert = TRUE)

```

**(ii)Calculating overall number of cells per Sample**
```
CountTable<-table(MyData_Subset@meta.data$PAGID)
write.table(CountTable,file="Filtered_CountTable_Per_Sample.txt",sep="\t")
```

**(iii)Calculating number of cells for each celltype per Sample**
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
**(iv)Calculating cell proportion for each celltype per Sample**
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
Idents(Data.filtered ) <- 'predicted.celltype.l1'
GData.SE<-as.SingleCellExperiment(Data.filtered)
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

# Session Info
The session info below is based on the composition analysis
```
R version 4.0.3 (2020-10-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] corrplot_0.90               Hmisc_4.5-0                 Formula_1.2-4               survival_3.2-13            
 [5] lattice_0.20-45             ggpubr_0.4.0                ggbiplot_0.55               scales_1.1.1               
 [9] plyr_1.8.6                  devtools_2.4.2              usethis_2.0.1               purrr_0.3.4                
[13] CATALYST_1.12.2             UpSetR_1.4.0                scater_1.16.2               SingleCellExperiment_1.10.1
[17] SummarizedExperiment_1.18.2 DelayedArray_0.14.1         matrixStats_0.61.0          Biobase_2.48.0             
[21] GenomicRanges_1.40.0        GenomeInfoDb_1.24.2         IRanges_2.22.2              S4Vectors_0.26.1           
[25] BiocGenerics_0.34.0         muscat_1.2.1                ggplot2_3.3.5               cowplot_1.1.1              
[29] RColorBrewer_1.1-2          patchwork_1.1.1             SeuratDisk_0.0.0.9013       SeuratObject_4.0.2         
[33] Seurat_4.0.4                dplyr_1.0.7                

loaded via a namespace (and not attached):
  [1] rsvd_1.0.5                  ica_1.0-2                   ps_1.6.0                    foreach_1.5.1              
  [5] lmtest_0.9-38               rprojroot_2.0.2             crayon_1.4.1                spatstat.core_2.3-0        
  [9] MASS_7.3-54                 nlme_3.1-153                backports_1.2.1             rlang_0.4.11               
 [13] XVector_0.28.0              ROCR_1.0-11                 readxl_1.3.1                irlba_2.3.3                
 [17] nloptr_1.2.2.2              callr_3.7.0                 limma_3.44.3                BiocParallel_1.22.0        
 [21] rjson_0.2.20                bit64_4.0.5                 glue_1.4.2                  sctransform_0.3.2          
 [25] pbkrtest_0.5.1              processx_3.5.2              vipor_0.4.5                 spatstat.sparse_2.0-0      
 [29] AnnotationDbi_1.50.3        spatstat.geom_2.3-0         haven_2.4.3                 tidyselect_1.1.1           
 [33] rio_0.5.27                  fitdistrplus_1.1-6          variancePartition_1.18.3    XML_3.99-0.8               
 [37] tidyr_1.1.4                 zoo_1.8-9                   nnls_1.4                    xtable_1.8-4               
 [41] magrittr_2.0.1              cli_3.0.1                   zlibbioc_1.34.0             rstudioapi_0.13            
 [45] miniUI_0.1.1.1              rpart_4.1-15                tinytex_0.34                shiny_1.7.1                
 [49] BiocSingular_1.4.0          xfun_0.26                   clue_0.3-59                 pkgbuild_1.2.0             
 [53] cluster_2.1.2               caTools_1.18.2              tibble_3.1.5                ggrepel_0.9.1              
 [57] listenv_0.8.0               png_0.1-7                   future_1.22.1               withr_2.4.2                
 [61] bitops_1.0-7                RBGL_1.64.0                 cellranger_1.1.0            ncdfFlow_2.34.0            
 [65] pillar_1.6.3                RcppParallel_5.1.4          gplots_3.1.1                GlobalOptions_0.1.2        
 [69] cachem_1.0.6                multcomp_1.4-17             fs_1.5.0                    CytoML_2.0.5               
 [73] hdf5r_1.3.4                 GetoptLong_1.0.5            DelayedMatrixStats_1.10.1   vctrs_0.3.8                
 [77] ellipsis_0.3.2              generics_0.1.0              tools_4.0.3                 foreign_0.8-81             
 [81] beeswarm_0.4.0              munsell_0.5.0               fastmap_1.1.0               compiler_4.0.3             
 [85] pkgload_1.2.2               abind_1.4-5                 httpuv_1.6.3                sessioninfo_1.1.1          
 [89] plotly_4.9.4.1              GenomeInfoDbData_1.2.3      gridExtra_2.3               glmmTMB_1.1.2.3            
 [93] edgeR_3.30.3                deldir_1.0-2                utf8_1.2.2                  later_1.3.0                
 [97] jsonlite_1.7.2              graph_1.66.0                pbapply_1.5-0               carData_3.0-4              
[101] genefilter_1.70.0           lazyeval_0.2.2              promises_1.2.0.1            car_3.0-11                 
[105] doParallel_1.0.16           latticeExtra_0.6-29         goftest_1.2-3               spatstat.utils_2.2-0       
[109] reticulate_1.22             checkmate_2.0.0             openxlsx_4.2.4              sandwich_3.0-1             
[113] blme_1.0-5                  Rtsne_0.15                  forcats_0.5.1               uwot_0.1.10                
[117] igraph_1.2.6                numDeriv_2016.8-1.1         yaml_2.2.1                  plotrix_3.8-2              
[121] cytolib_2.0.3               flowWorkspace_4.0.6         htmltools_0.5.2             memoise_2.0.0              
[125] locfit_1.5-9.4              viridisLite_0.4.0           digest_0.6.28               assertthat_0.2.1           
[129] mime_0.12                   RSQLite_2.2.8               future.apply_1.8.1          remotes_2.4.1              
[133] data.table_1.14.2           blob_1.2.2                  flowCore_2.0.1              drc_3.0-1                  
[137] splines_4.0.3               RCurl_1.98-1.5              broom_0.7.9                 hms_1.1.1                  
[141] colorspace_2.0-2            ConsensusClusterPlus_1.52.0 base64enc_0.1-3             ggbeeswarm_0.6.0           
[145] shape_1.4.6                 nnet_7.3-16                 Rcpp_1.0.7                  RANN_2.6.1                 
[149] mvtnorm_1.1-2               circlize_0.4.13             FlowSOM_1.20.0              RProtoBufLib_2.0.0         
[153] fansi_0.5.0                 parallelly_1.28.1           R6_2.5.1                    ggridges_0.5.3             
[157] lifecycle_1.0.1             zip_2.2.0                   curl_4.3.2                  ggsignif_0.6.3             
[161] minqa_1.2.4                 leiden_0.3.9                testthat_3.1.0              Matrix_1.3-4               
[165] desc_1.4.0                  RcppAnnoy_0.0.19            TH.data_1.1-0               iterators_1.0.13           
[169] TMB_1.7.22                  stringr_1.4.0               htmlwidgets_1.5.4           polyclip_1.10-0            
[173] ComplexHeatmap_2.4.3        mgcv_1.8-38                 globals_0.14.0              htmlTable_2.2.1            
[177] codetools_0.2-18            gtools_3.9.2                prettyunits_1.1.1           gtable_0.3.0               
[181] tsne_0.1-3                  DBI_1.1.1                   tensor_1.5                  httr_1.4.2                 
[185] KernSmooth_2.23-20          stringi_1.7.5               progress_1.2.2              reshape2_1.4.4             
[189] annotate_1.66.0             viridis_0.6.1               hexbin_1.28.2               Rgraphviz_2.32.0           
[193] xml2_1.3.2                  colorRamps_2.3              ggcyto_1.16.0               boot_1.3-28                
[197] BiocNeighbors_1.6.0         lme4_1.1-27.1               geneplotter_1.66.0          scattermore_0.7            
[201] DESeq2_1.28.1               bit_4.0.4                   jpeg_0.1-9                  spatstat.data_2.1-0        
[205] pkgconfig_2.0.3             lmerTest_3.1-3              rstatix_0.7.0               knitr_1.36      
```
