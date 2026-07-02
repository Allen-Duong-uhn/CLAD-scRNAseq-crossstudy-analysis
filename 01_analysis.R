#load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(scCustomize)

# File paths (EDIT THESE)
data_dir <- "path/to/data"

#load Toronto objects
objs <- list(
  DBALMac = readRDS(file.path(data_dir, "DBALMac.rds")),
  SBALMac = readRDS(file.path(data_dir, "SBALMac.rds")),
  CLADlungMac = readRDS(file.path(data_dir, "CLADlung_macs_integrated.rds"))
)

# Update objects and add metadata
for (nm in names(objs)) {
  objs[[nm]] <- UpdateSeuratObject(objs[[nm]])
  DefaultAssay(objs[[nm]]) <- "RNA"
  objs[[nm]]$StudyID <- "Moshkelgosha et al."
  objs[[nm]]$SampleID <- objs[[nm]]$orig.ident
}

# Object-specific metadata
objs$DBALMac$SampleType <- "BAL"
objs$SBALMac$SampleType <- "BAL"
objs$CLADlungMac$SampleType <- "Lung"

# unpack and remove list
DBALMac <- objs$DBALMac
SBALMac <- objs$SBALMac
CLADlungMac <- objs$CLADlungMac
rm(objs)

# helper functions
qc_filter <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 &
                 nCount_RNA > 1000 & nCount_RNA < 30000 &
                 percent.mt < 15)
}

process_seurat <- function(obj) {
  obj %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters(resolution = 0.5)
}

# Load and prepare Yan et al. datasets
yan_paths <- c("CLAD1", "CLAD2", "CLAD3", "CLAD5")
Yan_list <- lapply(yan_paths, function(sample_name) {
  counts <- Read10X(
    file.path(data_dir, sample_name)
  )
  obj <- CreateSeuratObject(
    counts,
    project = sample_name
  )
  obj$SampleID <- obj$orig.ident
  obj$SampleType <- "Lung"
  obj <- qc_filter(obj)
  obj <- process_seurat(obj)
  return(obj)
})
names(Yan_list) <- yan_paths

# Load and prepare Khatri et al. dataset
norm_mat <- read.csv(
  file.path(data_dir, "GSE224210_Normalized_Count_Matrix.csv.gz"),
  row.names = 1
)

Khatri_lung <- CreateSeuratObject(counts = as.matrix(norm_mat))

Khatri_lung <- SetAssayData(
  object = Khatri_lung,
  assay = "RNA",
  slot = "data",
  new.data = as.matrix(norm_mat)
)

# preprocessing
Khatri_lung <- Khatri_lung %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5, cluster.name = "unintegrated_clusters") %>%
  RunUMAP(dims = 1:30,reduction.name = "umap.unintegrated")

# Add metadata
IDs <- read.csv(
  file.path(data_dir, "GSE224210_Cell_Identities.csv.gz"),
  row.names = 1,
  header = TRUE,
  check.names = FALSE
)
rownames(IDs) <- gsub("-", ".", rownames(IDs))
common <- intersect(
  colnames(Khatri_lung),
  rownames(IDs)
)
Khatri_lung <- AddMetaData(Khatri_lung, IDs[common, ])
Khatri_lung$SampleID <- Khatri_lung$`sample ID`
Khatri_lung$`sample ID` <- NULL

# remove non-CLAD samples
Khatri_lung <- SetIdent(Khatri_lung, value = "SampleID")
Khatri_lung <- subset(Khatri_lung, ident = c("20-469_DT", "21-287_DT", "21-361_DT"), invert = TRUE)
table(Khatri_lung$SampleID)

# Merge and integrate Yan and Khatri data sets
public_CLAD <- merge(Khatri_lung, y = Yan_list)
public_CLAD <- JoinLayers(public_CLAD)
public_CLAD[["RNA"]] <- split(
  public_CLAD[["RNA"]],
  f = public_CLAD$SampleID
)

# preprocessing before integration
public_CLAD <- public_CLAD %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5, cluster.name = "unintegrated_clusters") %>%
  RunUMAP(dims = 1:30, reduction.name = "umap.unintegrated")

# Harmony integration
public_CLAD <- IntegrateLayers(public_CLAD, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony")

public_CLAD <- public_CLAD %>% 
  FindNeighbors(reduction = "integrated.harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5, cluster.name = "harmony_clusters") %>% 
  RunUMAP(reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony") %>% 
  JoinLayers()


public_CLAD$StudyID <- "Yan et al."

#identify and subset myeloid cells
public_CLAD.markers <- FindAllMarkers(public_CLAD, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
public_CLAD.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
public_CLAD.TOP10 <- public_CLAD.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
VlnPlot(public_CLAD, features = c("MARCO", "CD14", "FCN1"))
FeaturePlot(public_CLAD, reduction = "umap.harmony", features = c("MARCO", "CD14", "FCN1", "LGMN"), label = TRUE)
public_CLADmye <- subset(public_CLAD, ident = c("1", "5", "19", "22"))

# merge and integrate Yan, Khatri and Moshkelgosha_lungs
lungobj_mye <- merge(public_CLADmye, c(CLADlungMac))
lungobj_mye <- lungobj_mye %>% 
  JoinLayers()
lungobj_mye[["RNA"]] <- split(lungobj_mye[["RNA"]], f = lungobj_mye$SampleID)

# preprocessing before integration
lungobj_mye <- lungobj_mye %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.5, cluster.name = "unintegrated_clusters") %>% 
  RunUMAP(dims = 1:30, reduction.name = "umap.unintegrated")

# Harmony integration
lungobj_mye <- IntegrateLayers(lungobj_mye, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony")
lungobj_mye <- lungobj_mye %>% 
  FindNeighbors(reduction = "integrated.harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  RunUMAP(reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony") %>% 
  JoinLayers()

# verify integration
DimPlot(lungobj_mye, reduction = "umap.harmony", group.by = "StudyID", label = TRUE)
DimPlot(lungobj_mye, reduction = "umap.harmony", label = TRUE)
FeaturePlot_scCustom(lungobj_mye, reduction = "umap.harmony", label = TRUE, features = c("MARCO", "CD14", "FCN1", "S100A8", "S100A9", "LGMN"), num_columns = 5)

# annotation of lungobj_mye
DimPlot_scCustom(lungobj_mye, label = TRUE, reduction = "umap.harmony")
FeaturePlot_scCustom(lungobj_mye, label = TRUE, reduction = "umap.harmony", features = c("FABP4", "MARCO", "LGMN", "FCN1", "CD14", "MKI67"))
VlnPlot(lungobj_mye, features = c("FABP4", "MARCO", "LGMN", "FOLR2", "FCN1", "CD14", "S100A8", "S100A9"))
VlnPlot(lungobj_mye, features = c("SLAMF7", "IDO1", "CALHM6", "VAMP5", "GBP5", "CXCL9", "CXCL10", "CXCL11", "IL32", "IL4I1"))

lungobj_macs_ids <- c("Interstitial Macrophage" #0
                      ,"Alveolar Macrophage" #1 
                      ,"Monocyte" #2
                      ,"Alveolar Macrophage" #3
                      ,"Alveolar Macrophage" #4
                      ,"Interstitial Macrophage" #5
                      ,"Monocyte/Macrophage" #6
                      ,"Alveolar Macrophage" #7
                      ,"Alveolar Macrophage" #8
                      ,"Alveolar Macrophage" #9
                      ,"Monocyte/Macrophage" #10
                      ,"ISG/Super.Macro" #11 
                      ,"Interstitial Macrophage" #12
                      ,"Monocyte" #13
                      ,"Alveolar Macrophage" #14
                      ,"Monocyte" #15
                      ,"Alveolar Macrophage" #16
                      ,"ISG/Super.Macro" #17
                      ,"Monocyte/Macrophage" #18
                      ,"Alveolar Macrophage" #19
                      ,"Alveolar Macrophage" #20
                      ,"Monocyte/Macrophage" #21 
                      ,"Monocyte/Macrophage" #22
                      ,"Alveolar Macrophage" #23
)
names(lungobj_macs_ids) <- levels(lungobj_mye)
lungobj_mye <- RenameIdents(lungobj_mye, lungobj_macs_ids)
lungobj_mye$MacsLabels <- Idents(lungobj_mye)

DimPlot(lungobj_mye, reduction = "umap.harmony", label = TRUE)

# BAL integration
BALobj <- merge(SBALMac, c(DBALMac))
BALobj <- JoinLayers(BALobj)
BALobj[["RNA"]] <- split(BALobj[["RNA"]], f = BALobj$SampleID)

# preprocessing BEFORE integration
BALobj <- BALobj %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.5, cluster.name = "unintegrated_clusters") %>% 
  RunUMAP(dims = 1:30, reduction.name = "umap.unintegrated")

# CCA integration
BALobj <- IntegrateLayers(BALobj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca")
BALobj <- BALobj %>% 
  FindNeighbors(reduction = "integrated.cca", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  RunUMAP(reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.harmony") %>% 
  JoinLayers()

# annotation of BALobj
DimPlot_scCustom(BALobj, label = TRUE, reduction = "umap.cca")
FeaturePlot_scCustom(BALobj, label = TRUE, reduction = "umap.cca", features = c("FABP4", "MARCO", "LGMN", "FCN1", "CD14", "MKI67"))
VlnPlot(BALobj, features = c("FABP4", "MARCO", "LGMN", "FCN1", "CD14", "S100A8", "S100A9"))
VlnPlot(BALobj, features = c("SLAMF7", "IDO1", "CALHM6", "VAMP5", "GBP5", "CXCL9", "CXCL10", "CXCL11", "IL32", "IL4I1"))

BALobj_macs_ids <- c("Alveolar Macrophage" #0
                      ,"Alveolar Macrophage" #1 
                      ,"Alveolar Macrophage" #2
                      ,"Monocyte/Macrophage" #3
                      ,"Alveolar Macrophage" #4
                      ,"Alveolar Macrophage" #5
                      ,"Alveolar Macrophage" #6
                      ,"Cycling Macrophage" #7
                      ,"Alveolar Macrophage" #8
                      ,"ISG Macrophage" #9
)
names(BALobj_macs_ids) <- levels(BALobj)
BALobj <- RenameIdents(BALobj, BALobj_macs_ids)
BALobj$MacsLabels <- Idents(BALobj)

DimPlot(BALobj, reduction = "umap.cca", label = TRUE)


# save new .rds files
dir.create("results", showWarnings = FALSE)

saveRDS(lungobj_mye, "results/lungobj_mye_integrated.rds")
saveRDS(BALobj_mye, "results/BALobj.rds")

# save reproducibility information
capture.output(sessionInfo(), file = "results/sessionInfo.txt")
