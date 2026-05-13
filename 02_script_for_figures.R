# reproducibility
set.seed(1234)

# create output directory
dir.create("figures", showWarnings = FALSE)

# common UMAP theme
theme_umap <- theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title = element_blank(),
  plot.title = element_blank()
)

# marker genes for feature plots
features <- c(
  "SLAMF7",
  "IDO1",
  "CALHM6",
  "VAMP5",
  "GBP5",
  "CXCL9",
  "CXCL10",
  "CXCL11",
  "IL32",
  "IL4I1"
)

# Generate BAL macrophage cluster UMAP (Figure A, left)
p <- DimPlot(
  BALobj,
  group.by = "seurat_clusters",
  reduction = "umap.harmony",
  label = FALSE,
  raster = TRUE,
  pt.size = 3
) +
  theme_umap

ggsave(
  filename = file.path("figures", "BALobj.svg"),
  plot = p,
  dpi = 400,
  units = "px",
  width = 425,
  height = 400
)

# Generate BAL feature plots (Figure B)
for (gene in features) {

  p <- FeaturePlot_scCustom(
    BALobj,
    features = gene,
    reduction = "umap.harmony"
  ) +
    theme_umap +
    NoLegend()

  ggsave(
    filename = file.path(
      "figures",
      paste0("BAL_", gene, "_FeaturePlot.png")
    ),
    plot = p,
    dpi = 400,
    width = 2,
    height = 2
  )
}

# Generate integrated lung UMAP split by study (Figure A, right)
p <- DimPlot(
  lungobj_mye,
  reduction = "umap.harmony",
  split.by = "StudyID",
  raster = TRUE,
  pt.size = 3
) +
  theme_umap

ggsave(
  filename = file.path("figures", "lungobj_mye_split.svg"),
  plot = p,
  dpi = 400,
  units = "px",
  width = 950,
  height = 400
)

# Subset Moshkelgosha et al. samples 
lungobj_mye_moshkelgosha <- subset(
  lungobj_mye,
  subset = StudyID == "Moshkelgosha et al."
)

# Generate Moshkelgosha et al. feature plots (Figure C, left)
for (gene in features) {

  p <- FeaturePlot_scCustom(
    lungobj_mye_moshkelgosha,
    features = gene,
    reduction = "umap.harmony"
  ) +
    theme_umap +
    NoLegend()

  ggsave(
    filename = file.path(
      "figures",
      paste0("Moshkelgosha_", gene, "_FeaturePlot.png")
    ),
    plot = p,
    dpi = 400,
    width = 2,
    height = 2
  )
}

# Subset Yan et al. samples
lungobj_mye_yan <- subset(
  lungobj_mye,
  subset = StudyID == "Yan et al."
)

# Generate Yan et al. feature plots (Figure C, right)
for (gene in features) {

  p <- FeaturePlot_scCustom(
    lungobj_mye_yan,
    features = gene,
    reduction = "umap.harmony"
  ) +
    theme_umap +
    NoLegend()

  ggsave(
    filename = file.path(
      "figures",
      paste0("Yan_", gene, "_FeaturePlot.png")
    ),
    plot = p,
    dpi = 400,
    width = 2,
    height = 2
  )
}

# save reproducibility information
capture.output(
  sessionInfo(),
  file = file.path("figures", "sessionInfo_figures.txt")
)
