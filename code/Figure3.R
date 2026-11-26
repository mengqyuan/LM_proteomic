#!/usr/bin/env Rscript

# =============================================================================
# Figure 3: CSF 三套单细胞数据集的整合、分群和上皮功能分析
# =============================================================================

# =============================================================================

project_root <- normalizePath(
  Sys.getenv("PROTEO_PROJECT_ROOT", unset = getwd()),
  mustWork = FALSE
)
shared_data_dir <- Sys.getenv(
  "PROTEO_SHARED_DATA_DIR",
  unset = file.path(project_root, "shared_data")
)
workspace_dir <- Sys.getenv(
  "PROTEO_WORKSPACE_DIR",
  unset = file.path(project_root, "workspace")
)
results_dir <- Sys.getenv(
  "PROTEO_RESULTS_DIR",
  unset = file.path(project_root, "results")
)
upload_dir <- Sys.getenv(
  "PROTEO_UPLOAD_DIR",
  unset = file.path(project_root, "uploads")
)
single_cell_dir <- file.path(project_root, "single_cell", "new")
epi_dir <- file.path(single_cell_dir, "epi")
diff_dir <- file.path(project_root, "diff")

for (dir_path in c(
  workspace_dir,
  results_dir,
  single_cell_dir,
  epi_dir
)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(reshape2)
  library(svglite)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(GSVA)
  library(ggpubr)
  library(BiocParallel)
  library(ggnewscale)
  library(patchwork)
  library(grid)
})


# =============================================================================
# 1. Figure 3a-f, step 1: 读取三套 CSF scRNA-seq 数据集
# =============================================================================

message("Reading CSF scRNA-seq datasets")

# GSE134577: healthy controls, 10x matrices.
hc_root <- file.path(shared_data_dir, "GSE134577_HC", "extracted")
hc_samples <- sort(list.dirs(
  hc_root,
  full.names = TRUE,
  recursive = FALSE
))
hc_samples <- hc_samples[basename(hc_samples) != basename(hc_root)]

hc_list <- list()
for (sample_dir in hc_samples) {
  sample_id <- basename(sample_dir)
  counts <- Read10X(data.dir = sample_dir)
  rownames(counts) <- make.unique(rownames(counts))
  sample_object <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  sample_object$stim <- sample_id
  sample_object$dataset <- "GSE134577"
  sample_object$stimgroup <- "HC"
  sample_object <- RenameCells(
    sample_object,
    new.names = paste0(sample_id, "_", colnames(sample_object))
  )
  hc_list[[sample_id]] <- sample_object
}
hc_obj <- merge(hc_list[[1]], hc_list[-1])
rm(hc_list)
gc()

# GSE276139: LUAD leptomeningeal metastasis, 10x matrices.
lm_root <- file.path(
  shared_data_dir,
  "GSE276139_LUAD_LM",
  "extracted"
)
lm_samples <- sort(list.dirs(
  lm_root,
  full.names = TRUE,
  recursive = FALSE
))
lm_samples <- lm_samples[basename(lm_samples) != basename(lm_root)]

lm_list <- list()
for (sample_dir in lm_samples) {
  sample_id <- basename(sample_dir)
  counts <- Read10X(data.dir = sample_dir)
  rownames(counts) <- make.unique(rownames(counts))
  sample_object <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  sample_object$stim <- sample_id
  sample_object$dataset <- "GSE276139"
  sample_object$stimgroup <- "LM"
  sample_object <- RenameCells(
    sample_object,
    new.names = paste0(sample_id, "_", colnames(sample_object))
  )
  lm_list[[sample_id]] <- sample_object
}
lm_obj <- merge(lm_list[[1]], lm_list[-1])
rm(lm_list)
gc()

# GSE150660: mixed leptomeningeal metastasis, cells x genes CSV files.
gse150660_dir <- file.path(shared_data_dir, "GSE150660_LM")
csv_files <- sort(list.files(
  gse150660_dir,
  pattern = "_Patient_.*\\.csv\\.gz$",
  full.names = TRUE
))

lm2_list <- list()
for (csv_file in csv_files) {
  file_name <- basename(csv_file)
  patient_id <- sub(
    ".*_(Patient_[A-E])\\.csv\\.gz$",
    "\\1",
    file_name
  )
  patient_data <- read.csv(
    csv_file,
    row.names = 1,
    check.names = FALSE
  )
  patient_data <- patient_data[
    ,
    colnames(patient_data) != "CLUSTER",
    drop = FALSE
  ]
  counts <- as.matrix(t(patient_data))
  rownames(counts) <- make.unique(rownames(counts))
  sample_object <- CreateSeuratObject(
    counts = counts,
    project = patient_id,
    min.cells = 3,
    min.features = 200
  )
  sample_object$stim <- patient_id
  sample_object$dataset <- "GSE150660"
  sample_object$stimgroup <- "LM"
  sample_object <- RenameCells(
    sample_object,
    new.names = paste0(patient_id, "_", colnames(sample_object))
  )
  lm2_list[[patient_id]] <- sample_object
}
lm2_obj <- merge(lm2_list[[1]], lm2_list[-1])
rm(lm2_list, patient_data, counts)
gc()

saveRDS(hc_obj, file.path(workspace_dir, "hc_obj.rds"))
saveRDS(lm_obj, file.path(workspace_dir, "lm_obj.rds"))
saveRDS(lm2_obj, file.path(workspace_dir, "lm2_obj.rds"))


# =============================================================================
# 2. Figure 3a-f, step 2: 分数据集质量控制
# =============================================================================

message("Filtering cells by quality metrics")

hc_obj <- readRDS(file.path(workspace_dir, "hc_obj.rds"))
hc_obj <- PercentageFeatureSet(
  hc_obj,
  pattern = "^MT-",
  col.name = "percent.mt"
)
hc_obj <- subset(
  hc_obj,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    nCount_RNA > 500 &
    percent.mt < 20
)
saveRDS(hc_obj, file.path(workspace_dir, "hc_qc.rds"))
rm(hc_obj)
gc()

lm_obj <- readRDS(file.path(workspace_dir, "lm_obj.rds"))
lm_obj <- PercentageFeatureSet(
  lm_obj,
  pattern = "^MT-",
  col.name = "percent.mt"
)
lm_obj <- subset(
  lm_obj,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    nCount_RNA > 500 &
    percent.mt < 20
)
saveRDS(lm_obj, file.path(workspace_dir, "lm_qc.rds"))
rm(lm_obj)
gc()

lm2_obj <- readRDS(file.path(workspace_dir, "lm2_obj.rds"))
lm2_obj <- PercentageFeatureSet(
  lm2_obj,
  pattern = "^MT-",
  col.name = "percent.mt"
)
lm2_obj <- subset(
  lm2_obj,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    nCount_RNA > 500 &
    percent.mt < 20
)
saveRDS(lm2_obj, file.path(workspace_dir, "lm2_qc.rds"))
rm(lm2_obj)
gc()


# =============================================================================
# 3. Figure 3a-f, step 3: 共同基因、合并和低信息基因过滤
# =============================================================================

message("Merging quality-controlled datasets")

hc_obj <- readRDS(file.path(workspace_dir, "hc_qc.rds"))
lm_obj <- readRDS(file.path(workspace_dir, "lm_qc.rds"))
lm2_obj <- readRDS(file.path(workspace_dir, "lm2_qc.rds"))

common_genes <- Reduce(
  intersect,
  list(rownames(hc_obj), rownames(lm_obj), rownames(lm2_obj))
)
hc_obj <- subset(hc_obj, features = common_genes)
lm_obj <- subset(lm_obj, features = common_genes)
lm2_obj <- subset(lm2_obj, features = common_genes)

All <- merge(hc_obj, list(lm_obj, lm2_obj))
rm(hc_obj, lm_obj, lm2_obj)
gc()

remove_features <- unique(c(
  grep("^MT-", rownames(All)),
  grep("^RPL", rownames(All)),
  grep("^RPS", rownames(All))
))
if (length(remove_features) > 0L) {
  All <- subset(All, features = rownames(All)[-remove_features])
}

saveRDS(All, file.path(workspace_dir, "All_merged_raw.rds"))
rm(All)
gc()


# =============================================================================
# 4. Figure 3a-f, step 4: Harmony、降维和聚类
# =============================================================================
#
# 该步骤不能省略：前面只生成原始合并对象；本步骤生成后续注释和图 3
# 所依赖的 harmony、umap、tsne 和 cluster 信息。

message("Normalizing, reducing dimensions and clustering")

All <- readRDS(file.path(workspace_dir, "All_merged_raw.rds"))
All <- NormalizeData(
  All,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
All <- FindVariableFeatures(
  All,
  selection.method = "vst",
  nfeatures = 2000
)
All <- ScaleData(All, verbose = FALSE)
All <- RunPCA(
  All,
  features = VariableFeatures(All),
  verbose = FALSE
)

save_svg <- function(plot_object, file_name, width = 12, height = 5) {
  ggsave(
    file.path(results_dir, file_name),
    plot_object,
    device = svglite::svglite,
    width = width,
    height = height,
    limitsize = FALSE
  )
}

pre_harmony_stim <- DimPlot(
  All,
  reduction = "pca",
  pt.size = 0.1,
  group.by = "stim"
) +
  ggtitle("Pre-Harmony PCA by sample")
pre_harmony_pc1 <- VlnPlot(
  All,
  features = "PC_1",
  group.by = "stim",
  pt.size = 0.1
) +
  ggtitle("PC1 by sample")
save_svg(
  plot_grid(pre_harmony_stim, pre_harmony_pc1),
  "CSF_3dataset_pre_harmony_PCA.svg",
  width = 20,
  height = 5
)

pre_harmony_dataset <- DimPlot(
  All,
  reduction = "pca",
  pt.size = 0.1,
  group.by = "dataset"
) +
  ggtitle("Pre-Harmony PCA by dataset")
pre_harmony_dataset_pc1 <- VlnPlot(
  All,
  features = "PC_1",
  group.by = "dataset",
  pt.size = 0.1
) +
  ggtitle("PC1 by dataset")
save_svg(
  plot_grid(pre_harmony_dataset, pre_harmony_dataset_pc1),
  "CSF_3dataset_pre_harmony_PCA_dataset.svg",
  width = 20,
  height = 5
)

All <- RunHarmony(
  All,
  group.by.vars = c("dataset", "stim"),
  plot_convergence = FALSE
)

post_harmony_stim <- DimPlot(
  All,
  reduction = "harmony",
  pt.size = 0.1,
  group.by = "stim"
) +
  ggtitle("Post-Harmony by sample")
post_harmony_pc1 <- VlnPlot(
  All,
  features = "harmony_1",
  group.by = "stim",
  pt.size = 0.1
) +
  ggtitle("Harmony 1 by sample")
save_svg(
  plot_grid(post_harmony_stim, post_harmony_pc1),
  "CSF_3dataset_post_harmony_stim.svg",
  width = 20,
  height = 5
)

post_harmony_dataset <- DimPlot(
  All,
  reduction = "harmony",
  pt.size = 0.1,
  group.by = "dataset"
) +
  ggtitle("Post-Harmony by dataset")
post_harmony_dataset_pc1 <- VlnPlot(
  All,
  features = "harmony_1",
  group.by = "dataset",
  pt.size = 0.1
) +
  ggtitle("Harmony 1 by dataset")
save_svg(
  plot_grid(post_harmony_dataset, post_harmony_dataset_pc1),
  "CSF_3dataset_post_harmony_dataset.svg",
  width = 20,
  height = 5
)

All <- RunUMAP(All, reduction = "harmony", dims = 1:30)
All <- RunTSNE(All, reduction = "harmony", dims = 1:30)
All <- FindNeighbors(All, reduction = "harmony", dims = 1:30)
All <- FindClusters(All, resolution = 0.8)

saveRDS(
  All,
  file.path(workspace_dir, "All_harmony_final.rds")
)
saveRDS(
  All,
  file.path(results_dir, "CSF_3dataset_harmony_merged.rds")
)

save_svg(
  DimPlot(
    All,
    reduction = "umap",
    label = TRUE,
    pt.size = 0.1,
    raster = FALSE
  ) +
    ggtitle("UMAP - Clusters"),
  "CSF_3dataset_UMAP_cluster.svg",
  width = 10,
  height = 8
)
save_svg(
  DimPlot(
    All,
    reduction = "umap",
    group.by = "stimgroup",
    pt.size = 0.1,
    raster = FALSE
  ) +
    ggtitle("UMAP - HC vs LM"),
  "CSF_3dataset_UMAP_stimgroup.svg",
  width = 10,
  height = 8
)
save_svg(
  DimPlot(
    All,
    reduction = "umap",
    group.by = "dataset",
    pt.size = 0.1,
    raster = FALSE
  ) +
    ggtitle("UMAP - Dataset"),
  "CSF_3dataset_UMAP_dataset.svg",
  width = 10,
  height = 8
)
save_svg(
  DimPlot(
    All,
    reduction = "tsne",
    label = TRUE,
    pt.size = 0.1,
    raster = FALSE
  ) +
    ggtitle("tSNE - Clusters"),
  "CSF_3dataset_tSNE_cluster.svg",
  width = 10,
  height = 8
)

cluster_group_counts <- as.matrix(table(Idents(All), All$stimgroup))
cluster_group_prop <- prop.table(cluster_group_counts, margin = 2)
cluster_group_long <- melt(cluster_group_prop)
colnames(cluster_group_long) <- c("cluster", "group", "value")
save_svg(
  ggplot(
    cluster_group_long,
    aes(
      x = factor(cluster, levels = rownames(cluster_group_prop)),
      y = value,
      fill = group
    )
  ) +
    geom_bar(
      stat = "identity",
      position = "dodge",
      width = 0.6
    ) +
    theme(axis.text = element_text(size = 20)) +
    ggtitle("Cluster composition by group"),
  "CSF_3dataset_cluster_barplot.svg",
  width = 27,
  height = 12
)

cluster_stat_path <- file.path(
  workspace_dir,
  "CSF_3dataset_cluster_stat.txt"
)
sink(cluster_stat_path)
cat("=== Cluster x sample (counts) ===\n")
print(table(Idents(All), All$stim))
cat("\n=== Cluster x group (counts) ===\n")
print(table(Idents(All), All$stimgroup))
cat("\n=== Cluster x dataset (counts) ===\n")
print(table(Idents(All), All$dataset))
sink()
file.copy(
  cluster_stat_path,
  file.path(results_dir, basename(cluster_stat_path)),
  overwrite = TRUE
)

cluster_markers <- FindAllMarkers(
  All,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  only.pos = TRUE
)
write.table(
  cluster_markers,
  file.path(workspace_dir, "CSF_3dataset_all_markers.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
write.table(
  top_markers,
  file.path(workspace_dir, "CSF_3dataset_top20_markers.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)


# =============================================================================
# 5. Figure 3a-f, step 5: 细胞注释、比例比较和功能评分
# =============================================================================

message("Annotating cell types and comparing cell proportions")

All <- readRDS(
  file.path(results_dir, "CSF_3dataset_harmony_merged.rds")
)
All <- subset(
  All,
  idents = c("14", "20", "23"),
  invert = TRUE
)

cluster2type <- c(
  "0" = "Macrophage",
  "1" = "CD4 T",
  "2" = "Macrophage",
  "3" = "CD8 T",
  "4" = "Epithelial",
  "5" = "Monocyte",
  "6" = "Epithelial",
  "7" = "DC",
  "8" = "DC",
  "9" = "Monocyte",
  "10" = "CD4 Treg",
  "11" = "CD4 T",
  "12" = "Macrophage",
  "13" = "Macrophage",
  "15" = "NK",
  "16" = "Epithelial",
  "17" = "Macrophage",
  "18" = "Epithelial",
  "19" = "Plasma",
  "21" = "B",
  "22" = "Macrophage",
  "24" = "CD4 T",
  "25" = "Macrophage"
)
All$celltype <- unname(
  cluster2type[as.character(Idents(All))]
)

celltype_colors <- c(
  "CD4 T" = "#4E79A7",
  "CD4 Treg" = "#A0CBE8",
  "CD8 T" = "#76B7B2",
  "NK" = "#59A14F",
  "B" = "#FF9DA7",
  "Plasma" = "#EDC948",
  "Monocyte" = "#F28E2B",
  "Macrophage" = "#E15759",
  "DC" = "#B07AA1",
  "Epithelial" = "#9C755F"
)

celltype_umap <- DimPlot(
  All,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  pt.size = 0.1,
  raster = FALSE,
  cols = celltype_colors,
  repel = TRUE
) +
  ggtitle("UMAP - Cell types")
ggsave(
  file.path(single_cell_dir, "CSF_3dataset_UMAP_celltype.pdf"),
  celltype_umap,
  width = 6,
  height = 5
)

final_object_path <- file.path(
  single_cell_dir,
  "All_harmony_final.rds"
)
saveRDS(All, final_object_path)

case_object <- subset(All, stimgroup == "LM")
control_object <- subset(All, stimgroup == "HC")

case_umap <- DimPlot(
  case_object,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  pt.size = 0.1,
  raster = FALSE,
  cols = celltype_colors,
  repel = TRUE
) +
  ggtitle("UMAP - Cell types")
ggsave(
  file.path(single_cell_dir, "CSF_3dataset_UMAP_celltype_case.pdf"),
  case_umap,
  width = 10,
  height = 8
)

control_umap <- DimPlot(
  control_object,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  pt.size = 0.1,
  raster = FALSE,
  cols = celltype_colors,
  repel = TRUE
) +
  ggtitle("UMAP - Cell types")
ggsave(
  file.path(single_cell_dir, "CSF_3dataset_UMAP_celltype_control.pdf"),
  control_umap,
  width = 10,
  height = 8
)

metadata <- All@meta.data
cell_counts <- as.data.frame(table(metadata$stim, metadata$celltype))
colnames(cell_counts) <- c("stim", "celltype", "n")
sample_totals <- as.data.frame(table(metadata$stim))
colnames(sample_totals) <- c("stim", "total")
proportion_data <- merge(cell_counts, sample_totals, by = "stim")
proportion_data$prop <- proportion_data$n / proportion_data$total
sample_groups <- unique(metadata[, c("stim", "stimgroup", "dataset")])
proportion_data <- merge(proportion_data, sample_groups, by = "stim")
proportion_data$stimgroup <- factor(
  proportion_data$stimgroup,
  levels = c("HC", "LM")
)

proportion_stats <- proportion_data %>%
  group_by(celltype) %>%
  summarise(
    median_HC = median(prop[stimgroup == "HC"]),
    median_LM = median(prop[stimgroup == "LM"]),
    mean_HC = mean(prop[stimgroup == "HC"]),
    mean_LM = mean(prop[stimgroup == "LM"]),
    p_value = wilcox.test(
      prop[stimgroup == "LM"],
      prop[stimgroup == "HC"],
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    padj_BH = p.adjust(p_value, method = "BH"),
    log2FC = log2((mean_LM + 1e-4) / (mean_HC + 1e-4)),
    sig = case_when(
      padj_BH < 0.0001 ~ "****",
      padj_BH < 0.001 ~ "***",
      padj_BH < 0.01 ~ "**",
      padj_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

write.csv(
  proportion_stats,
  file.path(single_cell_dir, "CSF_3dataset_proportion_test.csv"),
  row.names = FALSE
)
write.csv(
  proportion_data[
    ,
    c("stim", "dataset", "stimgroup", "celltype", "n", "total", "prop")
  ],
  file.path(
    single_cell_dir,
    "CSF_3dataset_per_sample_proportions.csv"
  ),
  row.names = FALSE
)

celltype_order <- c(
  "Epithelial",
  "CD8 T",
  "CD4 T",
  "CD4 Treg",
  "NK",
  "Monocyte",
  "Macrophage",
  "DC",
  "B",
  "Plasma"
)
proportion_data$celltype <- factor(
  proportion_data$celltype,
  levels = celltype_order
)
max_proportions <- proportion_data %>%
  group_by(celltype) %>%
  summarise(y = max(prop) * 1.08, .groups = "drop")
stats_plot_data <- merge(
  proportion_stats,
  max_proportions,
  by = "celltype"
)
group_colors <- c("HC" = "#4E79A7", "LM" = "#E15759")

proportion_plot <- ggplot(
  proportion_data,
  aes(x = celltype, y = prop, fill = stimgroup)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.7,
    position = position_dodge(0.75)
  ) +
  geom_point(
    aes(color = stimgroup),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 1.8,
    alpha = 0.9
  ) +
  geom_text(
    data = stats_plot_data,
    aes(x = celltype, y = y, label = sig),
    inherit.aes = FALSE,
    size = 6,
    vjust = 0
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(
    x = NULL,
    y = "Proportion of cells",
    fill = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top"
  )
ggsave(
  file.path(single_cell_dir, "CSF_3dataset_proportion_boxplot.pdf"),
  proportion_plot,
  width = 6,
  height = 3
)

celltype_markers <- c(
  "KRT18",
  "KRT8",
  "CD3E",
  "CD8A",
  "CD8B",
  "CD4",
  "IL7R",
  "FOXP3",
  "TRDC",
  "KLRC1",
  "GNLY",
  "VCAN",
  "FCN1",
  "C1QA",
  "C1QB",
  "APOE",
  "CLEC10A",
  "FCER1A",
  "CD79A",
  "CD79B",
  "JCHAIN",
  "MZB1",
  "MS4A1",
  "IGHM"
)
celltype_markers <- celltype_markers[
  celltype_markers %in% rownames(All)
]
All$celltype <- factor(All$celltype, levels = rev(celltype_order))

marker_dotplot <- DotPlot(
  All,
  features = celltype_markers,
  group.by = "celltype"
) +
  scale_color_gradient(low = "olivedrab", high = "gold") +
  labs(
    x = NULL,
    y = NULL,
    color = "Average Expression",
    size = "Percent Expressed"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12)
  )
ggsave(
  file.path(
    single_cell_dir,
    "CSF_3dataset_celltype_marker_dotplot_ordered.pdf"
  ),
  marker_dotplot,
  width = 8,
  height = 4
)

All <- readRDS(final_object_path)
celltype_markers_table <- FindAllMarkers(
  All,
  group.by = "celltype",
  min.pct = 0.25,
  logfc.threshold = 0.5,
  only.pos = TRUE
)
top_celltype_markers <- celltype_markers_table %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)
write.table(
  top_celltype_markers,
  file.path(single_cell_dir, "celltype_top20_markers.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)

load(file.path(diff_dir, "protein_imp_gene.rds"))
protein_expression <- gene_exp
rownames(protein_expression) <- protein_expression[, 1]
protein_expression <- protein_expression[, -1]
protein_expression <- matrix(
  as.numeric(protein_expression),
  nrow = nrow(protein_expression),
  ncol = ncol(protein_expression)
)
rownames(protein_expression) <- gene_exp[, 1]
colnames(protein_expression) <- colnames(gene_exp)[-1]

filtered_samples <- read.table(
  file.path(diff_dir, "filter_sample.txt"),
  sep = "\t",
  quote = "",
  header = FALSE,
  check.names = FALSE
)
keep_samples <- setdiff(
  colnames(protein_expression),
  filtered_samples[, 1]
)
protein_expression <- protein_expression[, keep_samples]

marker_sets <- read.table(
  file.path(single_cell_dir, "celltype_top20_markers.txt"),
  sep = "\t",
  header = TRUE
)
marker_sets <- marker_sets[
  marker_sets[, "gene"] %in% rownames(protein_expression),
]
marker_gene_list <- split(
  as.character(marker_sets[, "gene"]),
  as.character(marker_sets[, "cluster"])
)
ssgsea_parameters <- ssgseaParam(
  protein_expression,
  marker_gene_list,
  minSize = 3
)
ssGSEA <- gsva(ssgsea_parameters)
save(
  ssGSEA,
  file = file.path(single_cell_dir, "ssGSEA_50.rds")
)


# =============================================================================
# 6. Figure 3g, step 1: 上皮细胞差异分析与功能富集
# =============================================================================

message("Analyzing epithelial-cell differential expression")

All <- readRDS(final_object_path)
Idents(All) <- All$celltype
epithelial_de <- FindMarkers(
  All,
  ident.1 = "Epithelial",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25
)
epithelial_de$gene <- rownames(epithelial_de)
epithelial_de <- epithelial_de %>%
  arrange(p_val_adj, desc(avg_log2FC))
write.csv(
  epithelial_de,
  file.path(results_dir, "CSF_3dataset_epithelial_vs_others_DE.csv"),
  row.names = FALSE
)

upregulated_genes <- epithelial_de %>%
  filter(p_val_adj < 0.01, avg_log2FC > 0.5)
if (nrow(upregulated_genes) > 3000L) {
  upregulated_genes <- epithelial_de %>%
    filter(p_val_adj < 0.01, avg_log2FC > 1)
}
write.csv(
  upregulated_genes,
  file.path(epi_dir, "CSF_3dataset_epithelial_up_genes.csv"),
  row.names = FALSE
)

object_universe <- rownames(All)
upregulated_entrez <- bitr(
  upregulated_genes$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
background_entrez <- bitr(
  object_universe,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

kegg_result <- enrichKEGG(
  gene = upregulated_entrez$ENTREZID,
  organism = "hsa",
  universe = background_entrez$ENTREZID,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
kegg_table <- setReadable(
  kegg_result,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
) %>%
  as.data.frame() %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust)
write.csv(
  kegg_table,
  file.path(epi_dir, "CSF_3dataset_epithelial_KEGG.csv"),
  row.names = FALSE
)

hallmark_collection <- msigdbr(
  species = "Homo sapiens",
  category = "H"
) %>%
  dplyr::select(gs_name, gene_symbol)
hallmark_result <- enricher(
  upregulated_genes$gene,
  TERM2GENE = hallmark_collection,
  universe = object_universe,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
hallmark_table <- as.data.frame(hallmark_result) %>%
  arrange(p.adjust)
write.csv(
  hallmark_table,
  file.path(epi_dir, "CSF_3dataset_epithelial_Hallmark.csv"),
  row.names = FALSE
)


# =============================================================================
# 7. Figure 3g, step 2: KEGG 和 Hallmark 功能条形图
# =============================================================================

plot_enrichment_bar <- function(
  data,
  output_path,
  fill_color,
  width = 5.5,
  height = NULL,
  bar_width = 0.55
) {
  data <- data %>%
    mutate(
      Description = gsub(
        "_",
        " ",
        gsub("^HALLMARK_", "", Description)
      ),
      Description = paste0(
        toupper(substr(Description, 1, 1)),
        tolower(substr(Description, 2, nchar(Description)))
      ),
      Description = dplyr::recode(
        Description,
        "E2f targets" = "E2F targets",
        "G2m checkpoint" = "G2M checkpoint",
        "Mtorc1 signaling" = "mTORC1 signaling",
        "Igsf cam signaling" = "IgSF CAM signaling"
      ),
      neglog10padj = -log10(p.adjust)
    ) %>%
    arrange(p.adjust) %>%
    mutate(Description = factor(Description, levels = rev(Description)))

  if (is.null(height)) {
    height <- max(2.5, 0.3 * nrow(data) + 1.3)
  }

  plot_object <- ggplot(
    data,
    aes(x = neglog10padj, y = Description)
  ) +
    geom_col(fill = fill_color, width = bar_width) +
    geom_text(
      aes(label = paste0(Count, " genes")),
      hjust = 1.15,
      size = 3.5,
      color = "grey20"
    ) +
    labs(
      x = expression(-log[10](adjusted~italic(P))),
      y = NULL
    ) +
    theme_classic(base_size = 12) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02)))

  ggsave(
    output_path,
    plot_object,
    width = width,
    height = height
  )
}

kegg_plot_data <- read.csv(
  file.path(epi_dir, "CSF_3dataset_epithelial_KEGG.csv")
)
hallmark_plot_data <- read.csv(
  file.path(epi_dir, "CSF_3dataset_epithelial_Hallmark.csv")
)

plot_enrichment_bar(
  kegg_plot_data,
  file.path(epi_dir, "CSF_3dataset_epithelial_KEGG_barplot.pdf"),
  fill_color = "#A9CFE4"
)
plot_enrichment_bar(
  hallmark_plot_data,
  file.path(epi_dir, "CSF_3dataset_epithelial_Hallmark_barplot.pdf"),
  fill_color = "#FB9B9A",
  width = 4.5
)

