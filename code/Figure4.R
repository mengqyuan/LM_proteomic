#!/usr/bin/env Rscript

# =============================================================================
# Figure 4: 蛋白、脂质和氨基酸的多组学功能与相关性
# =============================================================================
#
# 本文件已整合 Figure 4 的实际分析代码：
#
#   Figure 4A    蛋白/脂质/氨基酸差异表达热图
#   Figure 4B-E  关键通路蛋白 log2FC 热图
#   Figure 4F    hsa00071 蛋白-脂质相关性
#
# 路径通过环境变量配置：
#   PROTEO_PROJECT_ROOT
#   PROTEO_FIGURE4_OUTPUT_DIR
# =============================================================================

project_root <- normalizePath(
  Sys.getenv("PROTEO_PROJECT_ROOT", unset = getwd()),
  mustWork = FALSE
)
figure4_dir <- Sys.getenv(
  "PROTEO_FIGURE4_OUTPUT_DIR",
  unset = file.path(project_root, "Figure4")
)
figure4_dir <- normalizePath(figure4_dir, mustWork = FALSE)
dir.create(figure4_dir, recursive = TRUE, showWarnings = FALSE)

diff_dir <- file.path(project_root, "diff")
hsa00071_dir <- file.path(project_root, "hsa00071_pathway_rerun")
hsa00380_dir <- file.path(project_root, "hsa00380_pathway_rerun")
hsa04979_dir <- file.path(project_root, "hsa04979_pathway_rerun")
hsa00860_dir <- file.path(project_root, "hsa00860_pathway")

suppressPackageStartupMessages({
  library(pheatmap)
})


# =============================================================================
# 1. Figure 4A: 三条通路的多组学热图
# =============================================================================

mat_files <- c(
  FA = file.path(
    hsa00071_dir,
    "hsa00071_protein_lipid_heatmap_matrix_raw.tsv"
  ),
  TRP = file.path(
    hsa00380_dir,
    "hsa00380_protein_amino_heatmap_matrix_raw.tsv"
  ),
  CHOL = file.path(
    hsa04979_dir,
    "hsa04979_protein_lipid_heatmap_matrix_raw.tsv"
  )
)
pathway_name <- c(
  FA = "Fatty acid degradation",
  TRP = "Tryptophan metabolism",
  CHOL = "Cholesterol metabolism"
)

read_matrix <- function(file_path) {
  expression_matrix <- read.delim(
    file_path,
    row.names = 1,
    check.names = FALSE
  )
  expression_matrix <- as.matrix(expression_matrix)
  mode(expression_matrix) <- "numeric"
  expression_matrix
}

get_rows <- function(expression_matrix, feature_type) {
  row_ids <- rownames(expression_matrix)[
    grepl(
      paste0("^", feature_type, "\\|"),
      rownames(expression_matrix)
    )
  ]
  sub(
    paste0("^", feature_type, "\\|"),
    "",
    row_ids
  )
}

sort_samples <- function(samples) {
  prefix <- sub("\\.[0-9]+$", "", samples)
  number <- as.numeric(sub("^.*\\.", "", samples))
  samples[order(prefix, number)]
}

row_zscore <- function(expression_matrix) {
  z_score <- t(scale(t(expression_matrix)))
  z_score[is.na(z_score)] <- 0
  z_score[z_score > 2.5] <- 2.5
  z_score[z_score < -2.5] <- -2.5
  z_score
}

add_unique <- function(old, new) {
  c(old, setdiff(new, old))
}

pathway_matrices <- lapply(mat_files, read_matrix)
sample_order <- Reduce(
  intersect,
  lapply(pathway_matrices, colnames)
)
normal_samples <- sort_samples(
  sample_order[grepl("^Control", sample_order)]
)
cancer_samples <- sort_samples(
  sample_order[!grepl("^Control", sample_order)]
)
sample_order <- c(normal_samples, cancer_samples)

protein_order <- character(0)
protein_order <- add_unique(
  protein_order,
  get_rows(pathway_matrices$FA, "Protein")
)
protein_order <- add_unique(
  protein_order,
  get_rows(pathway_matrices$TRP, "Protein")
)
protein_order <- add_unique(
  protein_order,
  get_rows(pathway_matrices$CHOL, "Protein")
)

fa_lipid_order <- get_rows(pathway_matrices$FA, "Lipid")
amino_order <- get_rows(pathway_matrices$TRP, "Amino acid")
chol_lipid_order <- setdiff(
  get_rows(pathway_matrices$CHOL, "Lipid"),
  fa_lipid_order
)
lipid_order <- c(fa_lipid_order, chol_lipid_order)

find_expression <- function(feature, feature_type) {
  row_id <- paste0(feature_type, "|", feature)
  for (expression_matrix in pathway_matrices) {
    if (row_id %in% rownames(expression_matrix)) {
      return(expression_matrix[row_id, sample_order])
    }
  }
}

get_pathways <- function(feature, feature_type) {
  row_id <- paste0(feature_type, "|", feature)
  paste(
    pathway_name[
      vapply(
        pathway_matrices,
        function(expression_matrix) {
          row_id %in% rownames(expression_matrix)
        },
        logical(1)
      )
    ],
    collapse = "; "
  )
}

features <- data.frame(
  type = c(
    rep("Protein", length(protein_order)),
    rep("Lipid", length(fa_lipid_order)),
    rep("Amino acid", length(amino_order)),
    rep("Lipid", length(chol_lipid_order))
  ),
  feature = c(
    protein_order,
    fa_lipid_order,
    amino_order,
    chol_lipid_order
  ),
  stringsAsFactors = FALSE
)
features$pathway <- mapply(
  get_pathways,
  features$feature,
  features$type
)
features$label <- paste0(
  features$type,
  "|",
  features$feature
)

heatmap_raw <- t(mapply(
  find_expression,
  features$feature,
  features$type
))
rownames(heatmap_raw) <- features$label
colnames(heatmap_raw) <- sample_order
heatmap_z <- row_zscore(heatmap_raw)

sample_annotation <- data.frame(
  Group = factor(
    ifelse(
      grepl("^Control", sample_order),
      "Normal",
      "Cancer"
    ),
    levels = c("Normal", "Cancer")
  ),
  row.names = sample_order
)
row_annotation <- data.frame(
  Type = features$type,
  Pathway = features$pathway,
  row.names = features$label
)
annotation_colors <- list(
  Group = c(
    Normal = "#4DBBD5",
    Cancer = "#E64B35"
  ),
  Type = c(
    Protein = "#3C5488",
    Lipid = "#F39B7F",
    `Amino acid` = "#00A087"
  )
)

fa_protein_count <- length(
  get_rows(pathway_matrices$FA, "Protein")
)
trp_protein_count <- length(
  setdiff(
    get_rows(pathway_matrices$TRP, "Protein"),
    get_rows(pathway_matrices$FA, "Protein")
  )
)
chol_protein_count <- length(
  setdiff(
    get_rows(pathway_matrices$CHOL, "Protein"),
    c(
      get_rows(pathway_matrices$FA, "Protein"),
      get_rows(pathway_matrices$TRP, "Protein")
    )
  )
)
row_gaps <- c(
  fa_protein_count,
  fa_protein_count + trp_protein_count,
  length(protein_order),
  length(protein_order) + length(fa_lipid_order),
  length(protein_order) +
    length(fa_lipid_order) +
    length(amino_order)
)
row_gaps <- row_gaps[
  row_gaps > 0 & row_gaps < nrow(heatmap_z)
]

write.table(
  heatmap_raw,
  file.path(
    figure4_dir,
    "three_pathway_protein_lipid_amino_heatmap_matrix_raw.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  heatmap_z,
  file.path(
    figure4_dir,
    "three_pathway_protein_lipid_amino_heatmap_matrix_zscore.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  row_annotation,
  file.path(
    figure4_dir,
    "three_pathway_protein_lipid_amino_row_annotation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  sample_annotation,
  file.path(
    figure4_dir,
    "three_pathway_protein_lipid_amino_sample_annotation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

pheatmap(
  heatmap_z,
  color = colorRampPalette(
    c("#2166AC", "#F7F7F7", "#B2182B")
  )(101),
  breaks = seq(-2.5, 2.5, length.out = 102),
  annotation_col = sample_annotation,
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  gaps_col = length(normal_samples),
  gaps_row = row_gaps,
  show_colnames = FALSE,
  fontsize_row = 6,
  border_color = NA,
  main = paste(
    "hsa00071, hsa00380 and hsa04979",
    "proteins/lipids/amino acids"
  ),
  filename = file.path(
    figure4_dir,
    "three_pathway_protein_lipid_amino_pheatmap.pdf"
  ),
  width = 13,
  height = 12
)


# =============================================================================
# 2. Figure 4B-E: 关键通路蛋白 log2FC 热图
# =============================================================================

diff_file <- file.path(
  diff_dir,
  "case_control_filter_diff.txt"
)
pathway_files <- c(
  hsa00071 = file.path(
    hsa00071_dir,
    "hsa00071_protein_hits.tsv"
  ),
  hsa00380 = file.path(
    hsa00380_dir,
    "hsa00380_protein_hits.tsv"
  ),
  hsa04979 = file.path(
    hsa04979_dir,
    "hsa04979_protein_hits.tsv"
  ),
  hsa00860 = file.path(
    hsa00860_dir,
    "hsa00860_protein_hits.tsv"
  )
)
logfc_pdf <- file.path(
  figure4_dir,
  "three_pathway_protein_log2FC_heatmap.pdf"
)
logfc_table_path <- file.path(
  figure4_dir,
  "three_pathway_protein_log2FC_heatmap_table.tsv"
)
logfc_annotation_path <- file.path(
  figure4_dir,
  "three_pathway_protein_log2FC_heatmap_row_annotation.tsv"
)

pathway_priority <- names(pathway_files)
differential_result <- read.delim(
  diff_file,
  row.names = 1,
  check.names = FALSE
)

combined_pathways <- do.call(
  rbind,
  lapply(
    names(pathway_files),
    function(pathway_id) {
      pathway_table <- read.delim(
        pathway_files[[pathway_id]],
        check.names = FALSE
      )
      data.frame(
        Gene = pathway_table$symbol,
        Pathway = pathway_id,
        stringsAsFactors = FALSE
      )
    }
  )
)
combined_pathways <- combined_pathways[
  combined_pathways$Gene %in% rownames(differential_result),
]

pathways_by_gene <- aggregate(
  Pathway ~ Gene,
  combined_pathways,
  function(values) {
    paste(unique(values), collapse = "; ")
  }
)
pathways_by_gene$logFC <- differential_result[
  pathways_by_gene$Gene,
  "logFC"
]
pathways_by_gene$Primary <- vapply(
  strsplit(pathways_by_gene$Pathway, ";"),
  function(pathways) {
    pathways <- trimws(pathways)
    pathway_priority[
      which(pathway_priority %in% pathways)[1]
    ]
  },
  character(1)
)
pathways_by_gene$Primary <- factor(
  pathways_by_gene$Primary,
  levels = pathway_priority
)
pathways_by_gene <- pathways_by_gene[
  order(
    pathways_by_gene$Primary,
    -pathways_by_gene$logFC,
    pathways_by_gene$Gene
  ),
]
rownames(pathways_by_gene) <- pathways_by_gene$Gene

logfc_matrix <- matrix(
  pathways_by_gene$logFC,
  ncol = 1,
  dimnames = list(
    pathways_by_gene$Gene,
    "log2FC"
  )
)

pathway_levels <- c(
  "hsa00071",
  "hsa00071; hsa00380",
  "hsa00380",
  "hsa04979",
  "hsa00860"
)
pathways_by_gene$Pathway <- factor(
  pathways_by_gene$Pathway,
  levels = pathway_levels[
    pathway_levels %in% unique(pathways_by_gene$Pathway)
  ]
)
logfc_row_annotation <- data.frame(
  Pathway = pathways_by_gene$Pathway,
  row.names = pathways_by_gene$Gene,
  check.names = FALSE
)

minimum_logfc <- min(
  pathways_by_gene$logFC,
  na.rm = TRUE
)
maximum_logfc <- max(
  pathways_by_gene$logFC,
  na.rm = TRUE
)
logfc_breaks <- c(
  seq(minimum_logfc, -2, length.out = 41),
  seq(-2, -0.46, length.out = 181)[-1],
  seq(-0.46, 0, length.out = 31)[-1],
  seq(0, maximum_logfc, length.out = 41)[-1]
)
logfc_colors <- c(
  colorRampPalette(c("#0B5D1E", "#3E8E41"))(40),
  colorRampPalette(c("#3E8E41", "#A7D7A9"))(180),
  colorRampPalette(c("#A7D7A9", "#F7F7F7"))(30),
  colorRampPalette(c("#F7F7F7", "#F5A3A3", "#C62828"))(40)
)

row_group_sizes <- table(pathways_by_gene$Primary)
logfc_row_gaps <- cumsum(as.integer(row_group_sizes))
logfc_row_gaps <- logfc_row_gaps[
  logfc_row_gaps < nrow(logfc_matrix)
]

write.table(
  data.frame(
    Gene = rownames(logfc_matrix),
    logFC = pathways_by_gene$logFC,
    Primary = as.character(pathways_by_gene$Primary),
    Pathway = as.character(pathways_by_gene$Pathway),
    row.names = NULL
  ),
  logfc_table_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  logfc_row_annotation,
  logfc_annotation_path,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

pheatmap(
  logfc_matrix,
  color = logfc_colors,
  breaks = logfc_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = logfc_row_gaps,
  annotation_row = logfc_row_annotation,
  annotation_colors = list(
    Pathway = c(
      "hsa00071" = "#8ECFC9",
      "hsa00071; hsa00380" = "#E7C36A",
      "hsa00380" = "#F4A261",
      "hsa04979" = "#B9A5E8",
      "hsa00860" = "#64B5C3"
    )
  ),
  display_numbers = matrix(
    sprintf("%.2f", pathways_by_gene$logFC),
    ncol = 1
  ),
  number_color = "#2B2B2B",
  fontsize_number = 7,
  fontsize_row = 8,
  fontsize_col = 11,
  angle_col = 0,
  border_color = NA,
  main = paste(
    "Differential protein log2FC across",
    "hsa00071 / hsa00380 / hsa04979 / hsa00860"
  ),
  legend = TRUE,
  filename = logfc_pdf,
  width = 6.6,
  height = 7.8
)


# =============================================================================
# 3. Figure 4F: hsa00071 蛋白-脂质相关性
# =============================================================================

correlation_input <- file.path(
  hsa00071_dir,
  "hsa00071_protein_lipid_heatmap_matrix_raw.tsv"
)
correlation_output <- file.path(
  hsa00071_dir,
  "hsa00071_protein_lipid_pearson_spearman_cor.tsv"
)

correlation_expression <- read.delim(
  correlation_input,
  row.names = 1,
  check.names = FALSE
)
correlation_expression <- as.matrix(correlation_expression)
mode(correlation_expression) <- "numeric"

protein_expression <- correlation_expression[
  grepl("^Protein\\|", rownames(correlation_expression)),
  ,
  drop = FALSE
]
lipid_expression <- correlation_expression[
  grepl("^Lipid\\|", rownames(correlation_expression)),
  ,
  drop = FALSE
]
rownames(protein_expression) <- sub(
  "^Protein\\|",
  "",
  rownames(protein_expression)
)
rownames(lipid_expression) <- sub(
  "^Lipid\\|",
  "",
  rownames(lipid_expression)
)

correlation_results <- data.frame()
for (protein in rownames(protein_expression)) {
  for (lipid in rownames(lipid_expression)) {
    protein_values <- as.numeric(
      protein_expression[protein, ]
    )
    lipid_values <- as.numeric(
      lipid_expression[lipid, ]
    )
    valid_values <- is.finite(protein_values) &
      is.finite(lipid_values)

    pearson_result <- cor.test(
      protein_values[valid_values],
      lipid_values[valid_values],
      method = "pearson"
    )
    spearman_result <- suppressWarnings(
      cor.test(
        protein_values[valid_values],
        lipid_values[valid_values],
        method = "spearman",
        exact = FALSE
      )
    )

    correlation_results <- rbind(
      correlation_results,
      data.frame(
        protein = protein,
        lipid = lipid,
        pearson_cor = unname(pearson_result$estimate),
        pearson_p = pearson_result$p.value,
        spearman_cor = unname(spearman_result$estimate),
        spearman_p = spearman_result$p.value
      )
    )
  }
}

write.table(
  correlation_results,
  correlation_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
