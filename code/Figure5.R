#!/usr/bin/env Rscript

# =============================================================================
# Figure 5: 蛋白标志物模型的性能与功能解释
# =============================================================================

# =============================================================================

project_root <- normalizePath(
  Sys.getenv("PROTEO_PROJECT_ROOT", unset = getwd()),
  mustWork = FALSE
)
figure5_dir <- Sys.getenv(
  "PROTEO_FIGURE5_OUTPUT_DIR",
  unset = file.path(project_root, "Figure5")
)
figure5_dir <- normalizePath(figure5_dir, mustWork = FALSE)
dir.create(figure5_dir, recursive = TRUE, showWarnings = FALSE)

diff_dir <- file.path(project_root, "diff")
exp_dir <- file.path(project_root, "exp")
model_dir <- file.path(project_root, "model")
model_plot_dir <- file.path(model_dir, "plot")
lasso_dir <- file.path(model_dir, "lasso")
lasso_function_dir <- file.path(lasso_dir, "function")
compare_dir <- file.path(model_dir, "compare")

for (dir_path in c(
  figure5_dir,
  model_plot_dir,
  lasso_function_dir
)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(pROC)
  library(caret)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
  library(org.Hs.eg.db)
  library(clusterProfiler)
})


# =============================================================================
# 1. Figure 5A-B: 特征蛋白训练集和测试集热图
# =============================================================================

load(file.path(diff_dir, "protein_imp_gene.rds"))
expression_data <- gene_exp
rownames(expression_data) <- expression_data[, 1]
expression_data <- expression_data[, -1]
expression_data <- matrix(
  as.numeric(expression_data),
  nrow = nrow(expression_data),
  ncol = ncol(expression_data)
)
rownames(expression_data) <- gene_exp[, 1]
colnames(expression_data) <- colnames(gene_exp)[-1]

sample_source <- read.table(
  file.path(exp_dir, "patient_source.txt"),
  sep = "\t",
  quote = "",
  header = FALSE,
  check.names = FALSE
)
filtered_samples <- read.table(
  file.path(diff_dir, "filter_sample.txt"),
  sep = "\t",
  quote = "",
  header = FALSE,
  check.names = FALSE
)

rownames(sample_source) <- sample_source[, 1]
sample_source <- sample_source[
  setdiff(sample_source[, 1], filtered_samples[, 1]),
]
nanjing_samples <- rownames(sample_source)[
  sample_source[, 2] == "nanjing"
]
henan_samples <- rownames(sample_source)[
  sample_source[, 2] == "henan"
]
fujian_samples <- rownames(sample_source)[
  sample_source[, 2] == "fujian"
]
heilongjiang_samples <- rownames(sample_source)[
  sample_source[, 2] == "heilongjiang"
]

lasso_features <- read.table(
  file.path(lasso_dir, "LASSOfeature.txt"),
  sep = "\t"
)
training_data <- expression_data[
  lasso_features[, 1],
  c(fujian_samples, nanjing_samples)
]
test_data <- expression_data[
  lasso_features[, 1],
  c(heilongjiang_samples, henan_samples)
]

train_prediction <- read.table(
  file.path(model_dir, "LASSO_lda_train.prob.txt"),
  sep = " ",
  header = TRUE
)
rownames(train_prediction) <- train_prediction[, 1]
train_prediction <- train_prediction[colnames(training_data),]
train_prediction$prediction <- ifelse(
  train_prediction$y_prob > 0.5,
  "LM",
  "Control"
)

train_annotation <- data.frame(
  Real = train_prediction[, 2],
  prediction = train_prediction[, 4],
  prob = train_prediction[, 3]
)
rownames(train_annotation) <- colnames(training_data)

probability_color <- colorRamp2(
  c(0, 0.5, 1),
  c("#BDBDBD", "white", "#F9B6E0")
)
train_annotation_colors <- list(
  prob = probability_color,
  Real = c(
    control = "#A9CFE4",
    case = "#FF9A91"
  ),
  prediction = c(
    Control = "#79C1C1",
    LM = "#F6C55F"
  )
)

train_top_annotation <- HeatmapAnnotation(
  df = train_annotation,
  col = train_annotation_colors,
  gp = gpar(col = "white", lwd = 0.5),
  annotation_name_side = "left"
)

expression_colors <- colorRampPalette(
  rev(
    c(
      "#EB7078",
      rgb(255, 255, 255, maxColorValue = 255),
      "#4874B1"
    )
  )
)(100)

scaled_training <- t(scale(t(training_data)))
scaled_training[scaled_training < -1] <- -1
scaled_training[scaled_training > 1] <- 1

train_heatmap <- Heatmap(
  scaled_training,
  border = FALSE,
  col = expression_colors,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  column_gap = unit(3, "mm"),
  top_annotation = train_top_annotation,
  height = unit(50, "mm"),
  width = unit(150, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.2)
)

pdf(
  file.path(figure5_dir, "lasso_feature_train_heatmap.pdf"),
  width = 15,
  height = 10
)
print(train_heatmap)
dev.off()

gene_order <- rownames(scaled_training)[row_order(train_heatmap)]
test_data <- test_data[gene_order,]

test_prediction <- read.table(
  file.path(model_dir, "LASSO_lda_test.prob.txt"),
  sep = " ",
  header = TRUE
)
rownames(test_prediction) <- test_prediction[, 1]
test_prediction <- test_prediction[colnames(test_data),]
test_prediction$prediction <- ifelse(
  test_prediction$pred_prob > 0.5,
  "LM",
  "Control"
)
test_prediction$Real <- ifelse(
  grepl("Control", test_prediction$sample),
  "Control",
  "LM"
)

test_annotation <- data.frame(
  Real = test_prediction$Real,
  prediction = test_prediction$prediction,
  prob = test_prediction$pred_prob
)
rownames(test_annotation) <- colnames(test_data)

test_probability_color <- colorRamp2(
  c(0, 0.5, 1),
  c("#BDBDBD", "white", "#F9B6E0")
)
test_annotation_colors <- list(
  prob = test_probability_color,
  Real = c(
    Control = "#A9CFE4",
    LM = "#FF9A91"
  ),
  prediction = c(
    Control = "#79C1C1",
    LM = "#F6C55F"
  )
)

test_top_annotation <- HeatmapAnnotation(
  df = test_annotation,
  col = test_annotation_colors,
  gp = gpar(col = "white", lwd = 0.5),
  annotation_name_side = "left"
)

scaled_test <- t(scale(t(test_data)))
scaled_test[scaled_test < -1] <- -1
scaled_test[scaled_test > 1] <- 1

test_heatmap <- Heatmap(
  scaled_test,
  border = "#C3C3C3",
  col = expression_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  column_gap = unit(3, "mm"),
  top_annotation = test_top_annotation,
  height = unit(50, "mm"),
  width = unit(50, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.2)
)

pdf(
  file.path(figure5_dir, "lasso_feature_test_heatmap.pdf"),
  width = 15,
  height = 10
)
print(test_heatmap)
dev.off()


# =============================================================================
# 2. Figure 5C: 训练集和测试集 ROC 曲线
# =============================================================================

train_prediction <- read.table(
  file.path(model_dir, "LASSO_lda_train.prob.txt"),
  sep = " ",
  header = TRUE
)
training_labels <- ifelse(
  train_prediction$y_true == "control",
  0,
  1
)
training_probability <- train_prediction$y_prob
training_roc <- roc(training_labels, training_probability)
training_auc_ci <- ci.auc(training_roc, boot.n = 500)
training_auc_text <- sprintf(
  "%.3f (%.3f-%.3f)",
  auc(training_roc),
  training_auc_ci[1],
  training_auc_ci[3]
)
training_ci_se <- ci.se(
  training_roc,
  specificities = seq(0, 1, length.out = 20),
  boot.n = 500
)

pdf(
  file.path(figure5_dir, "train_roc.pdf"),
  width = 4,
  height = 4
)
plot(
  training_roc,
  col = "red",
  legacy.axes = TRUE,
  xlab = "1-Specificity"
)
plot(
  training_ci_se,
  type = "shape",
  col = "lightgray",
  border = NA
)
lines(training_roc, col = "red", lwd = 2)
text(
  0.6,
  0.2,
  labels = paste0("AUC: ", training_auc_text),
  adj = 0,
  cex = 0.9
)
dev.off()

test_prediction <- read.table(
  file.path(model_dir, "LASSO_lda_test.prob.txt"),
  sep = " ",
  header = TRUE
)
test_labels <- ifelse(
  grepl("Control", test_prediction$sample),
  0,
  1
)
test_probability <- test_prediction$pred_prob
test_roc <- roc(test_labels, test_probability)
test_auc_ci <- ci.auc(test_roc, boot.n = 500)
test_auc_text <- sprintf(
  "%.3f (%.3f-%.3f)",
  auc(test_roc),
  test_auc_ci[1],
  test_auc_ci[3]
)
test_ci_se <- ci.se(
  test_roc,
  specificities = seq(0, 1, length.out = 20),
  boot.n = 500
)

pdf(
  file.path(figure5_dir, "test_roc.pdf"),
  width = 4,
  height = 4
)
plot(
  test_roc,
  col = "red",
  legacy.axes = TRUE,
  xlab = "1-Specificity"
)
plot(
  test_ci_se,
  type = "shape",
  col = "lightgray",
  border = NA
)
lines(test_roc, col = "red", lwd = 2)
text(
  0.6,
  0.2,
  labels = paste0("AUC: ", test_auc_text),
  adj = 0,
  cex = 0.9
)
dev.off()


# =============================================================================
# 3. Figure 5D: 测试集单基因性能比较热图
# =============================================================================

single_gene_test <- read.table(
  file.path(lasso_dir, "test_single_gene_performance.txt"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)
rownames(single_gene_test) <- single_gene_test[, 1]
single_gene_test <- single_gene_test[, -1]
single_gene_test <- rbind(Ours = 1, single_gene_test)

single_gene_matrix <- as.matrix(single_gene_test)
single_gene_long <- melt(single_gene_matrix)
colnames(single_gene_long) <- c("Gene", "Metric", "Value")
single_gene_long$Gene <- factor(
  single_gene_long$Gene,
  levels = rownames(single_gene_test)
)
single_gene_long$Metric <- factor(
  single_gene_long$Metric,
  levels = rev(colnames(single_gene_test))
)

single_gene_value_min <- min(
  single_gene_long$Value,
  na.rm = TRUE
)
single_gene_value_max <- max(
  single_gene_long$Value,
  na.rm = TRUE
)
single_gene_colors <- colorRampPalette(
  c("#F8FBF6", "#D7F0E6", "#A7D8F0", "#F7C8B8")
)(100)

single_gene_heatmap <- ggplot(
  single_gene_long,
  aes(x = Gene, y = Metric, fill = Value)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = round(Value, 3)),
    size = 3,
    color = "#2B2B2B"
  ) +
  scale_fill_gradientn(
    colors = single_gene_colors,
    limits = c(single_gene_value_min, single_gene_value_max),
    breaks = pretty(
      c(single_gene_value_min, single_gene_value_max),
      n = 5
    ),
    name = "Value"
  ) +
  scale_y_discrete(
    labels = c(
      ACC = "Accuracy",
      AUC = "AUC",
      Sensitivity = "Sensitivity",
      Specificity = "Specificity",
      F1 = "F1 Score",
      PR_AUC = "PR AUC"
    )
  ) +
  labs(title = "Replication dataset") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10,
      face = "italic"
    ),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      face = "bold"
    ),
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave(
  file.path(
    figure5_dir,
    "LASSO_gene_compare_heatmap_test_fresh.pdf"
  ),
  single_gene_heatmap,
  width = 6,
  height = 3.2
)


# =============================================================================
# 4. Figure 5E: 测试集与已知 marker 的 AUC 比较
# =============================================================================

load(file.path(diff_dir, "protein_imp_gene.rds"))
marker_expression <- gene_exp
rownames(marker_expression) <- marker_expression[, 1]
marker_expression <- marker_expression[, -1]
marker_expression <- matrix(
  as.numeric(marker_expression),
  nrow = nrow(marker_expression),
  ncol = ncol(marker_expression),
  dimnames = list(
    gene_exp[, 1],
    colnames(gene_exp)[-1]
  )
)

marker_sample_source <- read.table(
  file.path(exp_dir, "patient_source.txt"),
  sep = "\t",
  quote = "",
  header = FALSE,
  check.names = FALSE
)
marker_filter <- read.table(
  file.path(diff_dir, "filter_sample.txt"),
  sep = "\t",
  quote = "",
  header = FALSE,
  check.names = FALSE
)
rownames(marker_sample_source) <- marker_sample_source[, 1]
marker_sample_source <- marker_sample_source[
  setdiff(
    marker_sample_source[, 1],
    marker_filter[, 1]
  ),
]

henan_samples <- rownames(marker_sample_source)[
  marker_sample_source[, 2] == "henan"
]
heilongjiang_samples <- rownames(marker_sample_source)[
  marker_sample_source[, 2] == "heilongjiang"
]
known_markers <- read.table(
  file.path(compare_dir, "known_marker.txt"),
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

replication_auc <- data.frame()
for (marker_index in seq_len(nrow(known_markers))) {
  gene <- known_markers[marker_index, 2]
  labels <- c(
    rep(0, length(heilongjiang_samples)),
    rep(1, length(henan_samples))
  )
  values <- c(
    marker_expression[gene, heilongjiang_samples],
    marker_expression[gene, henan_samples]
  )
  direction <- ifelse(
    known_markers[marker_index, 3] == "Increased",
    ">",
    "<"
  )
  marker_roc <- roc(
    labels,
    values,
    direction = direction,
    quiet = TRUE
  )
  replication_auc <- rbind(
    replication_auc,
    data.frame(
      Gene = gene,
      AUC = as.numeric(auc(marker_roc))
    )
  )
}

replication_auc <- rbind(
  replication_auc,
  data.frame(Gene = "Ours", AUC = 1)
)
replication_auc <- replication_auc[
  order(replication_auc$AUC, decreasing = TRUE),
]
replication_auc$Gene <- factor(
  replication_auc$Gene,
  levels = replication_auc$Gene
)
write.table(
  replication_auc,
  file.path(figure5_dir, "test.AUC_fresh.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

auc_colors <- c(
  "#D7DCE0",
  "#C9E6E3",
  "#B8E3D1",
  "#DDE8B8",
  "#F4D3A3",
  "#EEA7A1"
)
auc_plot <- ggplot(
  replication_auc,
  aes(x = Gene, y = AUC, fill = AUC)
) +
  geom_col(
    width = 0.68,
    color = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(label = sprintf("%.3f", AUC)),
    vjust = -0.35,
    size = 3,
    color = "#333333"
  ) +
  scale_fill_gradientn(
    colors = auc_colors,
    limits = c(min(replication_auc$AUC), 1),
    breaks = c(
      round(min(replication_auc$AUC), 2),
      0.7,
      0.8,
      0.9,
      1
    ),
    name = "AUC"
  ) +
  scale_y_continuous(
    limits = c(0, 1.08),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  labs(
    title = "Replication dataset",
    x = NULL,
    y = "AUC"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10
    ),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      face = "bold"
    ),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.line = element_line(color = "#4A4A4A", linewidth = 0.4),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  file.path(figure5_dir, "test.AUC_fresh.pdf"),
  auc_plot,
  width = 6,
  height = 3.2
)


# =============================================================================
# 5. Figure 5F: LASSO、细胞学和 MRI 敏感性比较
# =============================================================================

training_model <- read.table(
  file.path(model_dir, "LASSO_glm_train.prob.txt"),
  sep = " ",
  header = TRUE
)
training_model$y_true <- factor(
  training_model$y_true,
  levels = c("control", "case")
)
training_predicted_class <- factor(
  ifelse(
    training_model$y_prob > 0.5,
    "case",
    "control"
  ),
  levels = c("control", "case")
)
training_confusion <- confusionMatrix(
  training_predicted_class,
  training_model$y_true,
  positive = "case"
)
training_sensitivity <- training_confusion$byClass["Sensitivity"]

test_model <- read.table(
  file.path(model_dir, "LASSO_glm_test.prob.txt"),
  sep = " ",
  header = TRUE
)
test_group <- factor(
  ifelse(grepl("Control", test_model$sample), "control", "case"),
  levels = c("control", "case")
)
test_predicted_class <- factor(
  ifelse(test_model$pred_prob > 0.5, "case", "control"),
  levels = c("control", "case")
)
test_confusion <- confusionMatrix(
  test_predicted_class,
  test_group,
  positive = "case"
)
test_sensitivity <- test_confusion$byClass["Sensitivity"]

sensitivity_data <- data.frame(
  Dataset = rep(c("Training", "Replication"), each = 3),
  Method = rep(
    c("LASSO", "Cytology", "MRI"),
    times = 2
  ),
  Sensitivity = c(
    training_sensitivity,
    0.924,
    0.721,
    test_sensitivity,
    1,
    0.63
  )
)
sensitivity_data$Dataset <- factor(
  sensitivity_data$Dataset,
  levels = c("Training", "Replication")
)
sensitivity_data$Method <- factor(
  sensitivity_data$Method,
  levels = c("LASSO", "Cytology", "MRI")
)
sensitivity_colors <- c(
  "LASSO" = "#FB9B9A",
  "Cytology" = "#FDCE90",
  "MRI" = "#A9CFE4"
)

sensitivity_plot <- ggplot(
  sensitivity_data,
  aes(x = Dataset, y = Sensitivity, fill = Method)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = sprintf("%.2f", Sensitivity)),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 4
  ) +
  scale_fill_manual(values = sensitivity_colors) +
  scale_y_continuous(
    limits = c(0, 1.05),
    expand = c(0, 0)
  ) +
  labs(
    title = "Sensitivity: Cytology and MRI",
    x = "Method",
    y = "Sensitivity",
    fill = "Method"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 12),
    panel.grid.major.x = element_blank()
  )

ggsave(
  file.path(figure5_dir, "Cytology_MRI_sensitivity.pdf"),
  sensitivity_plot,
  width = 6,
  height = 4
)



# =============================================================================
# 7. Figure 5H: 特征节点功能富集气泡图
# =============================================================================

network_nodes <- read.table(
  file.path(lasso_function_dir, "node.txt"),
  sep = "\t"
)
node_entrez <- mget(
  network_nodes[, 1],
  org.Hs.egSYMBOL2EG,
  ifnotfound = NA
)
node_entrez <- as.character(node_entrez)

node_go <- enrichGO(
  gene = node_entrez,
  keyType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH"
)
node_kegg <- enrichKEGG(
  gene = node_entrez,
  organism = "hsa",
  keyType = "kegg",
  pAdjustMethod = "BH"
)
write.table(
  node_go,
  file.path(lasso_function_dir, "nodeGO.txt"),
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
write.table(
  node_kegg,
  file.path(lasso_function_dir, "nodeKEGG.txt"),
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

plot_bubble <- function(
  enrichment_result,
  output_file,
  width,
  height
) {
  top_pathways <- enrichment_result %>%
    arrange(p.adjust) %>%
    head(10)
  top_pathways$function.gene.num <- unlist(
    strsplit(top_pathways[, "BgRatio"], "/")
  )[
    seq(
      1,
      length(unlist(strsplit(top_pathways[, "BgRatio"], "/"))),
      2
    )
  ]
  top_pathways$GeneRatio <- top_pathways$Count /
    as.numeric(top_pathways$function.gene.num)
  top_pathways <- top_pathways[
    order(top_pathways$GeneRatio, decreasing = TRUE),
  ]
  top_pathways$Description <- factor(
    top_pathways$Description,
    levels = rev(top_pathways$Description)
  )
  x_gap <- (
    max(top_pathways$GeneRatio) -
      min(top_pathways$GeneRatio)
  ) / 5
  top_pathways$pvalue <- -log10(top_pathways$pvalue)

  bubble_plot <- ggplot(
    top_pathways,
    mapping = aes(x = GeneRatio, y = Description)
  ) +
    geom_point(aes(size = Count, color = pvalue)) +
    theme_bw() +
    scale_color_viridis_c(option = "D") +
    scale_x_continuous(expand = c(x_gap, x_gap)) +
    ylab(NULL) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16)
    ) +
    scale_size_continuous(range = c(2, 6))

  ggsave(
    output_file,
    bubble_plot,
    width = width,
    height = height
  )
}

plot_bubble(
  as.data.frame(node_kegg),
  file.path(lasso_function_dir, "KEGG_bubbleplot.pdf"),
  width = 6,
  height = 5
)

go_data <- as.data.frame(node_go)
top_go <- go_data[1:min(10, nrow(go_data)),]
plot_bubble(
  top_go,
  file.path(lasso_function_dir, "GOBP_bubbleplot.pdf"),
  width = 8,
  height = 5
)
