#!/usr/bin/env Rscript

# =============================================================================
# Figure 6: 蛋白组分子亚型的识别、生存和功能特征
# =============================================================================

# =============================================================================

project_root <- normalizePath(
  Sys.getenv("PROTEO_PROJECT_ROOT", unset = getwd()),
  mustWork = FALSE
)
figure6_dir <- Sys.getenv(
  "PROTEO_FIGURE6_OUTPUT_DIR",
  unset = file.path(project_root, "Figure6")
)
figure6_dir <- normalizePath(figure6_dir, mustWork = FALSE)

diff_dir <- file.path(project_root, "diff")
exp_dir <- file.path(project_root, "exp")
subtype_dir <- file.path(project_root, "subtype")
nmf_dir <- file.path(subtype_dir, "NMF")
nmf_plot_dir <- file.path(nmf_dir, "plot")
go_dir <- file.path(nmf_dir, "GO_KEGG")
single_cell_dir <- file.path(project_root, "single_cell", "new")

for (dir_path in c(
  figure6_dir,
  nmf_dir,
  nmf_plot_dir,
  go_dir
)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(openxlsx)
  library(NMF)
  library(CancerSubtypes)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(ggrepel)
  library(readxl)
  library(tidyr)
  library(svglite)
  library(pheatmap)
  library(ggnewscale)
  library(patchwork)
  library(rlang)
  library(reshape2)
  library(ggridges)
  library(ggpubr)
  library(gghalves)
  library(msigdbr)
})

# 在 pheatmap 中只保留指定基因标签，并绘制引线。
add.flag <- function(pheatmap_object, kept.labels, repel.degree) {
  heatmap_gtable <- pheatmap_object$gtable
  label_grob <- heatmap_gtable$grobs[[
    which(heatmap_gtable$layout$name == "row_names")
  ]]
  label_grob$label <- ifelse(
    label_grob$label %in% kept.labels,
    label_grob$label,
    ""
  )

  repel_positions <- function(distances, selected, spacing = repel.degree) {
    strip_npc <- function(distance) {
      if (!"unit.arithmetic" %in% class(distance)) {
        return(as.numeric(distance))
      }
      first_value <- strip_npc(distance$arg1)
      second_value <- strip_npc(distance$arg2)
      lazyeval::lazy_eval(
        paste(first_value, distance$fname, second_value)
      )
    }

    full_range <- vapply(
      seq_along(distances),
      function(index) strip_npc(distances[index]),
      numeric(1)
    )
    selected_distances <- distances[selected]
    selected_range <- vapply(
      seq_along(selected_distances),
      function(index) strip_npc(selected_distances[index]),
      numeric(1)
    )

    unit(
      seq(
        from = max(selected_range) +
          spacing * (max(full_range) - max(selected_range)),
        to = min(selected_range) -
          spacing * (min(selected_range) - min(full_range)),
        length.out = sum(selected)
      ),
      "npc"
    )
  }

  selected_labels <- label_grob$label != ""
  new_positions <- repel_positions(
    label_grob$y,
    selected_labels
  )
  label_flags <- segmentsGrob(
    x0 = label_grob$x,
    x1 = label_grob$x + unit(0.15, "npc"),
    y0 = label_grob$y[selected_labels],
    y1 = new_positions
  )

  label_grob$x <- label_grob$x + unit(0.2, "npc")
  label_grob$y[selected_labels] <- new_positions
  heatmap_gtable <- gtable::gtable_add_grob(
    x = heatmap_gtable,
    grobs = label_flags,
    t = 4,
    l = 4
  )
  heatmap_gtable$grobs[[
    which(heatmap_gtable$layout$name == "row_names")
  ]] <- label_grob

  grid.newpage()
  grid.draw(heatmap_gtable)
  invisible(heatmap_gtable)
}


# =============================================================================
# 1. Figure 6A-D: NMF 亚型识别和生存分析
# =============================================================================

survival_plot <- function(
  output_file,
  survival_data,
  survival_fit,
  plot_title
) {
  pdf(output_file, width = 4, height = 5)
  fit_plot <- ggsurvplot(
    data = survival_data,
    survival_fit,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    linetype = "strata",
    surv.median.line = "hv",
    ggtheme = theme_bw(),
    palette = c("#EF7B52", "#3C83C4", "#E7161A", "#9C4D9D"),
    title = plot_title,
    font.main = 10,
    font.x = 10,
    font.y = 10,
    font.tickslab = 10,
    risk.table.height = 0.3,
    pval.size = 3,
    risk.table.fontsize = 4,
    font.legend = 8
  )
  print(fit_plot)
  dev.off()
  fit_plot
}

run_nmf_subtype <- function(top_feature, number_of_subtypes) {
  top_feature <- as.numeric(top_feature)
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

  expression_data <- expression_data[
    ,
    -grep("Control", colnames(expression_data)),
    drop = FALSE
  ]
  expression_data <- expression_data[
    ,
    -grep("MM", colnames(expression_data)),
    drop = FALSE
  ]
  expression_data <- expression_data[
    ,
    -grep("EC", colnames(expression_data)),
    drop = FALSE
  ]

  filtered_samples <- read.table(
    file.path(diff_dir, "filter_sample.txt"),
    sep = "\t",
    quote = "",
    header = FALSE,
    check.names = FALSE
  )
  keep_samples <- setdiff(
    colnames(expression_data),
    filtered_samples[, 1]
  )
  expression_data <- expression_data[, keep_samples]

  selected_data <- FSbyMAD(
    expression_data,
    cut.type = "topk",
    value = top_feature
  )
  selected_data <- data.normalization(
    selected_data,
    type = "feature_zscore",
    log2 = FALSE
  )
  write.table(
    rownames(selected_data),
    file.path(nmf_dir, paste0("gene_", top_feature, ".txt")),
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
  )

  selected_data[selected_data < 0] <- 0
  rank_estimation <- nmf(
    selected_data,
    ranks = 2:6,
    nrun = 100,
    method = "lee"
  )
  pdf(
    file.path(
      nmf_dir,
      paste0("gene_2_6", top_feature, "estim.coad.pdf")
    )
  )
  plot(rank_estimation)
  dev.off()

  nmf_fit <- nmf(
    selected_data,
    rank = 2,
    nrun = 200,
    seed = 20250529,
    method = "lee"
  )
  subtype_group <- predict(nmf_fit)
  write.table(
    subtype_group,
    file.path(
      nmf_dir,
      paste0(
        "gene_",
        top_feature,
        "_group_",
        number_of_subtypes,
        "_NMF_group.txt"
      )
    ),
    sep = "\t",
    quote = FALSE,
    col.names = FALSE
  )

  clinical_data <- read.xlsx(
    file.path(exp_dir, "patients_clinical_information.xlsx"),
    startRow = 2
  )
  clinical_data[is.na(clinical_data)] <- "N"
  rownames(clinical_data) <- clinical_data[, 1]

  group_labels <- paste("c", subtype_group, sep = "")
  group_samples <- data.frame(
    sample = names(subtype_group),
    group = group_labels
  )

  os_data <- clinical_data[, c("Patient", "OS_time", "OS_status")]
  os_data <- na.omit(os_data)
  os_survival <- merge(
    os_data,
    group_samples,
    all = FALSE,
    by.x = "Patient",
    by.y = "sample"
  )
  os_fit <- survfit(
    Surv(OS_time, OS_status) ~ group,
    data = os_survival
  )
  bc_os_data <- os_survival[grep("BC", os_survival[, 1]),]
  la_os_data <- os_survival[grep("LA", os_survival[, 1]),]
  bc_os_fit <- survfit(
    Surv(OS_time, OS_status) ~ group,
    data = bc_os_data
  )
  la_os_fit <- survfit(
    Surv(OS_time, OS_status) ~ group,
    data = la_os_data
  )

  base_name <- paste("subtype", number_of_subtypes, top_feature, sep = "_")
  survival_plot(
    file.path(nmf_dir, paste0(base_name, "_OS_survival.pdf")),
    os_survival,
    os_fit,
    base_name
  )
  survival_plot(
    file.path(
      nmf_dir,
      paste0(
        "BC_subtype_",
        number_of_subtypes,
        "_",
        top_feature,
        "_OS_survival.pdf"
      )
    ),
    bc_os_data,
    bc_os_fit,
    paste("BC", base_name, sep = "_")
  )
  survival_plot(
    file.path(
      nmf_dir,
      paste0(
        "LA_subtype_",
        number_of_subtypes,
        "_",
        top_feature,
        "_OS_survival.pdf"
      )
    ),
    la_os_data,
    la_os_fit,
    paste("LA", base_name, sep = "_")
  )

  ipfs_data <- clinical_data[
    ,
    c("Patient", "iPFS_time", "iPFS_status")
  ]
  ipfs_data <- na.omit(ipfs_data)
  ipfs_survival <- merge(
    ipfs_data,
    group_samples,
    all = FALSE,
    by.x = "Patient",
    by.y = "sample"
  )
  ipfs_fit <- survfit(
    Surv(iPFS_time, iPFS_status) ~ group,
    data = ipfs_survival
  )
  bc_ipfs_data <- ipfs_survival[grep("BC", ipfs_survival[, 1]),]
  la_ipfs_data <- ipfs_survival[grep("LA", ipfs_survival[, 1]),]
  bc_ipfs_fit <- survfit(
    Surv(iPFS_time, iPFS_status) ~ group,
    data = bc_ipfs_data
  )
  la_ipfs_fit <- survfit(
    Surv(iPFS_time, iPFS_status) ~ group,
    data = la_ipfs_data
  )

  survival_plot(
    file.path(nmf_dir, paste0(base_name, "_iPFS_survival.pdf")),
    ipfs_survival,
    ipfs_fit,
    base_name
  )
  survival_plot(
    file.path(
      nmf_dir,
      paste0(
        "BC_subtype_",
        number_of_subtypes,
        "_",
        top_feature,
        "_iPFS_survival.pdf"
      )
    ),
    bc_ipfs_data,
    bc_ipfs_fit,
    paste("BC", base_name, sep = "_")
  )
  survival_plot(
    file.path(
      nmf_dir,
      paste0(
        "LA_subtype_",
        number_of_subtypes,
        "_",
        top_feature,
        "_iPFS_survival.pdf"
      )
    ),
    la_ipfs_data,
    la_ipfs_fit,
    paste("LA", base_name, sep = "_")
  )
}

run_nmf_subtype(400, 2)

# Figure 6A: NMF subtype landscape.
load(file.path(diff_dir, "protein_imp_gene.rds"))
expression_data <- gene_exp
rownames(expression_data) <- expression_data[, 1]
expression_data <- expression_data[, -1]
expression_data <- matrix(
  as.numeric(as.matrix(expression_data)),
  nrow = nrow(expression_data)
)
rownames(expression_data) <- gene_exp[, 1]
colnames(expression_data) <- colnames(gene_exp)[-1]

clinical_data <- read.xlsx(
  file.path(exp_dir, "patients_clinical_information.xlsx"),
  startRow = 2
)
clinical_data[is.na(clinical_data)] <- "N"
rownames(clinical_data) <- clinical_data[, 1]

nmf_group <- read.table(
  file.path(nmf_dir, "gene_400_group_2_NMF_group.txt"),
  sep = "\t"
)
rownames(nmf_group) <- nmf_group[, 1]
common_samples <- intersect(
  colnames(expression_data),
  nmf_group[, 1]
)
expression_data <- expression_data[
  ,
  common_samples,
  drop = FALSE
]

top_genes <- read.table(
  file.path(nmf_dir, "gene_400.txt"),
  sep = "\t"
)[, 1]
remaining_genes <- setdiff(rownames(expression_data), top_genes)

landscape_clinical <- clinical_data[common_samples,]
landscape_clinical$Group <- nmf_group[
  landscape_clinical$Patient,
  2
]
landscape_clinical$Group <- paste(
  "c",
  landscape_clinical$Group,
  sep = ""
)
landscape_clinical <- landscape_clinical %>%
  arrange(Group)
sample_order <- landscape_clinical$Patient

landscape_clinical <- landscape_clinical %>%
  mutate(
    Group = factor(Group),
    Response = factor(
      Response,
      levels = c("PR", "SD", "PD", "N")
    ),
    Age = as.numeric(Age),
    iPFS_time = as.numeric(ifelse(iPFS_time == "N", NA, iPFS_time)),
    OS_time = as.numeric(ifelse(OS_time == "N", NA, OS_time))
  )

clinical_variables <- c(
  "Group",
  "Source",
  "Age",
  "Sex",
  "Subtype",
  "iPFS_status",
  "iPFS_time",
  "OS_status",
  "OS_time",
  "Response",
  "MRI",
  "Cytology",
  "Cancer"
)
clinical_annotation <- landscape_clinical %>%
  select(all_of(clinical_variables))
rownames(clinical_annotation) <- landscape_clinical$Patient

missing_color <- "#EDEDED"
age_color <- colorRamp2(
  range(landscape_clinical$Age, na.rm = TRUE),
  c("#4575B4", "#D73027")
)
ipfs_color <- colorRamp2(
  c(0, max(landscape_clinical$iPFS_time, na.rm = TRUE)),
  c("white", "#762A83")
)
os_color <- colorRamp2(
  c(0, max(landscape_clinical$OS_time, na.rm = TRUE)),
  c("white", "#A50F15")
)
clinical_colors <- list(
  Group = c(c1 = "#FF9A91", c2 = "#CDE2F2"),
  Source = c(Henan = "#CFE2F3", Nanjing_Brain = "#FFF2CC"),
  Sex = c(Male = "#8FAADC", Female = "#EAC7CC"),
  Cancer = c(
    BC = "#8ECFC9",
    LA = "#FFBE7A",
    EC = "#FA7F6F",
    MM = "#82B0D2",
    N = missing_color
  ),
  Subtype = c(
    Adenocarcinoma = "#A3D1D1",
    TNBC = "#FADBD8",
    `HER2+` = "#D4EFDF",
    Other = "#B2B2B2",
    `Luminal A/B` = "#FDEBD0",
    N = missing_color
  ),
  Response = c(
    PR = "#00A087",
    SD = "#F39B7F",
    PD = "#DC0000",
    N = missing_color
  ),
  MRI = c(
    Positive = "#71CF8A",
    Negative = "#BFC9D4",
    N = missing_color
  ),
  Cytology = c(
    Positive = "#7C90C9",
    Negative = "#EAC7CC",
    N = missing_color
  ),
  iPFS_status = c("0" = "#C6D8E1", "1" = "#D9B4B5", N = missing_color),
  OS_status = c("0" = "#C6D8E1", "1" = "#D9B4B5", N = missing_color)
)

annotation_heatmap <- HeatmapAnnotation(
  df = clinical_annotation,
  col = c(
    list(
      Age = age_color,
      iPFS_time = ipfs_color,
      OS_time = os_color
    ),
    clinical_colors
  ),
  na_col = missing_color,
  gp = gpar(col = "white", lwd = 0.5),
  annotation_name_side = "left"
)

annotation_strip <- Heatmap(
  matrix(0, nrow = 1, ncol = nrow(clinical_annotation)),
  col = "white",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  column_gap = unit(2, "mm"),
  top_annotation = annotation_heatmap,
  height = unit(3, "mm"),
  width = unit(2, "mm") * nrow(clinical_annotation),
  rect_gp = gpar(col = NA)
)

top_matrix <- expression_data[top_genes, sample_order]
rest_matrix <- expression_data[remaining_genes, sample_order]
top_matrix <- t(scale(t(top_matrix)))
rest_matrix <- t(scale(t(rest_matrix)))
top_matrix[is.na(top_matrix)] <- 0
rest_matrix[is.na(rest_matrix)] <- 0

expression_color <- colorRamp2(
  c(-2, 0, 2),
  c("#A555EC", "white", "#FFD915")
)
top_heatmap <- Heatmap(
  top_matrix,
  col = expression_color,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  show_column_names = FALSE,
  column_split = landscape_clinical$Group,
  column_gap = unit(0, "mm"),
  row_title = "TOP400",
  height = unit(2.5, "cm"),
  width = unit(2, "mm") * ncol(top_matrix),
  border = NA,
  rect_gp = gpar(col = NA, lwd = 0)
)
rest_heatmap <- Heatmap(
  rest_matrix,
  col = expression_color,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  show_column_names = FALSE,
  column_split = landscape_clinical$Group,
  column_gap = unit(0, "mm"),
  row_title = paste0("Remaining (n=", nrow(rest_matrix), ")"),
  height = unit(5, "cm"),
  width = unit(2, "mm") * ncol(rest_matrix),
  border = NA,
  rect_gp = gpar(col = NA, lwd = 0)
)

pdf(
  file.path(figure6_dir, "NMF_final_landscape_heatmap.pdf"),
  width = 10,
  height = 14
)
draw(
  annotation_strip %v% top_heatmap %v% rest_heatmap,
  ht_gap = unit(c(2, 2), "mm"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE
)
dev.off()


# =============================================================================
# 2. Figure 6E: 亚型 Cox 回归森林图
# =============================================================================

nmf_df <- read.delim(
  file.path(nmf_dir, "gene_400_group_2_NMF_group.txt"),
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(nmf_df) <- c("Patient", "NMF_group")
nmf_df$NMF_group <- factor(
  nmf_df$NMF_group,
  levels = c(1, 2),
  labels = c("Group 1", "Group 2")
)

clinical_cox <- read_excel(
  file.path(exp_dir, "patients_clinical_information.xlsx"),
  sheet = 1,
  skip = 1
)
cox_data <- merge(
  nmf_df,
  clinical_cox,
  by = "Patient",
  all.x = TRUE
)
cox_data$OS_status <- as.integer(cox_data$OS_status)
cox_data$OS_time <- as.numeric(cox_data$OS_time)
cox_data$iPFS_status <- as.integer(cox_data$iPFS_status)
cox_data$iPFS_time <- as.numeric(cox_data$iPFS_time)
cox_data$Age <- as.numeric(cox_data$Age)
cox_data$Sex <- factor(cox_data$Sex, levels = c("Male", "Female"))
cox_data$Cancer <- factor(
  cox_data$Cancer,
  levels = c("LA", "BC")
)
cox_data$NMF_group_rev <- relevel(
  cox_data$NMF_group,
  ref = "Group 2"
)
cox_data$Response[cox_data$Response == "PD"] <- "SD"
cox_data$Response <- factor(
  cox_data$Response,
  levels = c("PR", "SD")
)

os_cox <- coxph(
  Surv(OS_time, OS_status) ~
    NMF_group_rev + Age + Cancer + Sex + Response,
  data = cox_data
)
ipfs_cox <- coxph(
  Surv(iPFS_time, iPFS_status) ~
    NMF_group_rev + Age + Cancer + Sex + Response,
  data = cox_data
)

extract_forest_data <- function(model, endpoint) {
  model_summary <- summary(model)
  data.frame(
    endpoint = endpoint,
    variable = rownames(model_summary$coefficients),
    HR = model_summary$coefficients[, "exp(coef)"],
    lower = model_summary$conf.int[, "lower .95"],
    upper = model_summary$conf.int[, "upper .95"],
    pval = model_summary$coefficients[, "Pr(>|z|)"]
  )
}

forest_data <- rbind(
  extract_forest_data(os_cox, "OS"),
  extract_forest_data(ipfs_cox, "iPFS")
)
forest_data$variable <- gsub(
  "NMF_group_revGroup 1",
  "NMF Group 1 vs 2",
  forest_data$variable
)
forest_data$variable <- gsub(
  "CancerBC",
  "Cancer: BC vs LA",
  forest_data$variable
)
forest_data$variable <- gsub(
  "SexFemale",
  "Sex: Female vs Male",
  forest_data$variable
)
forest_data$variable <- gsub(
  "Age",
  "Age (per year)",
  forest_data$variable
)
forest_data$label <- sprintf(
  "%.2f (%.2f-%.2f)",
  forest_data$HR,
  forest_data$lower,
  forest_data$upper
)
forest_data$pval_label <- ifelse(
  forest_data$pval < 0.001,
  "< 0.001",
  sprintf("%.3f", forest_data$pval)
)
forest_data$y_label <- paste0(
  forest_data$endpoint,
  " - ",
  forest_data$variable
)
forest_data$y <- nrow(forest_data):1

forest_plot <- ggplot(
  forest_data,
  aes(x = HR, y = y)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.2,
    linewidth = 0.6
  ) +
  geom_point(
    size = 3,
    shape = 21,
    fill = ifelse(
      forest_data$endpoint == "OS",
      "#0072B2",
      "#D55E00"
    ),
    color = "black",
    stroke = 0.5
  ) +
  geom_text(
    aes(
      x = max(forest_data$upper) * 1.15,
      label = label
    ),
    hjust = 0,
    size = 3.2
  ) +
  geom_text(
    aes(
      x = max(forest_data$upper) * 2.5,
      label = pval_label
    ),
    hjust = 0,
    size = 3.2,
    fontface = ifelse(
      forest_data$pval < 0.05,
      "bold",
      "plain"
    )
  ) +
  scale_y_continuous(
    breaks = forest_data$y,
    labels = forest_data$y_label,
    expand = c(0, 0.5)
  ) +
  scale_x_log10() +
  coord_cartesian(
    xlim = c(0.2, max(forest_data$upper) * 4)
  ) +
  labs(
    title = "Adjusted Cox Regression: Hazard Ratios",
    subtitle = "Reference: Group 2, LA, Male",
    x = "Hazard Ratio (log scale)",
    y = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    axis.text.y = element_text(size = 9, hjust = 1),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 120, 10, 10)
  )

ggsave(
  file.path(figure6_dir, "Figure6E_cox_forest.pdf"),
  forest_plot,
  width = 12,
  height = 5
)


# =============================================================================
# 3. Figure 6F: 亚型差异蛋白火山图
# =============================================================================

differential_data <- read.table(
  file.path(nmf_dir, "c1_c2_diff.txt"),
  header = TRUE,
  sep = "\t"
)
differential_data$Gene <- rownames(differential_data)
differential_data <- differential_data %>%
  mutate(
    abs_logFC = abs(logFC),
    group = case_when(
      adj.P.Val < 0.01 & logFC > 1 ~ "IPC1 up",
      adj.P.Val < 0.01 & logFC < -1 ~ "IPC2 up",
      adj.P.Val < 0.01 & abs(logFC) <= 1 ~
        "FDR < 0.01 & |log2 FC| < 1",
      adj.P.Val >= 0.01 & abs(logFC) > 1 ~
        "FDR > 0.01 & |log2 FC| > 1",
      TRUE ~ "No significant"
    )
  )

top_up <- differential_data %>%
  filter(group == "IPC1 up") %>%
  arrange(adj.P.Val) %>%
  head(5)
top_down <- differential_data %>%
  filter(group == "IPC2 up") %>%
  arrange(adj.P.Val) %>%
  head(5)
labeled_genes <- bind_rows(top_up, top_down)
number_up <- sum(differential_data$group == "IPC1 up")
number_down <- sum(differential_data$group == "IPC2 up")

volcano_plot <- ggplot(
  differential_data,
  aes(x = logFC, y = -log10(adj.P.Val))
) +
  annotate(
    "rect",
    xmin = 1,
    xmax = Inf,
    ymin = -log10(0.001),
    ymax = Inf,
    fill = "#FEE0D2",
    alpha = 0.2
  ) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = -1,
    ymin = -log10(0.001),
    ymax = Inf,
    fill = "#DEEBF7",
    alpha = 0.2
  ) +
  geom_point(
    aes(color = group, size = abs_logFC),
    alpha = 0.6
  ) +
  scale_color_manual(values = c(
    "IPC1 up" = "#FB6A4A",
    "IPC2 up" = "#6BAED6",
    "FDR < 0.01 & |log2 FC| < 1" = "#80CDC1",
    "FDR > 0.01 & |log2 FC| > 1" = "#FDB462",
    "No significant" = "lightgrey"
  )) +
  scale_size_continuous(
    range = c(0.5, 4),
    breaks = c(1, 1.5, 2, 2.5),
    name = "|log2 FC|"
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(0.001),
    linetype = "dashed",
    color = "grey40"
  ) +
  annotate(
    "segment",
    x = 1.2,
    xend = 2.5,
    y = 14,
    yend = 14,
    arrow = arrow(length = unit(0.2, "cm")),
    color = "#FB6A4A"
  ) +
  annotate(
    "text",
    x = 1.85,
    y = 15,
    label = "IPC1",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 2,
    y = 13,
    label = paste0("n=", number_up)
  ) +
  annotate(
    "segment",
    x = -1.2,
    xend = -2.5,
    y = 14,
    yend = 14,
    arrow = arrow(length = unit(0.2, "cm")),
    color = "#6BAED6"
  ) +
  annotate(
    "text",
    x = -1.85,
    y = 15,
    label = "IPC2",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = -2,
    y = 13,
    label = paste0("n=", number_down)
  ) +
  geom_text_repel(
    data = labeled_genes,
    aes(label = Gene),
    size = 3,
    fontface = "italic"
  ) +
  theme_classic() +
  labs(x = "log2 FC", y = "-log10 FDR") +
  guides(
    color = guide_legend(
      override.aes = list(size = 3)
    )
  )

ggsave(
  file.path(figure6_dir, "Figure6F_subtype_volcano.pdf"),
  volcano_plot,
  width = 8,
  height = 7
)


# =============================================================================
# 4. Figure 6G: 亚型差异蛋白热图
# =============================================================================

load(file.path(diff_dir, "protein_imp_gene.rds"))
subtype_expression <- gene_exp
rownames(subtype_expression) <- subtype_expression[, 1]
subtype_expression <- subtype_expression[, -1]
subtype_expression <- matrix(
  as.numeric(subtype_expression),
  nrow = nrow(subtype_expression),
  ncol = ncol(subtype_expression)
)
rownames(subtype_expression) <- gene_exp[, 1]
colnames(subtype_expression) <- colnames(gene_exp)[-1]

subtype_diff <- read.table(
  file.path(nmf_dir, "c1_c2_diff.txt"),
  header = TRUE,
  stringsAsFactors = FALSE
)
subtype_groups <- read.table(
  file.path(nmf_dir, "gene_400_group_2_NMF_group.txt"),
  header = FALSE,
  stringsAsFactors = FALSE
)
significant_diff <- subtype_diff[
  abs(subtype_diff[, 1]) > 1 & subtype_diff[, 5] < 0.01,
]
diff_expression <- subtype_expression[
  rownames(significant_diff),
]

subtype_groups <- subtype_groups[
  order(subtype_groups[, 2]),
]
subtype_groups <- subtype_groups[
  c(
    sample(which(subtype_groups[, 2] == 1), 35),
    sample(36:87, 51)
  ),
]
heatmap_expression <- diff_expression[
  ,
  subtype_groups[, 1]
]
group_counts <- as.numeric(table(subtype_groups[, 2]))
group_gap <- if (length(group_counts) > 1) {
  cumsum(group_counts)[1]
} else {
  NULL
}
sample_type <- ifelse(
  grepl("^BC", subtype_groups[, 1]),
  "BC",
  ifelse(grepl("^LA", subtype_groups[, 1]), "LA", "Unknown")
)

heatmap_annotation <- data.frame(
  Group = factor(subtype_groups[[2]], levels = c(1, 2)),
  Type = factor(sample_type, levels = c("BC", "LA"))
)
rownames(heatmap_annotation) <- subtype_groups[[1]]
heatmap_annotation_colors <- list(
  Group = c("1" = "#FB9B9A", "2" = "#A9CFE4"),
  Type = c("BC" = "#8A9EB5", "LA" = "#ED8A3F")
)

subtype_heatmap <- pheatmap(
  heatmap_expression,
  annotation_col = heatmap_annotation,
  annotation_colors = heatmap_annotation_colors,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  treeheight_row = 0,
  treeheight_col = 0,
  color = colorRampPalette(
    rev(
      c(
        "#DD0514",
        rgb(255, 255, 255, maxColorValue = 255),
        "#4874B1"
      )
    )
  )(100),
  breaks = seq(-2, 2, length.out = 101),
  main = "heatmap",
  gaps_col = group_gap,
  annotation_names_col = TRUE
)

labeled_heatmap_genes <- c(
  "PCNA",
  "EEF2",
  "RAN",
  "HSP90AA1",
  "HSP90AB1",
  "PSMA1",
  "PSMA4",
  "PSMB1",
  "PSMB2",
  "PSMB3",
  "PSMB5",
  "PSMB6",
  "PSME1",
  "FDPS",
  "ACAT2",
  "LDHA",
  "PGD",
  "AHCY"
)
pdf(
  file.path(figure6_dir, "Figure6G_subtype_diff_heatmap.pdf"),
  width = 6,
  height = 8
)
add.flag(
  subtype_heatmap,
  kept.labels = labeled_heatmap_genes,
  repel.degree = 0
)
dev.off()
# 5. Figure 6H-I: 亚型功能富集和 C1/C2 条形图
# =============================================================================

run_subtype_enrichment <- function(sample_group) {
  result <- read.table(
    file.path(nmf_dir, "c1_c2_diff.txt"),
    sep = "\t",
    header = TRUE
  )
  if (sample_group == "up") {
    selected_genes <- rownames(
      result[
        result[, 1] > 1 & result[, 5] < 0.01,
      ]
    )
  } else if (sample_group == "down") {
    selected_genes <- rownames(
      result[
        result[, 1] < -1 & result[, 5] < 0.01,
      ]
    )
  } else {
    stop("sample_group must be up or down", call. = FALSE)
  }

  gene_ids <- mget(
    selected_genes,
    org.Hs.egSYMBOL2EG,
    ifnotfound = NA
  )
  gene_ids <- as.character(gene_ids)
  mapped <- cbind(selected_genes, gene_ids)
  gene_ids <- gene_ids[gene_ids != "NA"]

  enrichment_scatter <- function(
    enrichment_result,
    output_file,
    label
  ) {
    enrichment_table <- as.data.frame(enrichment_result)
    enrichment_table$function.gene.num <- unlist(
      strsplit(enrichment_table[, "BgRatio"], "/")
    )[
      seq(
        1,
        length(
          unlist(
            strsplit(enrichment_table[, "BgRatio"], "/")
          )
        ),
        2
      )
    ]
    enrichment_table$GeneRatio <- enrichment_table$Count /
      as.numeric(enrichment_table$function.gene.num)
    top_terms <- head(enrichment_table, 10)
    top_terms <- top_terms[
      order(top_terms$GeneRatio, decreasing = TRUE),
    ]
    top_terms$Description <- factor(
      top_terms$Description,
      levels = rev(top_terms$Description)
    )
    x_gap <- (
      max(top_terms$GeneRatio) -
        min(top_terms$GeneRatio)
    ) / 5

    enrichment_plot <- ggplot(
      top_terms,
      aes(x = GeneRatio, y = Description)
    ) +
      geom_point(aes(size = Count, color = pvalue)) +
      theme_bw() +
      scale_color_gradient(
        low = "red",
        high = "blue"
      ) +
      scale_x_continuous(expand = c(x_gap, x_gap)) +
      ylab(NULL) +
      theme(
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16)
      ) +
      scale_size_continuous(range = c(6, 8))

    ggsave(
      output_file,
      enrichment_plot,
      width = 8,
      height = 4
    )
  }

  go_result <- enrichGO(
    gene = gene_ids,
    keyType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH"
  )
  if (!is.null(go_result)) {
    go_table <- as.data.frame(go_result)
    go_table <- go_table[go_table$p.adjust < 2,]
    write.table(
      go_table,
      file.path(go_dir, paste0(sample_group, "_enrichGO.txt")),
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    enrichment_scatter(
      go_result,
      file.path(go_dir, paste0(sample_group, "_enrichGO.pdf")),
      "GO"
    )

    go_by_ontology <- lapply(
      c("BP", "MF", "CC"),
      function(ontology) {
        ontology_table <- go_table[
          go_table$ONTOLOGY == ontology,
        ]
        head(ontology_table, 15)
      }
    )
    go_combined <- bind_rows(go_by_ontology)
    go_combined <- go_combined[
      order(go_combined$GeneRatio, decreasing = TRUE),
    ]
    go_combined$Description <- factor(
      go_combined$Description,
      levels = rev(go_combined$Description)
    )
    go_gap <- (
      max(go_combined$GeneRatio) -
        min(go_combined$GeneRatio)
    ) / 5
    go_plot <- ggplot(
      go_combined,
      aes(x = GeneRatio, y = Description)
    ) +
      geom_point(aes(size = Count, color = p.adjust)) +
      theme_bw() +
      scale_color_gradient(low = "red", high = "blue") +
      scale_x_continuous(expand = c(go_gap, go_gap)) +
      ylab(NULL) +
      theme(
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16)
      ) +
      facet_wrap(~ONTOLOGY, ncol = 1, scales = "free")
    ggsave(
      file.path(
        go_dir,
        paste0(sample_group, "_enrichGO_bp_mf_cc.pdf")
      ),
      go_plot,
      width = 15,
      height = 6
    )
  }

  kegg_result <- enrichKEGG(
    gene = gene_ids,
    organism = "hsa",
    keyType = "kegg",
    pAdjustMethod = "BH"
  )
  if (!is.null(kegg_result)) {
    kegg_table <- as.data.frame(kegg_result)
    kegg_table <- kegg_table[kegg_table$p.adjust < 2,]
    write.table(
      kegg_table,
      file.path(go_dir, paste0(sample_group, "_enrichKEGG.txt")),
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    enrichment_scatter(
      kegg_result,
      file.path(go_dir, paste0(sample_group, "_enrichKEGG.pdf")),
      "KEGG"
    )
  }

  hallmark_terms <- msigdbr(
    species = "Homo sapiens",
    category = "H"
  ) %>%
    dplyr::select(gs_name, entrez_gene)
  hallmark_result <- enricher(
    gene_ids,
    TERM2GENE = hallmark_terms
  )
  if (!is.null(hallmark_result)) {
    hallmark_table <- as.data.frame(hallmark_result)
    hallmark_table$function.gene.num <- unlist(
      strsplit(hallmark_table[, "BgRatio"], "/")
    )[
      seq(
        1,
        length(
          unlist(
            strsplit(hallmark_table[, "BgRatio"], "/")
          )
        ),
        2
      )
    ]
    hallmark_table$GeneRatio <- hallmark_table$Count /
      as.numeric(hallmark_table$function.gene.num)
    write.table(
      hallmark_table,
      file.path(
        go_dir,
        paste0(sample_group, "_enrichmsigdb.txt")
      ),
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    enrichment_scatter(
      hallmark_result,
      file.path(
        go_dir,
        paste0(sample_group, "_enrichmsigdb.pdf")
      ),
      "Hallmark"
    )
  }
}

run_subtype_enrichment("down")
run_subtype_enrichment("up")

highlight_label_grob <- function(label, keyword, fontsize = 10) {
  gTree(
    cl = "highlight_label",
    label = label,
    keyword = keyword,
    fontsize = fontsize
  )
}

makeContent.highlight_label <- function(x) {
  text_gp <- gpar(
    fontsize = x$fontsize,
    fontfamily = "sans"
  )
  keyword_position <- regexpr(
    x$keyword,
    x$label,
    fixed = TRUE
  )
  children <- list()
  if (keyword_position[1] > 0) {
    prefix <- substr(x$label, 1, keyword_position[1] - 1)
    measure_text <- function(text) {
      grobWidth(textGrob(text, gp = text_gp))
    }
    label_width <- measure_text(x$label)
    prefix_width <- measure_text(prefix)
    keyword_width <- measure_text(x$keyword)
    padding <- unit(0.6, "pt")
    children[[1]] <- rectGrob(
      x = unit(1, "npc") -
        label_width +
        prefix_width -
        padding,
      y = unit(0.5, "npc"),
      width = keyword_width + 2 * padding,
      height = unit(0.72, "npc"),
      just = c("left", "center"),
      gp = gpar(fill = "#FFE08A", col = NA)
    )
  }
  children[[length(children) + 1]] <- textGrob(
    x$label,
    x = unit(1, "npc"),
    y = unit(0.5, "npc"),
    just = c("right", "center"),
    gp = text_gp
  )
  setChildren(x, do.call(gList, children))
}

selected_terms <- tribble(
  ~group, ~category, ~source, ~ontology, ~description, ~label, ~keyword,
  "C1", "Proliferation", "msigdb", NA, "HALLMARK_MYC_TARGETS_V1", "MYC targets", "MYC",
  "C1", "Proliferation", "msigdb", NA, "HALLMARK_E2F_TARGETS", "E2F targets", "E2F",
  "C1", "Proliferation", "msigdb", NA, "HALLMARK_G2M_CHECKPOINT", "G2M checkpoint", "G2M",
  "C1", "Proliferation", "msigdb", NA, "HALLMARK_MTORC1_SIGNALING", "mTORC1 signaling", "mTORC1",
  "C1", "Proliferation", "msigdb", NA, "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "PI3K-AKT-mTOR signaling", "PI3K-AKT-mTOR",
  "C1", "Metabolism", "kegg", NA, "Glycolysis / Gluconeogenesis", "Glycolysis / gluconeogenesis", "Glycolysis",
  "C1", "Metabolism", "kegg", NA, "Pentose phosphate pathway", "Pentose phosphate pathway", "Pentose phosphate",
  "C1", "Metabolism", "kegg", NA, "Carbon metabolism", "Carbon metabolism", "Carbon",
  "C1", "Metabolism", "kegg", NA, "Biosynthesis of amino acids", "Biosynthesis of amino acids", "amino acids",
  "C1", "Proteostasis", "kegg", NA, "Proteasome", "Proteasome", "Proteasome",
  "C1", "Proteostasis", "go", "BP", "RNA splicing", "RNA splicing", "RNA splicing",
  "C1", "Proteostasis", "go", "BP", "cytoplasmic translation", "Cytoplasmic translation", "Cytoplasmic translation",
  "C1", "Proteostasis", "go", "BP", "chaperone-mediated protein folding", "Chaperone-mediated protein folding", "Chaperone-mediated",
  "C1", "Junction", "kegg", NA, "Tight junction", "Tight junction", "Tight junction",
  "C1", "Junction", "kegg", NA, "Adherens junction", "Adherens junction", "Adherens junction",
  "C1", "Junction", "go", "CC", "focal adhesion", "Focal adhesion", "Focal adhesion",
  "C1", "Junction", "go", "MF", "cadherin binding", "Cadherin binding", "Cadherin binding",
  "C1", "Junction", "kegg", NA, "IgSF CAM signaling", "IgSF CAM signaling", "IgSF CAM",
  "C1", "Junction", "kegg", NA, "Motor proteins", "Motor proteins", "Motor proteins",
  "C2", "Neural", "kegg", NA, "Neuroactive ligand-receptor interaction", "Neuroactive ligand-receptor interaction", "Neuroactive ligand-receptor",
  "C2", "Neural", "go", "BP", "neuropeptide signaling pathway", "Neuropeptide signaling pathway", "Neuropeptide",
  "C2", "Neural", "go", "MF", "G protein-coupled receptor binding", "G protein-coupled receptor binding", "G protein-coupled receptor",
  "C2", "Neural", "go", "CC", "synaptic membrane", "Synaptic membrane", "Synaptic",
  "C2", "Neural", "go", "CC", "distal axon", "Distal axon", "axon",
  "C2", "Immune", "go", "BP", "immunoglobulin mediated immune response", "Immunoglobulin mediated immune response", "Immunoglobulin",
  "C2", "Immune", "go", "BP", "B cell mediated immunity", "B cell mediated immunity", "B cell",
  "C2", "Immune", "go", "MF", "antigen binding", "Antigen binding", "Antigen",
  "C2", "Immune", "go", "BP", "antimicrobial humoral response", "Antimicrobial humoral response", "Antimicrobial",
  "C2", "Granule/ECM", "go", "CC", "secretory granule lumen", "Secretory granule lumen", "Secretory granule",
  "C2", "Granule/ECM", "go", "MF", "endopeptidase inhibitor activity", "Endopeptidase inhibitor activity", "Endopeptidase inhibitor",
  "C2", "Granule/ECM", "go", "MF", "extracellular matrix structural constituent", "Extracellular matrix structural constituent", "Extracellular matrix"
)

read_enrichment_file <- function(source_type, group_name) {
  file_path <- switch(
    source_type,
    msigdb = file.path(
      go_dir,
      paste0(group_name, "_enrichmsigdb.txt")
    ),
    kegg = file.path(
      go_dir,
      paste0(group_name, "_enrichKEGG.txt")
    ),
    go = file.path(
      go_dir,
      paste0(group_name, "_enrichGO.txt")
    )
  )
  read.delim(file_path, check.names = FALSE)
}

enrichment_rows <- list()
for (term_index in seq_len(nrow(selected_terms))) {
  selected_term <- selected_terms[term_index,]
  source_group <- ifelse(
    selected_term$group == "C1",
    "down",
    "up"
  )
  source_table <- read_enrichment_file(
    selected_term$source,
    source_group
  )
  term_hit <- source_table$Description ==
    selected_term$description
  if (
    !is.na(selected_term$ontology) &&
      "ONTOLOGY" %in% names(source_table)
  ) {
    term_hit <- term_hit &
      source_table$ONTOLOGY == selected_term$ontology
  }
  if (sum(term_hit) != 1) {
    stop(
      "无法唯一定位富集条目：",
      selected_term$description,
      call. = FALSE
    )
  }
  enrichment_rows[[term_index]] <- selected_term %>%
    mutate(
      Count = source_table$Count[term_hit],
      padj = source_table$p.adjust[term_hit],
      neglog10P = -log10(source_table$p.adjust[term_hit])
    )
}

enrichment_data <- bind_rows(enrichment_rows)
category_order <- c(
  "Proliferation",
  "Metabolism",
  "Proteostasis",
  "Junction",
  "Neural",
  "Immune",
  "Granule/ECM"
)
enrichment_data <- enrichment_data %>%
  mutate(
    category = factor(category, levels = category_order)
  ) %>%
  arrange(group, category, desc(neglog10P)) %>%
  mutate(label = factor(label, levels = rev(unique(label))))

category_colors <- c(
  "Proliferation" = "#E08214",
  "Metabolism" = "#8073AC",
  "Proteostasis" = "#1B9E77",
  "Junction" = "#737373",
  "Neural" = "#E08214",
  "Immune" = "#8073AC",
  "Granule/ECM" = "#1B9E77"
)
source_colors <- c(
  "msigdb" = "#66C2A5",
  "kegg" = "#FC8D62",
  "go" = "#8DA0CB"
)
source_labels <- c(
  "msigdb" = "Hallmark",
  "kegg" = "KEGG",
  "go" = "GO"
)

make_enrichment_panel <- function(
  panel_data,
  panel_title,
  fill_low,
  fill_high,
  x_maximum,
  show_source_legend = TRUE
) {
  panel_data <- panel_data %>%
    mutate(label = droplevels(label))
  panel_data$y_position <- as.integer(panel_data$label)
  category_x <- -x_maximum * 0.045
  source_x <- -x_maximum * 0.105

  enrichment_plot <- ggplot(panel_data) +
    geom_tile(
      aes(
        x = category_x,
        y = y_position,
        fill = category
      ),
      width = x_maximum * 0.028,
      height = 0.72,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = category_colors,
      guide = "none"
    ) +
    new_scale_fill() +
    geom_tile(
      aes(
        x = source_x,
        y = y_position,
        fill = source
      ),
      width = x_maximum * 0.028,
      height = 0.72,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = source_colors,
      labels = source_labels,
      name = NULL,
      limits = names(source_colors),
      drop = FALSE,
      guide = if (show_source_legend) "legend" else "none"
    ) +
    new_scale_fill() +
    geom_col(
      aes(
        x = Count,
        y = label,
        fill = neglog10P
      ),
      width = 0.62
    ) +
    scale_fill_gradient(
      low = fill_low,
      high = fill_high,
      name = expression(-log[10](italic(P) ~ value))
    ) +
    scale_x_continuous(
      position = "top",
      expand = expansion(mult = c(0, 0.03)),
      limits = c(-x_maximum * 0.14, x_maximum)
    ) +
    scale_y_discrete(expand = expansion(add = 0.6)) +
    coord_cartesian(clip = "off") +
    labs(
      title = panel_title,
      x = "Count",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y = element_text(
        size = 10,
        colour = NA,
        hjust = 1
      ),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 9),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black"),
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      ),
      plot.margin = margin(5, 5, 5, 5),
      legend.position = "right",
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.8, "cm"),
      legend.title = element_text(size = 9)
    )

  label_right <- -x_maximum * 0.155
  label_left <- -x_maximum * 1.3
  for (row_index in seq_len(nrow(panel_data))) {
    enrichment_plot <- enrichment_plot +
      annotation_custom(
        grob = highlight_label_grob(
          as.character(panel_data$label[row_index]),
          panel_data$keyword[row_index],
          fontsize = 10
        ),
        xmin = label_left,
        xmax = label_right,
        ymin = panel_data$y_position[row_index] - 0.45,
        ymax = panel_data$y_position[row_index] + 0.45
      )
  }

  enrichment_plot
}

c1_data <- enrichment_data %>%
  filter(group == "C1")
c2_data <- enrichment_data %>%
  filter(group == "C2")
c1_panel <- make_enrichment_panel(
  c1_data,
  "C1-enriched",
  "#FEDEC6",
  "#FE954B",
  55
)
c2_panel <- make_enrichment_panel(
  c2_data,
  "C2-enriched",
  "#C7DAEE",
  "#3182BD",
  max(c2_data$Count) * 1.25,
  show_source_legend = FALSE
)

combined_enrichment_plot <- c1_panel / c2_panel +
  plot_layout(
    heights = c(nrow(c1_data), nrow(c2_data)),
    guides = "collect"
  )

ggsave(
  file.path(figure6_dir, "C1C2_enrichment_barplot.pdf"),
  combined_enrichment_plot,
  width = 8,
  height = 0.32 * nrow(enrichment_data) + 1.6
)


# =============================================================================
# 6. Figure 6J: ssGSEA 热图和山脊图
# =============================================================================

plot_score_heatmap_ridge <- function(
  group_file,
  score_file,
  output_name,
  output_directory
) {
  subtype_group <- read.table(group_file)
  subtype_1 <- subtype_group[
    subtype_group[, 2] == 1,
    1
  ]
  subtype_2 <- subtype_group[
    subtype_group[, 2] == 2,
    1
  ]
  rownames(subtype_group) <- subtype_group[, 1]

  load(score_file)
  score_matrix <- ssGSEA[, subtype_group[, 1]]

  test_results <- data.frame()
  for (score_index in seq_len(nrow(score_matrix))) {
    wilcox_result <- wilcox.test(
      score_matrix[score_index, subtype_1],
      score_matrix[score_index, subtype_2]
    )
    test_results <- rbind(
      test_results,
      data.frame(
        score = rownames(score_matrix)[score_index],
        p_value = wilcox_result$p.value
      )
    )
  }
  selected_scores <- test_results$score[
    test_results$p_value < 0.05
  ]
  write.table(
    test_results,
    file.path(
      output_directory,
      paste0(output_name, "_ridge.txt")
    ),
    quote = FALSE,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE
  )

  selected_matrix <- score_matrix[
    selected_scores,
    c(subtype_1, subtype_2)
  ]
  subtype_labels <- paste(
    "c",
    subtype_group[c(subtype_1, subtype_2), 2],
    sep = ""
  )
  heatmap_annotation <- data.frame(
    Cancer_type = as.factor(
      substr(colnames(selected_matrix), 1, 2)
    ),
    Cluster = subtype_labels
  )
  rownames(heatmap_annotation) <- colnames(selected_matrix)
  annotation_colors <- list(
    Cluster = c(
      c1 = "#FB9B9A",
      c2 = "#A9CFE4"
    ),
    Cancer_type = c(
      BC = "#8A9EB5",
      LA = "#ED8A3F"
    )
  )

  scaled_scores <- t(apply(selected_matrix, 1, scale))
  rownames(scaled_scores) <- colnames(selected_matrix)
  scaled_scores <- t(scaled_scores)
  scaled_scores[scaled_scores < -2] <- -2
  scaled_scores[scaled_scores > 2] <- 2

  heatmap_colors <- colorRampPalette(
    c(
      "#4874B1",
      rgb(255, 255, 255, maxColorValue = 255),
      "#DD0514"
    )
  )(100)
  score_heatmap <- pheatmap(
    scaled_scores,
    fontsize = 6,
    scale = "none",
    treeheight_row = 0,
    cellwidth = 1.2,
    cellheight = 7,
    color = heatmap_colors,
    annotation_col = heatmap_annotation,
    annotation_colors = annotation_colors,
    clustering_method = "ward.D2",
    border_color = NA,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE
  )
  ggsave(
    file.path(
      output_directory,
      paste0(output_name, ".heatmap.pdf")
    ),
    score_heatmap$gtable,
    width = 6,
    height = max(4, 0.15 * length(selected_scores) + 2)
  )

  row_order <- score_heatmap$tree_row$order
  ordered_score_names <- rownames(scaled_scores)[row_order]
  score_long <- melt(
    selected_matrix,
    varnames = c("geneset", "sample")
  )
  score_long$subtype <- "NA"
  score_long$subtype[
    score_long$sample %in% subtype_1
  ] <- "c1"
  score_long$subtype[
    score_long$sample %in% subtype_2
  ] <- "c2"
  score_long$subtype <- factor(
    score_long$subtype,
    levels = c("c1", "c2")
  )

  ridge_plot <- ggplot(
    score_long,
    aes(
      x = value,
      y = geneset,
      fill = subtype
    )
  ) +
    geom_density_ridges() +
    scale_fill_manual(
      values = c("#FB9B9A", "#A9CFE4", "#9B8281")
    ) +
    theme_ridges() +
    xlab(output_name) +
    scale_y_discrete(limits = rev(ordered_score_names))

  ggsave(
    file.path(
      output_directory,
      paste0(output_name, "_ridge.pdf")
    ),
    ridge_plot,
    width = 6,
    height = max(4, 0.2 * length(selected_scores) + 2)
  )
}

plot_score_heatmap_ridge(
  file.path(nmf_dir, "gene_400_group_2_NMF_group.txt"),
  file.path(single_cell_dir, "ssGSEA_50.rds"),
  "single_cell",
  figure6_dir
)
plot_score_heatmap_ridge(
  file.path(nmf_dir, "gene_400_group_2_NMF_group.txt"),
  file.path(subtype_dir, "ssGSEA_29_Fges.rds"),
  "Fges",
  figure6_dir
)
plot_score_heatmap_ridge(
  file.path(nmf_dir, "gene_400_group_2_NMF_group.txt"),
  file.path(subtype_dir, "ssGSEA_28_immune.rds"),
  "immune_28",
  figure6_dir
)


# =============================================================================
