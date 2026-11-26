#!/usr/bin/env Rscript


project_root <- normalizePath(
  Sys.getenv("PROTEO_PROJECT_ROOT", unset = getwd()),
  mustWork = FALSE
)
exp_dir <- file.path(project_root, "exp")
diff_dir <- file.path(project_root, "diff")
figure1_dir <- file.path(project_root, "Figure1")
subtype_dir <- file.path(project_root, "subtype")
secretome_dir <- Sys.getenv(
  "PROTEO_SECRETOME_DIR",
  unset = file.path(figure1_dir, "SECRETOME")
)
add_flag_script <- Sys.getenv(
  "PROTEO_ADD_FLAG_SCRIPT",
  unset = file.path(project_root, "scripts", "add_flag.r")
)
extra_r_lib <- Sys.getenv("PROTEO_EXTRA_R_LIB", unset = "")


# =============================================================================
# 1. Figure 1b: 样本临床信息热图
# =============================================================================


library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(openxlsx)
#
setwd(figure1_dir)
data <- read.xlsx(file.path(exp_dir, "patients_clinical_information.xlsx"), startRow = 2)
data[is.na(data)] <- "N"
NA_COL <- "#EDEDED"
exp_raw=read.table(file.path(exp_dir, "protein_summary.txt"),sep="\t",quote="",header=T,check.names=F)
length(intersect(colnames(exp_raw),data[,1]))#107,olink有139
data=data[is.element(data[,1],colnames(exp_raw)),]
filter=read.table(file.path(diff_dir, "filter_sample.txt"),sep="\t",quote="",header=F,check.names=F)
data=data[!is.element(data[,1],filter[,1]),]

data2 <- data %>%
  mutate(
    # Cohort    = factor(Cohort, levels = c("Discovery", "Replication")),
    Group = factor(Group, levels = c("pwLM", "Control")),
    Response  = factor(Response, levels = c("PR","SD","PD","N"))
  ) %>%
  arrange( Group, Source, Cancer)
sample_order <- data2$Patient
data2 <- data2 %>%
  mutate(
    iPFS_time = as.numeric(ifelse(iPFS_time == "N", NA, iPFS_time)),
    OS_time   = as.numeric(ifelse(OS_time   == "N", NA, OS_time))
  )
clinical_vars <- c(
  "Group", "Source", 
  "Age", "Sex", "Cancer", "Subtype",
  "iPFS_status", "iPFS_time",
  "OS_status", "OS_time",
  "Response", "MRI", "Cytology"
)

anno_df <- data2 %>%
  select(all_of(clinical_vars))
col_fun_age <- colorRamp2(
  c(min(data2$Age, na.rm = TRUE), max(data2$Age, na.rm = TRUE)),
  c("#4575B4", "#D73027")
)

col_fun_iPFS <- colorRamp2(
  c(0, max(data2$iPFS_time, na.rm = TRUE)),
  c("white", "#762A83")
)

col_fun_OS <- colorRamp2(
  c(0, max(data2$OS_time, na.rm = TRUE)),
  c("white", "#A50F15")
)
col_list <- list(
  
  Group = c(
    pwLM = "#FF9A91",
    Control   = "#CDE2F2"
  ),
  
  Source = c(
    Fuzhou           = "#EAD1DC",  # 浅灰紫
    Heilongjiang     = "#D9EAD3",  # 浅灰绿
    Henan            = "#CFE2F3",  # 浅天蓝
    Nanjing_Brain    = "#FFF2CC",  # 浅香槟/淡黄
    Nanjing_province = "#F4CCCC"   # 浅粉红
  ),
  
  Sex = c(
    Male   = "#8FAADC",    # 莫兰迪蓝
    Female = "#EAC7CC"     
  ),

  Response = c(
    PR = "#00A087",
    SD = "#F39B7F",
    PD = "#DC0000",
    N  = NA_COL
  ),
  
  MRI = c(
    Positive = "#71cf8a",
    Negative = "#BFC9D4",
    N = NA_COL
  ),
  
  Cytology = c(
    Positive = "#7c90c9",
    Negative = "#EAC7CC",
    N = NA_COL
  ),
  
  Cancer = c(
    BC = "#8ECFC9",  # 浅湖绿
    LA = "#FFBE7A",  # 浅暖橙
    EC = "#FA7F6F",  # 浅珊瑚红
    MM = "#82B0D2",  # 浅灰蓝
    N  = NA_COL
  ),

  Subtype =c(
    Adenocarcinoma = "#A3D1D1",  # 浅青瓷色 (对应腺癌的纯净感)
    TNBC           = "#FADBD8",  # 极浅珊瑚粉 (替代刺眼的纯红)
    `HER2+`        = "#D4EFDF",  # 浅薄荷绿 (清新、专业)
    Other          = "#B2B2B2",  # 浅雾灰色 (让次要类别退到背景)
    `Luminal A/B`  = "#FDEBD0",  # 浅杏色 (柔和、温暖)
    N              = NA_COL      # 保持一致的缺失值色
  ),
  
  iPFS_status = c(
    "0" = "#C6D8E1",       # 亮冰蓝 (比 A9BDCE 更亮，更清透)
    "1" = "#D9B4B5",       # 浅烟粉 (比 C59395 更浅，温润感)
    "N" = NA_COL
  ),
  
  OS_status = c(
    "0" = "#C6D8E1",
    "1" = "#D9B4B5",
    "N" = NA_COL
  )
)

rownames(anno_df) <- data2$Patient

# ha <- HeatmapAnnotation(
#   df = anno_df,
#   col = c(
#     list(
#       Age       = col_fun_age,
#       iPFS_time = col_fun_iPFS,
#       OS_time   = col_fun_OS
#     ),
#     col_list
#   ),
#   na_col = NA_COL,
#   annotation_name_side = "left",
#   gp = gpar(col = "gray", lwd = 0.5)  
# )

ha <- HeatmapAnnotation(
  df = anno_df,
  col = c(
    list(Age = col_fun_age, iPFS_time = col_fun_iPFS, OS_time = col_fun_OS),
    col_list
  ),
  na_col = NA_COL,
  # 关键：lwd = 1 且 col = "white" 能让浅色块之间产生“呼吸感”
  gp = gpar(col = "white", lwd = 0.5), 
  annotation_name_side = "left"
)

#cohort_split <- anno_df$Cohort   
pdf("sample_info.pdf", width = 12)
p <- Heatmap(
  matrix(0, nrow = 1, ncol = nrow(anno_df)),
  
  col = c("0" = "white"),   # 关键1：把热图变成白色（不可见）
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  
  show_heatmap_legend = FALSE,  # 关键2：去掉 matrix_1 图例
  
  column_gap = unit(3, "mm"),
  
  top_annotation = ha,
  height = unit(3, "mm"),
  rect_gp = gpar(col = NA)
)
print(p)
dev.off()


# =============================================================================
# 2. Figure 1c: 样本 PCA 与 UMAP 降维图
# =============================================================================


info=read.table(file.path(exp_dir, "patient_source.txt"),sep="\t",quote="",header=F,check.names=F)
rownames(info)=info[,1]
nanjing=rownames(info)[which(is.element(info[,2],c("nanjing")))]
henan=rownames(info)[which(is.element(info[,2],c("henan")))]
fujian=rownames(info)[which(is.element(info[,2],c("fujian")))]
heilongjiang=rownames(info)[which(is.element(info[,2],c("heilongjiang")))]
con=c(fujian,heilongjiang)
filter=read.table(file.path(diff_dir, "filter_sample.txt"),sep="\t",quote="",header=F,check.names=F)
con=setdiff(con,filter[,1])
LA=info[grep("LA",info[,1]),1]
LA=setdiff(LA,filter[,1])
BC=info[grep("BC",info[,1]),1]
BC=setdiff(BC,filter[,1])
nanjing=setdiff(nanjing,filter[,1])
henan=setdiff(henan,filter[,1])
fujian=setdiff(fujian,filter[,1])
heilongjiang=setdiff(heilongjiang,filter[,1])
setwd(diff_dir)
load(file.path(diff_dir, "protein_imp.rds"))
exp_imp1=exp_imp1[,c(con,nanjing,henan)]

type=info[colnames(exp_imp1),2]
class <- ifelse(
  type %in% c("nanjing", "henan"), "case",
  ifelse(type %in% c("heilongjiang", "fujian"), "control", type)
)

getPCs <- function(exp){
	pcs <- prcomp(t(exp))
	pVar <- pcs$sdev^2/sum(pcs$sdev^2)
	pca.mat <- data.frame(pcs$x)
	pca.mat$type <- type
	res <- list(pcs = pca.mat, pVar=pVar)
	return(res)
}
library("ggplot2")


library(umap)
set.seed(1234) 
#pca_50 <- prcomp(t(exp_scaled), rank.=50)$x
pcs <- getPCs(t(scale(t(exp_imp1))))
umap_re=umap(pcs$pcs[,1:50])
#umap_re=umap(scale(t(exp_imp1)))
use=as.data.frame(cbind(umap_re$layout,type))

g3=ggplot(data=use, aes(x=as.numeric(V1), y=as.numeric(V2), color=type)) +
geom_point(size = 3) +
xlab("umap_1") +
ylab("umap_2") +
scale_colour_manual(values = c("#5AB4AC","#AF8DC3","#FCA689","#537E44")) +
coord_fixed() + theme_bw() +
theme(panel.grid = element_blank())
pdf("sample_umap_filter.pdf",width=4,height=4)
g3
dev.off()

g4=ggplot(data=use, aes(x=as.numeric(V1), y=as.numeric(V2), color=class)) +
geom_point(size = 3) +
xlab("umap_1") +
ylab("umap_2") +
scale_colour_manual(values = c("#F3A485","#527942")) +
coord_fixed() + theme_bw() +
theme(panel.grid = element_blank())
pdf("sample_umap_case_filter.pdf",width=4,height=4)
g4
dev.off()


# =============================================================================
# 3. Figure 1d-e: 蛋白丰度排序曲线
# =============================================================================


info=read.table(file.path(exp_dir, "patient_source.txt"),sep="\t",quote="",header=F,check.names=F)
rownames(info)=info[,1]
nanjing=rownames(info)[which(is.element(info[,2],c("nanjing")))]
henan=rownames(info)[which(is.element(info[,2],c("henan")))]
fujian=rownames(info)[which(is.element(info[,2],c("fujian")))]
heilongjiang=rownames(info)[which(is.element(info[,2],c("heilongjiang")))]
con=c(fujian,heilongjiang)
filter=read.table(file.path(diff_dir, "filter_sample.txt"),sep="\t",quote="",header=F,check.names=F)
con=setdiff(con,filter[,1])
LA=info[grep("LA",info[,1]),1]
LA=setdiff(LA,filter[,1])
BC=info[grep("BC",info[,1]),1]
BC=setdiff(BC,filter[,1])
nanjing=setdiff(nanjing,filter[,1])
henan=setdiff(henan,filter[,1])
fujian=setdiff(fujian,filter[,1])
heilongjiang=setdiff(heilongjiang,filter[,1])
setwd(diff_dir)
#load(file.path(diff_dir, "protein_imp.rds"))
load(file.path(diff_dir, "protein_imp_gene.rds"))
exp_imp1=gene_exp
rownames(exp_imp1)=exp_imp1[,1]
exp_imp1=exp_imp1[,-1]
exp_imp1=matrix(as.numeric(exp_imp1), nrow = nrow(exp_imp1), ncol = ncol(exp_imp1))
rownames(exp_imp1)=gene_exp[,1]
colnames(exp_imp1)=colnames(gene_exp)[-1]

protein_abundance=rowMeans(exp_imp1[,c(nanjing,henan)])
sorted_protein <- sort(protein_abundance, decreasing = TRUE)
x <- 1:length(protein_abundance)
original=2^sorted_protein
stat=sum(original[1:10])/sum(original)*100 #0.3158003
# 5. 绘制散点图（点点图）
pdf(file.path(diff_dir, "Protein_Abundance_case.pdf"),width=4,height=4)
plot(
  x = x,               # 横坐标：1到100
  y = sorted_protein,  # 纵坐标：排序后的表达值
  pch = 19,            # 点的形状（实心圆）
  
  cex = 0.2,           # 点的大小
  col = c(rep("red", 10), rep("#3498db", length(original)-10)),
  xlab = "Abundance rank",  # 横轴标签
  ylab = "Protein Abundance",                     # 纵轴标签
  main = "cancer",  # 标题
  las = 1,             # 坐标轴刻度文字方向
  lwd = 2,              # 点的边框粗细
  cex.axis = 0.5,   # 坐标轴刻度数字变小
  cex.lab = 0.7
)

text(3000, max(sorted_protein)*0.9, 
     paste0("Top10: ", round(stat,2), "%"), 
     col="red", cex=0.6)

# ===================== 4. 右侧标注前10个基因名 =====================
for (i in 1:10) {
  text(
    x = 1000,          # 放在图右侧
    y = 11-i*0.8,
    labels = names(sorted_protein)[i],
    col = "red",      # 基因名也红色
    cex = 0.6, 
  )
}
dev.off()

########control
protein_abundance=rowMeans(exp_imp1[,con])
sorted_protein <- sort(protein_abundance, decreasing = TRUE)
x <- 1:length(protein_abundance)
original=2^sorted_protein
stat=sum(original[1:10])/sum(original)*100 #0.4582172
# 5. 绘制散点图（点点图）
pdf(file.path(diff_dir, "Protein_Abundance_con.pdf"),width=4,height=4)
plot(
  x = x,               # 横坐标：1到100
  y = sorted_protein,  # 纵坐标：排序后的表达值
  pch = 19,            # 点的形状（实心圆）
  
  cex = 0.2,           # 点的大小
  col = c(rep("red", 10), rep("#3498db", length(original)-10)),
  xlab = "Abundance rank",  # 横轴标签
  ylab = "Protein Abundance",                     # 纵轴标签
  main = "control",  # 标题
  las = 1,             # 坐标轴刻度文字方向
  lwd = 2,              # 点的边框粗细
  cex.axis = 0.5,   # 坐标轴刻度数字变小
  cex.lab = 0.7
)

text(3000, max(sorted_protein)*0.9, 
     paste0("Top10: ", round(stat,2), "%"), 
     col="red", cex=0.6)

# ===================== 4. 右侧标注前10个基因名 =====================
for (i in 1:10) {
  text(
    x = 1000,          # 放在图右侧
    y = 11-i*0.8,
    labels = names(sorted_protein)[i],
    col = "red",      # 基因名也红色
    cex = 0.6, 
  )
}
dev.off()


##########
########LA
protein_abundance=rowMeans(exp_imp1[,LA])
sorted_protein <- sort(protein_abundance, decreasing = TRUE)
x <- 1:length(protein_abundance)
original=2^sorted_protein
stat=sum(original[1:10])/sum(original)*100 #0.4582172
# 5. 绘制散点图（点点图）
pdf(file.path(diff_dir, "Protein_Abundance_LA.pdf"),width=4,height=4)
plot(
  x = x,               # 横坐标：1到100
  y = sorted_protein,  # 纵坐标：排序后的表达值
  pch = 19,            # 点的形状（实心圆）
  
  cex = 0.2,           # 点的大小
  col = c(rep("red", 10), rep("#3498db", length(original)-10)),
  xlab = "Abundance rank",  # 横轴标签
  ylab = "Protein Abundance",                     # 纵轴标签
  main = "LA",  # 标题
  las = 1,             # 坐标轴刻度文字方向
  lwd = 2,              # 点的边框粗细
  cex.axis = 0.5,   # 坐标轴刻度数字变小
  cex.lab = 0.7
)

text(3000, max(sorted_protein)*0.9, 
     paste0("Top10: ", round(stat,2), "%"), 
     col="red", cex=0.6)

# ===================== 4. 右侧标注前10个基因名 =====================
for (i in 1:10) {
  text(
    x = 1000,          # 放在图右侧
    y = 11-i*0.8,
    labels = names(sorted_protein)[i],
    col = "red",      # 基因名也红色
    cex = 0.6, 
  )
}
dev.off()


##########
########BC
protein_abundance=rowMeans(exp_imp1[,BC])
sorted_protein <- sort(protein_abundance, decreasing = TRUE)
x <- 1:length(protein_abundance)
original=2^sorted_protein
stat=sum(original[1:10])/sum(original)*100 #0.4582172
# 5. 绘制散点图（点点图）
pdf(file.path(diff_dir, "Protein_Abundance_BC.pdf"),width=4,height=4)
plot(
  x = x,               # 横坐标：1到100
  y = sorted_protein,  # 纵坐标：排序后的表达值
  pch = 19,            # 点的形状（实心圆）
  
  cex = 0.2,           # 点的大小
  col = c(rep("red", 10), rep("#3498db", length(original)-10)),
  xlab = "Abundance rank",  # 横轴标签
  ylab = "Protein Abundance",                     # 纵轴标签
  main = "Breast cancer",  # 标题
  las = 1,             # 坐标轴刻度文字方向
  lwd = 2,              # 点的边框粗细
  cex.axis = 0.5,   # 坐标轴刻度数字变小
  cex.lab = 0.7
)

text(3000, max(sorted_protein)*0.9, 
     paste0("Top10: ", round(stat,2), "%"), 
     col="red", cex=0.6)

# ===================== 4. 右侧标注前10个基因名 =====================
for (i in 1:10) {
  text(
    x = 1000,          # 放在图右侧
    y = 11-i*0.8,
    labels = names(sorted_protein)[i],
    col = "red",      # 基因名也红色
    cex = 0.6, 
  )
}
dev.off()


# =============================================================================
# 4. Figure 1f: 蛋白亚细胞定位分布
# =============================================================================
#

# ==========================================
# 0. 载入核心库
# ==========================================
library(tidyverse)
library(ggplot2)
library(forcats)
library(patchwork)



# ==========================================
# 5. 可视化绘制
# ==========================================

plot_df <- final_counts %>%
  mutate(
    Final_Class = fct_reorder(Final_Class, Percentage),
    id = row_number()
  )

# count 分段（250一阶）
plot_df <- plot_df %>%
  mutate(
    Count_group = cut(
      n,
      breaks = c(0, 250, 500, 750, 1000, Inf),
      labels = c("0–250", "251–500", "501–750", "751–1000", ">1000"),
      include.lowest = TRUE,
      right = TRUE
    )
  )

# 径向线数据
spoke_df <- plot_df %>%
  transmute(
    x = Final_Class,
    xend = Final_Class,
    y = 0,
    yend = max(plot_df$Percentage)
  )

number_of_bar <- nrow(plot_df)
angle <- 90 - 360 * (plot_df$id - 0.5) / number_of_bar

# ==========================================
# 配色
# ==========================================

fill_colors <- c(
  "0–250" = "#7c7e81",
  "251–500" = "#9b8d86",
  "501–750" = "#b99d8b",
  "751–1000" = "#d9ad8f",
  ">1000" = "#f2a593"
)

# ==========================================
# A. 环形长条图
# ==========================================

p_circular_bar <- ggplot(
  plot_df,
  aes(x = Final_Class, y = Percentage, fill = Count_group)
) +

  geom_segment(
    data = spoke_df,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    color = "grey80",
    linetype = "dashed",
    linewidth = 0.35
  ) +

  geom_col(width = 0.85, color = NA) +

  coord_polar(start = 0, direction = -1) +

  scale_y_continuous(
    limits = c(-0.02, max(plot_df$Percentage)),
    expand = c(0, 0),
    breaks = seq(0, max(plot_df$Percentage), length.out = 4)
  ) +

  scale_fill_manual(values = fill_colors) +

  theme_void() +

  theme(
    axis.text.x = element_text(
      size = 9,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_blank(),

    panel.grid.major.y = element_line(
      color = "grey88",
      linewidth = 0.35
    ),
    panel.grid.major.x = element_blank(),

    legend.position = "none",

    plot.margin = margin(0, 0, 0, 0)
  )

# ==========================================
# A 图图例条带
# ==========================================

legend_bar <- data.frame(
  xmin = seq(0, 1000, by = 250),
  xmax = seq(250, 1250, by = 250),
  group = factor(1:5)
)

p_legend <- ggplot(legend_bar) +

  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = 0,
      ymax = 1,
      fill = group
    ),
    color = NA,
    linewidth = 0
  ) +

  scale_fill_manual(
    values = c(
      "#7c7e81",
      "#9b8d86",
      "#b99d8b",
      "#d9ad8f",
      "#f2a593"
    )
  ) +

  geom_text(
    data = data.frame(x = seq(0, 1250, 250)),
    aes(x = x, y = -0.25, label = x),
    inherit.aes = FALSE,
    size = 4
  ) +

  annotate(
    "text",
    x = 625,
    y = 1.6,
    label = "Numbers of protein",
    fontface = "bold",
    size = 5
  ) +

  coord_cartesian(
    xlim = c(-150, 1400),
    ylim = c(-0.4, 1.6),
    clip = "off"
  ) +

  theme_void() +

  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )

# 拼接 A 图
final_A <- p_circular_bar / p_legend +
  plot_layout(heights = c(10, 0.6))

ggsave(
  "Final_9class_circular.pdf",
  final_A,
  width = 4,
  height = 4
)


# =============================================================================
# 5. Figure 1g: 分泌蛋白分布
# =============================================================================

# =============================================================================

library(tidyverse)
library(ggplot2)
library(forcats)
library(patchwork)

# ==========================================
# 1. 数据读取
# ==========================================

load(file.path(diff_dir, "protein_imp.rds"))

exp_raw <- read.table(file.path(exp_dir, "protein_summary.txt"),
                      sep="\t", quote="", header=TRUE, check.names=FALSE)
rownames(exp_raw) <- exp_raw[,1]

df <- as.data.frame(exp_imp1)
df$gene <- toupper(trimws(as.character(exp_raw[rownames(exp_imp1), 3])))

secretome_files <- list.files(secretome_dir, pattern = "sa_.*\\.tsv$", full.names = TRUE)

secretome_list <- lapply(secretome_files, function(f) {
  tmp <- read.table(f, sep="\t", header=TRUE, quote="", check.names=FALSE)
  tmp$gene <- toupper(trimws(tmp$Gene))
  
  category <- gsub("sa_|_secreted\\.tsv", "", basename(f))
  tmp$Secretome_Class <- category
  
  tmp[, c("gene", "Secretome_Class")]
})

secretome_df <- bind_rows(secretome_list) %>%
  distinct(gene, Secretome_Class)

use=df$gene[which(df$gene!="--")]
final_all=as.data.frame(use)
colnames(final_all)[1]="gene"
final_all["Term"]="term"
final_all["Final_Class"]="Non-secreted"
final_all[is.element(final_all[,1],secretome_df[,1]),3]="Secreted"

left_df <- final_all %>%
  distinct(gene, Final_Class) %>%
  count(Final_Class) %>%
  mutate(Percentage = n / sum(n)) %>%
  arrange(desc(Final_Class))
  
left_df$ymax <- cumsum(left_df$Percentage)
left_df$ymin <- c(0, head(left_df$ymax, -1))


color_left <- c(
  "Secreted" = "#90cbfb",
  "Non-secreted" = "#b5b6b8"
)

# 右图颜色（自动生成 or 手动指定）
color_right <- c(
  "blood" = "#FF9A91",                     # muted red
  "brain" = "#5B7C99",                     # desaturated blue
  "digestive_system" = "#6A9F6B",          # muted green
  "extracellular_matrix" = "#d685b8",      # soft purple
  "female_reproductive_system" = "#dbb892",# muted orange
  "location_immunoglobulin" = "#f9cb1e",   # dull yellow
  "location_intracellular" = "#60bba2",    # brown-grey
  "male_reproductive_system" = "#f3c3db",  # soft pink-purple
  "other_tissues" = "#8f9497",             # neutral grey
  "unknown_location" = "#6FA3A3"           # muted teal
)

# ==========================================
# 8. 左 donut（无标注）
# ==========================================

left_df$label_pos <- (left_df$ymax + left_df$ymin) / 2

p_left <- ggplot(left_df) +
  geom_rect(aes(ymin=ymin, ymax=ymax, xmin=1, xmax=2, fill=Final_Class),
            color="white") +
  
  # ⭐ 环内标注
  geom_text(aes(
    x = 1.5,
    y = label_pos,
    label = paste0(Final_Class, "\n", scales::percent(Percentage, accuracy = 0.1))
  ), size = 4, fontface = "bold") +
  
  scale_fill_manual(values = color_left) +
  coord_polar(theta="y") +
  xlim(0,2) +
  theme_void() +
  labs(title="All Proteins") +
  theme(
    legend.position = "none",   # ❌ 去掉右侧图例
    plot.title = element_text(size=14, face="bold")
  )
  

secretome_annot=secretome_df[is.element(secretome_df[,1],use),]
# 右图（去掉.tsv等后缀）
right_df <- secretome_annot %>%
  filter(!is.na(Secretome_Class)) %>%
  mutate(Secretome_Class = gsub("\\..*", "", Secretome_Class)) %>%  # ⭐关键
  count(Secretome_Class) %>%
  mutate(Percentage = n / sum(n)) %>%
  arrange(desc(Percentage))

right_df$ymax <- cumsum(right_df$Percentage)
right_df$ymin <- c(0, head(right_df$ymax, -1))


right_df <- right_df %>%
  arrange(desc(Percentage))


right_df$Secretome_Class <- factor(
  right_df$Secretome_Class,
  levels = right_df$Secretome_Class
)

right_df$ymax <- cumsum(right_df$Percentage)
right_df$ymin <- c(0, head(right_df$ymax, -1))
right_df$label_pos <- (right_df$ymax + right_df$ymin) / 2

p_right <- ggplot(right_df) +
  geom_rect(aes(ymin=ymin, ymax=ymax, xmin=1, xmax=2, fill=Secretome_Class),
            color="white") +
  
  # ⭐ 图中只显示比例
  geom_text(aes(
    x = 1.5,
    y = label_pos,
    label = scales::percent(Percentage, accuracy = 0.1)
  ), size = 3.5) +
  
  scale_fill_manual(values = color_right) +
  
  coord_polar(theta="y", direction = -1) +  # ⭐ 顺时针
  xlim(0,2) +
  
  theme_void() +
  labs(title="Secretome Subtypes") +
  
  theme(
    legend.position = "right",   # ⭐ legend保留类别名
    plot.title = element_text(size=14, face="bold")
  )

# ==========================================
# 10. 拼接 + 图注
# ==========================================

p_final <- p_left + p_right +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Hierarchical Distribution of Secreted Proteins",
    caption = "Proteins were classified into secreted and non-secreted groups based on subcellular localization and HPA annotation. Secreted proteins were further subdivided into tissue-specific secretome categories."
  )

ggsave(file.path(figure1_dir, "secretome_clean.pdf"), p_final, width=7, height=4)


# =============================================================================
# 6. Figure 1h-i: 药物反应与生存比较
# =============================================================================

# =============================================================================

survial_plot=function(outname,surv_info,fit,title_name){
	if (nzchar(extra_r_lib)) .libPaths(c(.libPaths(), extra_r_lib))
	library(openxlsx)
	library(survival)
	library(survminer)
	pdf(outname, width = 2.8, height = 4)
	colour_man <- c("#427FB7","#BD7CA8")
	p1 <-  ggsurvplot(
		   data = surv_info,
		   fit,
		   pval = TRUE,
		   # conf.int = TRUE,
		   risk.table = TRUE,
		   risk.table.col = "strata",
		   linetype = "strata",
		   surv.median.line = "hv",
		   ggtheme = theme_bw(),
		   palette = colour_man,
		   title = title_name,
		   font.main = 10,     #标题字体大小
		   font.x = 10,        #x坐标字体大小
		   font.y = 10,        #y坐标字体大小
		   font.tickslab = 10,
		   risk.table.height = 0.3,
		   pval.size = 3,
		   risk.table.fontsize = 4,
		   font.legend = 8
		   # legend = "left"
		   # ggtheme = custom_theme()
		   )
	print(p1)
	dev.off()
	return(p1)
}

setwd(file.path(subtype_dir, "response_survival"))
if (nzchar(extra_r_lib)) .libPaths(c(.libPaths(), extra_r_lib))
	library(openxlsx)
	library(survival)
	library(survminer)
	#
	survival_info <- read.xlsx(file.path(exp_dir, "patients_clinical_information.xlsx"),na.strings = "NA",startRow = 2)
	surv_info <- survival_info[, c("Patient", "OS_time","OS_status", "Response", "iPFS_time", "iPFS_status")]
	surv_info=surv_info[grep("^Control",surv_info[,1],invert=T),]
	filter=read.table(file.path(diff_dir, "filter_sample.txt"),sep="\t",quote="",header=F,check.names=F)
	surv_info=surv_info[!is.element(surv_info[,1],filter[,1]),]
	is.element(surv_info[,"Response"],c("Response","SD","PD"))
	surv_info[which(surv_info[,4]=="PD"),4]="SD"
	#survival_os <- survival_os[!survival_os$Patient.Name%in%c("BC-2","LA-1-2"), ]#LA-1-2是用药后
	print(dim(surv_info))
	fit <- survfit(Surv(OS_time, OS_status) ~ Response, data = surv_info)
	fit_iPFS <- survfit(Surv(iPFS_time, iPFS_status) ~ Response, data = surv_info)
	surv_info_BC=surv_info[grep("BC",surv_info[,1]),]
	table(surv_info_BC[,4])
	surv_info_LA=surv_info[grep("LA",surv_info[,1]),]
	table(surv_info_LA[,4])
	fit_BC <- survfit(Surv(OS_time, OS_status) ~ Response, data = surv_info_BC)
	fit_iPFS_BC <- survfit(Surv(iPFS_time, iPFS_status) ~ Response, data = surv_info_BC)
	fit_LA <- survfit(Surv(OS_time, OS_status) ~ Response, data = surv_info_LA)
	fit_iPFS_LA <- survfit(Surv(iPFS_time, iPFS_status) ~ Response, data = surv_info_LA)
	title_name <- "response"
	outname=paste(title_name,"OS_survival.pdf",sep="_")
	p1=survial_plot(outname,surv_info,fit,title_name)

	title_name <- "response_iPFS"
	outname = paste(title_name, "OS_survival.pdf", sep = "_")
	p1 = survial_plot(outname, surv_info, fit_iPFS, title_name)
	#BC生存曲线

	title_name <- "BC_response"
	outname=paste(title_name,"OS_survival.pdf",sep="_")
	p1=survial_plot(outname,surv_info_BC,fit_BC,title_name)

	title_name <- "BC_response_iPFS"
	outname = paste(title_name, "OS_survival.pdf", sep = "_")
	p1 = survial_plot(outname, surv_info_BC, fit_iPFS_BC, title_name)
	#LA生存曲线
	title_name <- "LA_response"
	outname=paste(title_name,"OS_survival.pdf",sep="_")
	p1=survial_plot(outname,surv_info_LA,fit_LA,title_name)

	title_name <- "LA_response_iPFS"
	outname = paste(title_name, "OS_survival.pdf", sep = "_")
	p1 = survial_plot(outname, surv_info_LA, fit_iPFS_LA, title_name)


# =============================================================================
# 7. Figure 1j: 单因素/多因素 Cox 回归森林图
# =============================================================================

# =============================================================================

library(survival)
library(survminer)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(svglite)



# Build results table
extract_cox <- function(model, summary_obj, zph_obj, endpoint) {
  coef_df <- data.frame(
    Endpoint = endpoint,
    Variable = rownames(summary_obj$coefficients),
    HR = summary_obj$coefficients[, "exp(coef)"],
    HR_lower95 = summary_obj$conf.int[, "lower .95"],
    HR_upper95 = summary_obj$conf.int[, "upper .95"],
    z = summary_obj$coefficients[, "z"],
    p_value = summary_obj$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  # Map coefficient names to zph table rownames
  zph_p <- zph_obj$table[, "p"]  # named vector
  global_p <- zph_p["GLOBAL"]
  # Create mapping: strip factor level suffixes
  base_names <- gsub("Group 1$", "", coef_df$Variable)
  base_names <- gsub("BC$", "", base_names)
  base_names <- gsub("Female$", "", base_names)
  coef_df$PH_p_value <- zph_p[match(base_names, names(zph_p))]
  coef_df$PH_global_p <- global_p
  rownames(coef_df) <- NULL
  coef_df
}

info=read.table(file.path(exp_dir, "patient_source.txt"),sep="\t",quote="",header=F,check.names=F)
rownames(info)=info[,1]
nanjing=rownames(info)[which(is.element(info[,2],c("nanjing")))]
henan=rownames(info)[which(is.element(info[,2],c("henan")))]
fujian=rownames(info)[which(is.element(info[,2],c("fujian")))]
heilongjiang=rownames(info)[which(is.element(info[,2],c("heilongjiang")))]
con=c(fujian,heilongjiang)
filter=read.table(file.path(diff_dir, "filter_sample.txt"),sep="\t",quote="",header=F,check.names=F)
con=setdiff(con,filter[,1])
LA=info[grep("LA",info[,1]),1]
LA=setdiff(LA,filter[,1])
BC=info[grep("BC",info[,1]),1]
BC=setdiff(BC,filter[,1])
nanjing=setdiff(nanjing,filter[,1])
henan=setdiff(henan,filter[,1])

info=info[c(henan,nanjing),]
colnames(info) <- c("Patient", "NMF_group")



# Read clinical Excel (skip title row 1, header is row 2)
clinical <- read_excel(file.path(exp_dir, "patients_clinical_information.xlsx"),
                       sheet = 1, skip = 1)

# Merge: keep only NMF patients
dat <- merge(info, clinical, by = "Patient", all.x = TRUE)


# Type conversion
dat$OS_status   <- as.integer(dat$OS_status)
dat$OS_time     <- as.numeric(dat$OS_time)
dat$iPFS_status <- as.integer(dat$iPFS_status)
dat$iPFS_time   <- as.numeric(dat$iPFS_time)
dat$Age         <- as.numeric(dat$Age)
dat$Sex         <- factor(dat$Sex, levels = c("Male", "Female"))
dat$Cancer      <- factor(dat$Cancer, levels = c("LA", "BC"))  # LA as reference



dat[which(dat[,"Response"]=="PD"),"Response"]="SD"
dat$Response <- factor(dat$Response, levels = c("PR", "SD"))
dat$Source <- factor(dat$Source, levels = c("Nanjing_Brain", "Henan"))
dat$MRI <- factor(dat$MRI, levels = c("Positive", "Negative"))
dat$Cytology <- factor(dat$Cytology, levels = c("Positive", "Negative"))
# OS Cox model
os_cox <- coxph(Surv(OS_time, OS_status) ~  Age + Cancer + Sex+Response,
                data = dat)
os_summary <- summary(os_cox)
os_zph <- cox.zph(os_cox)
os_results   <- extract_cox(os_cox, os_summary, os_zph, "OS")
# iPFS Cox model
ipfs_cox <- coxph(Surv(iPFS_time, iPFS_status) ~ Age + Cancer + Sex+Response,
                  data = dat)
ipfs_summary <- summary(ipfs_cox)
ipfs_zph <- cox.zph(ipfs_cox)
ipfs_results <- extract_cox(ipfs_cox, ipfs_summary, ipfs_zph, "iPFS")
cox_results  <- rbind(os_results, ipfs_results)

## ---- 5. Forest Plot --------------------------------------------------------

# Prepare data for forest plot
forest_data <- rbind(
  data.frame(
    endpoint = "OS",
    variable = rownames(os_summary$coefficients),
    HR = os_summary$coefficients[, "exp(coef)"],
    lower = os_summary$conf.int[, "lower .95"],
    upper = os_summary$conf.int[, "upper .95"],
    pval = os_summary$coefficients[, "Pr(>|z|)"]
  ),
  data.frame(
    endpoint = "iPFS",
    variable = rownames(ipfs_summary$coefficients),
    HR = ipfs_summary$coefficients[, "exp(coef)"],
    lower = ipfs_summary$conf.int[, "lower .95"],
    upper = ipfs_summary$conf.int[, "upper .95"],
    pval = ipfs_summary$coefficients[, "Pr(>|z|)"]
  )
)

# Clean variable names for display

forest_data$variable <- gsub("CancerBC", "Cancer: BC vs LA", forest_data$variable)
forest_data$variable <- gsub("SexFemale", "Sex: Female vs Male", forest_data$variable)
forest_data$variable <- gsub("Age", "Age (per year)", forest_data$variable)

forest_data$label <- sprintf("%.2f (%.2f-%.2f)", forest_data$HR, forest_data$lower, forest_data$upper)
forest_data$pval_label <- ifelse(forest_data$pval < 0.001, "< 0.001",
                                  sprintf("%.3f", forest_data$pval))
forest_data$y_label <- paste0(forest_data$endpoint, " — ", forest_data$variable)
forest_data$y <- nrow(forest_data):1
svg_font="Arial"
# Forest plot
p_forest <- ggplot(forest_data, aes(x = HR, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 0.6) +
  geom_point(size = 3, shape = 21, fill = ifelse(forest_data$endpoint == "OS", "#0072B2", "#D55E00"),
             color = "black", stroke = 0.5) +
  geom_text(aes(x = max(forest_data$upper) * 1.15, label = label),
            hjust = 0, size = 3.2) +
  geom_text(aes(x = max(forest_data$upper) * 2.5, label = pval_label),
            hjust = 0, size = 3.2,
            fontface = ifelse(forest_data$pval < 0.05, "bold", "plain")) +
  scale_y_continuous(breaks = forest_data$y, labels = forest_data$y_label,
                     expand = c(0, 0.5)) +
  scale_x_log10() +
  coord_cartesian(xlim = c(0.2, max(forest_data$upper) * 4)) +
  labs(title = "Adjusted Cox Regression: Hazard Ratios",
       subtitle = "Reference: Group 2, LA, Male",
       x = "Hazard Ratio (log scale)",
       y = NULL) +
  annotate("text", x = max(forest_data$upper) * 1.15, y = max(forest_data$y) + 0.8,
           label = "HR (95% CI)", hjust = 0, size = 3.2, fontface = "bold") +
  annotate("text", x = max(forest_data$upper) * 2.5, y = max(forest_data$y) + 0.8,
           label = "P value", hjust = 0, size = 3.2, fontface = "bold") +
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
  
pdf(file.path(figure1_dir, "forest_plot.pdf"),width=12,height=5)
p_forest
dev.off()
