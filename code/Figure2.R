#!/usr/bin/env Rscript


project_root <- normalizePath(
  Sys.getenv("PROTEO_PROJECT_ROOT", unset = getwd()),
  mustWork = FALSE
)
exp_dir <- file.path(project_root, "exp")
diff_dir <- file.path(project_root, "diff")
add_flag_script <- Sys.getenv(
  "PROTEO_ADD_FLAG_SCRIPT",
  unset = file.path(project_root, "scripts", "add_flag.r")
)


# =============================================================================
# 1. Figure 2a: 差异基因数目图
# =============================================================================

# =============================================================================

# -----------------------------------------------------------------------------
# 1.1 生成各比较组的差异分析结果
# -----------------------------------------------------------------------------

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
dif_exp_limma <- function(exp_control, exp_case, logFC_cutoff, pValue_cutoff){
        library(limma)
        exp_mat <- cbind(exp_control, exp_case)
        group <- factor(c(rep("control", dim(exp_control)[2]), rep("case", dim(exp_case)[2])), levels = c("control", "case"))
        design <- model.matrix(~group)
        fit <- lmFit(exp_mat, design)
        fit <- eBayes(fit)
        DEGs <- topTable(fit, coef=2, adjust.method='fdr', n=Inf, sort.by="logFC")

        #设定方向
        DEGs$direction <- as.factor(ifelse(DEGs$adj.P.Val < pValue_cutoff &
        abs(DEGs$logFC) > logFC_cutoff, ifelse(DEGs$logFC > logFC_cutoff, "UP", "DOWN"), "NOT"))

        return(DEGs)
}
logFC_cutoff <- 1
pValue_cutoff <- 0.01

DEGs <- dif_exp_limma(exp_imp1[,heilongjiang],exp_imp1[,nanjing], logFC_cutoff, pValue_cutoff)
write.table(DEGs,file.path(diff_dir, "nanjing_heilongjiang_filter_diff.txt"),sep="\t",quote=F)
DEGs <- dif_exp_limma(exp_imp1[,fujian],exp_imp1[,nanjing], logFC_cutoff, pValue_cutoff)
write.table(DEGs,file.path(diff_dir, "nanjing_fujian_filter_diff.txt"),sep="\t",quote=F)
DEGs <- dif_exp_limma(exp_imp1[,heilongjiang],exp_imp1[,henan], logFC_cutoff, pValue_cutoff)
write.table(DEGs,file.path(diff_dir, "henan_fujian_filter_diff.txt"),sep="\t",quote=F)
DEGs <- dif_exp_limma(exp_imp1[,fujian],exp_imp1[,henan], logFC_cutoff, pValue_cutoff)
write.table(DEGs,file.path(diff_dir, "henan_heilongjiang_filter_diff.txt"),sep="\t",quote=F)

DEGs <- dif_exp_limma(exp_imp1[,con],exp_imp1[,BC], logFC_cutoff, pValue_cutoff)
write.table(DEGs,file.path(diff_dir, "BC_control_diff.txt"),sep="\t",quote=F)
DEGs <- dif_exp_limma(exp_imp1[,con],exp_imp1[,LA], logFC_cutoff, pValue_cutoff)
write.table(DEGs,file.path(diff_dir, "LA_control_diff.txt"),sep="\t",quote=F)
DEGs <- dif_exp_limma(exp_imp1[,con],exp_imp1[,c(nanjing,henan)], logFC_cutoff, pValue_cutoff)
write.table(DEGs,file.path(diff_dir, "case_control_filter_diff.txt"),sep="\t",quote=F)

DEGs <- dif_exp_limma(exp_imp1[,BC],exp_imp1[,LA], logFC_cutoff, pValue_cutoff)
write.table(DEGs,file.path(diff_dir, "BC_LA_filter_diff.txt"),sep="\t",quote=F)


# -----------------------------------------------------------------------------
# 1.2 绘制差异基因数目/倍数变化分布图
# -----------------------------------------------------------------------------

###  Figure 4A  ###
#install.packages("ggnewscale")
#install.packages("tidyverse")
library(ggnewscale)
library(tidyverse)
library(ggplot2)
library(ggrepel)
setwd(diff_dir)
BC <- read.table("BC_control_diff.txt",stringsAsFactors = F,
                            row.names = 1,check.names = F)
LA <- read.table("LA_control_diff.txt",stringsAsFactors = F,
                              row.names = 1,check.names = F)
ALL <- read.table("case_control_filter_diff.txt",stringsAsFactors = F,
                 row.names = 1,check.names = F)

BC_LA <- read.table("BC_LA_filter_diff.txt",stringsAsFactors = F,
                 row.names = 1,check.names = F)

### top 10 sig
BC.new <- BC[c(which(BC$direction == "UP"),
               which(BC$direction == "DOWN")),]
BC.new$direction <- ifelse(BC.new$direction == "UP", "BC-UP", "control-UP(BC)")
BC.new2 <- BC.new[order(BC.new$logFC),]
top10.direction.BC <- BC.new2[c(1:5,(dim(BC.new2)[1]-4):dim(BC.new2)[1]),]
top10.direction.BC$label <- rep("BC vs control",10)




ALL.new <- ALL[c(which(ALL$direction == "UP"),
               which(ALL$direction == "DOWN")),]
ALL.new$direction <- ifelse(ALL.new$direction == "UP", "pwLM-UP", "control-UP(pwLM)")
ALL.new2 <- ALL.new[order(ALL.new$logFC),]
top10.direction.ALL <- ALL.new2[c(1:5,(dim(ALL.new2)[1]-4):dim(ALL.new2)[1]),]
top10.direction.ALL$label <- rep("pwLM vs control",10)


LA.new <- LA[c(which(LA$direction == "UP"),
               which(LA$direction == "DOWN")),]
LA.new$direction <- ifelse(LA.new$direction == "UP", "LA-UP", "control-UP(LA)")
LA.new2 <- LA.new[order(LA.new$logFC),]
top10.direction.LA <- LA.new2[c(1:5,(dim(LA.new2)[1]-4):dim(LA.new2)[1]),]
top10.direction.LA$label <- rep("LA vs control",10)


BC_LA.new <- BC_LA[c(which(BC_LA$direction == "UP"),
               which(BC_LA$direction == "DOWN")),]
BC_LA.new$direction <- ifelse(BC_LA.new$direction == "UP", "LA-UP(BC)", "BC-UP(LA)")
BC_LA.new2 <- BC_LA.new[order(BC_LA.new$logFC),]
top10.direction.BC_LA <- BC_LA.new2[c(1:5,(dim(BC_LA.new2)[1]-4):dim(BC_LA.new2)[1]),]
top10.direction.BC_LA$label <- rep("LA vs BC",10)


BC.plot <- BC[which(BC$direction != "NOT"),] 
BC.plot$label <- rep("BC vs control",dim(BC.plot)[1])
BC.plot$ptname <- rownames(BC.plot)
BC.plot$direction <- ifelse(BC.plot$direction == "UP", "BC-UP", "control-UP(BC)")

LA.plot <- LA[which(LA$direction != "NOT"),] 
LA.plot$label <- rep("LA vs control",dim(LA.plot)[1])
LA.plot$ptname <- rownames(LA.plot)
LA.plot$direction <- ifelse(LA.plot$direction == "UP", "LA-UP", "control-UP(LA)")

ALL.plot <- ALL[which(ALL$direction != "NOT"),] 
ALL.plot$label <- rep("pwLM vs control",dim(ALL.plot)[1])
ALL.plot$ptname <- rownames(ALL.plot)
ALL.plot$direction <- ifelse(ALL.plot$direction == "UP", "pwLM-UP", "control-UP(pwLM)")


BC_LA.plot <- BC_LA[which(BC_LA$direction != "NOT"),] 
BC_LA.plot$label <- rep("LA vs BC",dim(BC_LA.plot)[1])
BC_LA.plot$ptname <- rownames(BC_LA.plot)
BC_LA.plot$direction <- ifelse(BC_LA.plot$direction == "UP", "LA-UP(BC)", "BC-UP(LA)")




merge.plot <- rbind(BC.plot,
                    LA.plot,
                    ALL.plot,
					BC_LA.plot)

merge.plot.top10 <- rbind(top10.direction.BC,
                          top10.direction.LA,
                          top10.direction.ALL,
						  top10.direction.BC_LA)

merge.plot.top10$gene <- rownames(merge.plot.top10)

merge.plot$detailed_group <- paste0(merge.plot$label, "_", merge.plot$direction)
merge.plot.top10$detailed_group <- paste0(merge.plot.top10$label, "_", merge.plot.top10$direction)

dbar <- 
  merge.plot.top10 %>% 
  group_by(label) %>% 
  summarise_all(list(min = min, max = max)) %>% 
  select(label, logFC_min, logFC_max, label) %>% 
  rename(label = label)
dbar$label=factor(dbar$label,levels=c("pwLM vs control","LA vs control","BC vs control","LA vs BC"))
pdf(file.path(diff_dir, "diff_num.pdf"),width=6,height=4)
ggplot()+
  geom_col(data = dbar,  # 绘制负向背景柱状图
           mapping = aes(x = label,y = logFC_min),
           fill = "#C3C3C3",alpha = 0.6, width = 0.9) +
  geom_col(data = dbar, # 绘制正向背景柱状图
           mapping = aes(x = label,y = logFC_max),
           fill = "#C3C3C3",alpha = 0.6, width = 0.9) +
  geom_jitter(data=merge.plot,
              aes(x=label,y=logFC,
                  color=direction),width =0.4,size=1) +
  geom_tile(data = merge.plot.top10, # 绘制中心分组标记图
            aes(x = label,
                y = 0,
                fill = label),
            height=0.5,
            color = "black",
            alpha = 0.6,
            show.legend = F) +
  ggsci::scale_fill_npg() + # 自定义颜色
  ggsci::scale_color_npg() + # 自定义颜色

  scale_color_manual(values = c('BC-UP' = '#64ACE8',           # 红色
								'control-UP(BC)' = '#8FB6CC',  # 绿色
                                'LA-UP' = '#E9B552',           # 蓝色
                                'control-UP(LA)' = '#F0CD8B' ,  
								'pwLM-UP' ='#E86E5A',
								'control-UP(pwLM)'='#F2B0A5', 
                                'LA-UP(BC)' = '#6C61CE',
                                'BC-UP(LA)' = '#777DCE'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("log (Fold-change)")+
  
# 去掉网格线的设置
  theme(
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.background = element_blank(),  # 去掉背景色
    axis.line = element_line(color = "black")  # 保留坐标轴线
  )
dev.off()

table(BC_LA.plot[,"direction"])
table(ALL.plot[,"direction"])
table(LA.plot[,"direction"])
table(BC.plot[,"direction"])


# =============================================================================
# 2. Figure 2b: 差异基因韦恩图
# =============================================================================

# =============================================================================

#conda activate r-env
library("ggVennDiagram")
library("ggplot2")
data1=read.table(file.path(diff_dir, "BC_control_diff.txt"), header = TRUE, sep = "\t", check.names = FALSE)
data2=read.table(file.path(diff_dir, "LA_control_diff.txt"), header = TRUE, sep = "\t", check.names = FALSE)
data3=read.table(file.path(diff_dir, "case_control_filter_diff.txt"), header = TRUE, sep = "\t", check.names = FALSE)
diff_1=rownames(data1[which(abs(data1[,1]) >1 & data1[,5]<0.01),])
diff_2=rownames(data2[which(abs(data2[,1]) >1 & data2[,5]<0.01),])
diff_3=rownames(data3[which(abs(data3[,1]) >1 & data3[,5]<0.01),])
write.table(diff_1,"BC_control_gene.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(diff_2,"LA_control_gene.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(diff_3,"case_control_gene.txt",sep="\t",quote=F,row.names=F,col.names=F)
# ========= 基础版：3种颜色 =========
p=ggVennDiagram(list(BC_control=diff_1,LA_control=diff_2,case_control=diff_3))+
  # 使用3种颜色的渐变填充
  scale_fill_gradient(low = "#B0DDAD", high = "#74C476") +
  # 设置3个分类的边框颜色
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))  +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )
pdf(file.path(diff_dir, "venn.pdf"),width=4,height=4)
p
dev.off()

###https://bioinfogp.cnb.csic.es/tools/venny/ 使用网址做的


# =============================================================================
# 3. Figure 2c: 差异基因火山图
# =============================================================================

# =============================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. 读入数据
data <- read.table("case_control_filter_diff.txt", header = TRUE, sep = "\t", check.names = FALSE)
# 假设第一列是基因名，如果不是，请用 data$Symbol <- rownames(data)
data$Symbol <- rownames(data) 

# 2. 设置阈值
PValue_cut <- 0.01
FC_cut <- 1

# 3. 定义目标基因列表
epithelial_genes <- c("WFDC2","SDC1","EFNA1","SLPI","SPINT1")
breast_genes <- c("KRT6B","PIP")
lung_genes <- c("SFTPB","SFTPA2","PIGR","SFTPD")

# 4. 核心逻辑：定义点的颜色(sig)和标签的颜色(label_group)
data <- data %>%
  mutate(
    neg_log10P = -log10(adj.P.Val),
    # 点的颜色分类：下调蓝，上调红，不显灰
    sig = case_when(
      adj.P.Val < PValue_cut & logFC <= -FC_cut ~ "Down",
      adj.P.Val < PValue_cut & logFC >= FC_cut  ~ "Up",
      TRUE ~ "NotSig"
    ),
    # 标签的分组逻辑（仅限指定基因）
    label_group = case_when(
      Symbol %in% epithelial_genes ~ "Epithelial_Label",
      Symbol %in% breast_genes     ~ "Breast_Label",
      Symbol %in% lung_genes       ~ "Lung_Label",
      TRUE ~ NA_character_
    ),
    # 提取要显示的文字
    label_text = ifelse(!is.na(label_group), Symbol, "")
  )

# 5. 绘图
p <- ggplot(data, aes(x = logFC, y = neg_log10P)) +
  # 绘制所有散点：下调蓝，不显灰，上调红
  geom_point(aes(color = sig), alpha = 0.8, size = 1.2) +
  scale_color_manual(values = c("Down" = "blue", "NotSig" = "gray80", "Up" = "red")) +
  
  # 辅助线
  geom_vline(xintercept = c(-FC_cut, FC_cut), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(PValue_cut), linetype = "dashed", color = "grey60") +
  
  # 标注指定基因名（使用不同于红蓝灰的特定颜色）
  geom_text_repel(
    aes(label = label_text, color = label_group),
    size = 4.5,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.8,
    point.padding = 0.3,
    segment.color = "black", # 连线用黑色方便识别
    force = 10,
    show.legend = FALSE # 取消标签颜色产生的图注
  ) +
  
  # 为标签颜色单独设置映射（不同于点的红蓝灰）
  # 使用 scale_color_manual 的扩展，注意：由于 ggplot 同一个 aes(color) 只能有一套 palette，
  # 这里的 values 需要包含上面 sig 的颜色，否则会报错。
  scale_color_manual(values = c(
    "Down" = "#0681D7", 
    "NotSig" = "gray80", 
    "Up" = "#D14424",
    "Epithelial_Label" = "#008B8B", # 上皮：深湖蓝 (DarkCyan)
    "Breast_Label"     = "#C71585", # 乳腺：中紫红 (MediumVioletRed)
    "Lung_Label"       = "#B8860B"  # 肺：暗金黄 (DarkGoldenrod)
  )) +
  
  # 样式美化
  theme_classic() +
  labs(
    title = "Volcano Plot: Highlighted Markers",
    x = "log2(Fold Change)",
    y = "-log10(Adjusted P-value)"
  ) +
  theme(
    legend.position = "none", # 取消右侧所有类别图注
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )


ggsave("volcano.pdf", width = 5, height = 5)


# =============================================================================
# 4. Figure 2d: 差异基因热图
# =============================================================================

# =============================================================================


setwd(diff_dir)
load(file.path(diff_dir, "protein_imp_gene.rds"))
gene_list <- read.table("case_control_filter_diff.txt", header=TRUE, stringsAsFactors=FALSE)
filter=read.table(file.path(diff_dir, "filter_sample.txt"),sep="\t",quote="",header=F,check.names=F)


gene_strik <- gene_list[which(abs(gene_list[,1]) >1 & gene_list[,5]<0.01),]
exp_imp1=gene_exp
rownames(exp_imp1)=exp_imp1[,1]
exp_imp1=exp_imp1[,-1]
exp_imp1=matrix(as.numeric(exp_imp1), nrow = nrow(exp_imp1), ncol = ncol(exp_imp1))
rownames(exp_imp1)=gene_exp[,1]
colnames(exp_imp1)=colnames(gene_exp)[-1]
gene_exp=exp_imp1
gene_exp=gene_exp[,setdiff(colnames(gene_exp),filter[,1])]

exp_strik <- gene_exp[rownames(gene_strik), c(91:110,1:90)]

group_list <- colnames(exp_strik)
new_names <- strsplit(group_list, "-")
sample_types <- sapply(new_names, function(x) x[1])
group_final=ifelse(grepl("^Control", group_list), "Control","LM")
heatmap_data=exp_strik


annotation_col <- data.frame(
  Group = factor(group_final, levels = c("Control", "LM")),
  Type = factor(sample_types, levels = c("Control","BC", "LA","MM","EC"))
)
rownames(annotation_col) <- group_list 

# 设置注释颜色
annotation_colors <- list(
  Type = c("Control" = "#CDE2F2","BC" = "#8ECFC9", "LA" = "#FFBE7A","EC" = "#FA7F6F","MM" = "#82B0D2"),
  Group = c("Control" = "#CDE2F2", "LM" = "#FF9A91")
)

#  绘制热图（按分组排序列，行聚类）
library(pheatmap)
library(grid)
source(add_flag_script)
gene_name=c("HAVCR2","LAYN", # T exhaustion基因高表达，预示LM中T细胞耗竭
			"CTSD","CTSC","CTSB","CTSA","LYVE1","MMP14", #M2 基因
			"WFDC2","SDC1","EFNA1","SLPI","SPINT1", #epithelial_genes
			"KRT6B","PIP",# breast_genes
			"SFTPB","SFTPA2","PIGR","SFTPD", #lung_genes 
			"APOC3","APOB","APOC4","APOC2","APOM" #正常里上调的
			)
heatmap_result <- pheatmap(heatmap_data,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE,      # 列不聚类，保持分组顺序
         cluster_rows = TRUE,       # 行聚类（默认就是TRUE）
         scale = "row",
         show_rownames = T,
         show_colnames = FALSE,
         treeheight_row = 0,       # 显示行聚类树的高度
         treeheight_col = 0,        # 列不聚类所以高度为0
         color = colorRampPalette(rev(c("#DD0514",rgb(255,255,255,maxColorValue = 255),"#4874B1")))(100),
         breaks = seq(-2,2,length.out = 101),
         main = "heatmap",
         #gaps_col = gaps_position,  # 在分组1和分组2之间添加分隔线
         annotation_names_col = TRUE)



pdf(file.path(diff_dir, "marker_heatmap.pdf"),width=5,height=6)
add.flag(heatmap_result,kept.labels = gene_name,repel.degree = 0)
dev.off()


# =============================================================================
# 5. Figure 2e: KEGG/GSEA 功能富集图
# =============================================================================

# =============================================================================

result=read.table(file.path(diff_dir, "case_control_filter_diff.txt"),sep="\t",header=T) 
library(clusterProfiler)
geneList <- result$logFC
names(geneList) <- rownames(result)  # 或ENTREZID
library(org.Hs.eg.db)
# 按值排序
geneList <- sort(geneList, decreasing = TRUE)
# 转换为 clusterProfiler 需要的格式
library(msigdbr)

gene_symbols <- names(geneList)

# 用org.Hs.eg.db转换
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = gene_symbols,
  column = "ENTREZID",
  keytype = "SYMBOL"
)

# 过滤掉未映射到的基因
valid_idx <- !is.na(entrez_ids)
geneList_entrez <- geneList[valid_idx]
names(geneList_entrez) <- entrez_ids[valid_idx]

# 重新排序（转换后顺序可能变化）
geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)

# 然后使用geneList_entrez进行gseKEGG
gsea_result <- gseKEGG(
  geneList = geneList_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)
result_df <- as.data.frame(gsea_result)
head(result_df)
sig_results <- subset(result_df, p.adjust < 0.05)
print(sig_results)
write.table(sig_results, file.path(diff_dir, "GSEA", "GSEA_KEGG_results.csv"), row.names = FALSE,quote=F,sep="\t")

library(enrichplot)
library(ggplot2)

p=ridgeplot(gsea_result,
          showCategory = 30,
          fill = "NES", #填充色 "pvalue", "p.adjust", "qvalue" 
          core_enrichment = TRUE,#是否只使用 core_enriched gene
          label_format = 30,
          orderBy = "NES",
          decreasing = FALSE
          )+
		   scale_fill_gradientn(
    colors = c("#74ADD1", "#ABD9E9", 
               "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61"))+
  theme(axis.text.y = element_text(size=8))
pdf(file.path(diff_dir, "GSEA", "ridgeplot_KEGG.pdf"),width=5,height=10)
p
dev.off()

gsego_bp <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  keyType      = "SYMBOL",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose      = TRUE,
  seed         = 123  # 设置随机种子确保结果可重复
)
result_df <- as.data.frame(gsego_bp)
sig_results <- subset(result_df, p.adjust < 0.01)
print(sig_results)
dim(sig_results)
write.table(sig_results, file.path(diff_dir, "GSEA", "GSEA_BP_results.csv"), row.names = FALSE,quote=F,sep="\t")


# =============================================================================
# 6. Figure 2f: Hallmark 功能富集图
# =============================================================================

# =============================================================================

result=read.table(file.path(diff_dir, "case_control_filter_diff.txt"),sep="\t",header=T) 
library(clusterProfiler)
library(msigdbr)
msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
# 只保留基因集和对应的基因
hallmark_gene_sets <- msigdb_hallmark %>%
    dplyr::select(gs_name, human_gene_symbol)
library(dplyr)

gene_sets <- hallmark_gene_sets %>%
    group_by(gs_name) %>%
    summarise(genes = list(unique(human_gene_symbol))) 
gene_sets_list <- setNames(
  lapply(gene_sets$genes, function(x) x),
  gene_sets$gs_name
)
term2gene <- data.frame(
  term = rep(names(gene_sets_list), times = sapply(gene_sets_list, length)),
  gene = unlist(gene_sets_list)
)

geneList <- result$logFC
names(geneList) <- rownames(result)  # 或ENTREZID

# 按值排序
geneList <- sort(geneList, decreasing = TRUE)
gsea_result <- GSEA(geneList,
                     TERM2GENE = data.frame(term=term2gene[,1], gene=term2gene[,2]),
                     pvalueCutoff=0.05)
result_df <- as.data.frame(gsea_result)
sig_results <- subset(result_df, p.adjust < 0.05)
print(sig_results)
write.table(sig_results, file.path(diff_dir, "GSEA", "GSEA_hallmark_results.csv"), row.names = FALSE,quote=F,sep="\t")

library(enrichplot)
library(ggplot2)

p=ridgeplot(gsea_result,
          showCategory = 20,
          fill = "NES", #填充色 "pvalue", "p.adjust", "qvalue"
          core_enrichment = TRUE,#是否只使用 core_enriched gene
          label_format = 30,
          orderBy = "NES",
          decreasing = FALSE
          )+
		   scale_fill_gradientn(
    colors = c("#74ADD1", "#ABD9E9",
               "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61"))+
  theme(axis.text.y = element_text(size=8))
## Picking joint bandwidth of 0.212

pdf(file.path(diff_dir, "GSEA", "ridgeplot_Hallmark.pdf"),width=5,height=4)
p
dev.off()


# =============================================================================
# 7. Figure 2g: 差异基因网络与模块识别
# =============================================================================

# =============================================================================

protein=read.table(file.path(diff_dir, "String", "9606.protein.links.addsymbol.v12.0.txt"),sep="\t",quote="",header=T)
dim(protein)
result=read.table(file.path(diff_dir, "case_control_filter_diff.txt"),sep="\t",quote="",header=T)

diff=rownames(result[which(abs(result[,1]) >1 & result[,5]<0.01),])
length(diff)
temp=which(is.element(protein[,1],diff))
temp1=which(is.element(protein[,2],diff))
use=protein[intersect(temp,temp1),]
length(unique(c(use[,1],use1[,2])))
dim(use)
length(which(use[,5]>900))

use1=use[which(use[,5]>900),]
library(dplyr)
edges_unique <- use1 %>%
  mutate(
    gene1 = pmin(symbol1, symbol2),
    gene2 = pmax(symbol1, symbol2),
  ) %>%
  distinct(gene1, gene2, .keep_all = TRUE) %>%
  select(from = gene1, to = gene2)
length(unique(c(edges_unique[,1],edges_unique[,2]))) #822

write.table(edges_unique,file.path(diff_dir, "String", "case_control_network_900.txt"),sep="\t",quote=F,row.names=F,col.names=F)

library(igraph)
g <- graph_from_data_frame(edges_unique[,1:2], directed = FALSE)
g <- decompose.graph(g)
g <- g[[which(sapply(g, vcount) == max(sapply(g, vcount)))]]
edges_matrix <- as_edgelist(g)
edges_largest <- data.frame(from = edges_matrix[, 1],
                            to = edges_matrix[, 2],
                            stringsAsFactors = FALSE)
write.table(edges_largest,file.path(diff_dir, "String", "case_control_network_900.txt"),sep="\t",quote=F,row.names=F,col.names=F)


# =============================================================================
# 8. Figure 2h: 模块功能富集
# =============================================================================

# =============================================================================

gene=read.table(file.path(diff_dir, "String", "mcode.txt"),sep="\t")
library("org.Hs.eg.db")
library("clusterProfiler")
setwd(file.path(diff_dir, "String", "clusterprofile"))
for(i in 1:10){
	name=unlist(strsplit(gene[i,5],", "))
	g_id <- mget(name, org.Hs.egSYMBOL2EG,ifnotfound=NA)
	g_id <- as.character(g_id)
	ego = enrichGO(gene=g_id, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH")
	ego.KEGG = enrichKEGG(gene=g_id, organism = "hsa",keyType = "kegg", pAdjustMethod = "BH")
	outFiletxt=paste("cluster",i,".txt",sep="")
	write.table(ego,outFiletxt,col.names=T,row.names=F,quote=F,sep="\t")
	outFiletxt=paste("cluster",i,"KEGG.txt",sep="")
	write.table(ego.KEGG,outFiletxt,col.names=T,row.names=F,quote=F,sep="\t")
}
