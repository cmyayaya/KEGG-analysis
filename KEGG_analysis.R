rm(list=ls())

setwd("D:/analysis")

genelist <- read.csv("genelist.csv",header =TRUE)
library(clusterProfiler)  
library(dplyr)  
library(ggplot2)  


# 筛选上调和下调基因  
up_genes <- genelist %>%   
  filter(log2FoldChange > 1.2, padj < 0.05) %>%   
  pull(entrez)  

down_genes <- genelist %>%   
  filter(log2FoldChange < -1.2, padj < 0.05) %>%   
  pull(entrez)  

# KEGG 富集分析 - 上调基因  
kegg_up <- enrichKEGG(  
  gene = up_genes,  
  organism = 'hsa',  
  pvalueCutoff = 0.05,  
  qvalueCutoff = 0.2  
)  

# KEGG 富集分析 - 下调基因  
kegg_down <- enrichKEGG(  
  gene = down_genes,  
  organism = 'hsa',  
  pvalueCutoff = 0.05,  
  qvalueCutoff = 0.2  
)  
# 检查 KEGG 结果  
if (!is.null(kegg_up) && nrow(kegg_up@result) > 0) {  
  # 准备数据  
  plot_data_up <- kegg_up@result %>%  
    mutate(neg_log10_padj = -log10(p.adjust)) %>%  
    arrange(desc(neg_log10_padj)) %>%  
    head(20)  # 选择前10个通路  
  
  # 绘制条形图  
  ggplot(plot_data_up, aes(x = neg_log10_padj, y = reorder(Description, neg_log10_padj), fill = p.adjust)) +  
    geom_bar(stat = "identity") +  
    labs(  
      title = "Top 10 Up-regulated Genes KEGG Enrichment",  
      x = "-LOG10(Adjusted P-value)",  
      y = "KEGG Pathway"  
    ) +  
    scale_fill_gradient(low = "#9575CD", high = "#311B92", name = "Adjusted P-value") +  
    theme_minimal() +  
    theme(axis.text.y = element_text(size = 8))  
  
  # 保存图形  
  ggsave("kegg_up_ggplot_barplot.pdf", width = 10, height = 6)  
}  
