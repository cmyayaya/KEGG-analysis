rm(list=ls())

setwd("D:/6-Mass specta/202406/20250102_analysis")

genelist <- read.csv("differential_expression_results.csv",header =TRUE)

read.csv("20250102_analysis/differential_expression_results.csv",header =TRUE)

############################去除第一列中编号的后缀


library(stringi)##加载包
genelist$id=stri_sub(genelist$id,1,15)##保留前15位

write.csv(data,"results.csv",quote = F)
#####################################

genelist <- read.table("GSE175815_BeWo_logCPM_exprTable.txt",header = T)
write.csv(genelist,"genelist.csv",quote=F,row.names = F)

###install.packages('R.utils')
getOption("clusterProfiler.download.method") #查看数据下载协议是什么
R.utils::setOption( "clusterProfiler.download.method",'auto' ) #更改数据下载的协议，改为“auto”
#R.utils::setOption( "clusterProfiler.download.method",'wininet' ) #更改数据下载的协议，改为“wininet”




##BiocManager::install("org.Mm.eg.db",force = TRUE)
options(connectionObserver = NULL)
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(org.Hs.eg.db))
library(org.Mm.eg.db)
library(org.Hs.eg.db)

head(genelist$symbol)
length(genelist$symbol)
genelist <- genelist[!is.na(genelist$symbol) & genelist$symbol != "", ]
genelist$entrez2 <- mapIds(org.Hs.eg.db,
                           keys = genelist$symbol,
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")

# ENSEMBL    SYMBOL  ENTREZID
genelist$SYMBOL <- mapIds(org.Hs.eg.db,
                          keys=genelist$Row.names,
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")




genelist$entrez2 <- mapIds(org.Hs.eg.db,
                           keys=genelist$SYMBOL,
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")


test <- data.frame(genelist)
print.csv(test,"genelist.csv",quote=F,row.names = F) #保存在本地

#genelist <- read.csv("genelist.csv",header =TRUE)

# 1. 更新 BiocManager 到最新版本  


if ("KEGG.db" %in% installed.packages()) {  
  remove.packages("KEGG.db")  
}  
# 1. 安装必要的包
options(repos = c(  
  CRAN = "https://cloud.r-project.org/",  
  BioCsoft = "https://bioconductor.org/packages/3.19/bioc",  
  BioCann = "https://bioconductor.org/packages/3.19/data/annotation"  
))  
if (!require("BiocManager", quietly = TRUE))  
  install.packages("BiocManager")  
install.packages("E:/software/R/KEGG.db_3.2.3.tar.gz", repos = NULL, type = "source")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "KEGG.db"))  

# 2. 加载需要的包  
library(clusterProfiler)  
library(org.Hs.eg.db)  
library(KEGG.db) 

if (!exists("genelist") || !"entrez2" %in% names(genelist)) {  
  stop("genelist or genelist$entrez2 does not exist")  
}  

head(genelist$entrez2)
length(genelist$entrez2)

tryCatch({  
  kegg <- enrichKEGG(  
    gene = genelist$entrez2,    # 你的 ENTREZ ID 列表  
    organism = "hsa",           # 人类是 "hsa"  
    pvalueCutoff = 0.05,        # p 值阈值  
    pAdjustMethod = "BH",       # p 值校正方法  
    use_internal_data = TRUE    # 使用内部数据  
  )  

# 3. 进行 KEGG 富集分析
kegg <- enrichKEGG(gene = genelist$entrez2,    # 你的 ENTREZ ID 列表
                   organism = "hsa",            # 人类是 "hsa"
                   pvalueCutoff = 0.05,         # p 值阈值
                   pAdjustMethod = "BH",        # p 值校正方法
                   use_internal_data = TRUE)    # 使用内部数据

# 4. 查看结果
head(kegg)

# 5. 如果需要，可以将结果保存为表格
write.csv(kegg, "kegg_enrichment_results.csv")








library(kegg.db)
library(clusterProfiler)

R.utils::setOption( "clusterProfiler.download.method",'auto' )
kegg <- enrichKEGG(genelist$entrez2, 
                   organism="hsa",   #### organism = "hsa", mmu
                   
                   pvalueCutoff = 0.05, 
                   pAdjustMethod="BH",
                   keyType="kegg",
                   use_internal_data =  T) #pvaluecutoff 是pvalue的阈值，显著富集性要<0.01

write.csv(kegg,file="up_kegg.csv",row.names = F)

# 得到kegg通路的表后直接进行排序表格的分析
dotplot(kegg, x = "GeneRatio",color = "p.adjust",size = "Count",showCategory = 15, font.size=10)


p  <- dotplot(kegg, x = "GeneRatio",color = "p.adjust",size = "Count",showCategory = 11, font.size=20, label_format = 100 )
p

p2<-p+ scale_color_continuous(low='red',high = 'green') 
p2
——————————————————————————————————————————————————————————————————————————————————————————————————————
#排序方式   首先以p.adjust<0.01或选择前20个通路，然后筛选得到的通路以counts进行排序，得到最终的前10通路
#以下正确，切莫修改


setwd("/home/admin888/project-2/rna-guo/09-25-result/Foldchage>1.2/Padj<0.05/")
x <-  read.csv("down_kegg1.csv",header =TRUE)
a <-x
a$type_order  <- factor(rev(as.integer(rownames(x))),labels=rev(x$Description))#
x <- a

x<- x[1:11,]
y <- ggplot(x, aes(x=GeneRatio, y=type_order,  size=Count, color=p.adjust)) +
  geom_point() +theme_bw() + 
  scale_colour_gradient(low="red",high="blue") +
  labs(color=expression(p.adjust),size="Gene number", x="GeneRatio",y="Pathway name",title="KEGG Pathway enrichment") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
y
y + theme(legend.position=c(0.85, 0.3),legend.title=element_text(size=14))
y + theme(axis.title.x =element_text(size=30), axis.title.y=element_text(size=30))
y + theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),axis.text=element_text(size=14,face = "bold"),axis.title.x=element_text(size=30),axis.title.y=element_text(size=30))

y + theme(axis.text=element_text(size=13),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),legend.title=element_text(size=14))+scale_size_continuous(range=c(3,8))+scale_x_continuous(limits = c(0.02, 0.045),breaks = seq(0.02, 0.04,0.01))

 #确定x轴参数，y轴参数，圆圈大小根据基因数改变，色卡的深浅依据padj+自定义色卡变换颜色+设置色卡名称+设置圆圈名称+x轴名称+y轴名称+图标题名称
y #显示气泡图，Rstudio会弹出一个新的对话框显示图片，图片可右键存为位图/图元文件等
————————————————
版权声明：本文为CSDN博主「twocanis」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
原文链接：https://blog.csdn.net/twocanis/article/details/121845958



————————————————————————————————————————————————————————————————————————————————
#kegg 柱状图绘制

ek.rt<-read.csv("up_kegg.csv",header =TRUE)

library(ggplot2)
a<-ek.rt
a$type_order  <- factor(rev(as.integer(rownames(ek.rt))),labels=rev(ek.rt$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱



ggplot(data=a, aes(x=type_order,y=LOG2padj, font.size = 100)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  geom_bar(stat="identity",fill="lightblue") + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") + 
  ylab("-LOG(padjust)") + 
  labs(title = "KEGG")+
  theme_bw()
                    
                    ————————————————
                    版权声明：本文为CSDN博主「Ryan_Yang_」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
                    原文链接：https://blog.csdn.net/Yunru_Yang/article/details/77327757


--------------------------------------------------------------------------------

library(tidyverse)


ek.rt <- read.csv("down_kegg1.csv")   #读取第1部分enrichKEGG分析输出的文件ek。
ek.rt <- separate(data=ek.rt, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
ek.rt <- separate(data=ek.rt, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
ek.rt <- mutate(ek.rt, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor
ek.rt <- mutate(ek.rt, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor





ek.rt10 <- ek.rt %>% filter(row_number() >= 1,row_number() <= 10)


library("ggplot2")

p <- ggplot(ek.rt10,aes(enrichment_factor, fct_reorder(factor(Description), enrichment_factor))) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue))) +
  scale_color_gradient(low="blue",high ="red") + 
  labs(color=expression(-log[10](p_value)),size="Count",
       x="Enrichment Factor",y="KEGG term",title="KEGG enrichment") + 
  theme_bw()
p

p <- ggplot(ek.rt10,aes(enrichment_factor,Description)) + 
  geom_point(aes(size=Count,color=p.adjust)) +
  scale_color_gradient(low="blue",high = "red") + 
  labs(color="pvalue",size="Count",
       x="Enrichment Factor",y="KEGG term",title="KEGG enrichment") + 
  theme_bw()
p

p <- ggplot(ek.rt10,aes(enrichment_factor,Description)) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue))) +
  scale_color_gradient(low="green",high ="red") + 
  labs(color=expression(-log[10](p_value)),size="Count",
       x="Enrichment Factor",y="KEGG term",title="KEGG enrichment") + 
  theme_bw()
p

p <- ggplot(ek.rt10,aes(enrichment_factor, fct_reorder(factor(Description), enrichment_factor))) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue))) +
  scale_color_gradient(low="green",high ="red") + 
  labs(color=expression(-log[10](p_value)),size="Count",
       x="Enrichment Factor",y="KEGG term",title="KEGG enrichment") + 
  theme_bw()
p


___________________________________________________________________________________________________


















library(enrichplot)#
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例

enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)#基因-通路关联网络图
enrichplot::cnetplot(kegg,circular=FALSE,colorEdge = TRUE)#circluar为指定是否环化，基因过多时建议设置为FALSE

enrichplot::heatplot(GO,showCategory = 50)#基因-通路关联热图
enrichplot::heatplot(kegg,showCategory = 50)

GO2 <- pairwise_termsim(GO)
KEGG2 <- pairwise_termsim(kegg)
enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk")#通路间关联网络图
enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk")


########################################################################################
library("clusterProfiler")
#读取数据
TGF-beta-signaling-pathway
Wnt-signaling-pathway



go_ythdf2 <- read.csv("Wnt-signaling-pathway")
go_ythdf2 <- t(go_ythdf2)

#ID转换


ensembl <- bitr(go_ythdf2,fromType = "ENTREZID",toType = "ENSEMBL",
                            OrgDb = "org.Mm.eg.db",drop = T)


ensembl$name  <- bitr(go_ythdf2,fromType = "ENTREZID",toType = "SYMBOL",
                      OrgDb = "org.Mm.eg.db",drop = T)

go_ythdf2_id_trance <- bitr(go_ythdf2,fromType = "ENTREZID",toType = "SYMBOL",
                            OrgDb = "org.Mm.eg.db",drop = T)


write.csv(ensembl,file="all-Wnt-signaling-pathway.csv",row.names = F)

#GO分析
result_go_ythdf2 <- enrichGO(go_ythdf2_id_trance$ENTREZID,OrgDb = "org.Hs.eg.db",
                             keyType = "ENTREZID",ont = "ALL",readable = T)

#绘图
dotplot(result_go_ythdf2,showCategory=50)
barplot(result_go_ythdf2,showCategory=20)

#保存
df_go_ythdf2 <- as.data.frame(result_go_ythdf2)
write.table(df_go_ythdf2,"df_go_ythdf2.xls",quote = F,sep = "\t")

################################################################################################
#安装加载VennDiagram包；
install.packages("grid")
install.packages("ggVennDiagram")
#加载VennDiagram包；
install.packages("VennDiagram")

install.packages("openxlsx")

library(VennDiagram)
library(ggVennDiagram)
library (VennDiagram)  

library(openxlsx)
#读入数据；
WNT <- read.csv("Wnt-signaling-pathway.csv")
TGF <- read.csv("TGF-beta-signaling-pathway.csv")

WNT1 <- as.vector(unlist(WNT[2]))
TGF1 <- as.vector(unlist(TGF[2]))

mydata<-list(WNT1,TGF1)
library(gplots)
venn(mydata, names("WNt", "TGF"))


BiocManager::install("grid")
BiocManager::install("futile.logger")

library(grid)
library(futile.logger)
library(VennDiagram)

dat <- read.csv("Diff_Groups.csv", header = T)

G8671_list <- dat$G8671_Groups[1:482]
G37364_list <- dat$G37364_Groups[1:318]
TCGA_list <- dat$TCGA_Groups[1:2984]

#韦恩图（VennDiagram 包，适用样本数 2-5）
library(VennDiagram)

#指定统计的分组列，并设置作图颜色、字体样式等
venn_list <- list(wntpath = WNT1, tgfpathway = TGF1)
venn.diagram(venn_list, filename = 'venn.png', imagetype = 'png', 
             fill = c('red', 'blue'), alpha = 0.50, 
             cat.col = c('red', 'blue'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('red', 'blue'), cex = 1.5, fontfamily = 'serif')

inter <- get.venn.partitions(venn_list)

write.csv(inter,file="venn",row.names = F)


help(venn)

