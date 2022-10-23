# 功能：热图展示
    # 1.


# 0工具包
rm(list = ls())
gc()
library(modifyTools)
library(writexl)
library(pheatmap)
library(ggplot2)
# 0函数


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()
dataFileName = fileNames[grepl(".csv$",fileNames)]
dataName = gsub(".csv$","",dataFileName)
# 数据导入

# 3处理
# 选取对应基因和样本
# 转换
expr = log(expr +1)
# 热图
# 分组注释信息+颜色
annotation_col = data.frame(Sample = factor(temp,levels = c("HC","NAG","AG","CA")),
                            row.names = names(exprTop20))
ann_colors = list(Sample = c(HC = colorRampPalette(c("white","firebrick"))(4)[1],
                             NAG = colorRampPalette(c("white","firebrick"))(4)[2],
                             AG = colorRampPalette(c("white","firebrick"))(4)[3],
                             CA = colorRampPalette(c("white","firebrick"))(4)[4]))
# 主函数
pheatmap(exprTop20,
         # 是否聚类
         cluster_rows = T,cluster_cols = F,
         # 标准化方式
         scale = "row",
         # 刻度调整
         breaks = seq(-2,2,4/100),
         # 添加间隔：没有聚类使用gaps指定间隔位置，聚类使用cutree指定类别层次。
         gaps_col = c(n1,n1+n2,n1+n2+n3),gaps_row = c(1,5,7),
         cutree_row = c(1),cutree_col = c(2),
         # 行列注释信息
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         # 是否显示行名或列名
         show_colnames = F,
         # 大小
         width = 10,height = 5,fontsize = 10,
         main = "Expression of Top20 genes of CA in samples of different stages",
         filename = paste(outputDir,"CA","Heatmap.png",sep = "")
         )

# 4保存
