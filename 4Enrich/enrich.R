# 功能：基因列表的功能富集
    # 主要内容：kegg、go
# 输入：
    # CA_HubGene.csv：基因列表
# 输出：


# 0工具包
rm(list = ls())
library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)
library(clusterProfiler)
library(modifyTools)
library(ggplot2)
library(stringr)
library(writexl)
# 0函数
enrichPlot = function(outputDir,dataName,fileName,enrichData){
    png(filename = paste(outputDir,dataName,"_",fileName,".png",sep = ""),
        width = 800,height = 600)
    print(
        dotplot(enrichData,font.size = 20) +
            theme(legend.text = element_text(size = 20),
                  legend.title = element_text(size = 20),
                  plot.title = element_text(hjust=0.5),
                  title = element_text(size = 30)) +
            scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
            labs(title = paste(fileName,"Enrichment Analysis"))
        )
    dev.off()
}


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileName = dir()
dataFileName = fileName[grep(".csv$",fileName)]
print("Input dataName：")
dataName = readline()
# 数据导入
data = read.csv(dataFileName,check.names = F,row.names = 1)


# 3处理
# 提取基因列表
geneList = data$GeneList
# 转换ID
geneListBitr = bitr(geneList,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
geneList = geneListBitr$ENTREZID
# KEGG、GO富集
kegg = enrichKEGG(geneList)
goBP = enrichGO(geneList,OrgDb = org.Hs.eg.db, ont = "BP")
goCC = enrichGO(geneList,OrgDb = org.Hs.eg.db, ont = "CC")
goMF = enrichGO(geneList,OrgDb = org.Hs.eg.db, ont = "MF")
# 绘图
try(enrichPlot(outputDir,dataName,"KEGG",kegg))
try(enrichPlot(outputDir,dataName,"BP",goBP))
try(enrichPlot(outputDir,dataName,"CC",goCC))
try(enrichPlot(outputDir,dataName,"MF",goMF))


# 4保存
try(dfSave(kegg,paste(dataName,"KEGG",sep = "_"),outputDir))
try(dfSave(goBP,paste(dataName,"BP",sep = "_"),outputDir))
try(dfSave(goCC,paste(dataName,"CC",sep = "_"),outputDir))
try(dfSave(goMF,paste(dataName,"MF",sep = "_"),outputDir))

