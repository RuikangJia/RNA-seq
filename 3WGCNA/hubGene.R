# 功能：
    # 1.进行核心基因的筛选：MM与变化基因的交集。
# 主要问题：
# 输入：
    # **_GeneAffectedByDrug：表达发生变化的基因。
    # **_MMSelected：MM标准得到的基因。
# 输出：
    # **_HubGene.csv、.xlsx、.RData：核心基因列表数据。
    # **_HubGene.png：结合运算VN图。


# 0工具包
rm(list = ls())
gc()
library(modifyTools)
library(writexl)
library(VennDiagram)
# 0函数


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()[grepl(".csv$",dir())]
MMFileName = fileNames[grepl("MM",fileNames)]
geneListFileName = fileNames[!grepl("MM",fileNames)]
print("Input dataName:")
dataName = readline()
# 数据导入
MM = read.csv(MMFileName,check.names = F,row.names = 1)
geneList = read.csv(geneListFileName,check.names = F,row.names = 1)


# 3处理
# MM基因
MMGene = rownames(MM)
# 差异基因(或药物影响基因)
DEGene = rownames(geneList)
# VN图展示S
venn.diagram(list("MM" = MMGene,"GeneList" = DEGene),
             fill=c("red","green"),
             paste(outputDir,dataName,"_HubGene.png",sep = ""),
             lty = 1,lwd = 1,col = "black",
             alpha = 0.5,cat.col = "black",
             cat.cex=1,cat.fontface = "bold",margin = 0.1,cex = 1,
             main.cex = 2,main.fontface = "bold",main.pos = c(0.5,0.95),
             main = paste("Hub Gene Screening（",dataName,"）")
)
# 核心基因列表
hubGene = data.frame( GeneList = intersect(MMGene,DEGene),check.names = F)
rownames(hubGene) = hubGene$GeneList

# 4保存
dfSave(hubGene,paste(dataName,"HubGene",sep = "_"),outputDir)