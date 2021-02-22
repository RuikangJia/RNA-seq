# 功能：获取两种方式差异基因的交集
    # 主要内容：癌旁样本和癌症阶段样本差异基因；癌症阶段样本和癌症其它阶段样本差异基因
# 输入：
    # I_CElseDEGSorted.csv：癌症阶段样本和癌症其它阶段样本差异基因
    # N_IDEGSorted.csv：癌旁样本和癌症阶段样本差异基因
# 输出：
    # I_CElseInter.RData：交集差异基因，使用的癌症阶段样本和癌症其它阶段样本差异基因数据
    # I_CElseInter.csv
# 工具包
library(VennDiagram)
library(modifyTools)

# 函数


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()[grep(".csv",dir())]
fileNameDENC = fileNames[!grepl("Else",fileNames)]
fileNameDECElse = fileNames[grepl("Else",fileNames)]
# 数据导入
DENC = read.csv(fileNameDENC,row.names = 1,check.names = F)
DECElse = read.csv(fileNameDECElse,row.names = 1,check.names = F)


# 3处理
DEInter = DECElse[rownames(DECElse) %in% rownames(DENC) ,]
venn.diagram(list(DENC = rownames(DENC),DECElse = rownames(DECElse)),
             filename = paste(outputDir,"Inter.png",sep = ""),
             fill = c("red","blue"),height = 2000,width = 2000,
             cat.dist = c(0.02, 0.02),cex = 1.45,cat.cex = 1.45,
             cat.pos = 0,ext.text = 2)

# 4保存
dfSave(DEInter,gsub("DEGSorted.csv","Inter",fileNameDECElse),outputDir)
