# 功能：
    # 1.ROC数据格式的整理：基因在列，样本在行，末尾样本类型列
# 主要问题：
    # 1.名称提取：stringr::str_extract
# 输入：
    # HB_AG_BD_TPM.csv:各种类型样本表达数据
    # HB_CA_TPM.csv:各种类型样本表达数据
    # HB_HC_TPM.csv:各种类型样本表达数据
    # HB_NAG_BD_TPM.csv:各种类型样本表达数据
    # hubGeneAll.csv：目标基因信息
# 输出：
    # ROCData.RData：完整ROC格式数据
    # ROCData.csv
    # ROCData.xlsx
    # hubGeneROCData.RData：目标基因ROC格式数据
    # hubGeneROCData.csv
    # hubGeneROCData.xlsx


# 0工具包
rm(list = ls())
gc()
library(modifyTools)
library(writexl)
# 0函数


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()
exprFileNames = fileNames[grepl("TPM.csv$",fileNames)]
geneFileNames = fileNames[grepl("HubGene",fileNames)]
dataName = stringr::str_extract(exprFileNames,"HC|NAG|AG|CA")
# 数据导入
exprList = lapply(exprFileNames,read.csv,check.names= F,row.names = 1)
hubGene = read.csv(geneFileNames,check.names = F,row.names = 1)


# 3处理
# 类型列
type = rep(dataName,sapply(exprList,ncol))
# 数据集整合
ROCData = data.frame(t(Reduce(cbind,exprList)),check.names = F)
# 目标基因集
hubGeneROCData = ROCData[,names(ROCData) %in% hubGene$GeneList]
ROCData$Type = type
hubGeneROCData$Type = type


# 4保存
dfSave(ROCData,"ROCData",outputDir)
dfSave(hubGeneROCData,"hubGeneROCData",outputDir)