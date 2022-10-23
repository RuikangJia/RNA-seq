# 功能：对总表数据进行整理
    # 1.筛选mRNA对应基因（行）。
    # 2.求和方式，合并对应多个探针的基因（行）。
    # 3.去除含有NA的基因（行）。
    # 4.分割基因表达信息和基因描述信息（列）。
# 主要问题：
# 输入：
    # rawCountAndGeneInfoH[B|T].csv：合成的总表
# 输出：
    # geneInfo.csv：基因描述信息
    # rawCount.csv：基因表达信息


# 0工具包
rm(list = ls())
library(modifyTools)
library(writexl)
# 0函数
# 根据列值筛选行
pickRowByColValue = function(data,colName,colValue){
    res = subset(data,data[,names(data) == colName] == colValue)
    return(res)
}
# 根据列名筛选列
pickColByColName = function(data,colPattern){
    colIndex = grepl(colPattern,colnames(data))
    res = data.frame(data[,colIndex])
    names(res) = names(data)[colIndex]
    rownames(res) = rownames(data)
    return(res)
}
# 合并重复行
mergeDupRow = function(data,colNameWithDup,fun){
    colIndex = (names(data) == colNameWithDup)
    res = aggregate(data[,!colIndex],by = list(data[,colIndex]),FUN = fun)
    names(res)[1] = names(data)[colIndex]
    names(res)[2:length(names(res))] = names(data)[!colIndex]
    return(res)
}


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()
dataFileName = fileNames[grepl(".csv$",fileNames)]
# 数据导入
dataInfo = read.csv(dataFileName,check.names = F,row.names = 1)
data = dataInfo


# 3处理
# 表达信息处理
#   提取mRNA对应基因
dataMRNA = pickRowByColValue(data,"gene_biotype","protein_coding")
#   去除NA基因
dataMRNADeNA  = dataMRNA[complete.cases(dataMRNA),]
#   提取表达信息列和gene_name列
exprInfo = pickColByColName(dataMRNADeNA,"^H|gene_name")
#   合并重复基因：表达值相加
exprInfoMergeDupRow = mergeDupRow(exprInfo,colNameWithDup = "gene_name",fun = sum)
#   将gene_name设置为行名
rownames(exprInfoMergeDupRow) = exprInfoMergeDupRow$gene_name
exprInfoMergeDupRow = exprInfoMergeDupRow[,-1]
# 描述信息整理
#   提取描述信息列
geneInfo = pickColByColName(dataMRNADeNA,"[gene|tf]")
#   基因信息直接去重复行
geneInfoMergeDupRow = geneInfo[!duplicated(geneInfo$gene_name),]
rownames(geneInfoMergeDupRow) = geneInfoMergeDupRow$gene_name
geneInfoMergeDupRow = geneInfoMergeDupRow[,-1]


# 4保存
dfSave(geneInfoMergeDupRow,"geneInfo",outputDir)
dfSave(exprInfoMergeDupRow,"rawCount",outputDir)
