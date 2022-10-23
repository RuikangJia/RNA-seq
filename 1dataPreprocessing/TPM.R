# 功能：
#   1.将rawCount数据标准化为TPM。 
# 主要问题：
# 输入：
#   rawCountFiltered.csv：待TPM转换的基因表达数据。
#   geneInfo.csv：基因的描述信息，含有基因长度。
# 输出：
#   TPM.csv：TPM转换后的基因表达数据。
#   TPM.xlsx：TPM转换后的基因表达数据。
#   TPM.RData：TPM转换后的基因表达数据。


# 0工具包
rm(list = ls())
library(modifyTools)
library(writexl)
# 0函数
TPM = function(rawCounts,geneLength){
    res = as.data.frame(
        apply(rawCounts / (geneLength/1000),2,function(x) (x / sum(x) * 10^6))
        )
    res = as.data.frame(res)
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
dataFileName = fileNames[grepl("^rawCount",fileNames)]
geneInfoFileName = fileNames[grepl("Info.csv$",fileNames)]
dataName = gsub(".csv$","",dataFileName)
# 数据导入
dataInfo = read.csv(dataFileName,check.names = F,row.names = 1)
geneInfo = read.csv(geneInfoFileName,check.names = F,row.names = 1)


# 3处理
# 获取基因长度
geneLength = geneInfo$gene_length[match(rownames(dataInfo),rownames(geneInfo))]


# TPM标准化
dataTPM = TPM(dataInfo,geneLength)


# 4保存
dfSave(dataTPM,"TPM",outputDir)
