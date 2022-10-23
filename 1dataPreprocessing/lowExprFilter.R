# 功能：
    # 1.去除整体表达量过低的基因。
# 主要问题：
    # 1.CPM的计算
# 输入：
    # rawCountBEAdj.csv：待过滤的基因表达数据
# 输出：
    # rawCountFiltered.RData:过滤后的基因表达数据
    # rawCountFiltered.csv:过滤后的基因表达数据
    # rawCountFiltered.xlsx:过滤后的基因表达数据


# 0工具包
rm(list = ls())
library(modifyTools)
library(writexl)
# 0函数
CPMFilter = function(data,sumCount){
    num = length(data)
    temp = (data / sumCount )* 1000000
    return(sum(temp >= 5) >= (num/10))
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
dataName = gsub(".csv$","",dataFileName)
# 数据导入
dataInfo = read.csv(dataFileName,check.names = F,row.names = 1)


# 3处理
# 过滤低表达基因：CPM
smapleReadSum = apply(dataInfo,2,sum)
index = apply(dataInfo,1,CPMFilter,smapleReadSum)
data = dataInfo[index,]


# 4保存
dfSave(data,"rawCountFiltered",outputDir)