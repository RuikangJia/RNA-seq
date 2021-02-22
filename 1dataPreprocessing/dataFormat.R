# 功能：统一TCGA的.txt数据文件格式，包括临床和表达值
    # 主要内容：转置；行列名调整
# 输入：
    # *.txt
# 输出：
    # *.Rdata:数据框
    # *.csv:excel文件

library(modifyTools)
# 0函数
dataFormat = function(data){
    # 转置
    data = data.frame(t(data),check.names = F)
    # 修改样本名称
    rownames(data) = toupper(rownames(data))
    return(data)
}


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()[grep(".txt$",dir())]
dataList = lapply(fileNames,read.table,sep = "\t",header = T,row.names = 1,check.names = F)


# 3整理
dataFormatedList = lapply(dataList,dataFormat)


# 4保存
dataNames = gsub(".txt","",fileNames)
for(i in 1:length(dataNames)){
    dfSave(dataFormatedList[[i]],dataNames[i],outputDir)
}
