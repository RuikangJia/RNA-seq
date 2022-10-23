# 功能：
    # 1.根据样本类别对表达数据进行划分
# 主要问题：
    # 1.划分不同类型样本的表达数据
    # 2.记录类别信息
# 输入：
    # rawCountFiltered.csv：所有类型样本的基因表达数据
    # TPM.csv：所有类型样本的基因表达数据
# 输出：
    # sampleTypeInfo.csv：样本类别信息
    # H*_**.csv：各类样本的表达数据


# 0工具包
rm(list = ls())
library(modifyTools)
library(writexl)
# 0函数
# 根据列名筛选列
pickColByColName = function(data,colPattern){
    colIndex = grepl(colPattern,colnames(data))
    res = data.frame(data[,colIndex])
    names(res) = names(data)[colIndex]
    rownames(res) = rownames(data)
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
dataName = gsub(".csv$","",dataFileName)
dataName = ifelse(grepl("TPM",dataName),"TPM","RC")
# 数据导入
dataInfo = read.csv(dataFileName,check.names = F,row.names = 1)
data = dataInfo


# 3处理
# 类别信息
typeInfo = c("^HB_HC_",
              "^HB_NAG_BD_","^HB_NAG_AD_","^HB_AG_BD_","^HB_AG_AD_",
              "^HB_CA_","^HT_CA_.+_01$","^HT_CA_.+_11$")
names(typeInfo) = c("HB_HC",
                     "HB_NAG_BD","HB_NAG_AD","HB_AG_BD","HB_AG_AD",
                     "HB_CA","HT_CA_01","HT_CA_11")
# 划分过程
sampleID = c()
type = c()
for(i in 1:length(typeInfo)){
    # 挑选列
    res = pickColByColName(data,typeInfo[i])
    if(dim(res)[2] != 0){
        sampleID = c(sampleID,names(res))
        type = c(type,rep(names(typeInfo)[i],ncol(res)))
        dfSave(res,paste(names(typeInfo)[i],dataName,sep = "_"),outputDir)
    }
    # 记录样本类别信息
}
# 整理记录信息
sampleTypeInfo = data.frame(SampleID = sampleID,Type = type)


# 4保存
dfSave(sampleTypeInfo,"sampleTypeInfo",outputDir)
