# 功能：
    # 1.合并诺禾返回的多个.csv格式的表达数据。
    # 2.记录样本批次。
# 主要问题：
    # 1.每张表的基因不完全一致，通过inner_join解决。
# 输入：
    # rawCountAndGeneInfo*.csv：多个需要合并的数据表
# 输出：
    # rawCountAndGeneInfo.csv：合并后的数据表
    # sampleBatch.csv：样本批次信息


# 0工具包
rm(list = ls())
gc()
library(modifyTools)
library(writexl)
library(dplyr)
library(stringr)
# 0函数



# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()
dataFileNames = fileNames[grep('.csv$',fileNames)]
# 数据导入
dataList = lapply(dataFileNames,read.csv,check.names = F,row.names = 1)


# 3处理
# 批次信息：Batch和文件名对应
batch = str_extract(dataFileNames,"\\d+")
# 拼接：合并表的同时记录批次信息
dataMerged = dataList[[1]]
batchInfo = data.frame(SampleID = names(dataList[[1]]),Batch = batch[1],check.names = F)
for(i in 2:length(dataList)){
    dataMerged = inner_join(dataMerged,dataList[[i]])
    batchInfo = rbind(batchInfo,data.frame(SampleID = names(dataList[[i]]),Batch = batch[i],check.names = F))
}
# 样本信息：去冗余
batchInfo = subset(batchInfo,grepl("^H",SampleID))


# 4保存
dfSave(batchInfo,"sampleBatchInfo",outputDir)
dfSave(dataMerged,"rawCountAndGeneInfo",outputDir)