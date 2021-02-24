# 功能：核心基因的筛选
    # 主要内容：差异基因、显著模块、GSMM
# 输入：
    # StageImoduleCorSig.csv：显著相关的模块
    # stageIMMGS.csv：MMGS信息
    # uniqueGeneStageI.csv：差异信息
# 输出：
    # stageIhubGeneSigMMGS.RData
    # stageIhubGeneSigMMGS.csv
    # stageIhubGeneSigMMGSDE.RData
    # stageIhubGeneSigMMGSDE.csv

# 工具包
library(modifyTools)
# 函数


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()
# 数据名称
fileNameMMGS = fileNames[grep("MMGS.csv$",fileNames)]
dataName = gsub("MMGS.+$","",fileNameMMGS)
# 数据导入
MMGS = read.csv(fileNameMMGS,check.names = F,row.names = 1)
geneDE = read.csv(fileNames[grep("^unique",fileNames)],check.names = F,row.names = 1)
modSig = read.csv(fileNames[grep("CorSig.csv$",fileNames)],check.names = F,row.names = 1)
modNameSig = rownames(modSig)


# 3处理
# 模块筛选
hubGeneSig = subset(MMGS,modName %in% modNameSig)
# GS、MM筛选
hubGeneSigMM = subset(hubGeneSig,abs(MM) > 0.8 & MMPvalue < 0.05)
hubGeneSigMMGS = subset(hubGeneSigMM, (abs(GST) > 0.2 & GSTPvalue < 0.05) | (abs(GSN) > 0.2 & GSNPvalue < 0.05))
# 差异基因筛选
hubGeneSigMMGSDE = subset(hubGeneSigMMGS, rownames(hubGeneSigMMGS) %in% getGeneID(geneDE$x))


# 4保存
dfSave(hubGeneSigMMGS,paste(dataName,"hubGeneSigMMGS",sep = ""),outputDir)
dfSave(hubGeneSigMMGSDE,paste(dataName,"hubGeneSigMMGSDE",sep = ""),outputDir)
