# 功能：对临床数据进行整理
    # 主要内容：分期缺失样本去除；划分不同分期样本
# 输入：
    # clin.csv：格式整理好的临床数据文件
    # RSEM.csv：格式整理好的基因表达数据文件
# 输出：
    # 阶段信息缺失
        # dataClinStageNA.RData：
        # dataClinStageNARemoved.RData
        # dataClinStageNA.csv
        # dataClinStageNARemoved.csv
    # RData格式临床分期数据
        # dataClinStageI.RData    
        # dataClinStageIV.RData  
        # dataClinStageII.RData   
        # dataClinStageIII.RData  
    # csv格式临床分期数据
        # dataClinStageI.csv
        # dataClinStageII.csv     
        # dataClinStageIII.csv   
        # dataClinStageIV.csv    
    # 统计信息
        # stageStasticInfo.csv
        # stageStasticInfo.RData
# 工具包
library(stringr)
library(modifyTools)

# 函数
sampleDivideNC = function(dataExpr){
    indexC = grepl("01$",rownames(dataExpr))
    dataExprC = dataExpr[indexC,]
    dataExprN = dataExpr[!indexC,]
    # dfSave(dataExprC,paste(dataName,"C",sep = ""),outputDir)
    # dfSave(dataExprN,paste(dataName,"N",sep = ""),outputDir)
    return(list(dataExprC,dataExprN))
}
attributeNaRemove = function(data,attribute){
    # 确定属性名称对应的位置
    indexAttribute = match(attribute,names(data))
    # 确定阶段信息缺失样本的索引:逻辑
    indexNA = is.na(data[,indexAttribute])
    return(indexNA)
}
stageDivide = function(dataClin,dataExpr,exprName){
    # 分期信息向量
    stageInfo = dataClin$pathologic_stage
    # 临床不同阶段样本索引
    indexStageI = (str_count(stageInfo,"i") == 1 & str_count(stageInfo,"v") == 0)
    indexStageII = (str_count(stageInfo,"i") == 2)
    indexStageIII = (str_count(stageInfo,"i") == 3)
    indexStageIV = (str_count(stageInfo,"i") == 1 & str_count(stageInfo,"v") == 1)
    # 临床信息不同阶段划分
    dataClinStageI = dataClin[indexStageI,]
    dataClinStageII = dataClin[indexStageII,]
    dataClinStageIII = dataClin[indexStageIII,]
    dataClinStageIV = dataClin[indexStageIV,]
    # 保存数据
    dfSave(dataClinStageI,"ClinStageI",outputDir)
    dfSave(dataClinStageII,"ClinStageII",outputDir)
    dfSave(dataClinStageIII,"ClinStageIII",outputDir)
    dfSave(dataClinStageIV,"ClinStageIV",outputDir)
    # 样本分期统计信息
    stageStasticInfo = data.frame(count = c(sum(indexStageI),sum(indexStageII),sum(indexStageIII),sum(indexStageIV)))
    rownames(stageStasticInfo) = c("stageI","stageII","stageIII","stageIV")
    # 临床不同阶段样本索引
    indexStageI = match(rownames(dataClinStageI),substr(rownames(dataExpr),1,12))
    indexStageII = match(rownames(dataClinStageII),substr(rownames(dataExpr),1,12))
    indexStageIII = match(rownames(dataClinStageIII),substr(rownames(dataExpr),1,12))
    indexStageIV = match(rownames(dataClinStageIV),substr(rownames(dataExpr),1,12))
    # 表达数据不同阶段划分
    dataExprStageI = dataExpr[indexStageI,]
    dataExprStageII = dataExpr[indexStageII,]
    dataExprStageIII = dataExpr[indexStageIII,]
    dataExprStageIV = dataExpr[indexStageIV,]
    # 数据保存
    dfSave(dataExprStageI,paste(exprName,"StageI",sep = ""),outputDir)
    dfSave(dataExprStageII,paste(exprName,"StageII",sep = ""),outputDir)
    dfSave(dataExprStageIII,paste(exprName,"StageIII",sep = ""),outputDir)
    dfSave(dataExprStageIV,paste(exprName,"StageIV",sep = ""),outputDir)
    return(stageStasticInfo)
}

# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()[grep(".csv$",dir())]
fileNameClin = fileNames[grepl("clin.csv$",fileNames)]
fileNameExpr = fileNames[!grepl("clin.csv$",fileNames)]
dataNameClin = gsub(".csv$","",fileNameClin)
dataNameExpr = gsub(".csv$","",fileNameExpr)
# 数据导入
dataClin = read.csv(fileNameClin,row.names = 1,check.names = F)
dataExpr = read.csv(fileNameExpr,row.names = 1,check.names = F)


# 3处理
# 表达：样本划分为癌旁和癌
dataExprCN = sampleDivideNC(dataExpr)
dataExprC = dataExprCN[[1]]
dataExprN = dataExprCN[[2]]
# 临床：去除阶段信息缺失的样本
indexNA = attributeNaRemove(dataClin,"pathologic_stage")
dataClinStageNARemoved = dataClin[!indexNA,]
dataClinStageNA = dataClin[indexNA,]
# 表达、临床：取临床和癌表达样本的交集
#   交集
sampleNameInter = intersect(substr(rownames(dataExprC),1,12),rownames(dataClinStageNARemoved))
#   临床数据
indexInterClin = rownames(dataClinStageNARemoved) %in% sampleNameInter
dataClinInter = dataClinStageNARemoved[indexInterClin,]
dataClinInterRemoved = dataClinStageNARemoved[!indexInterClin,]
#   表达数据
indexInterExpr = substr(rownames(dataExprC),1,12) %in% sampleNameInter
dataExprInter = dataExprC[indexInterExpr,]
dataExprInterRemoved = dataExprC[!indexInterExpr,]
# 表达：划分四个阶段
stageStasticInfo = stageDivide(dataClinInter,dataExprInter,dataNameExpr)


# 4保存
dfSave(dataExprC,paste(dataNameExpr,"StageC",sep = ""),outputDir)
dfSave(dataExprN,paste(dataNameExpr,"StageN",sep = ""),outputDir)
dfSave(dataClinStageNARemoved,paste(dataNameClin,"StageNARemoved",sep = ""),outputDir)
dfSave(dataClinStageNA,paste(dataNameClin,"StageNA",sep = ""),outputDir)
dfSave(dataClinInter,paste(dataNameClin,"Inter",sep = ""),outputDir)
dfSave(dataClinInterRemoved,paste(dataNameClin,"InterRemoved",sep = ""),outputDir)
dfSave(dataExprInter,paste(dataNameExpr,"Inter",sep = ""),outputDir)
dfSave(dataExprInterRemoved,paste(dataNameExpr,"InterRemoved",sep = ""),outputDir)
