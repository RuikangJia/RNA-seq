# 功能：计算GS和MM
    # 主要内容：GST、GSN、MM计算
# 输入：
    # ClinStageIFiltered.csv：临床数据
    # RSEMStageIFiltered2.csv：表达数据
    # stageIModuleSetList.RData：基因模块列表
    # stageINetwork.RData：网络信息
# 输出：
    # stageIGSMMFiltered.RData:三种标准过滤


# 工具包
library(modifyTools)
library(WGCNA)
# 函数
MMCompute = function(dataExpr,moduleSetList,mergedMEs){
    MM = data.frame()
    for(i in 1:length(moduleSetList)){
        modName = names(moduleSetList)[i]
        gene = moduleSetList[[i]]
        exprMod = dataExpr[,match(gene,colnames(dataExpr))]
        MEsMod = mergedMEs[,match(modName,colnames(mergedMEs))]
        MMMod = as.data.frame(cor(exprMod, MEsMod, use = "p"))
        MMPvalueMod = as.data.frame(corPvalueStudent(as.matrix(MMMod), dim(exprMod)[1]))
        modNameDF = data.frame(modeName = rep(modName,dim(MMMod)[1]))
        MMTemp = cbind(modNameDF,MMMod,MMPvalueMod)
        MM = rbind(MM,MMTemp)
    }
    names(MM) = c("MM","MMPvalue")
    return(MM)
}
factorizationT = function(data){
    data[data == "t1"] = 1
    data[data == "t1a"] = 1
    data[data == "t1b"] = 2
    data[data == "t2"] = 3
    data[data == "t2a"] = 3
    data[data == "t2b"] = 4
    data[data == "t3"] = 5
    data[data == "t4"] = 5
    data[data == "t4a"] = 6
    data[data == "t4b"] = 7
    return(data)
}
factorizationN = function(data){
    data[data == "n0"] = 1
    data[data == "n1"] = 2
    data[data == "n2"] = 3
    data[data == "n3"] = 4
    data[data == "n3a"] = 4
    data[data == "n3b"] = 5
    data[data == "nx"] = 2
    return(data)
}
GSCompute = function(dataExpr,moduleSetList,dataClin){
    GS = data.frame()
    for(i in 1:length(moduleSetList)){
        modName = names(moduleSetList)[i]
        gene = moduleSetList[[i]]
        exprMod = dataExpr[,match(gene,colnames(dataExpr))]
        GSMod = as.data.frame(cor(exprMod, dataClin, use = "p"))
        GSPvalueMod = as.data.frame(corPvalueStudent(as.matrix(GSMod),dim(exprMod)[1]))
        GSTemp = cbind(GSMod,GSPvalueMod)
        GS = rbind(GS,GSTemp)
    }
    names(GS) = c("GS","GSPvalue")
    return(GS)
}

# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()
dataName = fileNames[grep("Network.RData$",fileNames)]
dataName = gsub("Network.RData$","",dataName)
# 数据导入
dataExpr = read.csv(fileNames[grep("^RSEM.+csv$",fileNames)],check.names = F,row.names = 1)
dataClin = read.csv(fileNames[grep("^Clin.+csv$",fileNames)],check.names = F,row.names = 1)
load(fileNames[grep("Network.RData$",fileNames)])
load(fileNames[grep("SetList.RData$",fileNames)])


# 3处理
# 数据格式整理
colnames(dataExpr) = getGeneID(colnames(dataExpr))
dataClinT = data.frame(T = factorizationT(dataClin$pathology_T_stage))
rownames(dataClinT) = rownames(dataClin)
dataClinN = data.frame(N = factorizationN(dataClin$pathology_N_stage))
rownames(dataClinN) = rownames(dataClin)
# 计算
MM = MMCompute(dataExpr,moduleSetList,mergedMEs)
# GS条件筛选
GST = GSCompute(dataExpr,moduleSetList,dataClinT)
GSN = GSCompute(dataExpr,moduleSetList,dataClinN)
# 整合GS、MM计算结果
MMGS = cbind(MM,GST,GSN)
names(MMGS) = c("modName","MM","MMPvalue","GST","GSTPvalue","GSN","GSNPvalue")


# 4保存
dfSave(MMGS,paste(dataName,"MMGS",sep = ""),outputDir)
