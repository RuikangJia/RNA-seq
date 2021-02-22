# 功能：分析模块与分期性状之间的相关性;模块之间的相关性
# 输入：
    # ClinStageIVFiltered.csv：过滤后的临床数据
    # RSEMStageIVFiltered2.csv：过滤后的表达数据
    # stageIVNetwork.RData：共表达网络信息
# 输出：
    # StageIVMETraitsCor.png：相关性图
    # StageIVmoduleCorSig.RData：显著相关模块信息
    # StageIVmoduleCorSig.csv
# 工具包
library(WGCNA)
library(pheatmap)
library(modifyTools)
# 函数
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
factorizationDF = function(df){
    df[,1] = as.numeric(factorizationT(df[,1]))
    df[,2] = as.numeric(factorizationN(df[,2]))
    return(df)
}
METraitsCorPlot = function(corMETrait,pvalueCorMETrait,moduleCorSig,dataName,outputDir){
    png(filename = paste(outputDir,dataName,"METraitsCor.png",sep = ""),
        width = 800,height = 1500)
    par(mar = c(4, 7, 4, 4))
    # 标签文本
    textPvalue = paste(signif(corMETrait, 2), " (", signif(pvalueCorMETrait, 1), ")", sep = "")
    # 模块标签颜色
    color = rep("black",dim(corMETrait)[1])
    color[match(rownames(moduleCorSig),rownames(corMETrait))] = "red"
    # 相关性图
    labeledHeatmap(Matrix = corMETrait, 
                   xLabels = names(DFTrait), yLabels = names(mergedMEs), 
                   ySymbols = names(mergedMEs), 
                   colorLabels = T,colors.lab.y = color, 
                   colors = blueWhiteRed(50), 
                   textMatrix = textPvalue,
                   setStdMargins = FALSE, 
                   cex.text = 2, cex.lab = 2,cex.main = 2.5,
                   zlim = c(-1,1), xLabelsAngle = 0,xLabelsAdj = 0.5,
                   main = paste("Module-trait relationships")) 
    dev.off()
}
moduleSigCorToTrait = function(corMETrait,pvalueCorMETrait){
    # 数据整合
    cor = data.frame(corMETrait[,1],pvalueCorMETrait[,1],
                     corMETrait[,2],pvalueCorMETrait[,2],check.names = F)
    names(cor) = c("T","pvalueT","N","pvalueN")
    # 显著性筛选
    moduleCorSig = cor[cor$pvalueT < 0.05 | cor$pvalueN < 0.05,]
    return(moduleCorSig)
}
moduleCorPlot = function(MEs,outputDir,dataName){
    png(file = paste(outputDir,dataName,"MEs.png"),width = 1000,height = 800)
    pheatmap(cor(MEs),show_colnames = F,cex = 1.2)
    dev.off()
}

# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileName = dir()[grep("^ClinS",dir())]
dataName = gsub("Filtered.csv$","",fileName)
dataName = gsub("^Clin","",dataName)
# 数据导入
load(dir()[grep("Network.RData$",dir())])
clin = read.csv(fileName,row.names = 1,check.names = F)


# 3处理
# 绘图：模块相关性图
moduleCorPlot(mergedMEs,outputDir,dataName)
# 临床性状数据因子化
DFTrait = clin[,c(8:9)]
DFTraitToNum = factorizationDF(DFTrait)
# 计算模块与性状相关性、p值
corMETrait = cor(mergedMEs,DFTraitToNum,use = "p")
pvalueCorMETrait = corPvalueStudent(corMETrait,dim(clin)[1])
# 显著相关模块
moduleCorSig = moduleSigCorToTrait(corMETrait,pvalueCorMETrait)
# 绘图：模块性状相关性图
METraitsCorPlot(corMETrait,pvalueCorMETrait,moduleCorSig,dataName,outputDir)



# 4保存
dfSave(moduleCorSig,paste(dataName,"moduleCorSig",sep = ""),outputDir)

