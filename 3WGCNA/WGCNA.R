# 功能：进行WGCNA分析
    # 主要内容：软阈值选取、网络构建、动态剪切、模块合并、图形绘制
# 输入：
    # stage*.csv：RSEM后的基因表达矩阵文件
# 输出：
    #
    #
# 工具包
rm(list = ls())
gc()
library(WGCNA)
library(parallel)
library(modifyTools)
makeCluster(16)
# 函数
powerSelect = function(sft,dataName,outputDir,nSamples){
    # 绘图1：软阈值与R2
    png(paste(outputDir,dataName,"softThresholdPick(1).png",sep = ""),
        height = 500,width = 550)
    par(mar = c(5,5,3,3))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab = "Soft Threshold(power)",
         ylab = "Scale Free Topology Model Fit(R^2)",
         type = "n",main = paste("Scale independence"),
         cex.axis = 2,cex.lab = 2,cex.main = 2.3)
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels = powers,cex = 1.5,col = "red")
    abline(h = 0.9,col = "red")
    dev.off()
    # 绘图2：软阈值与平均连接度
    png(paste(outputDir,dataName,"softThresholdPick(2).png",sep = ""),
        height = 500,width = 550)
    par(mar = c(5,5,3,3))
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
         type = "n",main = paste("Mean connectivity"),
         cex.axis = 2,cex.lab = 2,cex.main = 2.3)
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers,
         cex = 1.5, col = "red")
    dev.off()
    if(is.na(sft$powerEstimate) | sft$powerEstimate == 1){
        if(nSamples < 20){
            power = 9
        }else if(nSamples < 30){
            power = 8
        }else if(nSamples < 40){
            power = 7
        }else{
            power = 6
        }
    }else{
        power = sft$powerEstimate
    }
    return(power)
}
networkConstruct = function(expr,softpower){
    # 计算邻接矩阵
    adjacency = adjacency(expr, power = softpower)
    # 计算拓扑重叠矩阵
    TOM = TOMsimilarity(adjacency)
    # 计算相异矩阵
    dissTOM = 1-TOM
    # 释放内存
    rm(adjacency,TOM)
    gc()
    return(dissTOM)
}
dynamicCut = function(geneTree){
    # 对基因树进行动态剪切
    dynamicMods = cutreeDynamic(dendro = geneTree,distM = dissTOM,deepSplit = 2,
                                pamRespectsDendro = FALSE,minClusterSize = 30)
    # 每个基因的颜色
    dynamicColors = labels2colors(dynamicMods)
    return(dynamicColors)
}
moduleMerge = function(MEs,dynamicColors){
    # 计算模块向量聚类树
    MEDiss = 1-cor(MEs);
    METree = hclust(as.dist(MEDiss), method = "average")
    # 绘图：模块合并
    png(file = paste(outputDir,dataName,"METree.png",sep = ""),
        width = 1000,height = 800)
    par(mar = c(5,5,3,3))
    plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "",
         cex.axis = 2.2,cex.lab = 2.2,cex.main = 2.5,cex = 1.3,)
    abline(h=0.25, col = "red")
    dev.off()
    # 合并距离小于0.25的模块
    merge = mergeCloseModules(expr, dynamicColors, cutHeight = 0.25, verbose = 3)
    return(merge)
}
geneTreeAndColorPlot = function(geneTree,dynamicColors,mergedColors,outputDir,dataName){
    png(file = paste(outputDir,dataName,"geneDendro.png"),
        width = 800,height = 600)
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                        c("Dynamic TreeCut", "Merged dynamic"),dendroLabels = F,
                        hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                        marAll = c(1, 9.5, 3, .01),cex.colorLabels = 1.5,
                        cex.dendroLabels = 2,cex.axis = 2,cex.lab = 2,
                        cex.main = 2.5)
    dev.off()
}
formModuleSetList = function(expr,mergedColors){
    colorset = names(table(mergedColors))
    temp = unlist(strsplit(colnames(expr),split = "|",fixed = T))
    temp = matrix(temp,ncol=2,byrow=T)
    geneSymbols = temp[,1]
    geneEntries = temp[,2]
    moduleSetList = list()
    for(i in colorset){
        colorgene = geneEntries[i== mergedColors]
        moduleSetList[i] = list(colorgene)
    }
    return(moduleSetList)
}


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileName = dir()[grep(".csv$",dir())]
dataName = gsub(".csv$","",fileName)
# 数据导入
expr = read.csv(fileName,row.names = 1,check.names = F)

# 
# # 3选择软阈值
# # 计算软阈值
# powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
# sft = pickSoftThreshold(expr, powerVector = powers, verbose = 5)
# # 绘制软阈值选取图
# softpower = powerSelect(sft,dataName,outputDir,dim(expr)[1])
# 
# 
# # 4网络构建与模块划分
# # 网络构建
# dissTOM = networkConstruct(expr,softpower)
# # 基因聚类树
# geneTree = hclust(as.dist(dissTOM), method = "average")
# # 动态剪切
# dynamicColors = dynamicCut(geneTree)
# 获得模块向量
load(file = paste(outputDir,dataName,"Network.RData",sep = ""))
MEList = moduleEigengenes(expr,colors = dynamicColors)
MEs = MEList$eigengenes
# # 模块合并
merge = moduleMerge(MEs,dynamicColors)
mergedColors = merge$colors
mergedMEs = merge$newMEs
# 绘图：基因聚类与模块颜色
geneTreeAndColorPlot(geneTree,dynamicColors,mergedColors,outputDir,dataName)
# 模块列表
moduleSetList = formModuleSetList(expr,mergedColors)
# 去除灰色模块
modGrey = moduleSetList$grey
moduleSetList$grey = NULL
moduleSetList["grey"]
mergedMEs$MEgrey = NULL
# X标记替代颜色
markMods = paste("X",1:dim(mergedMEs)[2],sep = "")
names(moduleSetList) = markMods[match(names(moduleSetList),gsub("^ME","",names(mergedMEs)))]
names(mergedMEs) = markMods
# 模块基因数目信息
moduleSetListDF = data.frame()
for(i in 1:length(moduleSetList)){
    name = names(moduleSetList[i])
    len = length(moduleSetList[[i]])
    moduleSetListDF = rbind(moduleSetListDF,c(name,len))
}
names(moduleSetListDF) = c("name","geneNum")

# 5保存
# save(sft,file = paste(outputDir,dataName,"SftSelect.RData",sep = ""))
save(moduleSetList,file = paste(outputDir,dataName,"ModuleSetList.RData",sep = ""))
save(geneTree,dynamicColors,MEs,mergedColors,mergedMEs,
     file = paste(outputDir,dataName,"Network.RData",sep = ""))
write.csv(moduleSetListDF,file = paste(outputDir,dataName,"ModuleSetList.csv",sep = ""),row.names = F)
