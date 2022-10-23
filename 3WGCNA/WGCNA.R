# 功能：
    # 软阈值选取、网络构建、动态剪切、模块合并、图形绘制
# 输入：
    # stage*.csv：RSEM后的基因表达矩阵文件
# 输出：
# 注意：
    # 本版本的基因为名称：moduleSetList程序部分不同


# 0工具包
rm(list = ls())
gc()
library(WGCNA)
library(writexl)
library(parallel)
library(modifyTools)
makeCluster(16)
# 0函数
WGCNAFilter = function(data){
    # WGCNA样本和基因检测：缺失值数目大于4，没有变化的基因
    gsg = goodSamplesGenes(data)
    if (!gsg$allOK){
        if (sum(!gsg$goodGenes)>0){
            printFlush(paste("Removing genes:", paste(colnames(data)[!gsg$goodGenes], collapse = ", ")))
        }
        if (sum(!gsg$goodSamples)>0)
            printFlush(paste("Removing samples:", paste(rownames(data)[!gsg$goodSamples], collapse = ", ")))
    }
    # 对应变量
    return(list(gsg$goodGenes,gsg$goodSamples))
}
powerSelect = function(sft,dataName,outputDir,nSamples){
    # 绘图1：软阈值与R2
    png(paste(outputDir,dataName,"_SoftThresholdPick(1).png",sep = ""),
        height = 500,width = 550)
    par(mar = c(5,5,3,3))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab = "Soft Threshold(power)",
         ylab = "Scale Free Topology Model Fit(R^2)",
         type = "n",main = paste("Scale independence（",dataName,"）",sep = ""),
         cex.axis = 2,cex.lab = 2,cex.main = 2.3)
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels = powers,cex = 1.5,col = "red")
    abline(h = 0.9,col = "red")
    dev.off()
    # 绘图2：软阈值与平均连接度
    png(paste(outputDir,dataName,"_SoftThresholdPick(2).png",sep = ""),
        height = 500,width = 550)
    par(mar = c(5,5,3,3))
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
         type = "n",main = paste("Mean connectivity（",dataName,"）",sep = ""),
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
moduleMerge = function(MEs,dynamicColors,expr,dataName){
    # 计算模块向量聚类树
    MEDiss = 1-cor(MEs);
    METree = hclust(as.dist(MEDiss), method = "average")
    # 绘图：模块合并
    png(file = paste(outputDir,dataName,"_METree.png",sep = ""),
        width = 550,height = 500)
    par(mar = c(5,5,3,3))
    plot(METree, main = paste("Clustering of module eigengenes（",dataName,"）",sep = ""),
         xlab = "", sub = "",labels = F,
         cex.axis = 1.5,cex.lab = 1.5,cex.main = 2,cex = 1.2,)
    abline(h=0.25, col = "red")
    dev.off()
    # 合并距离小于0.25的模块
    merge = mergeCloseModules(expr, dynamicColors, cutHeight = 0.25, verbose = 3)
    return(merge)
}
geneTreeAndColorPlot = function(geneTree,dynamicColors,mergedColors,outputDir,dataName){
    png(file = paste(outputDir,dataName,"_GeneDendro.png",sep = ""),
        width = 800,height = 600)
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                        c("Dynamic TreeCut", "Merged dynamic"),dendroLabels = F,
                        hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                        marAll = c(1, 9.5, 3, .01),cex.colorLabels = 1.5,
                        cex.dendroLabels = 1.5,cex.axis = 1.5,cex.lab = 1.5,
                        cex.main = 2,
                        main = paste("Cluster Dendrogram（",dataName,"）",sep ="")
                        )
    dev.off()
}
formModuleSetList = function(expr,mergedColors){
    colorset = names(table(mergedColors))
    moduleSetList = list()
    for(i in colorset){
        colorgene = names(expr)[i == mergedColors]
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
print("Input dataName:")
dataName = readline()
# 数据导入
expr = read.csv(fileName,row.names = 1,check.names = F)
expr = data.frame(t(expr),check.names = F)
expr = log(expr + 1)


# 3数据预处理
# Filter1:WGCNA标准
#   合格的基因和样本
goodSamplesGenes = WGCNAFilter(expr)
goodGenes = goodSamplesGenes[[1]]
goodSamples = goodSamplesGenes[[2]]
#   表达数据筛选
exprFiltered1 = expr[goodSamples,goodGenes]
# Filter2:离群样本
#   样本聚类数
sampleTree = hclust(dist(exprFiltered1),method = "average")
plot(sampleTree, main = paste("Sample clustering to detect outliers（",dataName,"）",sep = ""),
     sub = "",xlab = "",cex.lab = 1.5,cex.axis = 1.1)
#   根据图形结果选择及阈值
print("Input threshold:")
threshold = readline()
#   去除离群样本
clust = cutreeStatic(sampleTree, cutHeight = as.numeric(threshold), minSize = 1)
keepedSamples = (clust == 1)
#   表达数据筛选
exprFiltered2 = exprFiltered1[keepedSamples,]
# 样本聚类图最终绘制
png(filename = paste(outputDir,dataName,"_SampleTree.png",sep = ""),width = 550,height = 500)
plot(sampleTree, main = paste("Sample clustering to detect outliers（",dataName,"）",sep = ""),
     sub = "",xlab = "",labels = F,cex.lab = 1.5,cex.axis = 1.5,cex.main = 2)
abline(h = threshold,col = "red")
dev.off()


# 3选择软阈值
# 计算软阈值
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
sft = pickSoftThreshold(exprFiltered2, powerVector = powers, verbose = 5,RsquaredCut = 0.9)
# 绘制软阈值选取图
softpower = powerSelect(sft,dataName,outputDir,dim(exprFiltered2)[1])
write.table(softpower,file = paste(outputDir,dataName,"_Softpower.txt",sep = ""),
            quote = F,row.names = F,col.names = F)


# 4网络构建与模块划分
# 网络构建
dissTOM = networkConstruct(exprFiltered2,softpower)
# 基因聚类树
geneTree = hclust(as.dist(dissTOM), method = "average")
# 动态剪切
dynamicColors = dynamicCut(geneTree)
# 获得模块向量
# load(file = paste(outputDir,dataName,"Network.RData",sep = ""))
MEList = moduleEigengenes(exprFiltered2,colors = dynamicColors)
MEs = MEList$eigengenes
# 判断是否存在模块向量为NA
if(TRUE %in% is.na(MEs[1,])){
    NAME = MEs[which(is.na(MEs[1,]))]
    NoNAME = MEs[-which(is.na(MEs[1,]))]
}else{
    NoNAME = MEs
    NAME = NA
}
# 模块合并
merge = moduleMerge(NoNAME,dynamicColors,exprFiltered2,dataName)
mergedColors = merge$colors
mergedMEs = merge$newMEs
# 绘图：基因聚类与模块颜色
geneTreeAndColorPlot(geneTree,dynamicColors,mergedColors,outputDir,dataName)
# 模块列表
moduleSetList = formModuleSetList(exprFiltered2,mergedColors)
# 去除灰色模块
modGrey = moduleSetList$grey
moduleSetList$grey = NULL
mergedMEs$MEgrey = NULL
exprGreyRemoved = exprFiltered2[,dynamicColors != "grey"]
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
moduleSetListDF = moduleSetListDF[order(as.numeric(gsub("^X","",moduleSetListDF$name))),]
moduleSetListDF$geneNum = as.numeric(moduleSetListDF$geneNum)


# 5保存
# 汇总结果
sampleTotal = nrow(expr)
outlier = sum(clust != 1)
sampleReserved = sum(clust == 1)
modTotal = length(table(dynamicColors))
modMerged = length(table(mergedColors))
resultInfo = data.frame(SampleTotal = sampleTotal,Outlier = outlier,
                        SampleReserved = sampleReserved,Softpower = softpower,
                        ModTotal = modTotal,ModMerged = modMerged,
                        GeneNumGrey = length(modGrey),
                        GeneNumMax = max(moduleSetListDF$geneNum),
                        GeneNumMin = min(moduleSetListDF$geneNum),
                        GeneNumMean = round(mean(moduleSetListDF$geneNum),0))
# 保存
save(geneTree,dynamicColors,MEs,NAME,NoNAME,mergedColors,mergedMEs,sft,
     file = paste(outputDir,dataName,"_Network.RData",sep = ""))
save(sft,file = paste(outputDir,dataName,"_SftSelect.RData",sep = ""))
save(moduleSetList,file = paste(outputDir,dataName,"_ModuleSetList.RData",sep = ""))
write.csv(moduleSetListDF,file = paste(outputDir,dataName,"_ModuleSetList.csv",sep = ""),row.names = F)
write.csv(resultInfo,file = paste(outputDir,dataName,"_ResultInfo.csv",sep = ""),row.names = F)
dfSave(exprGreyRemoved,paste(dataName,"_GreyRemoved",sep = ""),outputDir)
dfSave(exprFiltered2,paste(dataName,"_Filtered2",sep = ""),outputDir)