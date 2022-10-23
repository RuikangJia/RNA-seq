# 功能：
    # 1.DEG中基因的提取。
    # 2.集合分析（逻辑表制作、绘图、提取共有特有基因）。
# 主要问题：
    # 1.图形参数设置。
    # 2.逻辑表制作和使用（分类特有基因）
# 输入：
    # AG_vs_HCDEGSorted.csv：基因运算的数据
    # CA_vs_HCDEGSorted.csv：基因运算的数据
    # NAG_vs_HCDEGSorted.csv：基因运算的数据
# 输出：
    # setAnalysis.png：集合运算图
    # interGene.RData：共有基因
    # interGene.csv：共有基因
    # interGene.xlsx：共有基因
    # uniqueGeneAG.RData：各阶段特有基因
    # uniqueGeneAG.csv
    # uniqueGeneAG.xlsx
    # uniqueGeneCA.RData
    # uniqueGeneCA.csv
    # uniqueGeneCA.xlsx
    # uniqueGeneNAG.RData
    # uniqueGeneNAG.csv
    # uniqueGeneNAG.xlsx


# 0工具包
rm(list = ls())
gc()
library(modifyTools)
library(writexl)
library(UpSetR)
library(ggplot2)
# 0函数
getGene = function(data){
    return(rownames(data))
}
getGroup = function(data){
    return(strsplit(data[1],split = "_")[[1]][1])
}
mutilUnion = function(list){
    res = list[[1]]
    for(i in 2:length(list)){
        res = union(res,list[[i]])
    }
    return(res)
}
mutilInter = function(list){
    res = list[[1]]
    for(i in 2:length(list)){
        res = intersect(res,list[[i]])
    }
    return(res)
}
upsetDraw = function(upsetData,setsOrder,plotName){
    print(upset(logicDF,sets = group,keep.order = T,
                main.bar.color = "blue",sets.bar.color = "red",
                matrix.color = "black",
                sets.x.label = plotName,
                order.by = "freq",point.size = 6,line.size = 1.5,
                text.scale = c(3,3,3,1.8,3,3))
    )
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
# 根据需求手动指定顺序：对应图中类别展示顺序（逆序）
dataFileName = c(
    dataFileName[grepl("^CA",dataFileName)],
    dataFileName[grepl("^AG",dataFileName)],
    dataFileName[grepl("^NAG",dataFileName)],
    dataFileName[grepl("^HC",dataFileName)]
)
dataName = gsub(".csv$","",dataFileName)
# 图片标题
if(grepl("DEG",dataName[1])){
    name = "Total DEGs"
}else if(grepl("up",dataName[1])){
    name = "Up DEGs"
}else{
    name = "Down DEGs"
}
plotName = name
type = unlist(strsplit(plotName,split = " "))[1]


# 3处理
# 3.1数据准备
#   提取分组名称
group = unlist(lapply(dataName,getGroup))
#   读入数据文件
dataList = lapply(dataFileName,read.csv,check.names = F,row.names = 1)
#   提取基因名称
geneList = lapply(dataList,getGene)
# 3.2制作01表
#   总基因
totalGene = mutilUnion(geneList)
#   逻辑表
logicList = list()
for(i in 1:length(geneList)){
    logic = totalGene %in% geneList[[i]]
    logicList = c(logicList,list(logic))
}
logicDF = data.frame(logicList)
names(logicDF) = group
rownames(logicDF) = totalGene
logicDF[logicDF == T] = 1
# 3.3绘制集合分析图
png(filename = paste(outputDir,paste("setAnalysis",type,".png",sep = ""),sep = ""),
    width = 800,height = 600)
upsetDraw(upsetData,group,plotName)
dev.off()
# 3.4提取基因
# 共有
interGene = mutilInter(geneList)
# 特有
uniqueLogicDF = logicDF[apply(logicDF,1,sum) == 1,]
uniqueGeneList = list()
for(g in group){
    uniqueGene = rownames(uniqueLogicDF)[uniqueLogicDF[,match(g,names(uniqueLogicDF))] == 1]
    uniqueGeneList = c(uniqueGeneList,list(uniqueGene))
}


# 4保存
interGeneDF = data.frame(GeneList = interGene,check.names = F)
dfSave(interGeneDF,paste("interGene",type,sep = ""),outputDir)
for(i in 1:length(group)){
    UniqueGeneDF = data.frame(UniqueGene = uniqueGeneList[[i]],check.names = F)
    dfSave(UniqueGeneDF,paste("uniqueGene",group[i],sep = ""),outputDir)
}
