# 功能：对基因和样本进行过滤
    # 主要内容：离群样本、低质量样本、低质量基因
# 输入： 
    # dataClinStageI.csv:临床数据
    # dataExprStageI.csv:表达数据
# 输出：
    # RSEMStageISampleTree.png:聚类图
    # ClinStageIFiltered.csv：去除了WGCNA过滤样本和离群样本
    # RSEMStageIFiltered1.csv：WGCNA过滤基因和样本
    # RSEMStageIFiltered2.csv：去除离群样本
    # RSEMStageIstaticInfo.csv
    # rawCountsStageIFiltered1.csv：WGCNA过滤样本
    # rawCountsStageIFiltered2.csv：去除离群样本
    # rawCountsStageIstaticInfo.csv
    # ClinStageIFiltered.RData
    # RSEMStageIFiltered1.RData
    # RSEMStageIFiltered2.RData
    # RSEMStageIstaticInfo.RData
    # rawCountsStageIFiltered1.RData
    # rawCountsStageIFiltered2.RData
    # rawCountsStageIstaticInfo.RData
#工具包
library(WGCNA)
library(modifyTools)

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

# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 文件名称
fileNameClin = dir()[grep("^Clin",dir())]
fileNameRawCounts = dir()[grep("^rawCounts",dir())]
fileNameRSEM = dir()[grep("^RSEM",dir())]
# 数据名称
dataNameClin = gsub(".csv","",fileNameClin)
dataNameRawCounts = gsub(".csv","",fileNameRawCounts)
dataNameRSEM = gsub(".csv","",fileNameRSEM)
# 数据导入
dataClin = read.csv(fileNameClin,check.names = F,row.names = 1)
dataRSEM = read.csv(fileNameRawCounts,check.names = F,row.names = 1)
dataRawCounts = read.csv(fileNameRSEM,check.names = F,row.names = 1)


# 3过滤
# Filter1:WGCNA标准
#   合格的基因和样本
goodSamplesGenes = WGCNAFilter(dataRSEM)
goodGenes = goodSamplesGenes[[1]]
goodSamples = goodSamplesGenes[[2]]
#   表达数据筛选
dataRSEMFiltered1 = dataRSEM[goodSamples,goodGenes]
# Filter2:离群样本
#   样本聚类数
sampleTree = hclust(dist(log(dataRSEMFiltered1 + 1)),method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers",sub = "",xlab = "",cex.lab = 1.5,cex.axis = 1.1)
#   根据图形结果选择及阈值
print("Input threshold:")
threshold = readline()
#   去除离群样本
clust = cutreeStatic(sampleTree, cutHeight = as.numeric(threshold), minSize = 1)
keepedSamples = (clust == 1)
#   表达数据筛选
dataRSEMFiltered2 = dataRSEMFiltered1[keepedSamples,]
# 样本聚类图最终绘制
png(filename = paste(outputDir,dataNameRSEM,"SampleTree.png",sep = ""),width = 400,height = 400)
plot(sampleTree, main = "Sample clustering to detect outliers",sub = "",xlab = "",labels = F,cex.lab = 1.5,cex.axis = 1.1,cex.main = 1.5)
abline(h = threshold,col = "red")
dev.off()


# 4临床数据和rawCounts数据调整
# rawCounts
dataRawCountsFiltered1 = dataRawCounts[goodSamples,]
dataRawCountsFiltered2 = dataRawCountsFiltered1[keepedSamples,]
# clin
dataClinFiltered = dataClin[rownames(dataClin) %in% substr(rownames(dataRSEMFiltered2),1,12),]


# 5统计信息
staticInfoRSEM = data.frame(data = dim(dataRSEM),filtered1 = dim(dataRSEMFiltered1),filtered2 = dim(dataRSEMFiltered2))
rownames(staticInfoRSEM) = c("sample","gene")
staticInfoRawCounts = data.frame(data = dim(dataRawCounts),filtered1 = dim(dataRawCountsFiltered1),filtered2 = dim(dataRawCountsFiltered2))
rownames(staticInfoRawCounts) = c("sample","gene")


# 6保存
dfSave(dataRSEMFiltered1,paste(dataNameRSEM,"Filtered1",sep = ""),outputDir)
dfSave(dataRSEMFiltered2,paste(dataNameRSEM,"Filtered2",sep = ""),outputDir)
dfSave(dataRawCountsFiltered1,paste(dataNameRawCounts,"Filtered1",sep = ""),outputDir)
dfSave(dataRawCountsFiltered2,paste(dataNameRawCounts,"Filtered2",sep = ""),outputDir)
dfSave(dataClinFiltered,paste(dataNameClin,"Filtered",sep = ""),outputDir)
dfSave(staticInfoRSEM,paste(dataNameRSEM,"staticInfo",sep = ""),outputDir)
dfSave(staticInfoRawCounts,paste(dataNameRawCounts,"staticInfo",sep = ""),outputDir)