# 功能：
    # 1.批次效应矫正。
    # 2.同时对样本类型进行排序
# 输入：
    # rawCount.csv：需要矫正批次效应的表达数据
    # sampleBatchInfo.csv：表达数据中对应样本批次信息
# 输出：
    # *AdjSampleTree.png：聚类树查看批次效应矫正的影响
    # *AdjSampleTree2.png：聚类树查看批次效应矫正的影响
    # *SampleTree.png：聚类树查看批次效应矫正的影响
    # *SampleTree2.png：聚类树查看批次效应矫正的影响
    # *NoBatchEffects.txt：提示对应样本类型没有批次效应
    # *OneSamplePerBatch.txt：提示对应样本类型存在一个批次只有一个样本的情况
    # rawCountBEAdj.RData：批次效应校正后表达数据
    # rawCountBEAdj.csv：批次效应校正后表达数据
    # rawCountBEAdj.xlsx：批次效应校正后表达数据


# 0工具包
rm(list = ls())
gc()
library(modifyTools)
library(writexl)
library(sva)
library(plyr)
# 0函数
batchEffect = function(data,sampleBatchInfo,dataName){
    # batch信息
    batch = sampleBatchInfo$Batch[match(names(data),sampleBatchInfo$SampleID)]
    # 判断batch情况
    #   去除单样本一个批次的样本
    index = duplicated(batch) | duplicated(batch,fromLast = T)
    if(sum(index) == length(index)){
        # 没有单样本一个批次的情况
        if(all(batch == batch[1])){
            # 没有批次效应
            dataAdj = data
            write.table(paste(dataName,"No Batch Effects",sep = ""),
                        row.names = F,col.names = F,quote = F,
                        file = paste(outputDir,dataName,"NoBatchEffects.txt",sep = ""))
        }else{
            # 有批次效应
            adj = ComBat_seq(as.matrix(data),batch = batch)
            dataAdj = data.frame(adj,check.names = F)
            tree1 = hclust(dist(t(data)),method = "average")
            # 检测批次效应：聚类树观察距离
            png(filename = paste(outputDir,dataName,"SampleTree2.png",sep = ""),
                width = 500,height = 400)
            plot(tree1,main = paste("Sample Clustering (",dataName,")"),
                 xlab = NA,sub = NA,cex.lab = 1.5,cex.axis = 1.1,label = batch)
            dev.off()
            tree2 = hclust(dist(t(dataAdj)),method = "average")
            png(filename = paste(outputDir,dataName,"AdjSampleTree2.png",sep = ""),
                width = 500,height = 400)
            plot(tree2,main = paste("Sample Clustering (",dataName,"Adj)"),
                 xlab = NA,sub = NA,cex.lab = 1.5,cex.axis = 1.1,label = batch)
            dev.off()
            png(filename = paste(outputDir,dataName,"SampleTree.png",sep = ""),
                width = 500,height = 400)
            plot(tree1,main = paste("Sample Clustering (",dataName,")"),
                 xlab = NA,sub = NA,cex.lab = 1.5,cex.axis = 1.1)
            dev.off()
            tree2 = hclust(dist(t(dataAdj)),method = "average")
            png(filename = paste(outputDir,dataName,"AdjSampleTree.png",sep = ""),
                width = 500,height = 400)
            plot(tree2,main = paste("Sample Clustering (",dataName,"Adj)"),
                 xlab = NA,sub = NA,cex.lab = 1.5,cex.axis = 1.1)
            dev.off()
        }
    }else{
        # 存在单样本一个批次的情况
        dataUnique = data.frame(data[,!index],check.names = F)
        names(dataUnique) = names(data)[!index]
        dataMuti =  data.frame(data[,index],check.names = F)
        names(dataMuti) = names(data)[index]
        batchUnique = batch[!index]
        batchMuti = batch[index]
        write.table(names(dataUnique),row.names = F,col.names = F,quote = F,
                    file = paste(outputDir,dataName,"OneSamplePerBatch.txt",sep = ""))
        if(length(batchMuti) == 0 | all(batchMuti == batchMuti[1])){
            # 没有批次效应
            dataAdj = data
            write.table("No Batch Effects",row.names = F,col.names = F,quote = F,
                        file = paste(outputDir,dataName,"NoBatchEffects.txt",sep = ""))
        }else{
            # 存在批次效应
            adj = ComBat_seq(as.matrix(dataMuti),batch = batchMuti)
            dataAdj = data.frame(adj,check.names = F)
            dataAdj = cbind(dataUnique,dataAdj)
            # 检测批次效应：聚类树观察距离
            tree1 = hclust(dist(t(data)),method = "average")
            png(filename = paste(outputDir,dataName,"SampleTree2.png",sep = ""),
                width = 500,height = 400)
            plot(tree1,main = paste("Sample Clustering (",dataName,")"),
                 xlab = NA,sub = NA,cex.lab = 1.5,cex.axis = 1.1,label = batch)
            dev.off()
            png(filename = paste(outputDir,dataName,"SampleTree.png",sep = ""),
                width = 500,height = 400)
            plot(tree1,main = paste("Sample Clustering (",dataName,")"),
                 xlab = NA,sub = NA,cex.lab = 1.5,cex.axis = 1.1)
            dev.off()
            # 批次效应效果：聚类树观察距离
            tree2 = hclust(dist(t(dataAdj)),method = "average")
            png(filename = paste(outputDir,dataName,"AdjSampleTree2.png",sep = ""),
                width = 500,height = 400)
            plot(tree2,main = paste("Sample Clustering (",dataName,"Adj)"),
                 xlab = NA,sub = NA,cex.lab = 1.5,cex.axis = 1.1,label = batch)
            dev.off()
            png(filename = paste(outputDir,dataName,"AdjSampleTree.png",sep = ""),
                width = 500,height = 400)
            plot(tree2,main = paste("Sample Clustering (",dataName,"Adj)"),
                 xlab = NA,sub = NA,cex.lab = 1.5,cex.axis = 1.1)
            dev.off()
        }
        
    }
    return(dataAdj)
}


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()
dataFileName = fileNames[grepl("rawCount",fileNames)]
sampleBatchInfoFileName = fileNames[grepl("Info.csv$",fileNames)]
# 数据导入
dataInfo = read.csv(dataFileName,check.names = F,row.names = 1)
sampleBatchInfo = read.csv(sampleBatchInfoFileName,check.names = F,row.names = 1)
data = dataInfo


# 3处理
# 样本分组信息
groups = c("HB_HC","HB_NAG_BD","HB_NAG_AD","HB_AG_BD","HB_AG_AD","HB_CA")
# 批次效应矫正：每组样本单独进行
resultList = list()
for(group in groups){
    index = grepl(group,names(data))
    dataTemp = data[,index]
    dataAdj = data.frame(batchEffect(dataTemp,sampleBatchInfo,group),check.names = F)
    resultList = c(resultList,list(dataAdj))
}
#   合并校正后的数据
result = data.frame(resultList[[1]])
for(i in 2:length(resultList)){
    result = cbind(result,resultList[[i]])
}


# 4保存
dfSave(result,"rawCountBEAdj",outputDir)