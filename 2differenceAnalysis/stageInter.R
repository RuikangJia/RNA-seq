# 功能：对不同阶段的差异基因进行集合分析
    # 主要内容：
# 输入：
    # III_CElseInter.csv  
    # IV_CElseInter.csv  
    # II_CElseInter.csv   
    # I_CElseInter.csv
# 输出：
    # InterPlot.png             
    # uniqueGeneStageIII.csv
    # uniqueGeneStageI.RData    
    # uniqueGeneStageIV.RData
    # uniqueGeneStageI.csv      
    # uniqueGeneStageIV.csv
    # uniqueGeneStageII.RData   
    # uniqueGeneStageI_II_III_IV.RData
    # uniqueGeneStageII.csv     
    # uniqueGeneStageI_II_III_IV.csv
    # uniqueGeneStageIII.RData
    # totalGeneList.csv
    # totalGeneList.RData
# 工具包
library(VennDiagram)
library(UpSetR)
library(ggplot2)
library(modifyTools)
# 函数
upsetDraw = function(upsetData,setsOrder){
    print(upset(upsetData,sets = setsOrder,
          keep.order = T,order.by = "freq",point.size = 3,line.size = 1,
          text.scale = c(3,3,3,3,3,3),
          queries = list(
              list(query = intersects,params = list("stageI"),color = 3,active = T,query.name = "stageI"),
              list(query = intersects,params = list("stageII"),color = 4,active = T,query.name = "stageII"),
              list(query = intersects,params = list("stageIII"),color = 5,active = T,query.name = "stageIII"),
              list(query = intersects,params = list("stageIV"),color = 6,active = T,query.name = "stageIV"),
              list(query = intersects,params = list("stageI","stageII","stageIII","stageIV"),color = 2,active = T,query.name = "allStage")
          )))
}
extractGene = function(upsetDataTF,... = ""){
    # 条件为真的列
    indexCol = match(list(...),colnames(upsetDataTF))
    # 目标基因行
    indexRes = rep(TRUE,dim(upsetDataTF)[1])
    for(i in 1:dim(upsetDataTF)[2]){
        if(i %in% indexCol){
            indexRes = indexRes & upsetData[,i]
        }else{
            indexRes = indexRes & (!upsetData[,i])
        }
    }
    return(rownames(upsetDataTF)[indexRes])
}


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()[grep(".csv$",dir())]
# 数据导入
dataList = lapply(fileNames, read.csv,row.names = 1,check.names = F)
geneListStageI = rownames(dataList[[1]])
geneListStageII = rownames(dataList[[2]])
geneListStageIII = rownames(dataList[[3]])
geneListStageIV = rownames(dataList[[4]])

            
# 3处理
# 所有基因
totalGeneList = union(
    union(geneListStageI,geneListStageII),
    union(geneListStageIII,geneListStageIV)
    )
# 基因阶段逻辑表
upsetDataTF = data.frame(
    row.names = totalGeneList, 
    stageI = totalGeneList %in% geneListStageI,
    stageII = totalGeneList %in% geneListStageII,
    stageIII = totalGeneList %in% geneListStageIII,
    stageIV = totalGeneList %in% geneListStageIV,
    check.names = F)
# 基因阶段01表
upsetData = upsetDataTF
upsetData[upsetData == T] = 1
# 绘制
png(filename = paste(outputDir,"InterPlot.png",sep = ""),
    width = 800,height = 600)
upsetDraw(upsetData,c("stageIV","stageIII","stageII","stageI"))
dev.off()
# 提取
uniqueGeneStageI = extractGene(upsetDataTF,"stageI")
uniqueGeneStageII = extractGene(upsetDataTF,"stageII")
uniqueGeneStageIII = extractGene(upsetDataTF,"stageIII")
uniqueGeneStageIV = extractGene(upsetDataTF,"stageIV")
uniqueGeneStageI_II_III_IV = extractGene(upsetDataTF,"stageI","stageII","stageIII","stageIV")


# 4保存
dfSave(uniqueGeneStageI,"uniqueGeneStageI",outputDir)
dfSave(uniqueGeneStageII,"uniqueGeneStageII",outputDir)
dfSave(uniqueGeneStageIII,"uniqueGeneStageIII",outputDir)
dfSave(uniqueGeneStageIV,"uniqueGeneStageIV",outputDir)
dfSave(uniqueGeneStageI_II_III_IV,"uniqueGeneStageI_II_III_IV",outputDir)
dfSave(totalGeneList,"totalGeneList",outputDir)
