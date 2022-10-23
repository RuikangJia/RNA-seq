# 功能：进行基因表达差异性分析
    # 主要内容：差异分析；筛选差异基因；绘制火山图
# 输入：
    # 1.csv:基因表达数据未处理 
    # 2.csv:基因表达数据处理
# 输出：
    # calcuRes.csv、calcuRes.RData:所有基因LFC计算结果
    # DEG.csv、DEG.RData:差异基因LFC计算结果
    # DEGSorted.csv、DEGSorted.csv：差异基因LFC排序
    # upRegulated.csv、upRegulated.csv：上调基因基因LFC排序
    # downRegulated.csv、downRegulated.csv：下调基因基因LFC排序
    # volanic.png:差异基因火山图
# 工具包
rm(list = ls())
library(DESeq2)
library(ggplot2)
library(modifyTools)
library(writexl)


# 0函数
groupInfo = function(data1,data2){
    group = as.factor(c(rep(1,dim(data1)[1]),rep(2,dim(data2)[1])))
    name = c(rownames(data1),rownames(data2))
    colData = data.frame(row.names = name, group_list=group)
    return(colData)
}
DECalculation = function(data,colData){
    # 方差估计
    dds = DESeqDataSetFromMatrix(round(data),colData = colData,design = ~ group_list)
    # 正常表达值调整为分母
    dds$group_list = relevel(dds$group_list, ref = "1")
    # 生成结果
    dds = DESeq(dds)
    res = results(dds)
    # 根据P值排序
    res  =  res[order(res$padj),]
    return(res)
}
DEGFilter = function(data,LFC = 1,p = 0.05){
    data = subset(data,abs(log2FoldChange) >= LFC & padj < p)
    return(data)
}
volcanicPlot = function(drawData,dataName,outputDir,title,LFC = 1,p = 0.05){
    values = c("blue","grey","red")
    levels = c("Down","NoSignifi","Up")
    # 去除NA
    drawData = subset(drawData,complete.cases(drawData))
    # p值区分不同基因，且因子化
    drawData$threshold = as.factor(ifelse(drawData$padj < p & abs(drawData$log2FoldChange) >= LFC, ifelse(drawData$log2FoldChange>= LFC ,'Up','Down'),'NoSignifi'))
    drawData$threshold = factor(drawData$threshold,levels = c("Down","NoSignifi","Up"),ordered = T)
    value = values[levels %in% drawData$threshold]
    #ggplot具体绘图
    ggplot(data = drawData, aes(x = log2FoldChange, y = -log10(padj), color=threshold)) +
        geom_point(alpha=0.5, size=2,) +
        scale_color_manual(values = value) +
        ylab("-log10(Padj)")+
        xlab("log2(FoldChange)")+
        labs(title = title) +
        geom_vline(xintercept = c(-1,1),lty = 4,col = "black",lwd = 0.8)+
        geom_hline(yintercept = -log10(0.05),lty = 4,col = "black",lwd = 0.8)+
        theme(panel.background = element_blank(),
              plot.title = element_text(size = 25,hjust = 0.5),
              axis.line = element_line(colour = "black"), 
              axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              legend.key = element_rect(fill = "white"),
              legend.title = element_blank(),
              legend.text = element_text(size = 20),
              legend.position = c(0.9,0.9),
              legend.key.size = unit(1,"cm"))
    ggsave(file = paste(outputDir,dataName,".png",sep = ""),width = 10,height = 8)
}


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileName1 = dir()[grep("1.csv$",dir())]
fileName2 = dir()[grep("2.csv$",dir())]
# 手动指定数据名称
print("Input tyepName(Blood or Tissue):")
typeName = readline()
print("Input case groupName:")
caseName = readline()
print("Input control groupName:")
controlName = readline()
# 半自动生成相关名称
dataName = paste(caseName,controlName,sep = "_vs_")
plotTitle = paste("Gene Expression Differences(",
                  typeName,": ",caseName," vs. ",controlName,")",sep = "")
# 导入数据
data1 = read.csv(fileName1,check.names = F,row.names = 1)
data2 = read.csv(fileName2,check.names = F,row.names = 1)
data1 = data.frame(t(data1),check.names = F)
data2 = data.frame(t(data2),check.names = F)


# 3差异分析
# 基因在行，样本在列
DEData = data.frame(t(rbind(data1,data2)),check.names = F)
# 分组信息
groupInfo = groupInfo(data1,data2)
# 计算LFC
calcuRes = DECalculation(DEData,groupInfo)
calcuRes = data.frame(calcuRes)
# 筛选差异基因
DEG = DEGFilter(calcuRes)
DEGGeneList = 
# LFC排序
DEGSorted = DEG[order(-abs(DEG$log2FoldChange)),]
# 上下调
upRegulated = subset(DEGSorted,DEGSorted$log2FoldChange > 0)
downRegulated = subset(DEGSorted,DEGSorted$log2FoldChange < 0)
# 绘制火山图
volcanicPlot(calcuRes,paste(dataName,"_Volcanic",sep = ""),outputDir,title = plotTitle)
#   数目
totalNum = nrow(DEData)
DENum = nrow(DEG)
upNum = nrow(upRegulated)
downNum = nrow(downRegulated)
#   比例
DERate = scales::percent(DENum/totalNum,0.01)
upRate = scales::percent(upNum/totalNum,0.01)
downRate = scales::percent(downNum/totalNum,0.01)
#   结果信息
resultInfo = data.frame(total = totalNum,
                        "DEG(%)" = paste(DENum,"(",DERate,")",sep = ""),
                        "up(%)" = paste(upNum,"(",upRate,")",sep = ""),
                        "down(%)" = paste(downNum,"(",downRate,")",sep = ""),
                        check.names = F)



# 4保存变量
dfSave(calcuRes,paste(dataName,"_CalcuRes",sep = ""),outputDir)
dfSave(DEG,paste(dataName,"_DEG",sep = ""),outputDir)
dfSave(DEGSorted,paste(dataName,"_DEGSorted",sep = ""),outputDir)
dfSave(upRegulated,paste(dataName,"_UpRegulated",sep = ""),outputDir)
dfSave(downRegulated,paste(dataName,"_DownRegulated",sep = ""),outputDir)
dfSave(resultInfo,paste(dataName,"_ResultInfo",sep = ""),outputDir)
