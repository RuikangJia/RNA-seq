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
library(DESeq2)
library(ggplot2)
library(modifyTools)

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
volcanicPlot = function(drawData,dataName,outputDir,LFC = 1,p = 0.05){
    # 去除NA
    drawData = subset(drawData,complete.cases(drawData))
    # p值区分不同基因，且因子化
    drawData$threshold = as.factor(ifelse(drawData$padj < p & abs(drawData$log2FoldChange) >= LFC, ifelse(drawData$log2FoldChange>= LFC ,'Up','Down'),'NoSignifi'))
    drawData$threshold = factor(drawData$threshold,levels = c("Down","NoSignifi","Up"),ordered = T)
    #ggplot具体绘图
    ggplot(data = drawData, aes(x = log2FoldChange, y = -log10(padj), color=threshold)) +
        geom_point(alpha=0.5, size=2,) +
        scale_color_manual(values=c("blue","grey","red")) +
        ylab("-log10(Padj)")+
        xlab("log2(FoldChange)")+
        geom_vline(xintercept = c(-1,1),lty = 4,col = "black",lwd = 0.8)+
        geom_hline(yintercept = -log10(0.05),lty = 4,col = "black",lwd = 0.8)+
        theme(panel.background = element_blank(),
              axis.line = element_line(colour = "black"), 
              axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              legend.key = element_rect(fill = "white"),
              legend.title = element_blank(),
              legend.text = element_text(size = 20),
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
# 数据名称
splitDir = strsplit(dirname(wholeName),split = "/")[[1]]
dataName = splitDir[length(splitDir)]

# 导入数据
data1 = read.csv(fileName1,check.names = F,row.names = 1)
data2 = read.csv(fileName2,check.names = F,row.names = 1)


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
# LFC排序
DEGSorted = DEG[order(-abs(DEG$log2FoldChange)),]
# 上下调
upRegulated = subset(DEGSorted,DEGSorted$log2FoldChange > 0)
downRegulated = subset(DEGSorted,DEGSorted$log2FoldChange < 0)
# 绘制火山图
volcanicPlot(calcuRes,paste(dataName,"volcanic",sep = ""),outputDir)


# 4保存变量
dfSave(calcuRes,paste(dataName,"calcuRes",sep = ""),outputDir)
dfSave(DEG,paste(dataName,"DEG",sep = ""),outputDir)
dfSave(DEGSorted,paste(dataName,"DEGSorted",sep = ""),outputDir)
dfSave(upRegulated,paste(dataName,"upRegulated",sep = ""),outputDir)
dfSave(downRegulated,paste(dataName,"downRegulated",sep = ""),outputDir)
