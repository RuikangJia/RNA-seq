# 功能：进行基因表达差异性分析
    # 主要内容：差异分析；筛选差异基因；绘制火山图
# 输入：
    # rawCount.RData:基因表达数据未标准化
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


# 0函数
envirPrepare = function(wholeName){
    # 文件名和路径名称
    pathName = dirname(wholeName)
    scriptName = basename(wholeName)
    # 设置路径
    setwd(pathName)
    # 文件夹名称
    dirName = gsub(".R$","",scriptName)
    # 创建脚本文件同名文件夹
    if(!dir.exists(dirName)){
        dir.create(dirName)
    }
    # 输出文件路径
    outputDir = gsub(".R","/",wholeName)
    return(outputDir)
}
dfSave = function(data,dataName,outputDir){
    # 保存excel
    write.csv(data,paste(outputDir,dataName,".csv",sep = ""))
    # 数据命名
    assign(dataName,data)
    # 保存RData
    save(list = dataName,file = paste(outputDir,dataName,".RData",sep = ""))
}
groupInfo = function(data){
    cName = names(data)
    group = as.factor(as.numeric(substr(cName,nchar(cName[1])-1,nchar(cName[1])-1)))
    colData = data.frame(row.names = cName, group_list=group)
    colData
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
volcanicPlot = function(drawData,dataName,outputDir){
    # 去除NA
    drawData = subset(drawData,complete.cases(drawData))
    # p值区分不同基因，且因子化
    drawData$threshold = as.factor(ifelse(drawData$padj < 0.05 & abs(drawData$log2FoldChange) >= 1, ifelse(drawData$log2FoldChange>= 1 ,'Up','Down'),'NoSignifi'))
    drawData$threshold = factor(drawData$threshold,levels = c("Down","NoSignifi","Up"),ordered = T)
    #ggplot具体绘图
    ggplot(data = drawData, aes(x = log2FoldChange, y = -log10(padj), color=threshold)) +
        geom_point(alpha=0.5, size=1) +
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
              legend.title = element_blank())
    ggsave(file = paste(outputDir,dataName,".png",sep = ""),width = 10,height = 8)
}


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileName = dir()[grep(".RData$",dir())]
# 导入RData数据
data = load(fileName)
data = get(data)


# 3差异分析
# 基因在行，样本在列
data = data.frame(t(data),check.names = F)
# 分组信息
groupInfo = groupInfo(data)
# 计算LFC
calcuRes = DECalculation(data,groupInfo)
calcuRes = data.frame(calcuRes)
# 筛选差异基因
DEG = DEGFilter(calcuRes)
# LFC排序
DEGSorted = DEG[order(-abs(DEG$log2FoldChange)),]
# 上下调
upRegulated = subset(DEGSorted,DEGSorted$log2FoldChange > 0)
downRegulated = subset(DEGSorted,DEGSorted$log2FoldChange < 0)
# 绘制火山图
volcanicPlot(calcuRes,"volcanic",outputDir)


# 4保存变量
dfSave(calcuRes,"calcuRes",outputDir)
dfSave(DEG,"DEG",outputDir)
dfSave(DEGSorted,"DEGSorted",outputDir)
dfSave(upRegulated,"upRegulated",outputDir)
dfSave(downRegulated,"downRegulated",outputDir)
