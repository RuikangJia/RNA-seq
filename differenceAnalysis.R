# 功能：进行基因表达差异性分析
    # 主要内容：差异分析；筛选差异基因；绘制火山图
# 输入：
    # rawCount.RData:基因表达数据未标准化
# 输出：
    # deGene.csv:差异基因信息
    # deGene.RData:差异基因信息
    # DERes.csv:所有基因LFC计算结果
    # DERes.RData:所有基因LFC计算结果
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
DEFilter = function(DERes,lfc = 1,p = 0.05){
    DERes = data.frame(DERes)
    deGenes = subset(DERes,abs(log2FoldChange) > lfc & padj < p)
    return(deGenes)
}
volcanicPlot = function(res,dataName,outputDir){
    # 去除NA
    res = subset(res,complete.cases(res))
    # 构建绘图使用的数据集
    drawData = data.frame(compound = res@rownames,res@listData)
    # p值区分不同基因，且因子化
    drawData$threshold = as.factor(ifelse(drawData$padj < 0.05 & abs(drawData$log2FoldChange) >= 1, ifelse(drawData$log2FoldChange>= 1 ,'Up','Down'),'NoSignifi'))
    drawData$threshold = factor(drawData$threshold,levels = c("Down","NoSignifi","Up"),ordered = T)
    #ggplot具体绘图
    ggplot(data = drawData, aes(x = log2FoldChange, y = -log10(pvalue), color=threshold)) +
        geom_point(alpha=0.5, size=1) +
        scale_color_manual(values=c("blue","grey","red")) +
        ylab("-log10(Pvalue)")+
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
DERes = DECalculation(data,groupInfo)
# 筛选差异基因
deGenes = DEFilter(DERes,1,.05)
# 绘制火山图
volcanicPlot(DERes,"volcanic",outputDir)


# 4保存变量
dfSave(DERes,"DERes",outputDir)
dfSave(deGenes,"deGenes",outputDir)
