# 功能：使用DESeq2进行差异分析
    # 输入：
        # rawCounts.RData:基因在行，样本在列；末尾01为病变样本，11为正常样本。
    # 输出：
        # rawCounts.csv:所有基因的LFC
        # deGenes.csv:差异基因的LFC

#函数
#0对应肿瘤，1对应正常，取自编号倒数第二位
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
DEFilter = function(res,lfc = 1,p = 0.05){
    res = data.frame(res)
    deGenes = subset(res,abs(log2FoldChange) > lfc & padj < p)
    deGenes
}
volcanicPlot = function(res,fileName,outputDir){
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
    ggsave(paste(outputDir,gsub(".RData","LCF",fileName),".png",sep = ""),width = 10,height = 8)
}


# 1环境
wholeName = parent.frame(2)$filename
# 文件名和路径名称
pathName = dirname(wholeName)
scriptName = basename(wholeName)
# 设置路径
setwd(pathName)
# 文件夹名称
dirName = gsub(".R","",scriptName)
# 创建脚本文件同名文件夹
if(!dir.exists(dirName)){
    dir.create(dirName)
}
# 输出文件路径
outputDir = gsub(".R$","/",wholeName)
library(DESeq2)
library(ggplot2)


# 2分析
# 数据导入
fileName = dir()[grep(".RData",dir())]
load(fileName)
# 样本分组信息
colData = groupInfo(rawCounts)
# DESeq计算LFC及对应P值
res = DECalculation(rawCounts,colData)
# 差异基因筛选
deGenes = DEFilter(res,1,.05)


# 3结果输出
write.csv(res,file = paste(outputDir,gsub(".RData","LCF",fileName),".csv",sep = ""))
write.csv(deGenes,file = paste(outputDir,"deGenesLCF.csv",sep = ""))
deGeneList = row.names(deGenes)
save(deGeneList,file = paste(outputDir,"deGeneList.RData",sep = ""))


# 4图片绘制
volcanicPlot(res,fileName,outputDir)

