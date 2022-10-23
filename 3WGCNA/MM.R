# 功能：
    # 1.计算MM
# 输入：
    # NAGGreyRemoved.csv：表达数据
    # NAGModuleSetList.RData：基因模块列表
    # NAGNetwork.RData：共表达网络
# 输出：
    # NAGMM.RData：MM值
    # NAGMM.csv：MM值
    # NAGMM.xlsx：MM值


# 0工具包
rm(list = ls())
library(modifyTools)
library(WGCNA)
library(ggplot2)
# 0函数
MMCompute = function(dataExpr,moduleSetList,mergedMEs){
    MM = data.frame()
    for(i in 1:length(moduleSetList)){
        modName = names(moduleSetList)[i]
        gene = moduleSetList[[i]]
        print(match(gene,colnames(dataExpr)))
        exprMod = dataExpr[,match(gene,colnames(dataExpr))]
        MEsMod = mergedMEs[,match(modName,colnames(mergedMEs))]
        MMMod = as.data.frame(cor(exprMod, MEsMod, use = "p"))
        MMPvalueMod = as.data.frame(corPvalueStudent(as.matrix(MMMod), dim(exprMod)[1]))
        modNameDF = data.frame(modeName = rep(modName,dim(MMMod)[1]))
        MMTemp = cbind(modNameDF,MMMod,MMPvalueMod)
        MM = rbind(MM,MMTemp)
    }
    names(MM) = c("ModName","MM","MMPvalue")
    return(MM)
}


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()
print("Input dataName:")
dataName = readline()
# 数据导入
dataExpr = read.csv(fileNames[grep("GreyRemoved.csv$",fileNames)],check.names = F,row.names = 1)
# dataClin = read.csv(fileNames[grep("^Clin.+csv$",fileNames)],check.names = F,row.names = 1)
load(fileNames[grep("Network.RData$",fileNames)])
load(fileNames[grep("SetList.RData$",fileNames)])


# 3处理
# 数据格式整理
# colnames(dataExpr) = getGeneID(colnames(dataExpr))
# 计算
MM = MMCompute(dataExpr,moduleSetList,mergedMEs)
# 筛选
MMSelected = subset(MM,MM>0.8 & MMPvalue < 0.05)
# 绘图
MM$threshold = 
    as.factor(ifelse(MM$MMPvalue < 0.05 & abs(MM$MM) > 0.8, 'Signifi','NoSignifi'))
MM$threshold = 
    factor(MM$threshold,levels = c("NoSignifi","Signifi"),ordered = T)
ggplot(data = MM, aes(y = -log(MMPvalue), x = MM,color = threshold)) +
           geom_point(alpha=0.5, size=2,) +
           scale_color_manual(values = c("grey","red")) +
           xlab("MM")+
           ylab("Pvalue")+
           labs(title = paste("Module Membership ( ",dataName," )",sep = "")) +
           geom_hline(yintercept = -log10(0.05),lty = 4,col = "black",lwd = 0.8)+
           geom_vline(xintercept = c(0.8),lty = 4,col = "black",lwd = 0.8)+
           geom_vline(xintercept = c(-0.8),lty = 4,col = "black",lwd = 0.8)+
           theme(panel.background = element_blank(),
                 plot.title = element_text(size = 35,hjust = 0.5),
                 axis.line = element_line(colour = "black"),
                 axis.text = element_text(size = 25),
                 axis.title = element_text(size = 25),
                 legend.key = element_rect(fill = "white"),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 25),
                 legend.position = c(0.9,0.9),
                 legend.key.size = unit(1,"cm"))
ggsave(file = paste(outputDir,dataName,"_MMPlot.png",sep = ""),width = 10,height = 8)


# 4保存
dfSave(MM,paste(dataName,"_MM",sep = ""),outputDir)
dfSave(MMSelected,paste(dataName,"_MMSelected",sep = ""),outputDir)