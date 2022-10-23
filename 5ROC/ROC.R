# 功能：进行ROC分析
    # 主要内容：曲线绘制；结果整理（带置信区间的AUC、阈值、特异性、敏感性）
# 输入：
    # hubGeneROCData_Else_HC.csv：ROC数据，名称设定control和case
# 输出：
    # AUCInfo_Else_HC.RData：ROC分类结果数据
    # AUCInfo_Else_HC.csv
    # AUCInfo_Else_HC.xlsx
    # ROC_Else_HC.png


# 0工具包
rm(list = ls())
library(modifyTools)
library(pROC)
library(ggplot2)
library(writexl)
options(digits=3)
# 0函数
getROC = function(input,output,data){
    roc = roc(data[,which(names(data) == output)],
              data[,which(names(data) == input)],
              auc = T,ci = T)
    return(roc)
}
getItem = function(data,itemName){
    # print(data[which(names(data) == itemName)])
    res = unlist(data[which(names(data) == itemName)])
    if(is.numeric(res)){
        res = round(res,3)
    }
    return(c(res))
}
getAUC = function(data){
    value = unlist(data)
    return(value[2])
    aucInterval = paste(value[2],"（",value[1],"-",value[3],"）",sep = "")
    return(aucInterval)
}
getAUCInterval = function(data){
    value = unlist(data)
    aucInterval = paste(value[1],"-",value[3],sep = "")
    return(aucInterval)
}


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()
dataFileName = fileNames[grep(".csv$",fileNames)]
dataName = gsub(".csv$","",dataFileName)
# 数据导入
data = read.csv(dataFileName,check.names = F,row.names = 1)
# 数据、预测因子、响应因子
predictor =  c(names(data)[-ncol(data)])
response = names(data)[ncol(data)]


# 3处理
# 分类信息
tempStr = strsplit(dataName,"_")
control = tempStr[[1]][2]
case = tempStr[[1]][3]
# 数据分类：分类值列可以含有NA，pROC会只选取有分类信息的样本
class = rep(0,nrow(data))
class[sapply(data$Type,grepl,case)] = 1
data$Class = class
# 进行ROC计算
listRoc = lapply(predictor,getROC,"Class",data)
names(listRoc) = predictor
# 绘制ROC曲线
ggroc(listRoc,legacy.axes = TRUE)+
    geom_line(size = 1) + labs(title = paste("ROC Curve ( Positive: ",case," )",sep = "")) +
    theme(panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 30,colour = "black"),
          axis.title = element_text(size = 30),
          plot.title = element_text(hjust = 0.5,size = 40),
          legend.position = "none"
          # legend.position =c(.79,.21),
          # legend.title = element_blank(),
          # legend.text = element_text(size = 10),
          # legend.key = element_rect(fill = "white")
          )+
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 colour = "grey", linetype = "dashed")
ggsave(width = 12,height = 10.5,filename = paste(outputDir,"ROC","_",control,"_",case,".png",sep = ""))
# 结果整理
Item = predictor
ci = lapply(listRoc,getItem,"ci")
AUC = sapply(ci, getAUC)
AUCInterval = sapply(ci, getAUCInterval)
sensitivities = lapply(listRoc,getItem,"sensitivities")
specificities = lapply(listRoc,getItem,"specificities")
thresholds = lapply(listRoc,getItem,"thresholds")
direction = unlist(lapply(listRoc,getItem,"direction"))
threshold = c()
sensitivitie = c()
specificitie = c()
for(i in 1:length(listRoc)){
    # 约登指数
    index = which.max(sensitivities[[i]] + specificities[[i]] -1)
    threshold = c(threshold,thresholds[[i]][index])
    sensitivitie = c(sensitivitie,sensitivities[[i]][index])
    specificitie = c(specificitie,specificities[[i]][index])
}
res = data.frame("Gene" = Item,"AUC" = AUC,
                 "AUC（95% CI）" = AUCInterval,
                 "Direction(Control)" = direction,
                 "Threshold" = threshold,"Sensitivity" = unlist(sensitivitie),
                 "Specificity" = unlist(specificitie),
                 check.names = F)


# 4保存
dfSave(res,paste("AUCInfo","_",control,"_",case,sep = ""),outputDir)

