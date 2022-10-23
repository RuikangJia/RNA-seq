# 功能：进行生存分析
    # 主要内容：COX模型构建；模型检验；生存曲线绘制
# 输入：

# 输出：

# 工具包
rm(list = ls())
library(modifyTools)
library(writexl)
library(survival)
library(survminer)
library(forestplot)
library(cowplot)
library(dplyr)


# 0函数
generateFormula = function(adjustVars){
    # 加号连接变量
    adjustVarsJoin = paste(adjustVars,collapse = "+")
    # 波浪号及映射关系的构建
    formulaString = paste("Surv(time,event = status)~ ",adjustVarsJoin)
    # 字符串转换为公式格式
    myFormula = as.formula(formulaString)
    return(myFormula)
}
formCoxMod = function(data,func,filename){
    # 生成模型
    coxFit = coxph(func,data)
    # 模型数据信息.RData
    save(coxFit,file = paste(outputDir,filename,".RData",sep = ""))
    return(coxFit)
}
coxGraph = function(coxFit,filename){
    # 可视化
    # 数据保存：森林图.png
    ggforest(coxFit,fontsize = 1,data = data,cpositions = c(0.02, 0.2, 0.38))
    ggsave(paste(outputDir,filename,".png",sep = ""),width = 10,height = 8,scale = 0.7)
}
generateFormula2 = function(adjustVar){
    # 波浪号及映射关系的构建
    formulaString = paste("Surv(time,event = status)~ ",adjustVar,sep = "")
    # 字符串转换为公式格式
    myFormula = as.formula(formulaString)
    return(myFormula)
}
coxModelCheck = function(data,cox,adjustVars,filename){
    # 1检查PH假设：变量产生的影响与时间无关
    PHTest = cox.zph(cox)
    multiP = c()
    for(var in adjustVars){
        p = ggcoxzph(PHTest,var = var,font.tickslab = c(16),font.x = c(16),font.y =c(16),font.main = c(16))
        multiP = c(multiP,p)
    }
    # 组合图
    plot_grid(plotlist = multiP,greedy = F)
    ggsave(paste(outputDir,filename,"PH.png",sep = ""),width = 14,height = 8)
    # 2检查异常值:dfbeta为删除观察值后对应的相关系数发生的变化
    ggcoxdiagnostics(cox,type = "dfbeta",font.tickslab = c(16),font.x = c(16),font.y =c(16),font.subtitle = c(16),font.legend = c(15)) + xlab("time") + 
        ggsave(paste(outputDir,filename,"DfBeta.png",sep = ""),width = 10,height = 8)
    # 3对数风险值与线性协变量之间的线性关系
    # multiP = c()
    # data2 = na.omit(data)
    # for(var in adjustVars[1]){
    #     cox = coxph(generateFormula2(var),data = data2)
    #     p = ggcoxfunctional(cox,data = data2)
    #     multiP = c(multiP,p)
    # }
    # plot_grid(plotlist = multiP)
    # ggsave(paste(outputDir,filename,"Linear.png",sep = ""),width = 10,height = 8)
}
survPlot = function(data,coxFit,varInfo,fileName){
    # 提取使用的变量
    data2 = data[,match(gsub("High","",names(coxFit$coefficients)),names(data))]
    # 构建两水平信息：连续变量平均值1，离散变量最低值0，目标变量分组-1
    for(i in 1:length(varInfo)){
        if(varInfo[i] == -1){
            # data2[,i] = ifelse(data2[,i] > mean(data2[,i],na.rm = T),"High","Low")
            # varName = names(data2)[i]
            # str(data2)
        }else if(varInfo[i] == 1){
            data2[,i] = mean(data2[,i],na.rm = T)
        }else if(varInfo[i] == 0){
            data2[,i] = min(data2[,i],na.rm = T)
        }
    }
    # 去除重复行
    dataTwolevel = distinct(data2)
    rownames(dataTwolevel) = paste(names(dataTwolevel)[1],"-",dataTwolevel[,1],sep = "")
    dataTwolevel = dataTwolevel[order(dataTwolevel[,1],decreasing = T),]
    # 提取HR和P值
    para = summary(coxFit)$coefficients
    HR = round(para[match(varName,gsub("High","",rownames(para))),match("exp(coef)",colnames(para))],5)
    PValue = round(para[match(varName,gsub("High","",rownames(para))),match("Pr(>|z|)",colnames(para))],5)
    # 构建生存分析
    survFit = survfit(coxFit, newdata=dataTwolevel)
    # 绘图
    ggsurvplot(survFit,data = data2,
               # conf.int = F,
               # legend.labs=c(paste(varName,"Low",sep = "-"),
               #               paste(varName,"High",sep ="-")),
               # palette = c("blue","red"),
               pval = paste("HR:",HR,"\nP:",PValue,sep = " "),
               legend.title = "group", font.legend = c(20),
               font.tickslab = c(20),
               font.x = c(20),font.y = c(20)
    )
    ggsave(paste(outputDir,varName,"_SurvKurve.png",sep = ""),width = 10,height = 8)
    return(survFit)
}


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()
dataFileName = fileNames[grepl(".csv$",fileNames)]
# 数据导入
data = read.csv(dataFileName,check.names = F,row.names = 1)
# 变量分类
varTime = c("time","status")
varFixed = c("sex","age")
varRest = setdiff(names(data),union(varTime,varFixed))
# 分组：高表达对应2，低表达对应1
for(var in varRest){
    expr = data[,var]
    match = (expr > mean(expr))
    expr[match] = "High"
    expr[!match] = "Low"
    data[,var] = factor(expr,levels = c("Low","High"))
}
# 3处理
# 汇总结果信息
resInfo = c()
for(varName in varRest){
    # 3.1模型构建
    #   生成公式
    f = generateFormula(c(varName,varFixed))
    #   构建模型
    coxFit = formCoxMod(data,f,paste(varName,"_CoxMod",sep = ""))
    if(summary(coxFit)[["coefficients"]][1,5] < 0.05){
        #   绘制图形
        coxGraph(coxFit,paste(varName,"_Forest",sep = ""))
        # 3.2模型检验
        coxModelCheck(data,coxFit,c(varName,varFixed),paste(varName,"_CoxCheck",sep = ""))
        # 3.3生存曲线绘制
        a = survPlot(data,coxFit,c(-1,0,1),varName)

    }
    # 汇总结果信息
    coxInfo = summary(coxFit)
    resInfo = c(resInfo,varName,
                coxInfo$coefficients[1,1:3],
                coxInfo$conf.int[1,3:4],
                coxInfo$coefficients[1,4:5])

}


# 4保存
# 汇总结果整理
result = data.frame(matrix(resInfo,ncol = 8,byrow = T))
names(result) = c("Gene",
                  colnames(coxInfo$coefficients)[1:3],
                  colnames(coxInfo$conf.int)[3:4],
                  colnames(coxInfo$coefficients)[4:5] )
resultSig = subset(result,result$`Pr(>|z|)` < 0.05)
dfSave(result,"resultInfo",outputDir)
dfSave(resultSig,"resultSigInfo",outputDir)
