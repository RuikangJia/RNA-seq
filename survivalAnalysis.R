# 功能：对数据进行CoxPH生存分析，绘制生存曲线
    # 输入：
        # 数据文件:带有time,status,mainVar和其它矫正变量
    # 输出：
        # cox.png：模型森林图
        # cox.RData：模型信息数据
        # coxdfbeta.png：离群值检测图
        # coxPH.png：PH假设检验图
        # coxLinear.png：线性检验图
        # survKurve.png：生存曲线
# 工具包
library(survival)
library(survminer)
library(forestplot)
library(cowplot)
library(dplyr)
# 函数
# 根据变量生成公式
generateFormula = function(adjustVars){
    # 加号连接变量
    adjustVarsJoin = paste(adjustVars,collapse = "+")
    # 波浪号及映射关系的构建
    formulaString = paste("Surv(time,event = status)~ ",adjustVarsJoin)
    # 字符串转换为公式格式
    myFormula = as.formula(formulaString)
    return(myFormula)
}
generateFormula2 = function(adjustVar){
    # 波浪号及映射关系的构建
    formulaString = paste("Surv(time,event = status)~ ",adjustVar,sep = "")
    # 字符串转换为公式格式
    myFormula = as.formula(formulaString)
    return(myFormula)
}
# Cox模型信息及森林图
coxModelGraphical = function(cox,filename){
    # 森林图.png
    ggforest(cox,fontsize = 1.1,data = data,cpositions = c(0.02, 0.2, 0.38))
    ggsave(paste(outputDir,filename,".png",sep = ""),width = 10,height = 8,scale = 0.7)
    # 模型数据信息.RData
    save(cox,file = paste(outputDir,filename,".RData",sep = ""))
}
# cox模型检测
coxModelCheck = function(cox,filename,data,adjustVars){
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
    ggsave(paste(outputDir,filename,"dfbeta.png",sep = ""),width = 10,height = 8)
    # 3对数风险值与线性协变量之间的线性关系
    multiP = c()
    data2 = na.omit(data)
    for(var in adjustVars){
        cox = coxph(generateFormula2(var),data = data2)
        p = ggcoxfunctional(cox,data = data2)
        multiP = c(multiP,p)
    }
    plot_grid(plotlist = multiP)
    ggsave(paste(outputDir,filename,"Linear.png",sep = ""),width = 10,height = 8)
}
coxAnalysis = function(data){
    # 获取除time和status的变量
    adjustVars = colnames(data)
    adjustVars = adjustVars[-which(adjustVars == "time" | adjustVars == "status")]
    # 最优模型选取
    num = 1
    isBestModel = FALSE
    while(!isBestModel){
        # 生成公式
        f = generateFormula(adjustVars)
        # 构建模型
        coxFit = coxph(f,data = data)
        # 模型信息及森林图
        coxModelGraphical(coxFit,paste("cox",num,sep = ""))
        # 模型检查
        coxModelCheck(coxFit,paste("cox",num,sep = ""),data,adjustVars)
        # 最优模型选取：判断模型变量的显著性
        coeffDF = summary(coxFit)[["coefficients"]]
        if(FALSE %in% (coeffDF[,5] < 0.05)){
            tempVars = rownames(coeffDF)[coeffDF[,5] < 0.05]
            select = c()
            for(var in adjustVars){
                select = c(select,TRUE %in% grepl(var,tempVars))
            }
            adjustVars = adjustVars[select]
            num = num +1
        }else{
            isBestModel = TRUE
        }
    }
    return(coxFit)
}
survPlot = function(data,coxFit,varInfo){
    # 提取使用的变量
    data2 = data[,match(names(coxFit$coefficients),names(data))]
    # 构建两水平信息：连续变量平均值1，离散变量最低值0，目标变量分组-1
    for(i in 1:length(varInfo)){
        if(varInfo[i] == -1){
            data2[,i] = ifelse(data2[,i] > mean(data2[,i],na.rm = T),"high","low")
            varName = names(data2)[i]
        }else if(varInfo[i] == 1){
            data2[,i] = mean(data2[,i],na.rm = T)
        }else if(varinfo[i] == 0){
            data2[,i] = min(data2[,i],na.rm = T)
        }
    }
    # 去除重复行
    dataTwolevel = distinct(data2)
    # 提取HR和P值
    para = summary(coxFit)$coefficients
    HR = round(para[match(varName,rownames(para)),match("exp(coef)",colnames(para))],5)
    PValue = round(para[match(varName,rownames(para)),match("Pr(>|z|)",colnames(para))],5)
    # 构建生存分析
    survFit = survfit(coxFit, newdata=dataTwolevel)
    # 绘图
    ggsurvplot(survFit,data = data,conf.int=T,
               legend.labs=c(paste(varName,"High",sep ="-"),paste(varName,"Low",sep = "-")),
               pval = paste("HR:",HR,"\nP:",PValue,sep = " "),
               legend.title = "group", font.legend = c(20),
               font.tickslab = c(20),
               font.x = c(20),font.y = c(20)
               )
    ggsave(paste(outputDir,"survKurve.png",sep = ""),width = 10,height = 8)
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


# 2数据导入
data(lung)
data = lung


# 3Cox分析
coxFit = coxAnalysis(data)


# 4矫正后生存曲线绘制
# survPlot(data,coxFit,c(-1,1))

