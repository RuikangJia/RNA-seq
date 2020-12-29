# ���ܣ������ݽ���CoxPH���������������������
    # ���룺
        # �����ļ�:����time,status,mainVar��������������
    # �����
        # cox.png��ģ��ɭ��ͼ
        # cox.RData��ģ����Ϣ����
        # coxdfbeta.png����Ⱥֵ���ͼ
        # coxPH.png��PH�������ͼ
        # coxLinear.png�����Լ���ͼ
        # survKurve.png����������
# ���߰�
library(survival)
library(survminer)
library(forestplot)
library(cowplot)
library(dplyr)
# ����
# ���ݱ������ɹ�ʽ
generateFormula = function(adjustVars){
    # �Ӻ����ӱ���
    adjustVarsJoin = paste(adjustVars,collapse = "+")
    # ���˺ż�ӳ���ϵ�Ĺ���
    formulaString = paste("Surv(time,event = status)~ ",adjustVarsJoin)
    # �ַ���ת��Ϊ��ʽ��ʽ
    myFormula = as.formula(formulaString)
    return(myFormula)
}
generateFormula2 = function(adjustVar){
    # ���˺ż�ӳ���ϵ�Ĺ���
    formulaString = paste("Surv(time,event = status)~ ",adjustVar,sep = "")
    # �ַ���ת��Ϊ��ʽ��ʽ
    myFormula = as.formula(formulaString)
    return(myFormula)
}
# Coxģ����Ϣ��ɭ��ͼ
coxModelGraphical = function(cox,filename){
    # ɭ��ͼ.png
    ggforest(cox,fontsize = 1.1,data = data,cpositions = c(0.02, 0.2, 0.38))
    ggsave(paste(outputDir,filename,".png",sep = ""),width = 10,height = 8,scale = 0.7)
    # ģ��������Ϣ.RData
    save(cox,file = paste(outputDir,filename,".RData",sep = ""))
}
# coxģ�ͼ��
coxModelCheck = function(cox,filename,data,adjustVars){
    # 1���PH���裺����������Ӱ����ʱ���޹�
    PHTest = cox.zph(cox)
    multiP = c()
    for(var in adjustVars){
        p = ggcoxzph(PHTest,var = var,font.tickslab = c(16),font.x = c(16),font.y =c(16),font.main = c(16))
        multiP = c(multiP,p)
    }
    # ���ͼ
    plot_grid(plotlist = multiP,greedy = F)
    ggsave(paste(outputDir,filename,"PH.png",sep = ""),width = 14,height = 8)
    # 2����쳣ֵ:dfbetaΪɾ���۲�ֵ���Ӧ�����ϵ�������ı仯
    ggcoxdiagnostics(cox,type = "dfbeta",font.tickslab = c(16),font.x = c(16),font.y =c(16),font.subtitle = c(16),font.legend = c(15)) + xlab("time") + 
    ggsave(paste(outputDir,filename,"dfbeta.png",sep = ""),width = 10,height = 8)
    # 3��������ֵ������Э����֮������Թ�ϵ
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
    # ��ȡ��time��status�ı���
    adjustVars = colnames(data)
    adjustVars = adjustVars[-which(adjustVars == "time" | adjustVars == "status")]
    # ����ģ��ѡȡ
    num = 1
    isBestModel = FALSE
    while(!isBestModel){
        # ���ɹ�ʽ
        f = generateFormula(adjustVars)
        # ����ģ��
        coxFit = coxph(f,data = data)
        # ģ����Ϣ��ɭ��ͼ
        coxModelGraphical(coxFit,paste("cox",num,sep = ""))
        # ģ�ͼ��
        coxModelCheck(coxFit,paste("cox",num,sep = ""),data,adjustVars)
        # ����ģ��ѡȡ���ж�ģ�ͱ�����������
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
    # ��ȡʹ�õı���
    data2 = data[,match(names(coxFit$coefficients),names(data))]
    # ������ˮƽ��Ϣ����������ƽ��ֵ1����ɢ�������ֵ0��Ŀ���������-1
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
    # ȥ���ظ���
    dataTwolevel = distinct(data2)
    # ��ȡHR��Pֵ
    para = summary(coxFit)$coefficients
    HR = round(para[match(varName,rownames(para)),match("exp(coef)",colnames(para))],5)
    PValue = round(para[match(varName,rownames(para)),match("Pr(>|z|)",colnames(para))],5)
    # �����������
    survFit = survfit(coxFit, newdata=dataTwolevel)
    # ��ͼ
    ggsurvplot(survFit,data = data,conf.int=T,
               legend.labs=c(paste(varName,"High",sep ="-"),paste(varName,"Low",sep = "-")),
               pval = paste("HR:",HR,"\nP:",PValue,sep = " "),
               legend.title = "group", font.legend = c(20),
               font.tickslab = c(20),
               font.x = c(20),font.y = c(20)
               )
    ggsave(paste(outputDir,"survKurve.png",sep = ""),width = 10,height = 8)
}


# 1����
wholeName = parent.frame(2)$filename
# �ļ�����·������
pathName = dirname(wholeName)
scriptName = basename(wholeName)
# ����·��
setwd(pathName)
# �ļ�������
dirName = gsub(".R","",scriptName)
# �����ű��ļ�ͬ���ļ���
if(!dir.exists(dirName)){
    dir.create(dirName)
}
# ����ļ�·��
outputDir = gsub(".R$","/",wholeName)


# 2���ݵ���
data(lung)
data = lung


# 3Cox����
coxFit = coxAnalysis(data)


# 4�������������߻���
# survPlot(data,coxFit,c(-1,1))
