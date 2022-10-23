# 功能：准备生存分析需要的数据，基因表达值和临床数据
    # 主要内容：
# 输入：

# 输出：

# 工具包
rm(list = ls())
library(modifyTools)
library(writexl)
library(clusterProfiler)
library(dplyr)
# 0函数


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()
geneListFileName = fileNames[grepl("List.csv$",fileNames)]
RSEMFileName = fileNames[grepl("^RSEM",fileNames)]
clinFileName = fileNames[grepl("^clin",fileNames)]
# 数据导入
geneListInfo = read.csv(geneListFileName,check.names = F,row.names = 1)
RSEM = read.csv(RSEMFileName,check.names = F,row.names = 1)
clinInfo = read.csv(clinFileName,check.names = F,row.names = 1)


# 3处理
# 临床信息整理
time = rep(0,nrow(clinInfo))
for(i in 1:nrow(clinInfo)){
    if(is.na(clinInfo$days_to_death[i]) & is.na(clinInfo$days_to_last_followup[i])){
        time[i] = NA
    }else if (is.na(clinInfo$days_to_death[i]) & !is.na(clinInfo$days_to_last_followup[i])){
        time[i] = clinInfo$days_to_last_followup[i]
    }else{
        time[i] = clinInfo$days_to_death[i]
    }
    
}
clin = data.frame(time = time,status = clinInfo$vital_status,
                  age = clinInfo$years_to_birth,sex = clinInfo$gender)
rownames(clin) = rownames(clinInfo)
# 表达数据整理
#   目的基因列表
geneList = geneListInfo$GeneList
geneListBitr = bitr(geneList,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
#   提取基因列表对应的表达数据
expr = RSEM
names(expr) = getGeneID(names(expr))
expr = expr[,match(geneListBitr$ENTREZID,names(expr))]

rownames(expr) = gsub("-01$","",rownames(expr))
names(expr) = geneListBitr$SYMBOL[match(names(expr),geneListBitr$ENTREZID)]
# 数据合并
clin$Name = rownames(clin)
expr$Name = rownames(expr)
data = inner_join(expr,clin,by = "Name")
#   整理
rownames(data) = data$Name
data$Name = NULL
data$sex[data$sex == "male"] = 1 
data$sex[data$sex != 1] = 2 
data$sex = as.numeric(data$sex)


# 4保存
dfSave(data,"sruvivalData",outputDir)