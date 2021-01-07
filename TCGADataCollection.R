# 功能：整理同一文件夹下的所有TCGA的.txt数据文件，包括临床和表达值
    # 主要内容：转置；行列名调整
# 输入：
    # *.txt 
# 输出：
    # *.Rdata:数据框
    # *.csv:excel文件


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
dataCollection = function(data){
    # 转置
    data = data.frame(t(data),check.names = F)
    # 修改样本名称
    rownames(data) = toupper(rownames(data))
    return(data)
}
dfSave = function(data,fileName,outputDir){
    dataName = gsub(".txt$","",fileName)
    # 保存excel
    write.csv(data,paste(outputDir,dataName,".csv",sep = ""))
    # 数据命名
    assign(dataName,data)
    # 保存RData
    save(list = dataName,file = paste(outputDir,dataName,".RData",sep = ""))
}


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()[grep(".txt$",dir())]
dataList = lapply(fileNames,read.table,sep = "\t",header = T,row.names = 1,check.names = F)


# 3整理
dataCollectedList = lapply(dataList,dataCollection)


# 4保存
for(i in 1:length(fileNames)){
    dfSave(dataCollectedList[[i]],fileNames[i],outputDir)
}
