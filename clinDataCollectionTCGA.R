# 功能：整理同一文件夹下的所有.txt临床数据文件
    # 输入：
        # *.txt 
    # 输出：
        # *.Rdata:数据框
        # *.csv:excel文件
# 函数
mySave = function(data,fileName,outputDir){
    # 去后缀，得名称
    dataName = gsub(".txt","",fileName)
    write.csv(data,paste(outputDir,dataName,".csv",sep = ""))
    # 变量值命名变量
    assign(dataName,data)
    save(list = dataName,file = paste(outputDir,dataName,".RData",sep = ""))
}
dataCollection = function(fileName,outputDir){
    # 读取数据
    data = read.table(fileName,sep = "\t",header = T,row.names = 1)
    data = data.frame(t(data))
    rownames(data) = toupper(gsub(".","-",rownames(data),fixed = T))
    mySave(data,fileName,outputDir)
}


# 1环境、目录
# 获取当前脚本文件路径全称
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
outputDir = gsub(".R","/",wholeName)


# 2数据整理
# 当前目录文件
fileNames = dir()[grep(".txt",dir())]
# 整理
lapply(fileNames,dataCollection,outputDir)
