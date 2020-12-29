# 功能：整理同一文件夹下的所有.txt基因表达数据
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
toNum = function(data){
    rnames = rownames(data)
    data = data.frame(lapply(data,as.numeric))
    rownames(data) = rnames
    return(data)
}
dataCollection = function(fileName,outputDir){
    # 读取数据
    data = read.table(fileName,sep = "\t",header = T,row.names = 1)
    data = toNum(data)
    mySave(data,fileName,outputDir)
    # RSEM表达值进行log转换
    if(fileName == "RSEM.txt"){
        data = log(data +1 )
        mySave(data,gsub(".txt","Log.txt",fileName),outputDir)
    }
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
outputDir = gsub(".R$","/",wholeName)


# 2数据整理
# 当前目录文件
fileNames = dir()[grep(".txt",dir())]
# 整理
lapply(fileNames,dataCollection,outputDir)


