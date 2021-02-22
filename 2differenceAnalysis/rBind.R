# 功能：
    # 主要内容：
# 输入：
    # 
# 输出：
    # 
    # 
# 工具包
library(modifyTools)
# 函数



# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()[grep(".csv$",dir())]
# 数据名
splitDir = strsplit(dirname(wholeName),split = "/")[[1]]
dataName = splitDir[length(splitDir)]
# 数据导入
dataList = lapply(fileNames,read.csv,check.names = F,row.names = 1)


# 3处理
rBindData = data.frame()
for(data in dataList){
    rBindData = rbind(rBindData,data)
}


# 4保存
dfSave(rBindData,dataName,outputDir)
