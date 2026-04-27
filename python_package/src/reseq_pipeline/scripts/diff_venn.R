#args1	diff_stat.csv
#args2	group_vs_group

library(ggplot2)
library(dplyr)
library(VennDiagram)

args <- commandArgs(trailingOnly=TRUE)


data <- read.table(args[1], header=T)
data2 <- read.table(args[2],header=T)

a = list(group1=data$genes,group2=data2$genes)
names(a)[1] = args[3]
names(a)[2] = args[4]
venn.diagram(
	x=a,
	scaled = F,
	filename = args[5],
	output = TRUE,

	# Output features
	imagetype="png",
	height = 700,
	width = 700,
	resolution = 300,
	compression = "lzw",# 压缩算法
	
	#
	alpha= 0.5, #透明度
	lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
	label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
	cex = 0.5, # 数字大小
	fontface = "bold",  # 字体粗细；加粗bold
	fill=c("#F8766D","#00BFC4"), # 填充色
	category.names = c(args[3], args[4]) , #标签名
	cat.dist = 0.02, # 标签距离圆圈的远近
	cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
	cat.cex = 0.5, #标签字体大小
	cat.fontface = "bold",  # 标签字体加粗
	cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
	cat.default.pos = "outer"  # 标签位置, outer内;text 外
)
