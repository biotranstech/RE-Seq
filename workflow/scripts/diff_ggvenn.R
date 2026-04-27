#args1	diff_stat.csv
#args2	group_vs_group

library(ggplot2)
library(dplyr)
library(ggvenn)
library(VennDiagram)

args <- commandArgs(trailingOnly=TRUE)

theme_style <- function(
                        )
{
######
      theme(
            text = element_text(family = "serif", face = "bold",size = 10),
            plot.title = element_text(family = "serif", size = 10, face = "bold",hjust = 0.5,colour = "black"),
            legend.title = element_text(family = "serif", size = 10,face = "bold",colour = "black"),
            legend.text = element_text(family = "serif", size = 10,face = "bold",colour = "black"),
            plot.subtitle = element_text(family = "serif",face = "bold", size = 10),
            plot.caption = element_text(family = "serif", face = "bold",size = 14),
      panel.grid.major = element_blank(),  # 移除大网格线
      panel.grid.minor = element_blank())   # 移除小网格线))
}

data <- read.table(args[1], header=T)
data2 <- read.table(args[2],header=T)

a = list(group1=data$genes,group2=data2$genes)
names(a)[1] = args[3]
names(a)[2] = args[4]
#p <- ggvenn(a, fill_color=c('#FFFFCC','#CCFFFF'), stroke_size = 0.5, set_name_size = 4) + theme_style()
p <- ggvenn(a, fill_color=c('#F8766D','#00BFC4'), stroke_size = 0.5, set_name_size = 8, text_size = 8) 
ggsave(p,file=args[5],units = "cm",height = 16,width = 20,dpi = 300)

