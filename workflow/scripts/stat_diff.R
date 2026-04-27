#args1	diff_stat.csv
#args2	group_vs_group

library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)

theme_style <- function(
                        )
{
######
      theme(axis.title.x = element_text(family = "serif", size = 16,face = "bold",colour = "black"),
            axis.text.x = element_text(family = "serif", size = 16,face = "bold",colour = "black"),
            axis.title.y = element_text(family = "serif", size = 16,face = "bold",colour = "black"),
            axis.text.y = element_text(family = "serif", size = 16,face = "bold",colour = "black"),
            text = element_text(family = "serif", face = "bold",size = 16),
            plot.title = element_text(family = "serif", size = 16, face = "bold",hjust = 0.5,colour = "black"),
            legend.title = element_text(family = "serif", size = 16,face = "bold",colour = "black"),
            legend.text = element_text(family = "serif", size = 16,face = "bold",colour = "black"),
            plot.subtitle = element_text(family = "serif",face = "bold", size = 16),
            plot.caption = element_text(family = "serif", face = "bold",size = 20),
      panel.grid.major = element_blank(),  # 移除大网格线
      panel.grid.minor = element_blank())   # 移除小网格线))
}


data <- read.csv(args[1], header=T)

data <- data[data$regulate %in% c('Down', 'Up'),]
a <- as.data.frame(table(data$regulate))
colnames(a) = c("Sig","num")
print(a)
aa = strsplit(args[2],'/')
group = aa[[1]][length(aa[[1]])]

p1 <- ggplot(a, aes(x=Sig,y=num,fill=Sig))+
	geom_bar(stat ="identity",width = 0.6,position = position_dodge(width=0.8))+
	theme_bw() +
	labs(title = "", x = group, y = "Number of Genes")+
	geom_text(aes(label = num),position=position_dodge(width = 0.9),size = 4,vjust = -0.25)+
	theme_style()


ggsave(p1,file=paste0(args[2],".pdf"),units = "cm",height = 16,width = 20,dpi = 300)
ggsave(p1,file=paste0(args[2],".png"),units = "cm",height = 16,width = 20,dpi = 300)
