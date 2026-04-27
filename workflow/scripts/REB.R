#Rscript - GATK_file group.xls outpath
#args1	GATK_information.csv
#args2	group.xls
#args3	out path 

library(ggplot2)
library(dplyr)
library(ggpubr)

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



if (!dir.exists(args[3])) {
	dir.create(args[3], recursive = TRUE)
}

data <- read.csv(args[1], header=T)
group <- read.table(args[2],header=T)
colnames(group) = c("Sample","group")

#df_data = data[data$ExonicFunc.refGene %in% "nonsynonymous SNV", ]
df_data = data

stat = as.data.frame(table(df_data$Sample))
colnames(stat) = c("Sample","number")

if (args[4] == "hg38"){
	bed_sum = 60507855/1000000
} else if (args[4] == "mm10"){
	bed_sum = 136649262/1000000
} else if (args[4] == "rn6"){
	bed_sum = 136649262/1000000
}
#hg38_bed_sum <- 60507855
#mm10_bed_sum <- 136649262
stat <- stat %>% mutate(REB = number / bed_sum )
write.table(stat, file = paste0(args[3], "/REB.xls"), col.name = T, row.names = F, quote = F, sep = "\t")

g_data = merge(stat,group,by="Sample")

print(g_data)
my_comparisons = list(c("Control","Rop"),c("Control","Rop_Dex"),c("Rop","Rop_Dex"))
p=ggplot(g_data,aes(x=group,y=REB)) + 
  stat_boxplot(geom="errorbar",width=0.1,size=0.8) + 
  geom_boxplot(aes(fill=group),outlier.colour = NA)+
  stat_compare_means(comparisons = my_comparisons,method = "t.test")+
  theme_style()

ggsave(p,file=paste0(args[3],"/REB_boxplot.pdf"), units = "cm",height = 16,width = 20,dpi = 300)
ggsave(p,file=paste0(args[3],"/REB_boxplot.png"), units = "cm",height = 16,width = 20,dpi = 300)
