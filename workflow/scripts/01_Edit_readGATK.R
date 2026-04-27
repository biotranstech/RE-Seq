#Rscript - out_path gatk_file
#args1	输出路径 04_Down_analysis/01_REseq_detection
#args2	GATK原始文件

library(stringr)
library(ggplot2)
library(dplyr)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)

theme_style <- function(
                        )
{
######
      theme(axis.title.x = element_text(family = "serif", size = 4,face = "bold",colour = "black"),
            axis.text.x = element_text(family = "serif", size = 4,face = "bold",colour = "black"),
            axis.title.y = element_text(family = "serif", size = 4,face = "bold",colour = "black"),
            axis.text.y = element_text(family = "serif", size = 4,face = "bold",colour = "black"),
            text = element_text(family = "serif", face = "bold",size = 4),
            plot.title = element_text(family = "serif", size = 4, face = "bold",hjust = 0.5,colour = "black"),
            legend.title = element_text(family = "serif", size = 4,face = "bold",colour = "black"),
            legend.text = element_text(family = "serif", size = 4,face = "bold",colour = "black"),
	    legend.key.size = unit(6, "pt"),
            plot.subtitle = element_text(family = "serif",face = "bold", size = 4),
            plot.caption = element_text(family = "serif", face = "bold",size = 6),
      panel.grid.major = element_blank(),  # 移除大网格线
      panel.grid.minor = element_blank())   # 移除小网格线))
}

Edit_Pile_maf_Statistics <- function(
                          out_path,
                          file_names)
{
#######初始统计样本的类型和频率信息
      print("Edit_Pile_maf_Statistics.start")
      #Edit_information <- file.path(out_path, "03_Edit", str_glue("01_Edit_type"))
      Edit_information <- out_path
      if (!dir.exists(Edit_information)) {
            dir.create(Edit_information, recursive = TRUE)
      }
      Initial_GATK_information <- str_glue("{file_names}")
      GATK_information <- read.csv(Initial_GATK_information,check.names = FALSE)
      GATK_information$maf <- GATK_information$MutationFrequency
      GATK_information$class <- paste0(GATK_information$Ref,"->",GATK_information$Alt)
      GATK_information$Loci_message <- paste0(GATK_information$Chr,"_",GATK_information$Pos,"_",GATK_information$Gene.refGene)
      write.csv(GATK_information, file = str_glue('{Edit_information}/GATK_information.csv'),row.names = F)
      write.table(GATK_information, file = str_glue('{Edit_information}/GATK_information.xls'),row.names = F)
      
      #colors1<-c("#934b43","#984EA3","cyan3","brown3","burlywood2","darkorange","chartreuse4","darkseagreen","dodgerblue4","gold","darkolivegreen3","grey50") 
      sample_num = length(unique(GATK_information$Sample))
      if (sample_num <= 12){
            colors1<-c("#f3a6a5","#FF7F00","#FFC125","#FFE4C4","#BDB76B","#4cb48a","#50b3bb","#1E90FF","#c6e2ff","#AB82FF","#ae62a4","#FFBBFF") 
      }else{
            colors1<-colorRampPalette(c("blue", "red"))(sample_num) 
      }
      colors2 <-c("#f3a6a5","#FF7F00","#FFC125","#FFE4C4","#BDB76B","#4cb48a","#50b3bb","#1E90FF","#c6e2ff","#AB82FF","#ae62a4","#FFBBFF","#FFD39B","#98F5FF","#BDB76B") 

      Frequency_Plot = ggplot(GATK_information, aes(x = maf, fill = Sample)) +
            geom_histogram(bins = 50) +
            facet_wrap(~ Sample)+theme_bw() +
            labs(title = "", x = "Editing frequency", y = "Frequency", fill = "Group")+
            theme_style()+
            scale_fill_manual(values=colors1)
      if (sample_num > 20){
            Frequency_Plot = Frequency_Plot + guides(fill = FALSE)
      }
      Frequency_Plot
      ggsave(Frequency_Plot,file=paste0(Edit_information,"/Frequency_Plot.pdf"),units = "cm",height = 6,width = 8,dpi = 300)
      ggsave(Frequency_Plot,file=paste0(Edit_information,"/Frequency_Plot.png"),units = "cm",height = 6,width = 8,dpi = 300)
      
      Type_Plot = ggplot(GATK_information, aes(x = Sample, fill = class, y = ..count..)) +
            geom_bar(position = "fill", stat = "count",width = 0.8)+
            scale_y_continuous(labels = scales::percent)+
            theme_bw() +
            labs(title = "", x = "Sample", y = "Percentage", fill = "Group")+ 
            theme_style()+
            theme(axis.text.x=element_text(angle = 45,hjust = 1))+
            scale_fill_manual(values=colors2)
      Type_Plot
      ggsave(Type_Plot,file=paste0(Edit_information,"/Type_Plot.pdf"),units = "cm",height = 6,width = 8,dpi = 300)
      ggsave(Type_Plot,file=paste0(Edit_information,"/Type_Plot.png"),units = "cm",height = 6,width = 8,dpi = 300)
      
      percentage_data <- GATK_information %>%
            dplyr::count(Sample, class) %>%  # 计算每个样本中每个类别的计数
            dplyr::group_by(Sample) %>%      # 按样本分组
            dplyr::mutate(percentage = n / sum(n) * 100) %>%  # 计算百分比
            ungroup()  # 取消分组
      matrix <- dcast(percentage_data, Sample ~ class, value.var = "percentage")
      write.csv(matrix,file=paste0(Edit_information,"/percentage_class.csv"),row.names=FALSE)
      write.table(matrix,file=paste0(Edit_information,"/percentage_class.xls"),row.names=FALSE)
      
      Function_Plot = ggplot(GATK_information, aes(x = Sample, fill = ExonicFunc.refGene, y = ..count..)) +
            geom_bar(position = "fill", stat = "count",width = 0.8)+
            scale_y_continuous(labels = scales::percent)+
            theme_bw() +
            labs(title = "", x = "Sample", y = "Percentage", fill = "Group")+ 
            theme_style()+
            theme(axis.text.x=element_text(angle = 45,hjust = 1))+
            scale_fill_manual(values=colors2)
      Function_Plot
      ggsave(Function_Plot,file=paste0(Edit_information,"/Function_Plot.pdf"),units = "cm",height = 6,width = 8,dpi = 300)
      ggsave(Function_Plot,file=paste0(Edit_information,"/Function_Plot.png"),units = "cm",height = 6,width = 8,dpi = 300)
      
      percentage_data <- GATK_information %>%
            dplyr::count(Sample, ExonicFunc.refGene) %>%  # 计算每个样本中每个类别的计数
            dplyr::group_by(Sample) %>%      # 按样本分组
            dplyr::mutate(percentage = n / sum(n) * 100) %>%  # 计算百分比
            ungroup()  # 取消分组
      matrix <- dcast(percentage_data, Sample ~ ExonicFunc.refGene, value.var = "percentage")
      matrix[is.na(matrix)] <- 0
      write.csv(matrix,file=paste0(Edit_information,"/percentage_ExonicFunc.refGene.csv"),row.names=FALSE)
      write.table(matrix,file=paste0(Edit_information,"/percentage_ExonicFunc.refGene.xls"),row.names=FALSE)

      Exonic_Plot = ggplot(GATK_information, aes(x = Sample, fill = Func.refGene, y = ..count..)) +
            geom_bar(position = "fill", stat = "count",width = 0.8)+
            scale_y_continuous(labels = scales::percent)+
            theme_bw() +
            labs(title = "", x = "Sample", y = "Percentage", fill = "Group")+ 
            theme_style()+
            theme(axis.text.x=element_text(angle = 45,hjust = 1))+
            scale_fill_manual(values=colors2)
      Exonic_Plot
      ggsave(Exonic_Plot,file=paste0(Edit_information,"/Exonic_Plot.pdf"),units = "cm",height = 6,width = 8,dpi = 300)
      ggsave(Exonic_Plot,file=paste0(Edit_information,"/Exonic_Plot.png"),units = "cm",height = 6,width = 8,dpi = 300)
      
      table_ExonicFunc.refGene <- table(GATK_information$Sample, GATK_information$ExonicFunc.refGene)
      write.csv(table_ExonicFunc.refGene, file = str_glue('{Edit_information}/table_ExonicFunc.refGene.csv'))
      write.table(table_ExonicFunc.refGene, file = str_glue('{Edit_information}/table_ExonicFunc.refGene.xls'), sep = "\t")
      table_Func.refGene <- table(GATK_information$Sample, GATK_information$Func.refGene)
      write.csv(table_Func.refGene, file = str_glue('{Edit_information}/table_Func.refGene.csv'))
      write.table(table_Func.refGene, file = str_glue('{Edit_information}/table_Func.refGene.xls'), sep = "\t")
      
      table_ALL_Func.refGene <- as.data.frame(table(GATK_information$Sample, GATK_information$Func.refGene, GATK_information$class))
      colnames(table_ALL_Func.refGene) <- c("Sample","Func.refGene","Type","Freq")
      write.csv(table_ALL_Func.refGene, file = str_glue('{Edit_information}/table_ALL_Func.refGene.csv'),row.names = F)
      write.table(table_ALL_Func.refGene, file = str_glue('{Edit_information}/table_ALL_Func.refGene.xls'),row.names = F, sep = "\t")
      
      table_ALL_ExonicFunc.refGene <- as.data.frame(table(GATK_information$Sample, GATK_information$ExonicFunc.refGene, GATK_information$class))
      colnames(table_ALL_ExonicFunc.refGene) <- c("Sample","ExonicFunc.refGene","Type","Freq")
      write.csv(table_ALL_ExonicFunc.refGene, file = str_glue('{Edit_information}/table_ALL_ExonicFunc.refGene.csv'),row.names = F)
      write.table(table_ALL_ExonicFunc.refGene, file = str_glue('{Edit_information}/table_ALL_ExonicFunc.refGene.xls'),row.names = F, sep = "\t")
      
      #file.create(paste0(Edit_information,"/Edit_Pile_maf_Statistics.done"))
      print("Edit_Pile_maf_Statistics.done")
      return(GATK_information)
}

GATK_information <- Edit_Pile_maf_Statistics(
      out_path = args[1],
      file_names = args[2]
      )
