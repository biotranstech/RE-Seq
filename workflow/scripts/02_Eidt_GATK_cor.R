#Rscript - output_path GATKfile groupfile
#args1	输出路径	04_Down_analysis/02_REseq_PCA
#args2	第一步输出的GATK file
#args3	分组文件group file

library(stringr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(corrplot)
library(factoextra)
library(ggcorrplot)

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


Edit_correlation_Statistics<- function(
                            out_path,
                            GATK_information,
                            group_file)
{
######
      print("Edit_correlation_Statistics.start")
      #Edit_correlation <- file.path(out_path, "03_Edit", str_glue("02_Edit_Loci"))
      Edit_correlation <- out_path
      if (!dir.exists(Edit_correlation)) {
            dir.create(Edit_correlation, recursive = TRUE)
      }
      group <- read.table(group_file, header=T)
      group_num = length(unique(group$Species))
      #colors1<-c("#E87D72","#FF7F00","#FFC125","#BDB76B","#","#56BD8F","#56BCC2","#1E90FF","#C6E2FF","#AB82FF","#E057E0","#FFBBFF") 
      colors1<-c("#FFC125","#56BD8F","#AB82FF","#ABFF00","#33e6cc","#ff69b4") 
      colors2<-colors1[1:group_num]
      
      df <- GATK_information
      data <- df %>%
            dplyr::select(Sample, Pos, maf) %>%
            dplyr::collect()
      express_matrix <- dcast(data, Sample ~ Pos, value.var = "maf", fun.aggregate = mean)  ##三列信息构建表达数据库
      express_matrix[is.na(express_matrix)] <- 0
      rownames(express_matrix) <- express_matrix[, 1]
      species <- as.matrix(express_matrix$Sample)
      express_matrix <- express_matrix[, -1]
      express_matrix_scaled <- scale(express_matrix)
      df_pca <- prcomp(express_matrix_scaled)
      df_pcs <-data.frame(df_pca$x, species = species)
      df_pcs <- merge(df_pcs, group, by="species")
      df_pcs <- relocate(df_pcs, species, .before=Species)
      #df_pcs$Species <- gsub("-\\d+", "", df_pcs$species)
      percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
      percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
      PCA <- ggplot(df_pcs, aes(x=PC1, y=PC2, color=Species)) +  
            geom_point(size=6) +  
            xlab(percentage[1]) +   
            ylab(percentage[2]) +  
            stat_ellipse(level = 0.95, show.legend = FALSE) +  
            labs(title="", subtitle="") +  
            labs(color="PCA Group") +  
            theme_bw() +  
            geom_text(aes(label = species), family = "serif", size = 6, show.legend = FALSE) + # 确保这个图层不创建图例  
            scale_color_manual(values = colors2) + # 使用 scale_color_manual 来设置颜色  
            theme(legend.position = "right")+
            theme_style()
      PCA
      ggsave(PCA,file=paste0(Edit_correlation,"/PCA_Plot.pdf"),units = "cm",height = 16,width = 20,dpi = 300)      
      ggsave(PCA,file=paste0(Edit_correlation,"/PCA_Plot.png"),units = "cm",height = 16,width = 20,dpi = 300)      
      
      t_express_matrix <- t(express_matrix_scaled)
      M <- cor(t_express_matrix)
      P <-round(cor_pmat(t_express_matrix,method = "pearson"),3)
      
      m=par(no.readonly = TRUE) 
      col1=colorRampPalette(colors =c("red","white","darkgreen"),space="Lab")
      pdf(file = paste0(Edit_correlation, "/Cor_Plot.pdf"), width = 8, height = 8,family="Times") # 打开PDF设备  
      par(mfrow=c(1,1)) 
      corrplot(corr =M,p.mat = P,type="upper",
               tl.pos="lt",tl.col="black", 
               insig = "label_sig", sig.level = c(.01, .05),
               pch.cex=2,pch.col = "black",order = "AOE",col = col1(10))
      corrplot(corr = M,type="lower",add=TRUE,method="number",
               tl.pos="n",tl.col="black",tl.cex=2.0,
               col="black",diag=FALSE, cl.pos="n",pch.col = "black",
               number.cex = 0.8,order = "AOE")

      dev.off() 

      system(paste("/usr/bin/pdftoppm -png ", Edit_correlation, "/Cor_Plot.pdf ", Edit_correlation, "/Cor_Plot", sep = ""))
      system(paste("mv ", Edit_correlation, "/Cor_Plot-1.png ", Edit_correlation, "/Cor_Plot.png", sep = ""))
      #png
      #m=par(no.readonly = TRUE) 
      #col1=colorRampPalette(colors =c("red","white","darkgreen"),space="Lab")
      #pdf(file = paste0(Edit_correlation, "/Cor_Plot.png"), width = 16, height = 20,family="Times",dpi = 300) # 打开PDF设备  
      #par(mfrow=c(2,2)) 
      #corrplot(corr =M,p.mat = P,type="upper",
      #         tl.pos="lt",tl.col="black", 
      #         insig = "label_sig", sig.level = c(.01, .05),
      #         pch.cex=2,pch.col = "black",order = "AOE",col = col1(10))
      #corrplot(corr = M,type="lower",add=TRUE,method="number",
      #         tl.pos="n",tl.col="black",tl.cex=2.0,
      #         col="black",diag=FALSE, cl.pos="n",pch.col = "black",
      #         number.cex = 0.8,order = "AOE")

      #dev.off() 

      d<-dist(express_matrix_scaled)
      fit1<-hclust(d,method = "average")
      fviz_dend(fit1)
      Tree <- fviz_dend(fit1, k = 3, rect = TRUE, rect_fill = TRUE, type = "circular", rect_border = c("#00AFBB", "#E7B800", "#FC4E07"))
      Tree
      ggsave(Tree,file=paste0(Edit_correlation,"/Tree_Plot.pdf"),units = "cm",height = 16,width = 20,dpi = 300)
      ggsave(Tree,file=paste0(Edit_correlation,"/Tree_Plot.png"),units = "cm",height = 16,width = 20,dpi = 300)
      #file.create(paste0(Edit_correlation,"/Edit_correlation.done"))
      print("Edit_correlation_Statistics.done")
      return()
}

read_GATK <- function(file_names) {
      Initial_GATK_information <- str_glue("{file_names}")

      GATK_information <- read.csv(Initial_GATK_information,header = T,check.names = FALSE)
      return(GATK_information)
}

GATK_information <- read_GATK(file_names = args[2])

Edit_correlation_Statistics(
      out_path = args[1],
      GATK_information = GATK_information,
      group_file = args[3]
      )

