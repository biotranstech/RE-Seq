#Rscript - output_path TPM_file group_file
#args1  输出路径,初始Result路径output_path
#args2  计算得到的基因层面TPM file
#args3  分组文件group file

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

Edit_TMB_cor_Statistics<- function(
                            out_path,
                            Edit_TMB,
                            group_file)
{
      ######
      print("Edit_TMB_cor_Statistics.start")
      #Edit_TMB_correlation <- file.path(out_path, "03_Edit", str_glue("04_Edit_QC"))
      Edit_TMB_correlation <- out_path
      if (!dir.exists(Edit_TMB_correlation)) {
            dir.create(Edit_TMB_correlation, recursive = TRUE)
      }
      #colors1<-c("#E87D72","#FF7F00","#FFC125","#BDB76B","#","#56BD8F","#56BCC2","#1E90FF","#C6E2FF","#AB82FF","#E057E0","#FFBBFF") 
      #colors2<-c("#FFC125","#56BD8F","#AB82FF")  
      group <- read.table(group_file, header=T)
      group_num = length(unique(group$Species))
      colors1<-c("#FFC125","#56BD8F","#AB82FF","#ABFF00","#33e6cc","#ff69b4") 
      colors2<-colors1[1:group_num]

      express_matrix <- Edit_TMB
      species <- as.matrix(colnames(express_matrix))
      express_matrix <- express_matrix[which(rowSums(express_matrix)>0),]
      express_matrix <- t(express_matrix)
      express_matrix_scaled <- scale(express_matrix)

      express_matrix_scaled[is.na(express_matrix_scaled)] <- 0
      #express_matrix_scaled <- express_matrix_scaled[which(rowSums(express_matrix_scaled)>0),]
      #print(head(express_matrix_scaled))
      df_pca <- prcomp(express_matrix_scaled)
      df_pcs <-data.frame(df_pca$x, species = species)
      #df_pcs$Species <- gsub("-\\d+", "", df_pcs$species)
      df_pcs <- merge(df_pcs, group, by="species")
      df_pcs <- relocate(df_pcs, species, .before=Species)
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
      ggsave(PCA,file=paste0(Edit_TMB_correlation,"/PCA_Plot.pdf"),units = "cm",height = 16,width = 20,dpi = 300)
      ggsave(PCA,file=paste0(Edit_TMB_correlation,"/PCA_Plot.png"),units = "cm",height = 16,width = 20,dpi = 300)
      
      t_express_matrix <- t(express_matrix_scaled)
      M <- cor(t_express_matrix)
      P <-round(cor_pmat(t_express_matrix,method = "pearson"),3)
      
      m=par(no.readonly = TRUE) 
      col1=colorRampPalette(colors =c("red","white","darkgreen"),space="Lab")
      pdf(file = paste0(Edit_TMB_correlation, "/Cor_Plot.pdf"),width = 8, height = 8,family="Times")
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

      system(paste("/usr/bin/pdftoppm -png ", Edit_TMB_correlation, "/Cor_Plot.pdf ",Edit_TMB_correlation, "/Cor_Plot", sep = ""))
      system(paste("mv ", Edit_TMB_correlation, "/Cor_Plot-1.png ", Edit_TMB_correlation, "/Cor_Plot.png", sep = ""))

      d<-dist(express_matrix_scaled)
      fit1<-hclust(d,method = "average")
      fviz_dend(fit1)
      Tree <- fviz_dend(fit1, k = 3, rect = TRUE, rect_fill = TRUE, type = "circular", rect_border = c("#00AFBB", "#E7B800", "#FC4E07"))

      ggsave(Tree,file=paste0(Edit_TMB_correlation,"/Tree_Plot.pdf"),units = "cm",height = 16,width = 20)
      ggsave(Tree,file=paste0(Edit_TMB_correlation,"/Tree_Plot.png"),units = "cm",height = 16,width = 20)
      print("Edit_TMB_cor_Statistics.done")
}

read_TPM <- function(file_names) {
      TPM_info <- str_glue("{file_names}")

      TPM <- read.csv(TPM_info,header = T,check.names = FALSE,row.names = 1)
      return(TPM)
}

TPM <- read_TPM(file_names = args[2])


Edit_TMB_cor_Statistics(
      out_path = args[1],
      Edit_TMB = TPM,
      group_file = args[3])

