#Rscript - diff_work_path TPM_file output_path method group_file
#args1      差异分析路径
#args2      编辑统计结果路径
#args3      输出路径
#args4      差异分析方法
#args5      分组文件

library(stringr)
library(ggplot2)
library(glue)
library(pheatmap)
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

Complex_pheat <- function(DEG,
                          filepath,
                          i,
                          express_matrix,
                          method)
{
      ######
      filepath = filepath
      DEG = DEG
      i = i
      method = method
      express_matrix = express_matrix
      Sig_genes <- as.data.frame(DEG$genes[DEG$regulate == "Up" | DEG$regulate == "Down"])
      colnames(Sig_genes)[1] <- "genes"
      Sig_TMB <- merge(Sig_genes,express_matrix,by="genes")
      row.names(Sig_TMB) <- Sig_TMB$genes
      Sig_TMB <- select(Sig_TMB,-genes)
      Sig_TMB <- Sig_TMB[rowSums(Sig_TMB != 0) > 0, ]
      #print(head(Sig_TMB))
      pheat <- pheatmap(Sig_TMB, scale="row", border="white",
                        angle_col = 45,cellwidth = 15,
                        clustering_distance_rows = "minkowski",
                        clustering_method="complete",
                        show_rownames = F,
                        color = colorRampPalette(colors = c("#1E90FF", "grey", "#FF4500"))(10),
                        cluster_cols = F,treeheight_col = 20,
                        cluster_rows = T,treeheight_row = 20,
                        fontfamily= "serif",fontsize = 15)
      
      ggsave(plot = pheat,file=paste0(filepath,"/",i,"_",method,"_Pheatmap.pdf"),  units = "cm", height = 20, width = 20,dpi = 300)
      ggsave(plot = pheat,file=paste0(filepath,"/",i,"_",method,"_Pheatmap.png"),  units = "cm", height = 20, width = 20,dpi = 300)

}

Sig_Statistics <- function(
                             out_path,
                             Edit_TMB,
                             Edit_DEG_list,
                             method,
                             comparisons)
{
#####
      #Edit_Sig_information <- file.path(out_path, "03_Edit", str_glue("08_Edit_sig_Pheatmap"))
      Edit_Sig_information <- out_path
      if (!dir.exists(Edit_Sig_information)) {
            dir.create(Edit_Sig_information, recursive = TRUE)
      }

      #express_matrix <- as.data.frame(Edit_TMB)
      #express_matrix$genes <- rownames(express_matrix)
      
      method <- method
      Volcano <- Edit_DEG_list
      Stat <- Edit_TMB
      comparisons <- comparisons
      for ( i in comparisons) {
	    i <- gsub("\\+", "_", i)
            # 使用 get() 函数获取名为 i 的数据框
            Limma_frame <- get("i")
            DEG <- Volcano[[i]]
	    express_matrix <- Stat[[i]]
	    express_matrix$genes <- rownames(express_matrix)
            Complex_pheat(DEG = DEG,
                          filepath = Edit_Sig_information,
                          i = i,
                          express_matrix = express_matrix,
                          method = method)
      }
}

read_group <- function(
      group_file)
{
      group_te = read.table(group_file, header=TRUE)
      a = t(group_te)
      colnames(a) = group_te$species
      group_names = a[-1,]
      unique_groups <- sort(unique(group_names))
      comparisons <- combn(unique_groups, 2, FUN = function(x) paste(x[2], x[1], sep = "_vs_") )

      return(comparisons)
}

read_list <- function(
      method,
      work_path,
      comparisons)
{
      comparisons = comparisons
      DEG_list = list()
      method = method
      work_path = work_path
      print(comparisons)
      for (cmp in comparisons){
	    cmp <- gsub("\\+", "_", cmp)
            DEG <- read.csv(paste0(work_path,"/",cmp,"_",method,".csv"), header = TRUE)
            DEG_list[[cmp]] <- DEG
      }
      return(DEG_list)
}

read_stat <- function(
      stat_path,
      comparisons)
{
      comparisons = comparisons
      stat_list = list()
      stat_path = stat_path
      for (cmp in comparisons){
	    cmp2 <- gsub("\\+", "_", cmp)
            stat_l <- read.csv(paste0(stat_path,"/",cmp,".csv"), header = TRUE,check.names = FALSE, row.names = 1)
            stat_list[[cmp2]] <- stat_l
      }
      return(stat_list)
}

comparisons <- read_group(group_file = args[5])

DEG_list <- read_list(
      work_path = args[1],
      method = args[4],
      comparisons = comparisons)
stat_list <- read_stat(
      stat_path = args[2],
      comparisons = comparisons
)

Sig_Statistics(
      out_path = args[3],
      Edit_DEG_list = DEG_list,
      Edit_TMB = stat_list,
      method = args[4],
      comparisons = comparisons
      )
