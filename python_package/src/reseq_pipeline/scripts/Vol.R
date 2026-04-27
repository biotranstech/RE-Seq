#Rscript - diff_work_path output_path method group_file
#args1      差异分析路径
#args2      输出路径:./03_Edit/06_Edit_Volcanogram/
#args3      差异分析方法:Limma; DESeq2; Permutation
#args4      分组文件

library(stringr)
library(ggplot2)
library(glue)

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

Volcano_Statistics <- function(
                                    out_path,
                                    DEG_list,
                                    method,
                                    comparisons)
{
      #######初始统计样本的类型和频率信息
      print("Edit_Volcano_Statistics.start")
      #Edit_Volcano_information <- file.path(out_path, "03_Edit", str_glue("06_Edit_Volcanogram"))
      Edit_Volcano_information <- out_path
      
      if (!dir.exists(Edit_Volcano_information)) {
            dir.create(Edit_Volcano_information, recursive = TRUE)
      }
      
      method <- method
      Volcano <- DEG_list
      comparisons <- comparisons
      #DEG_frame_names <- sapply(Volcano, attr, "name")
      for ( i in comparisons) {
	    i <- gsub("\\+", "_", i)
            # 使用 get() 函数获取名为 i 的数据框
            DEG_frame <- get("i")
            DEG <- Volcano[[i]]
	    DEG <- DEG[!is.na(DEG$regulate),]
            x_breaks <- seq(-40, 40, by = 5)
            y_breaks <- c(seq(0, 100, by = 2), 1.3)
            p0=ggplot(DEG,aes(x=logFC,y=-log10(P.Value)))+ #x轴logFC,y轴adj.p.value
                  geom_point(alpha=0.7,size=2,aes(color=regulate))+ #点的透明度，大小
                  scale_color_manual(values = c("blue", "grey", "red"))+ #点的颜色
                  geom_vline(xintercept = c(-1,1),lty=6,col ="black",lwd=1.3)+ #logFC分界线
                  geom_hline(yintercept=-log10(0.05),lty=6,col = "black",lwd=1.3)+ #adj.p.val分界线
                  theme_bw() + #火山图绘制
                  labs(title = "", x = "Log2(Fold Change)", y = "-Log10(Pvalue)", fill = "Sig",color="Sig")+theme_style()+
                  scale_x_continuous(breaks = x_breaks) +
                  scale_y_continuous(breaks = y_breaks)
            x_breaks <- sort(unique(c(ggplot_build(p0)$layout$panel_params[[1]]$x$breaks, -1, 1)))
            y_breaks <- sort(unique(c(ggplot_build(p0)$layout$panel_params[[1]]$y$breaks, 1.3)))
            
            p1 <- p0 + scale_x_continuous(breaks = x_breaks)+
                  scale_y_continuous(breaks = y_breaks)+
                  theme_style()
            output_file <- glue("{Edit_Volcano_information}/{i}_{method}_Vol.pdf")
            ggsave(p1,file=output_file,units = "cm",height = 16,width = 20,dpi = 300)
            output_file2 <- glue("{Edit_Volcano_information}/{i}_{method}_Vol.png")
            ggsave(p1,file=output_file2,units = "cm",height = 16,width = 20,dpi = 300)}

      print("Edit_Volcano_Statistics.done")
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
      for (cmp in comparisons){
	    cmp <- gsub("\\+", "_", cmp)
            DEG <- read.csv(paste0(work_path,"/",cmp,"_",method,".csv"), header = TRUE, row.names = 1)
            DEG_list[[cmp]] <- DEG
      }
      return(DEG_list)
}

comparisons <- read_group(group_file = args[4])

DEG_list <- read_list(
      work_path = args[1],
      method = args[3],
      comparisons = comparisons)


Volcano_Statistics(
      out_path = args[2],
      DEG_list = DEG_list,
      method = args[3],
      comparisons = comparisons)
