#Rscript - RNA_diff_path Edit_diff_path out_path method group_file
#args1      RNA表达差异分析路径
#args2      编辑差异分析结果路径
#args3      输出路径
#args4      差异分析方法
#args5      分组方案
#args6      loci:y |gene:n

args <- commandArgs(trailingOnly=TRUE)


library(data.table)
library(stringr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggExtra)
library(ggrepel)


smart_ggsave <- function(plot, 
                         file, 
                         units = "cm", 
                         max_elements = 100, 
                         max_columns = 40) 
{
#######
      library(ggplot2)
      # 计算图形中的行数和列数
      row_count <- nrow(plot$data)
      col_count <- ncol(plot$data)
      
      # 计算元素总数
      element_count <- row_count * col_count
      
      # 计算基本宽度和高度
      base_width <- sqrt(element_count / max_elements) * 24
      base_height <- sqrt(element_count / max_elements) * 20
      
      # 动态调整宽度和高度
      width <- base_width
      height <- base_height
      
      # 调整字体大小
      font_size <- 12 / sqrt(element_count / max_elements)
      plot <- plot +
            theme(axis.title.x = element_text(size = font_size + 8),
                  axis.text.x = element_text(size = font_size + 8,hjust = 0),
                  axis.title.y = element_text(size = font_size + 8),
                  axis.text.y = element_text(size = font_size + 8,hjust = 0),
                  axis.text = element_text(size = font_size +4),
                  axis.title = element_text(size = font_size + 4),
                  plot.title = element_text(size = font_size + 4),
                  legend.text = element_text(size = font_size +8),
                  legend.title = element_text(size = font_size + 8)
            )
      
      # 保存图像
      ggsave(plot, file = file, units = units, width = width, height = height, dpi = 300)
}


#######初始统计样本的类型和频率信息
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


Edit_RNA_association_Quadrant <- function(RNA_DEG,
                                          Edit_DEG,
                                          filepath,
                                          i,
                                          method)
{
      ######
      print("Edit_RNA_association_Quadrant.start")
      RNA_DEG = RNA_DEG
      Edit_DEG = Edit_DEG
      i = i
      folder_path = filepath
      method = method

      RNA_DEG <- RNA_DEG[ !is.na(RNA_DEG$P.Value),]
      Edit_DEG <- Edit_DEG[ !is.na(Edit_DEG$P.Value),]
      print(head(RNA_DEG))
      print(head(Edit_DEG))
      combine= merge(RNA_DEG,Edit_DEG,
                     by.x="genes",
                     by.y="genes",
                     suffixes = c("_RNA","_Edit") ,
                     all.x=FALSE,
                     all.y=FALSE)
      write.csv(combine,file=paste0(folder_path,"/",i,"_",method,"_Nine_Quadrant.csv"),row.names=FALSE)
      write.table(combine,file=paste0(folder_path,"/",i,"_",method,"_Nine_Quadrant.xls"),row.names=FALSE,sep="\t")
      data <- data.frame(combine[ , c("genes","logFC_RNA","P.Value_RNA","logFC_Edit","P.Value_Edit")])
      data <- data %>% drop_na()
      #计算两个组学差异倍数的相关性，并取2位小数;
      cor = round(cor(data$logFC_RNA,data$logFC_Edit),2)
      #准备作为图形的标题;
      lab = paste("correlation=",cor,sep="")
      data$part <- case_when(abs(data$logFC_RNA) >= 1 & abs(data$logFC_Edit) >= 1 ~ "part1379",
                             abs(data$logFC_RNA) < 1 & abs(data$logFC_Edit) > 1 ~ "part28",
                             abs(data$logFC_RNA) > 1 & abs(data$logFC_Edit) < 1 ~ "part46",
                             abs(data$logFC_RNA) < 1 & abs(data$logFC_Edit) < 1 ~ "part5")
      
      mycolor <- c("#FF7F00","#FFC125","#E057E0","gray90","#AB82FF")
      x_breaks_cor <- c(seq(-100, 100, by = 10),1.3,-1.3)
      y_breaks_cor <- c(seq(-100, 100, by = 10), 1.3,-1.3)
      
      p0 <-ggplot(data,aes(logFC_RNA,logFC_Edit,color=part))+
            geom_point(size=3)+guides(color="none")+
            scale_colour_manual(name="",values=alpha(mycolor,0.7))+
            geom_hline(yintercept = c(-1,1),
                       size = 1,
                       color = "grey40",
                       lty = "dashed")+
            geom_vline(xintercept = c(-1,1),
                       size = 1,
                       color = "grey40",
                       lty = "dashed") + 
            theme_bw()+
            labs(title="",fill = "Sig",color="Sig",x = "Log2(Fold Change) of RNA_expression",y = "Log2(Fold Change) of Edit_expression")+
            theme_style()
      p0
      x_breaks <- sort(unique(c(ggplot_build(p0)$layout$panel_params[[1]]$x$breaks, -1, 1)))
      y_breaks <- sort(unique(c(ggplot_build(p0)$layout$panel_params[[1]]$y$breaks, -1, 1)))
      p0 <- p0 + scale_x_continuous(breaks = x_breaks)+
            scale_y_continuous(breaks = y_breaks)+
            theme_style()
      p1 <- ggMarginal(p0, type="densigram",fill="cyan3")
      p1
      ggsave(p1,file=paste0(folder_path,"/",i,"_",method,"_Quadrant1.pdf"),units = "cm",height = 16,width = 20,dpi = 300)
      ggsave(p1,file=paste0(folder_path,"/",i,"_",method,"_Quadrant1.png"),units = "cm",height = 16,width = 20,dpi = 300)
      
      #如果考虑差异的显著性，则需要进一步分组；
      #生成至少在一个组学显著上下调的数据标签；
      data$sig <- case_when(data$logFC_RNA < -1 & data$logFC_Edit > 1 & data$P.Value_RNA < 0.05 & data$P.Value_Edit < 0.05 ~ "sig1",
                            data$logFC_RNA > 1 & data$logFC_Edit > 1 & data$P.Value_RNA < 0.05 & data$P.Value_Edit < 0.05 ~ "sig3",
                            data$logFC_RNA < -1 & data$logFC_Edit < -1 & data$P.Value_RNA < 0.05 & data$P.Value_Edit < 0.05 ~ "sig7",
                            data$logFC_RNA > 1 & data$logFC_Edit < -1 & data$P.Value_RNA < 0.05 & data$P.Value_Edit < 0.05 ~ "sig9",
                            data$P.Value_RNA >= 0.05 | data$P.Value_Edit >=0.05 | data$logFC_RNA < 1 | data$logFC_Edit< 1 ~ "sig24568")
      sig <- filter(data, sig %in% c("sig1", "sig3","sig7","sig9"))
      non <- filter(data,sig == "sig24568")
      
      mycolor <- c("#FF7F00","gray90","#1E90FF","#E057E0","#56BD8F")
      p2 <- ggplot(data,aes(logFC_RNA,logFC_Edit,color=sig))+
            geom_point(size=3)+guides(color="none")+
            scale_colour_manual(name="",values=alpha(mycolor,0.7))+
            geom_hline(yintercept = c(-1,1),
                       size = 1,
                       color = "grey3",
                       lty = "dashed")+
            geom_vline(xintercept = c(-1,1),
                       size = 1,
                       color = "grey3",
                       lty = "dashed") + 
            theme_bw()+
            labs(title="",fill = "Sig",color="Sig",x = "Log2(Fold Change) of RNA_expression",y = "Log2(Fold Change) of Edit_expression")+
            geom_label_repel(data = sig, aes(label = genes), family = "serif", size = 6,nudge_y = 0.5)+
            geom_point(data=sig,size=6)+
            theme_style()
      p2
      x_breaks <- sort(unique(c(ggplot_build(p2)$layout$panel_params[[1]]$x$breaks, -1, 1)))
      y_breaks <- sort(unique(c(ggplot_build(p2)$layout$panel_params[[1]]$y$breaks, -1, 1)))
      p2 <- p2 + scale_x_continuous(breaks = x_breaks)+
            scale_y_continuous(breaks = y_breaks)+
            theme_style()
      p3 <- ggMarginal(p2, type="densigram",fill="cyan3")
      p3
      ggsave(p3,file=paste0(folder_path,"/",i,"_",method,"_Quadrant2.pdf"),units = "cm",height = 16,width = 20,dpi = 300)
      ggsave(p3,file=paste0(folder_path,"/",i,"_",method,"_Quadrant2.png"),units = "cm",height = 16,width = 20,dpi = 300)
      
      data$RNA_Rank <- data$logFC_RNA * -log10(data$P.Value_RNA)
      data$Edit_Rank <- data$logFC_Edit * -log10(data$P.Value_Edit)
      data$PP_Rank <- data$P.Value_RNA * data$P.Value_Edit
      sig$RNA_Rank <- sig$logFC_RNA * -log10(sig$P.Value_RNA)
      sig$Edit_Rank <- sig$logFC_Edit * -log10(sig$P.Value_Edit)
      sig$PP_Rank <- sig$P.Value_RNA * sig$P.Value_Edit
      nosig <- data[!data$genes %in% sig$genes, ]
      
      
      p4 <- ggplot(data,aes(RNA_Rank,Edit_Rank,color=sig))+
            geom_point(size=3)+guides(color="none")+
            scale_colour_manual(name="",values=alpha(mycolor,0.7))+
            geom_hline(yintercept = c(-1.3,1.3),
                       size = 1,
                       color = "grey3",
                       lty = "dashed")+
            geom_vline(xintercept = c(-1.3,1.3),
                       size = 1,
                       color = "grey3",
                       lty = "dashed") + 
            theme_bw()+
            labs(title="",fill = "Sig",color="Sig",x = "RNA_expression Rank",y = "Edit_expression Rank")+
            geom_label_repel(data = sig, aes(label = genes), family = "serif", size = 6,nudge_y = 0.5)+
            geom_point(data=sig,aes(size=-log10(PP_Rank)))+
            theme_style()+
            scale_x_continuous(breaks = x_breaks_cor) +
            scale_y_continuous(breaks = y_breaks_cor)
      p4
      x_breaks <- sort(unique(c(ggplot_build(p4)$layout$panel_params[[1]]$x$breaks, -1.3, 1.3)))
      y_breaks <- sort(unique(c(ggplot_build(p4)$layout$panel_params[[1]]$y$breaks, -1.3, 1.3)))
      p4 <- p4 + scale_x_continuous(breaks = x_breaks)+
            scale_y_continuous(breaks = y_breaks)+
            theme_style()
      p4
      ggsave(p4,file=paste0(folder_path,"/",i,"_",method,"_Quadrant3.pdf"),units = "cm",height = 16,width = 20,dpi = 300)
      ggsave(p4,file=paste0(folder_path,"/",i,"_",method,"_Quadrant3.png"),units = "cm",height = 16,width = 20,dpi = 300)
      #ggMarginal(p, type = "densigram", groupColour = TRUE, groupFill = TRUE, alpha = 0.4)
      
      p6 <- ggplot(data,aes(RNA_Rank,Edit_Rank))+
            geom_point(aes(size=3,color=-log10(PP_Rank)))+
            scale_color_gradient(low="dodgerblue4",high = "brown3")+
            geom_hline(yintercept = c(-1.3,1.3),
                       size = 1,
                       color = "grey3",
                       lty = "dashed")+
            geom_vline(xintercept = c(-1.3,1.3),
                       size = 1,
                       color = "grey3",
                       lty = "dashed") + 
            theme_bw()+
            labs(title="",fill = "Sig",color="-Lg(FDR)",x = "RNA_expression Rank",y = "Edit_expression Rank")+
            geom_label_repel(data = sig, aes(label = genes), family = "serif", size = 6,nudge_y = 0.5)+
            geom_point(data=sig,aes(size=4,color=-log10(PP_Rank)))+
            theme_style()+guides(size = "none")+
            scale_x_continuous(breaks = x_breaks_cor) +
            scale_y_continuous(breaks = y_breaks_cor)
      p6
      x_breaks <- sort(unique(c(ggplot_build(p6)$layout$panel_params[[1]]$x$breaks, -1.3, 1.3)))
      y_breaks <- sort(unique(c(ggplot_build(p6)$layout$panel_params[[1]]$y$breaks, -1.3, 1.3)))
      p7 <- p6 + scale_x_continuous(breaks = x_breaks)+
            scale_y_continuous(breaks = y_breaks)+
            theme_style()
      p7
      ggsave(p7,file=paste0(folder_path,"/",i,"_",method,"_Quadrant4.pdf"),units = "cm",height = 16,width = 20,dpi = 300)
      ggsave(p7,file=paste0(folder_path,"/",i,"_",method,"_Quadrant4.png"),units = "cm",height = 16,width = 20,dpi = 300)
      
      p8 <- ggplot(data, aes(RNA_Rank, Edit_Rank, color = sig)) +
            guides(color = "none") +
            geom_point(data = nosig, aes(size = 3), shape = 19,stroke = 0) +
            scale_colour_manual(name = "", values = alpha(mycolor, 0.7)) +
            geom_hline(yintercept = c(-1.3, 1.3), size = 1, color = "grey3", linetype = "dashed") +
            geom_vline(xintercept = c(-1.3, 1.3), size = 1, color = "grey3", linetype = "dashed") +
            theme_bw() +
            labs(title = "", fill = "Sig" , x = "RNA_expression Rank", y = "Edit_expression Rank",size="-Log10(PP_Rank)") +
            geom_label_repel(data = sig, aes(label = genes), family = "serif", nudge_y = 0.8) +
            geom_point(data = sig, aes(size = -log10(PP_Rank)), shape = 19,stroke = 0) +
            theme_style() +
            theme(
                  panel.grid.major = element_blank(),  # 移除大网格线
                  panel.grid.minor = element_blank()   # 移除小网格线
            )
      
      p8
      x_breaks <- sort(unique(c(ggplot_build(p8)$layout$panel_params[[1]]$x$breaks, -1.3, 1.3)))
      y_breaks <- sort(unique(c(ggplot_build(p8)$layout$panel_params[[1]]$y$breaks, -1.3, 1.3)))
      x_breaks <- x_breaks[x_breaks != 0]
      y_breaks <- y_breaks[y_breaks != 0]
      p9 <- p8 + scale_x_continuous(breaks = x_breaks)+
            scale_y_continuous(breaks = y_breaks)+
            theme_style()
      p9
      ggsave(p9,file=paste0(folder_path,"/",i,"_",method,"_Quadrant5.pdf"),units = "cm",height = 16,width = 20,dpi = 300)
      ggsave(p9,file=paste0(folder_path,"/",i,"_",method,"_Quadrant5.png"),units = "cm",height = 16,width = 20,dpi = 300)
      print("Edit_RNA_association_Quadrant.done")
}

Edit_RNA_association_Statistics <- function(
                                            out_path,
                                            RNA_DEG_list,
                                            Edit_DEG_list,
					    chr_work,
                                            method)
{
      ######
      print("Edit_RNA_association_Statistics.start")
      #Edit_RNA_association_DEG <- file.path(out_path, "09_association", str_glue("01_Gene_Edit_RNA"))
      Edit_RNA_association_DEG <- out_path

      if (!dir.exists(Edit_RNA_association_DEG)) {
            dir.create(Edit_RNA_association_DEG, recursive = TRUE)
      }

      
      method <- method
      
      RNA_Volcano_Limma <- RNA_DEG_list
      Edit_Volcano_Limma <- Edit_DEG_list
      comparisons <- comparisons
      for ( i in comparisons) {
            # 使用 get() 函数获取名为 i 的数据框
	    i <- gsub("\\+", "_", i)
            Limma_frame <- get("i")
            RNA_DEG <- RNA_Volcano_Limma[[i]]
            DEG <- Edit_Volcano_Limma[[i]]
	    if(chr_work == "y"){
	        DEG$Gene <- ifelse(grepl("chr", DEG$genes), sub(".*_(.*)", "\\1", DEG$genes), NA)
	        DEG <- select(DEG,-genes)
		DEG <- dplyr::rename(DEG, genes = Gene )
	    }
            Edit_RNA_association_Quadrant(RNA_DEG = RNA_DEG,
                                          Edit_DEG = DEG,
                                          filepath = Edit_RNA_association_DEG,
                                          i = i,
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
      for (cmp in comparisons){
	    cmp <- gsub("\\+", "_", cmp)
            DEG <- read.csv(paste0(work_path,"/",cmp,"_",method,".csv"), header = TRUE)
            DEG_list[[cmp]] <- DEG
      }
      return(DEG_list)
}

comparisons <- read_group(group_file = args[5])

RNA_DEG_list <- read_list(
      work_path = args[1],
      method = args[4],
      comparisons = comparisons)

Edit_DEG_list <- read_list(
      work_path = args[2],
      method = args[4],
      comparisons = comparisons)

Edit_RNA_association_Statistics(
      out_path = args[3],
      RNA_DEG_list = RNA_DEG_list,
      Edit_DEG_list = Edit_DEG_list,
      chr_work=args[6],
      method = args[4])
