#Rscript - diff_work_path output_path method group_file
#args1      富集结果路径
#args2      差异分析方法
#args3      分组文件

library(stringr)
library(ggplot2)
library(glue)
library(dplyr)
library(pathview)


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

smart_ggsave <- function(plot,
                         file,
                         units = "cm",
                         max_elements = 100,
                         size_num = 12,
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
      font_size <- size_num / sqrt(element_count / max_elements)
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

plot_GO <- function(
    method,
    work_path,
    go_enrich_df,
    up_down,
    cmp)
{
    go_enrich_df = go_enrich_df
    for(ii in 1:nrow(go_enrich_df)){
        description_splite=strsplit(go_enrich_df$Description[ii],split = " ")
        description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
        go_enrich_df$Description[ii]=description_collapse
        go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
    }
    go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(str_to_title(go_enrich_df$Description))) #这一步是必须的，为了让柱子按顺序显示
    COLS <- c("burlywood2", "chartreuse4", "#984EA3")#设定颜色
    P1 = ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
        geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
        scale_fill_manual(values = COLS) + ###颜色
        coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
        xlab("GO term") +
        ylab("_Number") +
        labs(title = "",fill = "Type",color="Type")+
        theme_bw() + theme_style()
    max1 <- max(go_enrich_df$GeneNumber)
    if (max1 <= 3 & nrow(go_enrich_df) <= 10){
        smart_ggsave(P1,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_GO1_padj.pdf"),size_num=6)
        smart_ggsave(P1,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_GO1_padj.png"),size_num=6)
    }else{
        smart_ggsave(P1,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_GO1_padj.pdf"))
        smart_ggsave(P1,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_GO1_padj.png"))
    }

    P2 = ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) +
        geom_bar(stat="identity", width=0.8) +
        scale_fill_manual(values = COLS) +
        theme_bw() +
        xlab("GO term") +
        ylab("Num of Genes") +
        labs(title = "",fill = "Type",color="Type")+
        theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 90,vjust = 0, hjust = 0))  + #angle是坐标轴字体倾斜的角度，可以自己设置
        theme_style()
    if (max1 <= 3 & nrow(go_enrich_df) <= 10){
        smart_ggsave(P2,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_GO2_padj.pdf"),size_num=6)
        smart_ggsave(P2,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_GO2_padj.png"),size_num=6)
    }else{
        smart_ggsave(P2,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_GO2_padj.pdf"))
        smart_ggsave(P2,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_GO2_padj.png"))
    }
}

plot_KEGG <- function(
    method,
    work_path,
    hh,
    up_down,
    cmp)
{
    hh = hh
    P3 = ggplot(hh,aes(y=Description,x=Count,fill=p.adjust))+
        geom_bar(stat = "identity",width=0.7)+####柱子宽度
        #coord_flip()+##颠倒横纵轴
        scale_fill_gradient(low = "brown3",high ="dodgerblue4" )+#颜色自己可以换
        labs(title = "",
            x = "Gene numbers",
            y = "Pathways",fill = "p.adjust",color="p.adjust") +
        theme_bw() +theme_style()
    max1 <- max(hh$Count)
    if (max1 <= 3 & nrow(hh) <= 10){
        smart_ggsave(P3,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_KK1_padj.pdf"),size_num=6)
        smart_ggsave(P3,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_KK1_padj.png"),size_num=6)
    }else{
        smart_ggsave(P3,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_KK1_padj.pdf"))
        smart_ggsave(P3,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_KK1_padj.png"))
    }
    P4 = ggplot(hh,aes(y=Description,x=Count))+
        geom_point(aes(size=Count,color=p.adjust))+# 修改点的大小
        scale_color_gradient(low="dodgerblue4",high = "brown3")+
        labs(color=expression(p.adjust,size="Count"),
            x="Gene Count",y="Pathways",title="",fill = "p.adjust",color="p.adjust")+
        theme_bw()+ theme_style()
    if (max1 <= 3 & nrow(hh) <= 10){
        smart_ggsave(P4,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_KK2_padj.pdf"),size_num=6)
        smart_ggsave(P4,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_KK2_padj.png"),size_num=6)
    }else{
        smart_ggsave(P4,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_KK2_padj.pdf"))
        smart_ggsave(P4,file=paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_",up_down,"_KK2_padj.png"))
    }
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
        print(cmp)
        ego_result_BP <- read.table(paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_Up_BP.xls"), header = TRUE, sep = "\t")
        ego_result_CC <- read.table(paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_Up_CC.xls"), header = TRUE, sep = "\t")
        ego_result_MF <- read.table(paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_Up_MF.xls"), header = TRUE, sep = "\t")
        ego_result_BP <- ego_result_BP[ego_result_BP$p.adjust<0.05,]
        ego_result_CC <- ego_result_CC[ego_result_CC$p.adjust<0.05,]
        ego_result_MF <- ego_result_MF[ego_result_MF$p.adjust<0.05,]
	if (nrow(ego_result_BP) > 100){
	    ego_result_BP <- ego_result_BP[1:100,]
	}
	if (nrow(ego_result_CC) > 100){
            ego_result_CC <- ego_result_CC[1:100,]
        }
	if (nrow(ego_result_MF) > 100){
            ego_result_MF <- ego_result_MF[1:100,]
        }
        go_enrich_df <- data.frame(
            ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
            Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
            GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
            type=factor(c(rep("biological process", nrow(ego_result_BP)),
                rep("cellular component", nrow(ego_result_CC)),
                rep("molecular function", nrow(ego_result_MF))),
                levels=c("biological process", "cellular component","molecular function" )))
        go_enrich_df <- go_enrich_df[!is.na(go_enrich_df$Description), ]
	print(nrow(go_enrich_df))
        if (nrow(go_enrich_df) > 0){
            plot_GO(
                method = method,
                work_path = work_path,
                go_enrich_df = go_enrich_df,
                up_down = "Up",
                cmp = cmp
            )
        }

        ego_result_BP2 <- read.table(paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_Down_BP.xls"), header = TRUE, sep = "\t")
        ego_result_CC2 <- read.table(paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_Down_CC.xls"), header = TRUE, sep = "\t")
        ego_result_MF2 <- read.table(paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_Down_MF.xls"), header = TRUE, sep = "\t")
        ego_result_BP2 <- ego_result_BP2[ego_result_BP2$p.adjust<0.05,]
        ego_result_CC2 <- ego_result_CC2[ego_result_CC2$p.adjust<0.05,]
        ego_result_MF2 <- ego_result_MF2[ego_result_MF2$p.adjust<0.05,]
	if (nrow(ego_result_BP2) > 100){
            ego_result_BP2 <- ego_result_BP2[1:100,]
        }
        if (nrow(ego_result_CC2) > 100){
            ego_result_CC2 <- ego_result_CC2[1:100,]
        }
        if (nrow(ego_result_MF2) > 100){
            ego_result_MF2 <- ego_result_MF2[1:100,]
        }
        go_enrich_df2 <- data.frame(
            ID=c(ego_result_BP2$ID, ego_result_CC2$ID, ego_result_MF2$ID),
            Description=c(ego_result_BP2$Description, ego_result_CC2$Description, ego_result_MF2$Description),
            GeneNumber=c(ego_result_BP2$Count, ego_result_CC2$Count, ego_result_MF2$Count),
            type=factor(c(rep("biological process", nrow(ego_result_BP2)),
                rep("cellular component", nrow(ego_result_CC2)),
                rep("molecular function", nrow(ego_result_MF2))),
                levels=c("biological process", "cellular component","molecular function" )))
        go_enrich_df2 <- go_enrich_df2[!is.na(go_enrich_df2$Description), ]
	print(nrow(go_enrich_df2))
        if(nrow(go_enrich_df2) > 0){
            plot_GO(
                method = method,
                work_path = work_path,
                go_enrich_df = go_enrich_df2,
                up_down = "Down",
                cmp = cmp
            )
        }
        #print(go_enrich_df) 
        #print(ego_result_BP)
        #print(ego_result_CC)
        #print(ego_result_MF)
        hh <- read.table(paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_Up_KK.xls"), header = TRUE, sep = "\t")
        hh <- hh[hh$p.adjust<0.05,]
        hh <- hh[!is.na(hh$Description), ]
	if (nrow(hh) > 300){
	    hh <- hh[1:300,]
	}
	print(nrow(hh))
        if(nrow(hh) > 0){
            plot_KEGG(
                method = method,
                work_path = work_path,
                hh = hh,
                up_down = "Up",
                cmp = cmp
            )
        }

        hh2 <- read.table(paste0(work_path,"/",cmp,"/R_Enrichment/",cmp,"_",method,"_Down_KK.xls"), header = TRUE, sep = "\t")
        hh2 <- hh2[hh2$p.adjust<0.05,]
        hh2 <- hh2[!is.na(hh2$Description), ]
	if (nrow(hh2) > 300){
            hh2 <- hh2[1:300,]
        }
	print(nrow(hh2))
        if(nrow(hh2) > 0){
            plot_KEGG(
                method = method,
                work_path = work_path,
                hh = hh2,
                up_down = "Down",
                cmp = cmp
            )
        }
    }
    return()
}

comparisons <- read_group(group_file = args[3])

read_list(
    work_path = args[1],
    method = args[2],
    comparisons = comparisons)
