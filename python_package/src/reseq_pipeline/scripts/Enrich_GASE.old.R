#Rscript - diff_work_path output_path method group_file
#args1      差异分析路径
#args2      输出路径
#args3      差异分析方法
#args4      分组文件
#args5      ref: hsa | mmu | rno
#args6      loci:y |gene:n

library(stringr)
library(ggplot2)
library(glue)
library(dplyr)
library(org.Hs.eg.db)#基因注释包
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(clusterProfiler)
library(pathview)
library(data.table)
library(enrichplot)
library(patchwork)

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
      base_width <- sqrt(element_count / max_elements) * 20
      base_height <- sqrt(element_count / max_elements) * 16

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

Enrichment_GASE <- function(DEG,
                            filepath,
                            i,
			    database,
			    database2,
                            method)
{
      ######
      #loadfonts(device = "win") 
      filepath = filepath
      DEG = DEG
      i = i
      method = method
      folder_path = paste0(filepath,'/',i)
      GSEA_Enrichment_path = file.path(folder_path,"/GSEA_Enrichment")
      if (!dir.exists(folder_path)) {
            dir.create(folder_path, recursive = TRUE)
      }
      if (!dir.exists(GSEA_Enrichment_path)) {
            dir.create(GSEA_Enrichment_path, recursive = TRUE)
      }
      DEG <- DEG
      DEG <- DEG[is.finite(DEG$logFC), ]
      gene_aninly <- DEG[, c("genes", "logFC")]
      gene.df <- bitr(gene_aninly$genes,fromType="SYMBOL",toType="ENTREZID", OrgDb = database,drop=TRUE)
      final <- merge(gene_aninly,gene.df,gene_aninly,by.x="genes",by.y="SYMBOL")
      entrezid_logFC <- dplyr::select(final,ENTREZID,logFC)
      final.sort <- arrange(entrezid_logFC,desc(logFC))
      genelist <- final.sort$logFC
      names(genelist)=final.sort$ENTREZID
      genelist=sort(genelist,decreasing = T)
      gsemf <- gseGO(genelist,
                     OrgDb = database,
                     keyType = "ENTREZID",
                     ont="BP",
                     pvalueCutoff = 0.8)
      convert_geneID_to_symbol <- function(geneID_string) {
            # 分割 geneID 字符串
            entrez_ids <- unlist(strsplit(geneID_string, "/"))
            
            # 将 ENTREZID 转换为 Symbol 号
            gene_symbols <- bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = database)
            
            # 检查是否所有的 ENTREZID 都成功转换
            if (nrow(gene_symbols) > 0) {
                  # 返回转换后的 Symbol 号字符串
                  return(paste(gene_symbols$SYMBOL, collapse = "/"))
            } else {
                  # 如果没有成功转换，返回原始 ENTREZID 字符串
                  return(geneID_string)
            }
      }
      # 将函数应用到富集结果的 geneID 列
      gsemf@result$SYMBOL <- sapply(gsemf@result$core_enrichment, convert_geneID_to_symbol)
      gsemf_1 = data.frame(gsemf@result)
      gsemf_1 = gsemf_1 %>% select(ID,Description,setSize,enrichmentScore,NES,pvalue,p.adjust,SYMBOL)
      write.csv(gsemf@result,file=paste0(GSEA_Enrichment_path,'/',i,'_',method,'_GASE_GO.csv'),quote=FALSE,row.names=F)
      write.table(gsemf_1,file=paste0(GSEA_Enrichment_path,'/',i,'_',method,'_GASE_GO.xls'),quote=FALSE,row.names=F,sep="\t")
      if (nrow(gsemf@result) > 0) {  
      gsemf <- gsemf %>%
            dplyr::filter(abs(NES) > 1.5)
      } else {
            message("No enriched GO terms found.")  
      } 
      GSEA_conbin <- function(Plot){
            gseaPlotList <- lapply(Plot, function(p) {
                  p + theme(axis.title.x = element_text(family = "serif", size = 20,colour = "black"),
                            axis.text.x = element_text(family = "serif", size = 20,colour = "black"),
                            axis.title.y = element_text(family = "serif", size = 20,colour = "black"),
                            axis.text.y = element_text(family = "serif", size = 20,colour = "black"),
                            plot.title = element_text(family = "serif", size = 20, face = "bold",hjust = 0.5,colour = "black"),
                            legend.title = element_text(family = "serif", size = 20,colour = "black"),
                            legend.text = element_text(family = "serif", size = 20,colour = "black"),# 设置p-value表格字体大小
                            plot.subtitle = element_text(family = "serif", size = 20),
                            plot.caption = element_text(family = "serif", size = 20) #,  aspect.ratio = 9/16 # + theme(text = element_text(family = "serif", size = 20))  # geom_text(aes(x = 4, y = 10, label = "This is a point"), color = "red") 
                  )
            })
            # 使用gridExtra包的arrangeGrob函数合并多个图
            plot1 <- gseaPlotList[[1]]
            plot2 <- gseaPlotList[[2]] + annotate("text", x = -Inf, y = Inf, label = "Case", hjust = 0, vjust = 1.5, size = 6, family = "serif",color = "brown3",fontface = "bold") +annotate("text", x = Inf, y = Inf, label = "Control", hjust = 1, vjust = 1.5, size = 6, family = "serif",color = "dodgerblue4",fontface = "bold")
            plot3 <- gseaPlotList[[3]]
            
            # 使用patchwork将三个图形按照X轴对齐
            GSEA_GO_Plot <- plot1 / plot2 / plot3 + plot_layout(heights = c(2, 1, 1))
            return(GSEA_GO_Plot)
      }
      
      if (nrow(gsemf@result) > 0) {
      for(j in 1:nrow(gsemf@result)){
            gseaPlot = gseaplot2(gsemf, 
                                 title = str_to_title(paste(gsemf$Description[j], "\nP-value:", formatC(gsemf$pvalue[j], format = "e", digits = 4),"\nNES:",round(gsemf$NES[j], 4))), 
                                 geneSetID = j,pvalue_table=F,
                                 base_size = 20)
            GSEA_GO_Plot <- GSEA_conbin(gseaPlot)
            ggsave(plot = GSEA_GO_Plot,file=paste0(GSEA_Enrichment_path,"/",i,"_",method,"_GSEA_GO",j,".pdf"),  units = "cm", height = 24, width = 28,dpi = 300)
            ggsave(plot = GSEA_GO_Plot,file=paste0(GSEA_Enrichment_path,"/",i,"_",method,"_GSEA_GO",j,".png"),  units = "cm", height = 24, width = 28,dpi = 300)}
      } else { 
      message("No enriched GO terms found")  }  
      
      ##### gsea富集 ####
      tryCatch({
            gsekk <- gseKEGG(geneList     = genelist,
                       organism = database2,  ##大鼠rno 小鼠mmu
                       keyType = "kegg",
                       pvalueCutoff = 0.5,
                       verbose      = FALSE)
      gsekk@result$SYMBOL <- sapply(gsekk@result$core_enrichment, convert_geneID_to_symbol)
      gsekk_1 = data.frame(gsekk@result)
      gsekk_1 = gsekk_1 %>% select(ID,Description,setSize,enrichmentScore,NES,pvalue,p.adjust,SYMBOL)
      write.csv(gsekk,file=paste0(GSEA_Enrichment_path,"/",i,"_",method,'_GASE_KK.csv'),quote=FALSE,row.names=F)
      write.table(gsekk,file=paste0(GSEA_Enrichment_path,"/",i,"_",method,'_GASE_KK.xls'),quote=FALSE,row.names=F,sep="\t")
      # 检查是否有富集结果并且结果不为空
      if (is.null(gsekk) || length(gsekk@result) == 0) {  # 确保gsekk对象的正确检查
            cat("在", i, "中没有富集到任何通路。\n")
            return(FALSE)  # 返回 FALSE 表示未富集到通路
      }
      
      # 如果数据框没有行，返回 FALSE
      if (nrow(gsekk@result) == 0) {
            cat("在", i, "中没有富集到任何通路。\n")
            return(FALSE)
      }
      # 检查并处理缺失值，然后筛选
      gsekk <- gsekk %>%
            filter(!is.na(NES) & abs(NES) > 1.5)  # 过滤掉 NES 中的缺失值并筛选绝对值大于 1.5 的行
      
      # 如果没有符合条件的通路，返回 FALSE
      if (nrow(gsekk) == 0) {
            cat("在", i, "中没有符合条件的通路。\n")
            return(FALSE)
      }
      
      cat("在", i, "中富集到符合条件的通路，继续进行后续分析。\n")
      gsekk <- gsekk %>%
            dplyr::mutate(Description = str_to_title(gsub(" - Mus musculus \\(house mouse\\)", "", Description)))
      for(k in 1:nrow(gsekk)){
            gseaPlot = gseaplot2(gsekk, 
                                 title = str_to_title(paste(gsekk$Description[k], "\nP-value:", formatC(gsekk$pvalue[k], format = "e", digits = 4),"\nNES:",round(gsekk$NES[k], 4))), 
                                 geneSetID = k,pvalue_table=F,
                                 base_size = 20)
            GSEA_KK_Plot <- GSEA_conbin(gseaPlot)
            ggsave(plot = GSEA_KK_Plot,file=paste0(GSEA_Enrichment_path,"/",i,"_",method,"_GSEA_KK",k,".pdf"),  units = "cm", height = 24, width = 28,dpi = 300)# 返回 TRUE 表示富集到通路
            ggsave(plot = GSEA_KK_Plot,file=paste0(GSEA_Enrichment_path,"/",i,"_",method,"_GSEA_KK",k,".png"),  units = "cm", height = 24, width = 28,dpi = 300)}# 返回 TRUE 表示富集到通路
      return(TRUE)
}, error = function(e) {
      message("Error in Enrichment_GASE: ", e)
      return(FALSE)
})
}


Edit_Functional_Statistics <- function(
                                    out_path,
                                    Edit_DEG_list,
                                    method,
				    database,
				    database2,
				    chr_work,
                                    comparisons)
{
      #######初始统计样本的类型和频率信息
      #Edit_Functional_information <- file.path(out_path, "03_Edit", str_glue("07_Edit_Func_enrich"))
      #Edit_Functional_path <- file.path(out_path, "03_Edit", str_glue("07_Edit_Func_enrich"),method)
      Edit_Functional_information <- out_path
      #Edit_Functional_path <- file.path(out_path, method)
      Edit_Functional_path <- out_path
      outpath = getwd()
      
      if (!dir.exists(Edit_Functional_information)) {
            dir.create(Edit_Functional_information, recursive = TRUE)
      }
      if (!dir.exists(Edit_Functional_path)) {
            dir.create(Edit_Functional_path, recursive = TRUE)
      }

      Volcano <- Edit_DEG_list
      comparisons <- comparisons
      for ( i2 in comparisons) {
	    i2 <- gsub("\\+", "_", i2)
            # 使用 get() 函数获取名为 i 的数据框
            Limma_frame <- get("i2")
            DEG <- Volcano[[i2]]
	    if (chr_work == "y"){
	    	DEG$Gene <- ifelse(grepl("chr", DEG$genes), sub(".*_(.*)", "\\1", DEG$genes), NA)
            	DEG <- select(DEG,-genes)
            	DEG <- dplyr::rename(DEG, genes = Gene )
	    }
            setwd(outpath)
            if (! Enrichment_GASE(DEG = DEG,
                            filepath = Edit_Functional_path,
                            i = i2,
			    database = database,
			    database2 = database2,
                            method = method)){next}
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

comparisons <- read_group(group_file = args[4])

DEG_list <- read_list(
      work_path = args[1],
      method = args[3],
      comparisons = comparisons)

if(args[5] == "hg38"){
        database="org.Hs.eg.db"
        database2="hsa"
}else if(args[5] == "mm10"){
        database="org.Mm.eg.db"
        database2="mmu"
}else if(args[5] == "rno"){
        database="org.Rn.eg.db"
        database2="rno"
}

Edit_Functional_Statistics(
      out_path = args[2],
      Edit_DEG_list = DEG_list,
      method = args[3],
      database=database,
      database2=database2,
      chr_work=args[6],
      comparisons = comparisons)
