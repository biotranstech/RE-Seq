#Rscript - diff_work_path output_path method group_file
#args1      差异分析路径
#args2      输出路径
#args3      差异分析方法
#args4      分组文件
#args5      ref: hsa | mmu | rno
#args6      loci:y |gene:n
#args7      up_down; Up | Down

library(stringr)
library(ggplot2)
library(glue)
library(dplyr)
library(org.Hs.eg.db)#基因注释包
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(clusterProfiler)
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

Enrichment_R <- function(DEG,
                         filepath,
                         i,
			 database,
			 database2,
			 up_down,
                         method)
{
      ######
      filepath = filepath
      DEG = DEG
      i = i
      print(i)
      method = method
      folder_path = paste0(filepath,'/',i)
      if (!dir.exists(folder_path)) {
            dir.create(folder_path, recursive = TRUE)
      }
      gene_aninly <- DEG[, c("genes", "regulate","logFC")]
      #diff <- gene_aninly$genes[DEG$regulate == "Up" | DEG$regulate == "Down"]
      diff <- gene_aninly$genes[DEG$regulate == up_down]
      diffed <- unique(as.matrix(diff))
      #gene.df <- bitr(diffed,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Mm.eg.db) #TCGA数据框如果没有进行基因注释
      #gene.df <- bitr(diffed,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)#TCGA数据框如果没有进行基因注释
      gene.df <- bitr(diffed,fromType="SYMBOL",toType="ENTREZID", OrgDb = database)
      gene <- gene.df$ENTREZID
      Symol <- gene.df$SYMBOL
      
      logFC_information <- merge(gene_aninly,gene.df,by.x="genes",by.y="SYMBOL")
      logFC_numbers <- logFC_information[, c("ENTREZID","logFC")]
      ego_ALL <- enrichGO(gene = gene,#我们上面定义了
                          #OrgDb=org.Hs.eg.db,
			  OrgDb=database,
                          keyType = "ENTREZID",
                          ont = "ALL",#富集的GO类型
                          pAdjustMethod = "BH",#这个不用管，一般都用的BH
                          minGSSize = 1,
                          pvalueCutoff = 0.5,#P值可以取0.05
                          #qvalueCutoff = 0.05,
                          readable = TRUE)
      
      ego_CC <- enrichGO(gene = gene,
                         #OrgDb=org.Hs.eg.db,
			 OrgDb=database,
                         keyType = "ENTREZID",
                         ont = "CC",
                         pAdjustMethod = "BH",
                         minGSSize = 1,
                         pvalueCutoff = 0.1,
                         #qvalueCutoff = 0.05,
                         readable = TRUE)
      
      ego_BP <- enrichGO(gene = gene,
                         #OrgDb=org.Hs.eg.db,
			 OrgDb=database,
                         keyType = "ENTREZID",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         minGSSize = 1,
                         pvalueCutoff = 0.1,
                         #qvalueCutoff = 0.05,
                         readable = TRUE)
      
      ego_MF <- enrichGO(gene = gene,
                         #OrgDb=org.Hs.eg.db,
			 OrgDb=database,
                         keyType = "ENTREZID",
                         ont = "MF",
                         pAdjustMethod = "BH",
                         minGSSize = 1,
                         pvalueCutoff = 0.1,
                         #qvalueCutoff = 0.05,
                         readable = TRUE)
      
      #4、将结果保存到当前路径
      print(head(ego_ALL))
      ego_ALL <- as.data.frame(ego_ALL)
      ego_result_BP <- as.data.frame(ego_BP)
      ego_result_CC <- as.data.frame(ego_CC)
      ego_result_MF <- as.data.frame(ego_MF)
      ego <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)#或者这样也能得到ego_ALL一样的结果
      R_Enrichment_path = file.path(folder_path,"/R_Enrichment")
      if (!dir.exists(R_Enrichment_path)) {
            dir.create(R_Enrichment_path, recursive = TRUE)
      }
      print(R_Enrichment_path)
      #change yjc
      ego_ALL_1 <- select(ego_ALL,ID,Description,ONTOLOGY,GeneRatio,BgRatio,pvalue,p.adjust,Count,geneID)
      ego_result_BP_1 <- select(ego_result_BP,ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,Count,geneID)
      ego_result_CC_1 <- select(ego_result_CC,ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,Count,geneID)
      ego_result_MF_1 <- select(ego_result_MF,ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,Count,geneID)

      write.csv(ego_ALL,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_ALL.csv"),quote=FALSE,row.names=F)
      write.csv(ego_result_BP,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_BP.csv"),quote=FALSE,row.names=F)
      write.csv(ego_result_CC,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_CC.csv"),quote=FALSE,row.names=F)
      write.csv(ego_result_MF,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_MF.csv"),quote=FALSE,row.names=F)
      
      write.table(ego_ALL_1,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_ALL.xls"),quote=FALSE,row.names=F,sep="\t")
      write.table(ego_result_BP_1,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_BP.xls"),quote=FALSE,row.names=F,sep="\t")
      write.table(ego_result_CC_1,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_CC.xls"),quote=FALSE,row.names=F,sep="\t")
      write.table(ego_result_MF_1,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_MF.xls"),quote=FALSE,row.names=F,sep="\t")
      
      #5、但是有时候我们富集到的通路太多了，不方便绘制图片，可以选择部分绘制，这时候跳过第4步，来到第5步
      display_number = c(20, 20, 20)  #这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
      ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
      ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
      ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]
      
      ##将以上我们摘取的部分通路重新组合成数据框
      go_enrich_df <- data.frame(
            ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
            GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
            type=factor(c(rep("biological process", display_number[1]), 
                          rep("cellular component", display_number[2]),
                          rep("molecular function", display_number[3])), 
                        levels=c("biological process", "cellular component","molecular function" )))
      
      ##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
      for(ii in 1:nrow(go_enrich_df)){
            description_splite=strsplit(go_enrich_df$Description[ii],split = " ")
            description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
            go_enrich_df$Description[ii]=description_collapse
            go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
      }
      
      ##开始绘制GO柱状图
      ###横着的柱状图
      go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(str_to_title(go_enrich_df$Description))) #这一步是必须的，为了让柱子按顺序显示
      COLS <- c("burlywood2", "chartreuse4", "#984EA3")#设定颜色
      P1 = ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
            geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
            scale_fill_manual(values = COLS) + ###颜色
            coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
            xlab("GO term") + 
            ylab("Gene_Number") + 
            labs(title = "",fill = "Type",color="Type")+
            theme_bw() + theme_style()
      smart_ggsave(P1,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_GO1.pdf"))
      smart_ggsave(P1,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_GO1.png"))
      
      ###竖着的柱状图 
      P2 = ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
            geom_bar(stat="identity", width=0.8) + 
            scale_fill_manual(values = COLS) + 
            theme_bw() + 
            xlab("GO term") + 
            ylab("Num of Genes") + 
            labs(title = "",fill = "Type",color="Type")+ 
            theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 90,vjust = 0, hjust = 0))  + #angle是坐标轴字体倾斜的角度，可以自己设置
            theme_style()
      smart_ggsave(P2,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_GO2.pdf"))
      smart_ggsave(P2,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_GO2.png"))
      
      #1、KEGG富集
      kk <- enrichKEGG(gene = gene,keyType = "kegg",organism= database2, qvalueCutoff = 0.05, pvalueCutoff=0.05)
      print(head(gene))
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
      kk@result$SYMBOL <- sapply(kk@result$geneID, convert_geneID_to_symbol)
      kk_1 = data.frame(kk@result)
      kk_1 = kk_1 %>% select(ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,Count,SYMBOL)
      write.csv(kk@result,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_KK.csv"),quote=FALSE,row.names=F)
      write.table(kk_1,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_KK.xls"),quote=FALSE,row.names=F,sep="\t")
      #2、可视化
      ###柱状图
      hh <- as.data.frame(kk@result)#自己记得保存结果哈！
      rownames(hh) <- 1:nrow(hh)
      #display_number = c(20)  #这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
      hh <- as.data.frame(hh)[1:20, ]
      print(head(hh))
      #hh <- hh %>%
      #      dplyr::mutate(
      #            order = factor(rev(as.integer(rownames(hh))), labels = rev(hh$Description))
      #      )
      #hh <- hh[!is.na(hh$order), ]
      #P3 = ggplot(hh,aes(y=order,x=Count,fill=p.adjust))+
      P3 = ggplot(hh,aes(y=Description,x=Count,fill=p.adjust))+
            geom_bar(stat = "identity",width=0.7)+####柱子宽度
            #coord_flip()+##颠倒横纵轴
            scale_fill_gradient(low = "brown3",high ="dodgerblue4" )+#颜色自己可以换
            labs(title = "",
                 x = "Gene numbers", 
                 y = "Pathways",fill = "p.adjust",color="p.adjust") +
            theme_bw() +theme_style()
      smart_ggsave(P3,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_KK1.pdf"))
      smart_ggsave(P3,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_KK1.png"))
      ###气泡图
      
      #P4 = ggplot(hh,aes(y=order,x=Count))+
      P4 = ggplot(hh,aes(y=Description,x=Count))+
            geom_point(aes(size=Count,color=-1*p.adjust))+# 修改点的大小
            scale_color_gradient(low="dodgerblue4",high = "brown3")+
            labs(color=expression(p.adjust,size="Count"), 
                 x="Gene Count",y="Pathways",title="",fill = "p.adjust",color="p.adjust")+
            theme_bw()+ theme_style()
      smart_ggsave(P4,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_KK2.pdf"))
      smart_ggsave(P4,file=paste0(R_Enrichment_path,"/",i,"_",method,"_",up_down,"_KK2.png"))
      R_Enrichment_path_pathway = paste0(R_Enrichment_path,"/",method,"_",up_down,"_pathway")
      if (!dir.exists(R_Enrichment_path_pathway)) {
            dir.create(R_Enrichment_path_pathway, recursive = TRUE)
      }
      setwd(R_Enrichment_path_pathway)
      #hh <- hh[hh$ID != "mmu00510" & hh$ID != "mmu00512" & hh$ID != "mmu01230" & hh$ID != "mmu01232" & hh$ID != "mmu01200" & hh$ID != "mmu00533" & hh$ID !="mmu00513", ]
      for (j in 1:nrow(hh)) {
            pathway_id <- hh$ID[j]
            print(pathway_id)
            pathway_genes <- hh$geneID[j]
            #gene_list <- unlist(strsplit(pathway_genes, "/"))
            gene_list <- logFC_numbers$ENTREZID
            logFC_values <- logFC_numbers$logFC
            gene.data <- setNames(as.numeric(logFC_values), gene_list)
            result <- try({
                  pathview(gene.data = gene.data, pathway.id = pathway_id, gene.idtype = "entrez", species = database2, out.suffix = database2, kegg.native = TRUE)
            }, silent = TRUE)
            
            # 检查 result 是否为错误对象
            if (inherits(result, "try-error")) {
                  cat("无法绘制路径图:", pathway_id, "\n")
                  next  # 跳过当前条目，继续处理下一个条目
            }
      }
      
}

Edit_Functional_Statistics <- function(
				    i2,
                                    out_path,
                                    Edit_DEG_list,
                                    method,
				    database,
				    database2,
				    up_down,
				    chr_work)
{
      #######初始统计样本的类型和频率信息
      #Edit_Functional_information <- file.path(out_path, "03_Edit", str_glue("07_Edit_Func_enrich"))
      #Edit_Functional_path <- file.path(out_path, "03_Edit", str_glue("07_Edit_Func_enrich"),method)
      Edit_Functional_information <- out_path
      Edit_Functional_path <- Edit_Functional_information
      outpath = getwd()
      
      if (!dir.exists(Edit_Functional_information)) {
            dir.create(Edit_Functional_information, recursive = TRUE)
      }
      if (!dir.exists(Edit_Functional_path)) {
            dir.create(Edit_Functional_path, recursive = TRUE)
      }

      Volcano <- Edit_DEG_list
      i2 <- gsub("\\+", "_", i2)
      # 使用 get() 函数获取名为 i 的数据框
      DEG <- Volcano
      if (chr_work == "y"){
            DEG$Gene <- ifelse(grepl("chr", DEG$genes), sub(".*_(.*)", "\\1", DEG$genes), NA)
	    DEG <- select(DEG,-genes)
	    DEG <- dplyr::rename(DEG, genes = Gene )
      }
      setwd(outpath)
      Enrichment_R(DEG = DEG,
                   filepath = Edit_Functional_path,
                   i = i2,
                   database = database,
                   database2 = database2,
                   up_down = up_down,
                   method = method)
}


if(args[5] == "hsa"){
	database="org.Hs.eg.db"
	database2="hsa"
}else if(args[5] == "mmu"){
	database="org.Mm.eg.db"
	database2="mmu"
}else if(args[5] == "rno"){
	database="org.Rn.eg.db"
	database2="rno"
}

Edit_DEG <- read.csv(args[1], header = TRUE)

Edit_Functional_Statistics(
      i2=args[4],
      out_path = args[2],
      Edit_DEG_list = Edit_DEG,
      method = args[3],
      database=database,
      database2=database2,
      chr_work=args[6],
      up_down=args[7])
