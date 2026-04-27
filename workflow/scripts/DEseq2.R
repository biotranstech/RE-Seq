#Rscript - output_path matrix_file  TF_file group_file 
#args1  output_path:./03_Edit/05_Edit_dif/
#args2  matrix file:count file
#args3  TF file
#args4  group file

library(stringr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(DESeq2)
library(limma)
library(tidyverse)
library(glue)
library(coin)

args <- commandArgs(trailingOnly=TRUE)

Edit_Difference_analysis_Statistics <- function(
                                     out_path,
                                     Edit_TMB,
                                     TF_file_names,
                                     group_file,
                                     stat_file)
{
######
      #Edit_DEG <- file.path(out_path, "03_Edit", str_glue("05_Edit_dif"))
      Edit_DEG <- out_path
       
      if (!dir.exists(Edit_DEG)) {
            dir.create(Edit_DEG, recursive = TRUE)
      }

      TF_information <- read.csv(TF_file_names,sep = "\t",check.names = FALSE)
      tf_list <- as.character(TF_information$Symbol)
      group_te = read.table(group_file, header=TRUE)
      dir_path <- dirname(stat_file)
      #number_stat <- read.csv(stat_file, header=TRUE,check.names = FALSE)
      #colnames(number_stat)[1] <- "genes"
      
      ##DEseq2
      expr_data <- Edit_TMB
      matrix_reads <- expr_data[,-1]
      row.names(matrix_reads) <- Edit_TMB[,1]


      colData <- read.table(group_file, header=TRUE)
      colnames(colData) <- c("a","condition")

      a = t(group_te)
      colnames(a) = group_te$species
      group_names = a[-1,]
      group_names <- gsub("\\+", "_", group_names)
      unique_groups <- sort(unique(group_names))
      #group_counts <- table(group_names)
      #condition <- factor(rep(names(group_counts), times = group_counts))
      #colData <- data.frame(row.names=colnames(matrix_reads), condition) 
      
      aa = data.frame(colnames(matrix_reads))
      colnames(aa) = "a"
      colData2 = colData[match(aa$a,colData$a),]
      rownames(colData2) = colData2$a
      colData <- data.frame(row.names=colData2$a, colData2$condition)
      colnames(colData) = "condition"
      print(colData)
      #构建比对矩阵
      comparisons <- combn(unique_groups, 2, FUN = function(x) paste(x[2], x[1], sep = "_vs_") ) 
      dds <- DESeqDataSetFromMatrix(countData = matrix_reads, colData = colData, design = ~ condition)
      dds <- DESeq(dds)
      # 初始化一个空列表来存储结果
      DEG_DEseq2_list <- list()
      for (cmp in comparisons) {
            groups <- unlist(strsplit(cmp, "_vs_"))
            group1 <- groups[1]
            group2 <- groups[2]
	    s_file <- glue("{dir_path}/{cmp}.csv")
	    number_stat <- read.csv(s_file, header=TRUE,check.names = FALSE)
	    colnames(number_stat)[1] <- "genes"
            contrast <- c("condition", group1, group2)
            res <- results(dds, contrast = contrast, cooksCutoff = FALSE)
            res <- na.omit(res)
            res$genes <- row.names(res)
            colnames(res)[colnames(res) == "log2FoldChange"] <- "logFC"
            colnames(res)[colnames(res) == "pvalue"] <- "P.Value"
            res$regulate <- ifelse(res$P.Value > 0.05, "Undiff",
                                   ifelse(res$logFC > 1, "Up",
                                          ifelse(res$logFC < -1, "Down", "Undiff")))
            res$TF <- sapply(res$genes, function(gene) {
                  if (gene %in% tf_list) {
                        return('TF')
                  } else {
                        return('-')
                  }
            })
            # 构建比较名称
            comparison_name <- paste(cmp)
            DEG_DEseq2_list[[cmp]] <- as.data.frame(res)
            output_file <- glue("{Edit_DEG}/{cmp}_DESeq2.csv")
            output_file2 <- glue("{Edit_DEG}/{cmp}_DESeq2.xls")
            df = as.data.frame(res)
            df = df %>% relocate(genes)
            write.csv(df, file = output_file, row.names = FALSE)
            DEG_1 <- df[,c("genes", "logFC", "P.Value", "padj", "regulate", "TF")]
            r_data <- merge(number_stat,DEG_1,by.x="genes",by.y="genes")
            r_data <- r_data[which(r_data$regulate %in% c("Down" ,"Up")),]
            write.table(r_data, file = output_file2, row.names = FALSE, sep = "\t")
            #write.csv(as.data.frame(res), file = output_file, row.names = FALSE)
            #write.table(as.data.frame(res), file = output_file2, row.names = FALSE)
      }    
}

read_matrix <- function(file_names) {
      matrix_info <- str_glue("{file_names}")

      Matrix <- read.csv(matrix_info,header = T,check.names = FALSE)
      return(Matrix)
}


Matrix <- read_matrix(file_names = args[2])


Edit_Difference_analysis_Statistics(
      out_path = args[1],
      Edit_TMB = Matrix,
      TF_file_names = args[3],
      group_file = args[4],
      stat_file = args[2])
