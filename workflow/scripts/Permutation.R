#Rscript - output_path matrix_file TPM_file TF_file group_file
#args1  完整输出路径output_path:./03_Edit/05_Edit_dif/
#args2  数值数据框 matrix file:count file
#args3  计算的 TPM file
#args4  参考基因组TF文件,TF file
#args5  分组文件group file

library(stringr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyverse)
library(glue)
library(coin)

args <- commandArgs(trailingOnly=TRUE)

Edit_Difference_analysis_Statistics <- function(
                                     out_path,
                                     Edit_TMB2,
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
      dir_path <- dirname(stat_file)
      group_te = read.table(group_file, header=TRUE)
      
      ##Permutation   

      Data <- as.data.frame(Edit_TMB2)
      print(head(Data))

      group_te2 = group_te[match(colnames(Data),group_te$species),]
      colData = as.data.frame(group_te2[,-1])
      rownames(colData) = group_te2$species
      colnames(colData) = "condition"
      print(colData)
      a = t(group_te)
      colnames(a) = group_te$species
      group_names = a[-1,]
      group_names <- gsub("\\+", "_", group_names)
      unique_groups <- sort(unique(group_names))
      #group_counts <- table(group_names)
      #condition <- factor(rep(names(group_counts), times = group_counts))
      #colData <- data.frame(row.names=colnames(Data), condition)
      #print(colData)
      
      Data$genes <- row.names(Data)
      grades_long <- melt(Data, id.vars = "genes", measure.vars = setdiff(colnames(Data), "genes"),variable.name = "Individual", value.name = "TPM")
      #output_file <- glue("{Edit_DEG}/test2_Permutation.csv")
      #write.csv(as.data.frame(grades_long), file = output_file, row.names = FALSE)
      grades_long <- grades_long %>%
            dplyr::mutate(Group = colData[Individual, "condition"])
      grades_long$Group <- as.factor(grades_long$Group)
      genes_numbers <- unique(grades_long$genes)
      groups <- unique(grades_long$Group)
      #print(head(grades_long))
      #output_file <- glue("{Edit_DEG}/test_Permutation.csv")
      #write.csv(as.data.frame(grades_long), file = output_file, row.names = FALSE)
      
      comparisons <- combn(unique_groups, 2, FUN = function(x) paste(x[2], x[1], sep = "_vs_") ) 
      # 存储每个组别比较结果的列表
      DEG_Permutation_list <- list()
      for (cmp in comparisons) {
            groups <- unlist(strsplit(cmp, "_vs_")) 
            group1 <- groups[1]
            group2 <- groups[2]

	    s_file <- glue("{dir_path}/{cmp}.csv")
	    print(s_file)
	    number_stat <- read.csv(s_file, header=TRUE,check.names = FALSE)
	    colnames(number_stat)[1] <- "genes"
            
            # 创建一个空的数据框，用于存储当前组别比较的结果
            comparison_results <- data.frame()
            
            for (gene in genes_numbers) {
                  gene_data <- subset(grades_long, genes == gene)
                  
                  # 子集选择并去除多余因子水平
                  comparison_data <- droplevels(subset(gene_data, Group %in% c(group1, group2)))
                  
                  # 计算每组的平均值
                  group_means <- comparison_data %>%
                        dplyr::group_by(Group) %>%
                        dplyr::summarise(mean_TPM = mean(TPM))
                  
                  # 提取各组的平均值
                  group1_means <- as.numeric(group_means[group_means$Group == group1, "mean_TPM"])
                  group2_means <- as.numeric(group_means[group_means$Group == group2, "mean_TPM"])
                  epsilon = 0.00
                  # 计算logFC，避免除以0的情况
                  if (group2_means == 0) {
                        logfc <- Inf 
                  } else if (group1_means == 0){
		        logfc <- -Inf
		  } else {
                        logfc <- log2((group1_means + epsilon) / (group2_means + epsilon))
                  }
                  
                  # 进行独立性检验
                  test_result <- independence_test(TPM ~ Group, data = comparison_data)
                  p_value <- as.numeric(pvalue(test_result))
                  
                  # 保存当前基因的结果到数据框中
                  comparison_results <- rbind(comparison_results, data.frame(
                        genes = gene,
                        Group1 = group1,
                        Group2 = group2,
                        P.Value = p_value,
                        Group1_means = group1_means,
                        Group2_means = group2_means,
                        logFC = logfc,
                        stringsAsFactors = FALSE
                  ))
            }
            
            # 将当前组别比较的结果添加到列表中
            comparison_results$regulate <- ifelse(comparison_results$P.Value > 0.05, "Undiff",
                                                  ifelse(comparison_results$logFC > 1, "Up",
                                                         ifelse(comparison_results$logFC < -1, "Down", "Undiff")))
            comparison_results$TF <- sapply(comparison_results$genes, function(gene) {
                  if (gene %in% tf_list) {
                        return('TF')
                  } else {
                        return('-')
                  }
            })
            
            output_file <- glue("{Edit_DEG}/{cmp}_Permutation.csv")
            output_file2 <- glue("{Edit_DEG}/{cmp}_Permutation.xls")
            DEG_Permutation_list[[cmp]] <- comparison_results
            write.csv(as.data.frame(comparison_results), file = output_file, row.names = FALSE)
            test_data <- as.data.frame(comparison_results)
            test_data <- test_data[,c("genes", "logFC", "P.Value","regulate","TF")]
            r_data <- merge(number_stat,test_data,by.x="genes",by.y="genes")
            r_data <- r_data[which(r_data$regulate %in% c("Down" ,"Up")),]
            write.table(r_data, file = output_file2, row.names = FALSE,sep="\t")
      }
}

read_matrix <- function(file_names) {
      matrix_info <- str_glue("{file_names}")

      Matrix <- read.csv(matrix_info,header = T,check.names = FALSE, row.names = 1)
      return(Matrix)
}


Matrix <- read_matrix(file_names = args[2])


Edit_DEG_list <- Edit_Difference_analysis_Statistics(
      out_path = args[1],
      Edit_TMB2 = Matrix,
      TF_file_names = args[3],
      group_file = args[4],
      stat_file = args[2])
