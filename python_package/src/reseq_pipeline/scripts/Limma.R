#Rscript - output_path matrix_file  TF_file group_file locr.csv
#args1  output_path:./03_Edit/05_Edit_dif/
#args2  matrix file:count file
#args3  TF file
#args4  group file
#args5  number stat file

library(stringr)
library(ggplot2)
library(dplyr)
library(reshape2)
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
				     log2,
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
      
      ##Limma
      expr_data <- Edit_TMB
      expr_data <- expr_data[,-1]
      row.names(expr_data) <- Edit_TMB[,1]

      test=t(expr_data)
      test=data.frame(species=rownames(test))
      group_te = group_te[match(group_te$species,test$species),]

      a = t(group_te)
      colnames(a) = group_te$species
      group_names = a[-1,]
      group_names <- gsub("\\+", "_", group_names)
      unique_groups <- sort(unique(group_names))
      group_counts <- table(group_names)
      #group <- rep(names(group_counts), times = group_counts)
      group <- rep(unique(group_names), times = group_counts)
      design <- model.matrix(~0+factor(group))
      colnames(design) <- levels(factor(group))
      rownames(design) <- colnames(expr_data)
      expr_data <- expr_data[which(rowSums(expr_data)!=0),] 
      if (log2 == "y"){
          expr_data = log2(expr_data) #log化处理
          expr_data[expr_data == -Inf] = 0
      }else{
          expr_data = scale(expr_data) #标准化处理
      }
     
      #构建比对矩阵
      comparisons <- combn(unique_groups, 2, FUN = function(x) paste(x[2], x[1], sep = "_vs_") )  
      print(comparisons)
      DEG_limma_list <- list() 
      for (cmp in comparisons) {
            groups <- unlist(strsplit(cmp, "_vs_"))  
            s_file <- glue("{dir_path}/{cmp}.csv")
	    number_stat <- read.csv(s_file, header=TRUE,check.names = FALSE)
	    colnames(number_stat)[1] <- "genes"
            contrast_string <- paste(groups, collapse = "-") 
            contrast <- makeContrasts(contrasts = contrast_string, levels = design)
	    print(cmp)
	    print(contrast)
            fit <- lmFit(expr_data, design)
            fit2 <- contrasts.fit(fit, contrast)
            fit2 <- eBayes(fit2)
            DEG <- topTable(fit2, coef = 1, n = Inf, sort.by = "logFC")
            DEG <- na.omit(DEG)
            DEG<- DEG %>% rownames_to_column("genes")
            DEG$regulate <- ifelse(DEG$P.Value > 0.05, "Undiff",
                                   ifelse(DEG$logFC > 1, "Up",
                                          ifelse(DEG$logFC < -1, "Down", "Undiff")))
            DEG$TF <- sapply(DEG$genes, function(gene) {
                  if (gene %in% tf_list) {
                        return('TF')
                  } else {
                        return('-')
                  }
            })
            
            DEG_limma_list[[cmp]] <- DEG
            write.csv(DEG, file = paste0(Edit_DEG,"/",cmp, "_Limma.csv"), row.names = FALSE)
	    DEG_1 <- DEG[,c("genes", "logFC", "P.Value", "adj.P.Val", "regulate", "TF")]
	    r_data <- merge(number_stat,DEG_1,by.x="genes",by.y="genes")
	    r_data <- r_data[which(r_data$regulate %in% c("Down" ,"Up")),]
            write.table(r_data, file = paste0(Edit_DEG,"/",cmp, "_Limma.xls"), row.names = FALSE, sep = "\t")
      }
}

read_matrix <- function(file_names) {
      matrix_info <- str_glue("{file_names}")

      Matrix <- read.csv(matrix_info,header = T,check.names = FALSE)
      return(Matrix)
}


Matrix <- read_matrix(file_names = args[2])

file1 = unlist(strsplit(args[2], "/"))
file2 = file1[length(file1)]
if (file2 == "locr_f_stat.csv" | file2 == "gener_f_stat.csv"){
      log2="n"
}else{
      log2="y"
}

Edit_Difference_analysis_Statistics(
      out_path = args[1],
      Edit_TMB = Matrix,
      log2=log2,
      TF_file_names = args[3],
      group_file = args[4],
      stat_file = args[2])
