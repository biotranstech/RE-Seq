#Rscript - ../hg38_gene.length.txt RNA_count.csv output_path  
#args1      gene长度文件
#args2      RNA 表达数据框
#args3      输出路径

library(stringr)


args <- commandArgs(trailingOnly=TRUE)


RNA_TPM_Statistics<- function(
                               out_path,
                               gene_length,
                               RNA_Count)
{
      ######
      print("RNA_TPM_Statistics.start")
      RNA_TPM_path <- out_path
      RNA_TPM_path2 <- str_glue('{RNA_TPM_path}/RNA_Count')
      if (!dir.exists(RNA_TPM_path)) {
            dir.create(RNA_TPM_path, recursive = TRUE)
      }
      if (!dir.exists(RNA_TPM_path2)) {
            dir.create(RNA_TPM_path2, recursive = TRUE)
      }
      RNA_Count = read.csv(RNA_Count,check.names = FALSE)
      count_df = RNA_Count
      colnames(RNA_Count)[1] <- "Gene.refGene"
      write.csv(RNA_Count,file = str_glue('{RNA_TPM_path2}/RNA_Count.csv'), quote = F,row.names=F)
      write.table(RNA_Count,file = str_glue('{RNA_TPM_path2}/RNA_Count.xls'), quote = F,row.names=F,sep="\t")
      gene_length = read.csv(gene_length,sep = '\t', check.names = FALSE)
      names(gene_length) <- c("Gene", "efflen")
      data <- merge(count_df,gene_length,by.x="genes",by.y="Gene")
      expr_df <- data[,2:(ncol(data)-1)]
      gene_length_kb <- data$efflen / 1000
      data_rpk <- expr_df /gene_length_kb
      RNA_TPM <- t(t(data_rpk) / colSums(data_rpk) * 10^6)
      rownames(RNA_TPM)  <- data$genes
      b = data.frame(RNA_TPM,check.names = F)
      b = b %>% tibble::rownames_to_column(var = "Gene.refGene")
      write.csv(b,file = str_glue('{RNA_TPM_path}/RNA_TPM.csv'), quote = F,row.names=F)
      write.table(b,file = str_glue('{RNA_TPM_path}/RNA_TPM.xls'), quote = F,row.names=F,sep="\t")
      RNA_FPKM <- t(t(data_rpk) / colSums(expr_df) * 10^6)
      rownames(RNA_FPKM)  <- data$genes
      b = data.frame(RNA_FPKM)
      b = b %>% tibble::rownames_to_column(var = "Gene.refGene")
      write.csv(b,file = str_glue('{RNA_TPM_path}/RNA_FPKM.csv'),quote = F,row.names=F)
      write.table(b,file = str_glue('{RNA_TPM_path}/RNA_FPKM.xls'),quote = F,row.names=F,sep="\t")

      print("RNA_TPM_Statistics.done")
}

RNA_TPM_Statistics(
      out_path = args[3],
      gene_length = args[1],
      RNA_Count = args[2])
