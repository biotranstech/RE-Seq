#Rscript - hg38_gene.length.txt GATKfile output_path
#args1  参考基因组基因长度文件hg38_gene.length.txt
#args2  第一步GATK文件GATK file
#args3  输出路径 04_Down_analysis/03_locr/01_locr_stat
#args4	RNA表达数据文件

library(stringr)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)

Edit_Loci_TMB_Statistics<- function(
                               out_path,
                               gene_length,
                               GATK_information,
                               RNA_file)
{
      ######
      #Edit_Loci_TMB <- file.path(out_path, "05_Edit_loci", str_glue("01_Edit_Loci_option"))
      Edit_Loci_TMB <- out_path
      if (!dir.exists(Edit_Loci_TMB)) {
            dir.create(Edit_Loci_TMB, recursive = TRUE)
      }
      RNA_information <- read.csv(RNA_file,header = T,check.names = FALSE, row.names = 1)
      rawdata <- GATK_information
      matrix_reads <- dcast(data = rawdata, formula = Loci_message ~ Sample, value.var = "maf")
      matrix_reads[is.na(matrix_reads)] <- 0
      total_reads <- RNA_information
      #位点编辑负荷
      count_df = matrix_reads
      count_df$Gene <- ifelse(grepl("chr", count_df$Loci_message), 
                          sub(".*_(.*)", "\\1", count_df$Loci_message), 
                          NA)
      gene_length = read.csv(gene_length,sep = '\t', check.names = FALSE)
      names(gene_length) <- c("Gene", "efflen")
      data <- merge(count_df,gene_length,by.x="Gene",by.y="Gene")
      #print(head(data))
      expr_df <- data[,3:(ncol(data)-1)]
      #print(head(expr_df))
      gene_length_kb <- data$efflen / 1000
      data_rpk <- expr_df /gene_length_kb
      #print(head(data_rpk))
      Edit_TMB_TPM <- t(t(data_rpk) / colSums(total_reads) * 10^6)
      rownames(Edit_TMB_TPM)  <- data$Loci_message
      #print(colSums(expr_df))
      b = data.frame(Edit_TMB_TPM,check.names = F)
      b = b %>% tibble::rownames_to_column(var = "Loci_message")
      write.csv(b,file = str_glue('{Edit_Loci_TMB}/RELB_stat.csv'), quote = F,row.names=F)
      write.table(b,file = str_glue('{Edit_Loci_TMB}/RELB_stat.xls'), quote = F,row.names=F,sep="\t")
      #Edit_TMB_FPKM <- t(t(data_rpk) / colSums(expr_df) * 10^6)
      #rownames(Edit_TMB_FPKM)  <- data$Loci_message
      #b = data.frame(Edit_TMB_FPKM)
      #b = b %>% tibble::rownames_to_column(var = "Loci_message")
      #write.csv(b,file = str_glue('{Edit_Loci_TMB}/Edit_Loci_TMB_FPKM.csv'),quote = F,row.names=F)
      #write.table(b,file = str_glue('{Edit_Loci_TMB}/Edit_Loci_TMB_FPKM.xls'),quote = F,row.names=F,sep="\t")
}

read_GATK <- function(file_names) {
      Initial_GATK_information <- str_glue("{file_names}")

      GATK_information <- read.csv(Initial_GATK_information,header = T,check.names = FALSE)
      return(GATK_information)
}

GATK_information = read_GATK(file_names = args[2])


Edit_Loci_TMB_Statistics(
      out_path = args[3],
      gene_length = args[1],
      RNA_file = args[4],
      GATK_information = GATK_information)
