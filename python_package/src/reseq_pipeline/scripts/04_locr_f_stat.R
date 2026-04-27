#Rscript - output_path  TF_file GATK_file group_file outpath
#args1  TF文件
#args2  处理后GATK文件
#args3  输出路径

library(stringr)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)

Maf_Difference_analysis_Statistics <- function(
                                               out_path,
                                               GATK_information,
                                               TF_file_names)
{
      ######
      Maf_Loci_DEG <- out_path
      if (!dir.exists(Maf_Loci_DEG)) {
            dir.create(Maf_Loci_DEG, recursive = TRUE)
      }

      TF_information <- str_glue("{TF_file_names}")
      TF_information <- read.csv(TF_information,sep = "\t",check.names = FALSE)
      tf_list <- as.character(TF_information$Symbol)
      
      GATK_information <- GATK_information
      rawdata <- GATK_information
      matrix_reads <- dcast(data = rawdata, formula = Loci_message ~ Sample, value.var = "maf")
      matrix_reads[is.na(matrix_reads)] <- 0
      write.csv(matrix_reads,file = str_glue('{Maf_Loci_DEG}/locr_f_stat.csv'),row.names = F, quote = F)
      write.table(matrix_reads,file = str_glue('{Maf_Loci_DEG}/locr_f_stat.xls'),row.names = F, quote = F,sep="\t")
}

read_GATK <- function(file_names) {
      Initial_GATK_information <- str_glue("{file_names}")

      GATK_information <- read.csv(Initial_GATK_information,header = T,check.names = FALSE)
      return(GATK_information)
}

GATK_information = read_GATK(file_names = args[2])

Maf_Loci_DEG_list <- Maf_Difference_analysis_Statistics(
      out_path = args[3],
      GATK_information = GATK_information,
      TF_file_names = args[1]
      )
