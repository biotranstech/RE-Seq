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
                                               TF_file_names
                                               )
{
      ######
      #Gene_Maf_Loci_DEG <- file.path(out_path, "07_Gene_Maf", str_glue("01_Gene_Maf_Dif"))
      Gene_Maf_Loci_DEG <- out_path
      if (!dir.exists(Gene_Maf_Loci_DEG)) {
            dir.create(Gene_Maf_Loci_DEG, recursive = TRUE)
      }

      TF_information <- str_glue("{TF_file_names}")
      TF_information <- read.csv(TF_information,sep = "\t",check.names = FALSE)
      tf_list <- as.character(TF_information$Symbol)
      
      GATK_information <- GATK_information
      rawdata <- GATK_information
      reads_result <- rawdata %>%
            dplyr::group_by(Sample,Gene.refGene) %>%
            dplyr::summarise(total_genes_mutation_reads = sum(MutationReadsCount), total_genes_reads = sum(ReadsCount))
      
      reads_result <- transform(reads_result, gener_f = total_genes_mutation_reads/total_genes_reads)
      reads_result <- reads_result[,c(1,2,5)]
      matrix_reads <- dcast(data = reads_result, formula = Gene.refGene ~ Sample, value.var = "gener_f")
      matrix_reads[is.na(matrix_reads)] <- 0
      write.csv(matrix_reads,file = str_glue('{Gene_Maf_Loci_DEG}/gener_f_stat.csv'),row.names = F, quote = F)
      write.table(matrix_reads,file = str_glue('{Gene_Maf_Loci_DEG}/gener_f_stat.xls'),row.names = F, quote = F,sep="\t")
}

read_GATK <- function(file_names) {
      Initial_GATK_information <- str_glue("{file_names}")

      GATK_information <- read.csv(Initial_GATK_information,header = T,check.names = FALSE)
      return(GATK_information)
}

GATK_information = read_GATK(file_names = args[2])

Gene_Maf_Loci_DEG_list <- Maf_Difference_analysis_Statistics(
      out_path = args[3],
      GATK_information = GATK_information,
      TF_file_names = args[1])
