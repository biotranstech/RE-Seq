#Rscript - hg38_gene.length.txt GATKfile outpath
#args1  第一步生成GATK文件
#args2  输出路径 

library(stringr)
library(reshape2)


args <- commandArgs(trailingOnly=TRUE)

Edit_TMB_Statistics<- function(
                    out_path,
                    gene_length,
                    GATK_information)
{
######

      Edit_TMB <- out_path
      if (!dir.exists(Edit_TMB)) {
            dir.create(Edit_TMB, recursive = TRUE)
      }
      rawdata <- GATK_information
      reads_result <- rawdata %>%
            dplyr::group_by(Sample,Gene.refGene) %>%
            dplyr::summarise(total_MutationReadsCount = sum(MutationReadsCount))
      matrix_reads <- dcast(data = reads_result, formula = Gene.refGene ~ Sample, value.var = "total_MutationReadsCount")
      matrix_reads[is.na(matrix_reads)] <- 0
      write.csv(matrix_reads,file = str_glue('{Edit_TMB}/gener_stat.csv'),row.names = F, quote = F)
      write.table(matrix_reads,file = str_glue('{Edit_TMB}/gener_stat.xls'),row.names = F, quote = F, sep = "\t")
}

read_GATK <- function(file_names) {
      Initial_GATK_information <- str_glue("{file_names}")

      GATK_information <- read.csv(Initial_GATK_information,header = T,check.names = FALSE)
      return(GATK_information)
}

GATK_information <- read_GATK(file_names = args[1])

Edit_TMB <- Edit_TMB_Statistics(
      out_path = args[2],
      GATK_information = GATK_information)

