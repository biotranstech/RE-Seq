#Rscript - hg38_gene.length.txt GATKfile output_path
#args1  第一步GATK文件GATK file
#args2  输出路径 04_Down_analysis/03_locr/01_locr_stat

library(stringr)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)

Edit_Loci_TMB_Statistics<- function(
                               out_path,
                               gene_length,
                               GATK_information)
{
      ######
      #Edit_Loci_TMB <- file.path(out_path, "05_Edit_loci", str_glue("01_Edit_Loci_option"))
      Edit_Loci_TMB <- out_path
      if (!dir.exists(Edit_Loci_TMB)) {
            dir.create(Edit_Loci_TMB, recursive = TRUE)
      }
      rawdata <- GATK_information
      reads_result <- rawdata %>%
            dplyr::group_by(Sample,Loci_message) %>%
            dplyr::summarise(total_MutationReadsCount = sum(MutationReadsCount))
      matrix_reads <- dcast(data = reads_result, formula = Loci_message ~ Sample, value.var = "total_MutationReadsCount")
      matrix_reads[is.na(matrix_reads)] <- 0
      write.csv(matrix_reads,file = str_glue('{Edit_Loci_TMB}/locr_stat.csv'),row.names = F, quote = F)
      write.table(matrix_reads,file = str_glue('{Edit_Loci_TMB}/locr_stat.xls'),row.names = F, quote = F,sep="\t")
      
}

read_GATK <- function(file_names) {
      Initial_GATK_information <- str_glue("{file_names}")

      GATK_information <- read.csv(Initial_GATK_information,header = T,check.names = FALSE)
      return(GATK_information)
}

GATK_information = read_GATK(file_names = args[1])


Loci_Reads_information = Edit_Loci_TMB_Statistics(
      out_path = args[2],
      GATK_information = GATK_information)
