#Rscript - hg38_gene.length.txt GATKfile outpath
#args1	参考基因组基因长度文件: hg38_gene.length.txt
#args2  第一步生成GATK文件
#args3  输出路径 
#args4	RNA表达数据

library(stringr)
library(reshape2)


args <- commandArgs(trailingOnly=TRUE)

Edit_TMB_Statistics<- function(
                    out_path,
                    gene_length,
		    RNA_file,
                    GATK_information)
{
######
      #Edit_TMB <- file.path(out_path, "03_Edit", str_glue("03_Edit_TMB"))
      Edit_TMB <- out_path
      if (!dir.exists(Edit_TMB)) {
            dir.create(Edit_TMB, recursive = TRUE)
      }
      RNA_information <- read.csv(RNA_file,header = T,check.names = FALSE, row.names = 1)
      rawdata <- GATK_information
      reads_result <- rawdata %>%
            dplyr::group_by(Sample,Gene.refGene) %>%
            dplyr::summarise(total_genes_mutation_reads = sum(MutationReadsCount), total_genes_reads = sum(ReadsCount))
      reads_result <- transform(reads_result, gener_f = total_genes_mutation_reads/total_genes_reads)
      reads_result <- reads_result[,c(1,2,5)]
      matrix_reads <- dcast(data = reads_result, formula = Gene.refGene ~ Sample, value.var = "gener_f")
      matrix_reads[is.na(matrix_reads)] <- 0
      #print(head(matrix_reads))
      total_reads <- RNA_information

      count_df = matrix_reads
      gene_length = read.csv(gene_length,sep = '\t', check.names = FALSE)
      names(gene_length) <- c("Gene", "efflen")
      data <- merge(count_df,gene_length,by.x="Gene.refGene",by.y="Gene")
      expr_df <- data[,2:(ncol(data)-1)]
      gene_length_kb <- data$efflen / 1000
      data_rpk <- expr_df /gene_length_kb
      Edit_TMB_TPM <- t(t(data_rpk) / colSums(total_reads) * 10^6)
      rownames(Edit_TMB_TPM)  <- data$Gene.refGene
      
      b = data.frame(Edit_TMB_TPM,check.names = F)
      b = b %>% tibble::rownames_to_column(var = "Gene.refGene")
      write.csv(b,file = str_glue('{Edit_TMB}/REGB_stat.csv'), quote = F, row.names=F)
      write.table(b,file = str_glue('{Edit_TMB}/REGB_stat.xls'), quote = F, row.names=F, sep = "\t")
      
      #Edit_TMB_FPKM <- t(t(data_rpk) / colSums(expr_df) * 10^6)
      #rownames(Edit_TMB_FPKM)  <- data$Gene.refGene
      #b = data.frame(Edit_TMB_FPKM)
      #b = b %>% tibble::rownames_to_column(var = "Gene.refGene")
      #write.csv(b,file = str_glue('{Edit_TMB}/Edit_TMB_FPKM.csv'),quote = F, row.names=F)
      #write.table(b,file = str_glue('{Edit_TMB}/Edit_TMB_FPKM.xls'),quote = F, row.names=F, sep = "\t")
}

read_GATK <- function(file_names) {
      Initial_GATK_information <- str_glue("{file_names}")

      GATK_information <- read.csv(Initial_GATK_information,header = T,check.names = FALSE)
      return(GATK_information)
}

GATK_information <- read_GATK(file_names = args[2])

Edit_TMB_Statistics(
      out_path = args[3],
      gene_length = args[1],
      RNA_file = args[4],
      GATK_information = GATK_information)
