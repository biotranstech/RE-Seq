import os.path
from collections import OrderedDict

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")
BIN = os.path.join(ROOT, "bin")
TEMPLATES = os.path.join(ROOT, "template")
SCRIPTS = os.path.join(ROOT, "scripts")

ref_path = {"hg38": "/home/data/jc1/pipeline/01_RE-seq/Refer/hg38/hg38",
			"mm10": "/home/data/jc1/pipeline/01_RE-seq/Refer/mm10/mm10",
			"rn6":  "/home/data/jc1/pipeline/01_RE-seq/Refer/rn6/rn6"
}

ref_fasta = {"hg38": "/home/data/jc1/pipeline/01_RE-seq/Refer/hg38/hg38.fa",
			"mm10": "/home/data/jc1/pipeline/01_RE-seq/Refer/mm10/mm10.fa",
			"rn6":  "/home/data/jc1/pipeline/01_RE-seq/Refer/rn6/rn6.fa"
}

ref_bed = {"hg38": "/home/data/jc1/pipeline/01_RE-seq/Refer/hg38/hg38.bed",
			"mm10": "/home/data/jc1/pipeline/01_RE-seq/Refer/mm10/mm10.bed",
			"rn6":  "/home/data/jc1/pipeline/01_RE-seq/Refer/rn6/rn6.bed"
}

ref_gtf = {"hg38": "/home/data/jc1/pipeline/01_RE-seq/Refer/hg38/hg38.refGene.gtf",
			"mm10": "/home/data/jc1/pipeline/01_RE-seq/Refer/mm10/mm10.refGene.gtf",
			"rn6":  "/home/data/jc1/pipeline/01_RE-seq/Refer/rn6/rn6.refGene.gtf"
}

ref_repeat = {"hg38": "/home/data/jc1/pipeline/01_RE-seq/Refer/hg38/hg38.repeat.txt",
			"mm10": "/home/data/jc1/pipeline/01_RE-seq/Refer/mm10/mm10.repeat.txt",
			"rn6":  "/home/data/jc1/pipeline/01_RE-seq/Refer/rn6/rn6.repeat.txt"
}

ref_dbsnp = {"hg38": "/home/data/jc1/pipeline/01_RE-seq/Refer/hg38/hg38.dbsnp.txt",
			"mm10": "/home/data/jc1/pipeline/01_RE-seq/Refer/mm10/mm10.dbsnp.txt",
			"rn6":  "/home/data/jc1/pipeline/01_RE-seq/Refer/rn6/rn6.dbsnp.txt"
}

ref_alu = {"hg38": "/home/data/jc1/pipeline/01_RE-seq/Refer/hg38/hg38.alu.bed",
			"mm10": "/home/data/jc1/pipeline/01_RE-seq/Refer/mm10/mm10.alu.bed",
			"rn6":  "/home/data/jc1/pipeline/01_RE-seq/Refer/rn6/rn6.alu.bed"
}

ref_snp_d = {"hg38": "/home/data/jc1/pipeline/01_RE-seq/Refer/hg38/snp151CodingDbSnp.txt",
			"mm10": "/home/data/jc1/pipeline/01_RE-seq/Refer/mm10/snp151CodingDbSnp.txt",
			"rn6":  "/home/data/jc1/pipeline/01_RE-seq/Refer/rn6/snp151CodingDbSnp.txt"
}

ref_tf = {"hg38": "/home/data/jc1/pipeline/01_RE-seq/Refer/hg38/Homo_sapiens_TF.txt",
			"mm10": "/home/data/jc1/pipeline/01_RE-seq/Refer/mm10/Mus_musculus_TF.txt",
			"rn6":  "/home/data/jc1/pipeline/01_RE-seq/Refer/rn6/Rattus_norvegicus_TF.txt"
}

ref_gene_length = {"hg38": "/home/data/jc1/pipeline/01_RE-seq/Refer/hg38/hg38_gene.length.txt",
			"mm10": "/home/data/jc1/pipeline/01_RE-seq/Refer/mm10/mm10_lengths.txt",
			"rn6":  "/home/data/jc1/pipeline/01_RE-seq/Refer/rn6/rn6_gene.length.txt"
}

ref_name = {"hg38": "GRCh38(hg38)",
                "mm10": "GRCm38(mm10)",
                "rn6": "RN6"
}

FASTP_BIN = "/home/data/jc1/miniconda3/bin/fastp"
FASTQC_BIN = "/home/data/jc1/miniconda3/bin/fastqc"
PYTHON_BIN = "/home/data/jc1/miniconda3/bin/python"
PERL_BIN = "/home/data/jc1/miniconda3/bin/perl"
JAVA_BIN = "/home/data/jc1/miniconda3/bin"
HISAT2_BIN = "/home/data/jc1/miniconda3/bin/hisat2"
SAMTOOLS_BIN = "/home/data/jc1/miniconda3/bin/samtools"
SPRINT_BIN = "/home/data/jc1/pipeline/01_RE-seq/scripts/SPRINT-0.1.8/bin/sprint"
SPRINT_BAM = "/home/data/jc1/pipeline/01_RE-seq/scripts/SPRINT-0.1.8/bin/sprint_from_bam"
REDITOOL_BIN = "/home/data/jc1/pipeline/01_RE-seq/scripts/reditools2.0-master/src/cineca/reditools.py"
GATK_BIN = "/home/data/jc1/miniconda3/bin/gatk"
ANNOVAR_BIN = "/home/data/jc1/software/annovar"
BGZIP = "/home/data/jc1/miniconda3/bin/bgzip"
TABIX = "/home/data/jc1/miniconda3/bin/tabix"
HTSEQ = "/home/data/jc1/miniconda3/bin/htseq-count"
RSCRIPT = "/home/data/jc1/miniconda3/envs/r_environment/bin/Rscript"
