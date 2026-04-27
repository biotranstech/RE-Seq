#!/usr/bin/env python
# -*- coding:utf-8 -*-

__AUTHOR__ = "2469823707@qq.com"
__VERSION__ = "v1.0_2024"
#change pca boxplot cor_heatmap name

import os
import sys
import re
import importlib
import glob
import shutil
import jinja2
import logging
import argparse
from time import localtime
from docxtpl import DocxTemplate

#编码设置为utf8
importlib.reload(sys)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
#BIN_DIR = "/Work/user/tangdong/pipline/backup_pipline/RNASeq/Pipeline/noRef_Isoseq/Report"
#BIN_DIR = "/Work/user/lvzhaopeng/pipeline/RNASeq/Report2.0"
#sys.path.insert(1, BIN_DIR) 

from report_utils import safecopy,readtbl,schfile,loopcopy,readraw,loopcopy20,readtbl2,readtbl3,readtbl_go

logging.basicConfig(level=logging.DEBUG, format="[%(levelname)s] %(asctime)s: %(message)s")
logger = logging.getLogger('Bioyigene RNA-Seq Report')


def copyfigs(project_dir, report_dir):
    if os.path.exists(report_dir):
        shutil.rmtree(report_dir)
    shutil.copytree(os.path.join(BASE_DIR, "template"), report_dir)
    if ASAN:
        safecopy("%s/report.html" % BASE_DIR, report_dir)
    else:
        safecopy("%s/report_noAS.html" % BASE_DIR, report_dir)
    imgdir = os.path.join(report_dir, "Result")
    loopcopy("%s/01_QC/QC_result/*.base_quality.png" % project_dir, os.path.join(imgdir, "01_QC"))
    loopcopy("%s/01_QC/QC_result/*.base_content.png" % project_dir, os.path.join(imgdir, "01_QC"))
    #loopcopy("%s/02_Hisat/mapping_rate.txt" % project_dir, os.path.join(imgdir, "02_Hisat"))
    #loopcopy("%s/02_Hisat/*.hisat_content.png" % project_dir, os.path.join(imgdir, "02_Hisat"))
    
    loopcopy("%s/04_Down_analysis/01_REseq_detection/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "01_REseq_detection"))
    loopcopy("%s/04_Down_analysis/02_REseq_PCA/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "02_REseq_PCA"))
    loopcopy("%s/04_Down_analysis/03_locr/01_locr_stat/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "03_locr", "01_locr_stat"))
    loopcopy("%s/04_Down_analysis/03_locr/02_locr_diff/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "03_locr", "02_locr_diff"))
    loopcopy("%s/04_Down_analysis/03_locr/02_locr_diff/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "03_locr", "02_locr_diff"))
    loopcopy("%s/04_Down_analysis/03_locr/02_locr_diff/veen/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "03_locr", "02_locr_diff", "veen"))
    loopcopy("%s/04_Down_analysis/03_locr/03_locr_vol/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "03_locr", "03_locr_vol"))
    loopcopy("%s/04_Down_analysis/03_locr/04_locr_pheatmap/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "03_locr", "04_locr_pheatmap"))
    loopcopy("%s/04_Down_analysis/03_locr/05_Enrich_R/*/R_Enrichment/*GO*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "03_locr", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/03_locr/05_Enrich_R/*/R_Enrichment/*KK*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "03_locr", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/03_locr/05_Enrich_R/*/R_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "03_locr", "05_Enrich_R", "R_Enrichment"))
    loopcopy20("%s/04_Down_analysis/03_locr/05_Enrich_R/*/GSEA_Enrichment/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "03_locr", "05_Enrich_R", "GSEA_Enrichment"))
    loopcopy("%s/04_Down_analysis/03_locr/05_Enrich_R/*/GSEA_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "03_locr", "05_Enrich_R", "GSEA_Enrichment"))
    """
    loopcopy("%s/04_Down_analysis/04_locr_f/01_locr_f_stat/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "01_locr_f_stat"))
    loopcopy("%s/04_Down_analysis/04_locr_f/02_locr_f_diff/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "02_locr_f_diff"))
    loopcopy("%s/04_Down_analysis/04_locr_f/02_locr_f_diff/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "02_locr_f_diff"))
    loopcopy("%s/04_Down_analysis/04_locr_f/02_locr_f_diff/veen/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "02_locr_f_diff", "veen"))
    loopcopy("%s/04_Down_analysis/04_locr_f/03_locr_f_vol/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "03_locr_f_vol"))
    loopcopy("%s/04_Down_analysis/04_locr_f/04_locr_f_pheatmap/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "04_locr_f_pheatmap"))
    loopcopy("%s/04_Down_analysis/04_locr_f/05_Enrich_R/*/R_Enrichment/*GO*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/04_locr_f/05_Enrich_R/*/R_Enrichment/*KK*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/04_locr_f/05_Enrich_R/*/R_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "05_Enrich_R", "R_Enrichment"))
    loopcopy20("%s/04_Down_analysis/04_locr_f/05_Enrich_R/*/GSEA_Enrichment/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "05_Enrich_R", "GSEA_Enrichment"))
    loopcopy("%s/04_Down_analysis/04_locr_f/05_Enrich_R/*/GSEA_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "05_Enrich_R", "GSEA_Enrichment"))
    """
    loopcopy("%s/04_Down_analysis/05_RELB/01_RELB_stat/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "05_RELB", "01_RELB_stat"))
    loopcopy("%s/04_Down_analysis/05_RELB/02_RELB_diff/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "05_RELB", "02_RELB_diff"))
    loopcopy("%s/04_Down_analysis/05_RELB/02_RELB_diff/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "05_RELB", "02_RELB_diff"))
    loopcopy("%s/04_Down_analysis/05_RELB/03_RELB_vol/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "05_RELB", "03_RELB_vol"))
    loopcopy("%s/04_Down_analysis/05_RELB/04_RELB_pheatmap/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "05_RELB", "04_RELB_pheatmap"))
    loopcopy("%s/04_Down_analysis/05_RELB/05_Enrich_R/*/R_Enrichment/*GO*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "05_RELB", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/05_RELB/05_Enrich_R/*/R_Enrichment/*KK*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "05_RELB", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/05_RELB/05_Enrich_R/*/R_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "05_RELB", "05_Enrich_R", "R_Enrichment"))
    loopcopy20("%s/04_Down_analysis/05_RELB/05_Enrich_R/*/GSEA_Enrichment/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "05_RELB", "05_Enrich_R", "GSEA_Enrichment"))
    loopcopy("%s/04_Down_analysis/05_RELB/05_Enrich_R/*/GSEA_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "05_RELB", "05_Enrich_R", "GSEA_Enrichment"))

    loopcopy("%s/04_Down_analysis/06_gener/01_gener_stat/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "06_gener", "01_gener_stat"))
    loopcopy("%s/04_Down_analysis/06_gener/02_gener_diff/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "06_gener", "02_gener_diff"))
    loopcopy("%s/04_Down_analysis/06_gener/02_gener_diff/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "06_gener", "02_gener_diff"))
    loopcopy("%s/04_Down_analysis/06_gener/02_gener_diff/veen/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "06_gener", "02_gener_diff", "veen"))
    loopcopy("%s/04_Down_analysis/06_gener/03_gener_vol/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "06_gener", "03_gener_vol"))
    loopcopy("%s/04_Down_analysis/06_gener/04_gener_pheatmap/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "06_gener", "04_gener_pheatmap"))
    loopcopy("%s/04_Down_analysis/06_gener/05_Enrich_R/*/R_Enrichment/*GO*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "06_gener", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/06_gener/05_Enrich_R/*/R_Enrichment/*KK*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "06_gener", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/06_gener/05_Enrich_R/*/R_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "06_gener", "05_Enrich_R", "R_Enrichment"))
    loopcopy20("%s/04_Down_analysis/06_gener/05_Enrich_R/*/GSEA_Enrichment/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "06_gener", "05_Enrich_R", "GSEA_Enrichment"))
    loopcopy("%s/04_Down_analysis/06_gener/05_Enrich_R/*/GSEA_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "06_gener", "05_Enrich_R", "GSEA_Enrichment"))
    """
    loopcopy("%s/04_Down_analysis/07_gener_f/01_gener_f_stat/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "01_gener_f_stat"))
    loopcopy("%s/04_Down_analysis/07_gener_f/02_gener_f_diff/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "02_gener_f_diff"))
    loopcopy("%s/04_Down_analysis/07_gener_f/02_gener_f_diff/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "02_gener_f_diff"))
    loopcopy("%s/04_Down_analysis/07_gener_f/02_gener_f_diff/veen/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "02_gener_f_diff", "veen"))
    loopcopy("%s/04_Down_analysis/07_gener_f/03_gener_f_vol/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "03_gener_f_vol"))
    loopcopy("%s/04_Down_analysis/07_gener_f/04_gener_f_pheatmap/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "04_gener_f_pheatmap"))
    loopcopy("%s/04_Down_analysis/07_gener_f/05_Enrich_R/*/R_Enrichment/*GO*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/07_gener_f/05_Enrich_R/*/R_Enrichment/*KK*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/07_gener_f/05_Enrich_R/*/R_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "05_Enrich_R", "R_Enrichment"))
    loopcopy20("%s/04_Down_analysis/07_gener_f/05_Enrich_R/*/GSEA_Enrichment/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "05_Enrich_R", "GSEA_Enrichment"))
    loopcopy("%s/04_Down_analysis/07_gener_f/05_Enrich_R/*/GSEA_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "05_Enrich_R", "GSEA_Enrichment"))
    """
    loopcopy("%s/04_Down_analysis/08_REGB/01_REGB_stat/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "01_REGB_stat"))
    loopcopy("%s/04_Down_analysis/08_REGB/02_REGB_diff/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "02_REGB_diff"))
    loopcopy("%s/04_Down_analysis/08_REGB/02_REGB_diff/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "02_REGB_diff"))
    loopcopy("%s/04_Down_analysis/08_REGB/03_REGB_vol/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "03_REGB_vol"))
    loopcopy("%s/04_Down_analysis/08_REGB/04_REGB_pheatmap/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "04_REGB_pheatmap"))
    loopcopy("%s/04_Down_analysis/08_REGB/05_Enrich_R/*/R_Enrichment/*GO*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/08_REGB/05_Enrich_R/*/R_Enrichment/*KK*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/08_REGB/05_Enrich_R/*/R_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"))
    loopcopy20("%s/04_Down_analysis/08_REGB/05_Enrich_R/*/GSEA_Enrichment/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "GSEA_Enrichment"))
    loopcopy("%s/04_Down_analysis/08_REGB/05_Enrich_R/*/GSEA_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "GSEA_Enrichment"))

    loopcopy("%s/04_Down_analysis/09_REB/*" % project_dir, os.path.join(imgdir, "04_Down_analysis", "09_REB"))

    loopcopy("%s/04_Down_analysis/10_RNA/01_RNA_stat/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "10_RNA", "01_RNA_stat"))
    loopcopy("%s/04_Down_analysis/10_RNA/01_RNA_stat/RNA_PCA/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "10_RNA", "01_RNA_stat"))
    loopcopy("%s/04_Down_analysis/10_RNA/02_RNA_diff/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "10_RNA", "02_RNA_diff"))
    loopcopy("%s/04_Down_analysis/10_RNA/02_RNA_diff/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "10_RNA", "02_RNA_diff"))
    loopcopy("%s/04_Down_analysis/10_RNA/02_RNA_diff/veen/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "10_RNA", "02_RNA_diff", "veen"))
    loopcopy("%s/04_Down_analysis/10_RNA/03_RNA_vol/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "10_RNA", "03_RNA_vol"))
    loopcopy("%s/04_Down_analysis/10_RNA/04_RNA_pheatmap/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "10_RNA", "04_RNA_pheatmap"))
    loopcopy("%s/04_Down_analysis/10_RNA/05_Enrich_R/*/R_Enrichment/*GO*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "10_RNA", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/10_RNA/05_Enrich_R/*/R_Enrichment/*KK*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "10_RNA", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/10_RNA/05_Enrich_R/*/R_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "10_RNA", "05_Enrich_R", "R_Enrichment"))
    loopcopy20("%s/04_Down_analysis/10_RNA/05_Enrich_R/*/GSEA_Enrichment/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "10_RNA", "05_Enrich_R", "GSEA_Enrichment"))
    loopcopy("%s/04_Down_analysis/10_RNA/05_Enrich_R/*/GSEA_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "10_RNA", "05_Enrich_R", "GSEA_Enrichment"))

    if ASAN:
        loopcopy("%s/04_Down_analysis/11_AS/01_AS_analysis_base/*" % project_dir, os.path.join(imgdir, "08_AS", "01_AS_analysis_base"))
        loopcopy("%s/04_Down_analysis/12_association/01_RNA_locr/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "12_association", "01_RNA_locr"))
        #loopcopy("%s/04_Down_analysis/12_association/02_RNA_locr_f/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "12_association", "02_RNA_locr_f"))
        loopcopy("%s/04_Down_analysis/12_association/02_RNA_RELB/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "12_association", "02_RNA_RELB"))
        loopcopy("%s/04_Down_analysis/12_association/03_RNA_gener/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "12_association", "03_RNA_gener"))
        #loopcopy("%s/04_Down_analysis/12_association/05_RNA_gener_f/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "12_association", "05_RNA_gener_f"))
        loopcopy("%s/04_Down_analysis/12_association/04_RNA_REGB/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "12_association", "04_RNA_REGB"))
    else:
        loopcopy("%s/04_Down_analysis/11_association/01_RNA_locr/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "11_association", "01_RNA_locr"))
        #loopcopy("%s/04_Down_analysis/11_association/02_RNA_locr_f/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "11_association", "02_RNA_locr_f"))
        loopcopy("%s/04_Down_analysis/11_association/02_RNA_RELB/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "11_association", "02_RNA_RELB"))
        loopcopy("%s/04_Down_analysis/11_association/03_RNA_gener/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "11_association", "03_RNA_gener"))
        #loopcopy("%s/04_Down_analysis/11_association/05_RNA_gener_f/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "11_association", "05_RNA_gener_f"))
        loopcopy("%s/04_Down_analysis/11_association/04_RNA_REGB/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "11_association", "04_RNA_REGB"))

def readtbls(project_dir,diff):
    d = {}
    d["table_qc"] = readtbl("%s/01_QC/QC_stat.xls" % project_dir)
    #d["table_mapping"] = readtbl("%s/02_Hisat/mapping_rate.txt" % project_dir)
    
    d["table_REB"] = readtbl3("%s/04_Down_analysis/09_REB/REB.xls" % project_dir)

    d["table_gene_count_header"],d["table_gene_count"] = readtbl2("%s/04_Down_analysis/06_gener/01_gener_stat/gener_stat.xls" % project_dir, only5=True)
    d["table_gene_TPM_header"],d["table_gene_TPM"] = readtbl2("%s/04_Down_analysis/08_REGB/01_REGB_stat/REGB_stat.xls" % project_dir, only5=True)
    

    d["table_loci_count_header"],d["table_loci_count"] = readtbl2("%s/04_Down_analysis/03_locr/01_locr_stat/locr_stat.xls" % project_dir, only5=True)
    d["table_loci_TPM_header"],d["table_loci_TPM"] = readtbl2("%s/04_Down_analysis/05_RELB/01_RELB_stat/RELB_stat.xls" % project_dir, only5=True)
    #d["table_gene_maf_header"],d["table_gene_maf"] = readtbl2("%s/04_Down_analysis/07_gener_f/01_gener_f_stat/gener_f_stat.xls" % project_dir, only5=True)
    #d["table_loc_maf_header"],d["table_loc_maf"] = readtbl2("%s/04_Down_analysis/04_locr_f/01_locr_f_stat/locr_f_stat.xls" % project_dir, only5=True)

    d["table_loci_DESeq_header"],d["table_loci_DESeq"] = readtbl2("%s/04_Down_analysis/03_locr/02_locr_diff/%s_DESeq2.xls" % (project_dir, diff),  only5=True)
    #d["table_loci_maf_Per_header"],d["table_loci_maf_Per"] = readtbl2("%s/04_Down_analysis/04_locr_f/02_locr_f_diff/%s_Permutation.xls" % (project_dir, diff),  only5=True)
    d["table_RELB_Per_header"],d["table_RELB_Per"] = readtbl2("%s/04_Down_analysis/05_RELB/02_RELB_diff/%s_Permutation.xls" % (project_dir, diff),  only5=True)
    d["table_gene_DESeq_header"],d["table_gene_DESeq"] = readtbl2("%s/04_Down_analysis/06_gener/02_gener_diff/%s_DESeq2.xls" % (project_dir, diff),  only5=True)
    #d["table_gene_maf_Per_header"],d["table_gene_maf_Per"] = readtbl2("%s/04_Down_analysis/07_gener_f/02_gener_f_diff/%s_Permutation.xls" % (project_dir, diff),  only5=True)
    d["table_REGB_Per_header"],d["table_REGB_Per"] = readtbl2("%s/04_Down_analysis/08_REGB/02_REGB_diff/%s_Permutation.xls" % (project_dir, diff),  only5=True)
    d["table_RNA_Per_header"],d["table_RNA_Per"] = readtbl2("%s/04_Down_analysis/10_RNA/02_RNA_diff/%s_Permutation.xls" % (project_dir, diff),  only5=True)

    d["table_loci_go"] = readtbl_go(glob.glob("%s/04_Down_analysis/03_locr/05_Enrich_R/*/R_Enrichment/*_ALL.xls" % project_dir)[0], only5=True)
    d["table_loci_kegg"] = readtbl_go(glob.glob("%s/04_Down_analysis/03_locr/05_Enrich_R/*/R_Enrichment/*_KK.xls" % project_dir)[0], only5=True)
    d["table_loci_gsea"] = readtbl(glob.glob("%s/04_Down_analysis/03_locr/05_Enrich_R/*/GSEA_Enrichment/*_GO.xls" % project_dir)[0], only5=True)
    d["table_gene_go"] = readtbl_go(glob.glob("%s/04_Down_analysis/06_gener/05_Enrich_R/*/R_Enrichment/*_ALL.xls" % project_dir)[0], only5=True)
    d["table_gene_kegg"] = readtbl_go(glob.glob("%s/04_Down_analysis/06_gener/05_Enrich_R/*/R_Enrichment/*_KK.xls" % project_dir)[0], only5=True)
    d["table_gene_gsea"] = readtbl(glob.glob("%s/04_Down_analysis/06_gener/05_Enrich_R/*/GSEA_Enrichment/*_GO.xls" % project_dir)[0], only5=True)
    """
    d["table_loci_maf_go"] = readtbl_go(glob.glob("%s/04_Down_analysis/04_locr_f/05_Enrich_R/*/R_Enrichment/*_ALL.xls" % project_dir)[0], only5=True)
    d["table_loci_maf_kegg"] = readtbl_go(glob.glob("%s/04_Down_analysis/04_locr_f/05_Enrich_R/*/R_Enrichment/*_KK.xls" % project_dir)[0], only5=True)
    d["table_loci_maf_gsea"] = readtbl(glob.glob("%s/04_Down_analysis/04_locr_f/05_Enrich_R/*/GSEA_Enrichment/*_GO.xls" % project_dir)[0], only5=True)
    d["table_gene_maf_go"] = readtbl_go(glob.glob("%s/04_Down_analysis/07_gener_f/05_Enrich_R/*/R_Enrichment/*_ALL.xls" % project_dir)[0], only5=True)
    d["table_gene_maf_kegg"] = readtbl_go(glob.glob("%s/04_Down_analysis/07_gener_f/05_Enrich_R/*/R_Enrichment/*_KK.xls" % project_dir)[0], only5=True)
    d["table_gene_maf_gsea"] = readtbl(glob.glob("%s/04_Down_analysis/07_gener_f/05_Enrich_R/*/GSEA_Enrichment/*_GO.xls" % project_dir)[0], only5=True)
    """
    d["table_RELB_go"] = readtbl_go(glob.glob("%s/04_Down_analysis/05_RELB/05_Enrich_R/*/R_Enrichment/*_ALL.xls" % project_dir)[0], only5=True)
    d["table_RELB_kegg"] = readtbl_go(glob.glob("%s/04_Down_analysis/05_RELB/05_Enrich_R/*/R_Enrichment/*_KK.xls" % project_dir)[0], only5=True)
    d["table_RELB_gsea"] = readtbl(glob.glob("%s/04_Down_analysis/05_RELB/05_Enrich_R/*/GSEA_Enrichment/*_GO.xls" % project_dir)[0], only5=True)
    d["table_REGB_go"] = readtbl_go(glob.glob("%s/04_Down_analysis/08_REGB/05_Enrich_R/*/R_Enrichment/*_ALL.xls" % project_dir)[0], only5=True)
    d["table_REGB_kegg"] = readtbl_go(glob.glob("%s/04_Down_analysis/08_REGB/05_Enrich_R/*/R_Enrichment/*_KK.xls" % project_dir)[0], only5=True)
    d["table_REGB_gsea"] = readtbl(glob.glob("%s/04_Down_analysis/08_REGB/05_Enrich_R/*/GSEA_Enrichment/*_GO.xls" % project_dir)[0], only5=True)

    d["table_rna_go"] = readtbl_go(glob.glob("%s/04_Down_analysis/10_RNA/05_Enrich_R/*/R_Enrichment/*_ALL.xls" % project_dir)[0], only5=True)
    d["table_rna_kegg"] = readtbl_go(glob.glob("%s/04_Down_analysis/10_RNA/05_Enrich_R/*/R_Enrichment/*_KK.xls" % project_dir)[0], only5=True)
    d["table_rna_gsea"] = readtbl(glob.glob("%s/04_Down_analysis/10_RNA/05_Enrich_R/*/GSEA_Enrichment/*_GO.xls" % project_dir)[0], only5=True)

    d["table_RNA_TPM_header"],d["table_RNA_TPM"] = readtbl2("%s/04_Down_analysis/10_RNA/01_RNA_stat/RNA_TPM.xls" % project_dir, only5=True)
    #d["table_expression"] = readtbl("%s/04_RNA/01_RNA_TPM/Expression_stat.xls" % project_dir)

    return d

def listfigs(report_dir):
    d = {}
    imgdir = os.path.join(report_dir, "Result")
    d["figure_base_quality"] = schfile(os.path.join(imgdir, "01_QC"), "*.base_quality.png")
    d["figure_base_content"] = schfile(os.path.join(imgdir, "01_QC"), "*.base_content.png")
    #d["figure_hisat_content"] = schfile(os.path.join(imgdir, "02_Hisat"))
    d["figure_Annovar"] = schfile(os.path.join(imgdir, "04_Down_analysis", "01_REseq_detection"), "*png")
    d["figure_Edit_dif_stat1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "02_locr_diff"), "*png")
    d["figure_Edit_dif_Volcanogram1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "03_locr_vol"), "*png")
    d["figure_Edit_dif_heatmap1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "04_locr_pheatmap"), "*png")
    #d["figure_Edit_dif_stat2"] = schfile(os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "02_locr_f_diff"), "*png")
    #d["figure_Edit_dif_Volcanogram2"] = schfile(os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "03_locr_f_vol"), "*png")
    #d["figure_Edit_dif_heatmap2"] = schfile(os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "04_locr_f_pheatmap"), "*png")
    d["figure_Edit_dif_stat3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "02_RELB_diff"), "*png")
    d["figure_Edit_dif_Volcanogram3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "03_RELB_vol"), "*png")
    d["figure_Edit_dif_heatmap3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "04_RELB_pheatmap"), "*png")
    d["figure_Edit_dif_stat4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "02_gener_diff"), "*png")
    d["figure_Edit_dif_Volcanogram4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "03_gener_vol"), "*png")
    d["figure_Edit_dif_heatmap4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "04_gener_pheatmap"), "*png")
    #d["figure_Edit_dif_stat5"] = schfile(os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "02_gener_f_diff"), "*png")
    #d["figure_Edit_dif_Volcanogram5"] = schfile(os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "03_gener_f_vol"), "*png")
    #d["figure_Edit_dif_heatmap5"] = schfile(os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "04_gener_f_pheatmap"), "*png")
    d["figure_Edit_dif_stat6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "02_REGB_diff"), "*png")
    d["figure_Edit_dif_Volcanogram6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "03_REGB_vol"), "*png")
    d["figure_Edit_dif_heatmap6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "04_REGB_pheatmap"), "*png")
    d["figure_Edit_dif_stat7"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "02_RNA_diff"), "*png")
    d["figure_Edit_dif_Volcanogram7"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "03_RNA_vol"), "*png")
    d["figure_Edit_dif_heatmap7"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "04_RNA_pheatmap"), "*png")

    d["figure_veen1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "02_locr_diff", "veen"), "*png")
    #d["figure_veen2"] = schfile(os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "02_locr_f_diff", "veen"), "*png")
    d["figure_veen3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "02_gener_diff", "veen"), "*png")
    #d["figure_veen4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "02_gener_f_diff", "veen"), "*png")
    d["figure_veen5"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "02_RNA_diff", "veen"), "*png")

    d["figure_Edit_dif_GO1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "05_Enrich_R", "R_Enrichment"), "*GO1.png")
    d["figure_Edit_dif_KEGG1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "05_Enrich_R", "R_Enrichment"), "*KK1.png")
    d["figure_Edit_dif_KEGG2"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "05_Enrich_R", "R_Enrichment"), "*KK2.png")
    d["figure_Edit_dif_GSEA1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "05_Enrich_R", "GSEA_Enrichment"), "*png")
    #d["figure_Edit_dif_GO2"] = schfile(os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "05_Enrich_R", "R_Enrichment"), "*GO1.png")
    #d["figure_Edit_dif_KEGG3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "05_Enrich_R", "R_Enrichment"), "*KK1.png")
    #d["figure_Edit_dif_KEGG4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "05_Enrich_R", "R_Enrichment"), "*KK2.png")
    #d["figure_Edit_dif_GSEA2"] = schfile(os.path.join(imgdir, "04_Down_analysis", "04_locr_f", "05_Enrich_R", "GSEA_Enrichment"), "*png")
    d["figure_Edit_dif_GO3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "05_Enrich_R", "R_Enrichment"), "*GO1.png")
    d["figure_Edit_dif_KEGG5"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "05_Enrich_R", "R_Enrichment"), "*KK1.png")
    d["figure_Edit_dif_KEGG6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "05_Enrich_R", "R_Enrichment"), "*KK2.png")
    d["figure_Edit_dif_GSEA3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "05_Enrich_R", "GSEA_Enrichment"), "*png")
    d["figure_Edit_dif_GO4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "05_Enrich_R", "R_Enrichment"), "*GO1.png")
    d["figure_Edit_dif_KEGG7"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "05_Enrich_R", "R_Enrichment"), "*KK1.png")
    d["figure_Edit_dif_KEGG8"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "05_Enrich_R", "R_Enrichment"), "*KK2.png")
    d["figure_Edit_dif_GSEA4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "05_Enrich_R", "GSEA_Enrichment"), "*png")
    #d["figure_Edit_dif_GO5"] = schfile(os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "05_Enrich_R", "R_Enrichment"), "*GO1.png")
    #d["figure_Edit_dif_KEGG9"] = schfile(os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "05_Enrich_R", "R_Enrichment"), "*KK1.png")
    #d["figure_Edit_dif_KEGG10"] = schfile(os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "05_Enrich_R", "R_Enrichment"), "*KK2.png")
    #d["figure_Edit_dif_GSEA5"] = schfile(os.path.join(imgdir, "04_Down_analysis", "07_gener_f", "05_Enrich_R", "GSEA_Enrichment"), "*png")
    d["figure_Edit_dif_GO6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"), "*GO1.png")
    d["figure_Edit_dif_KEGG11"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"), "*KK1.png")
    d["figure_Edit_dif_KEGG12"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"), "*KK2.png")
    d["figure_Edit_dif_GSEA6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "GSEA_Enrichment"), "*png")
    d["figure_Edit_dif_GO7"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "05_Enrich_R", "R_Enrichment"), "*GO1.png")
    d["figure_Edit_dif_KEGG13"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "05_Enrich_R", "R_Enrichment"), "*KK1.png")
    d["figure_Edit_dif_KEGG14"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "05_Enrich_R", "R_Enrichment"), "*KK2.png")
    d["figure_Edit_dif_GSEA7"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "05_Enrich_R", "GSEA_Enrichment"), "*png")

    if ASAN:
        d["figure_association1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "12_association", "01_RNA_locr"), "*png")
        #d["figure_association2"] = schfile(os.path.join(imgdir, "04_Down_analysis", "12_association", "02_RNA_locr_f"), "*png")
        d["figure_association3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "12_association", "02_RNA_RELB"), "*png")
        d["figure_association4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "12_association", "03_RNA_gener"), "*png")
        #d["figure_association5"] = schfile(os.path.join(imgdir, "04_Down_analysis", "12_association", "05_RNA_gener_f"), "*png")
        d["figure_association6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "12_association", "04_RNA_REGB"), "*png")
    else:
        d["figure_association1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "11_association", "01_RNA_locr"), "*png")
        #d["figure_association2"] = schfile(os.path.join(imgdir, "04_Down_analysis", "11_association", "02_RNA_locr_f"), "*png")
        d["figure_association3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "11_association", "02_RNA_RELB"), "*png")
        d["figure_association4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "11_association", "03_RNA_gener"), "*png")
        #d["figure_association5"] = schfile(os.path.join(imgdir, "04_Down_analysis", "11_association", "05_RNA_gener_f"), "*png")
        d["figure_association6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "11_association", "04_RNA_REGB"), "*png")
    return d

def main():
    parser = argparse.ArgumentParser(description='''
        create HTML report for REseq pipeline.''')
    parser.add_argument('-id', '--project_dir', help="input project dir", required=True)
    parser.add_argument('-od', '--report_dir', help="output report dir", default="./")
    parser.add_argument('-cn', '--contract', help="contract number", required=True)
    #parser.add_argument('-cu', '--customer_unit', help="customer unit", required=True)
    parser.add_argument('-pn', '--name', help="project name", required=True)
    parser.add_argument('-sp', '--species', help="species name", required=True, nargs='+')
    parser.add_argument('-rf', '--reference', help="reference fasta", default='')
    parser.add_argument('-fz', '--group', help="sample group project; LZ:LZ1/LZ2;CN:CN1/CN2", required=True, nargs='+')
    parser.add_argument('-de', '--diffrent', help="differential analysis;LZ_vs_CN", required=True, nargs='+')
    parser.add_argument('-type', '--type', help="padj or pval", default="padj")
    parser.add_argument('-type_num', '--type_num', help="padj or pval threshold value", default="0.01")
    parser.add_argument('-as', '--asan', help="AS analysis or not, yes or no, default no",
                        choices=['yes', 'no'], default='no')
    #parser.add_argument('-as', '--hbase', help="report.html or report_noAS.html", default="report.html")
    #parser.add_argument('-f',  help="Foldchange", default="2")

    args = parser.parse_args()
    global ASAN
    ASAN = False if args.asan == 'no' else True
    html_base = "report_noAS_v1.2.html" if args.asan == 'no' else "report.html"

    report_dir = os.path.join(args.report_dir, "%s_Report" % args.contract)
    img_dir = os.path.join(report_dir, "Result")
    copyfigs(args.project_dir, report_dir)
    rnd = readtbls(args.project_dir, args.diffrent[0])
    figs = listfigs(report_dir)
    rnd.update(figs)
    rnd["project_name"] = args.name
    rnd["project_number"] = args.contract
    #rnd["customer_unit"] = args.customer_unit
    rnd["report_time"] = "%s年%s月%s日" % (localtime()[0], localtime()[1], localtime()[2])
    g1 = "&nbsp;&nbsp;&nbsp;&nbsp;".join(args.species)
    g2 = ""
    for i in args.group:
        new1 = re.sub(r"^", "<b>", i)
        new2 = re.sub(r"/", "&nbsp;", new1)
        new3 = re.sub(r":", ":&nbsp;</b>", new2)
        g2 = g2 + new3 + "&nbsp;&nbsp;&nbsp;&nbsp;"
    n_g2 = re.sub(r"&nbsp;&nbsp;&nbsp;&nbsp;", "&nbsp;&nbsp;", g2)
    g2 = "&nbsp;&nbsp;&nbsp;&nbsp;".join(args.group)
    g3 = "&nbsp;&nbsp;&nbsp;&nbsp;".join(args.diffrent)
    rnd["species_name"] = g1
    rnd["reference"] = args.reference
    rnd["group"] = n_g2
    rnd["diffrent"] = g3
    rnd["diff_project"] = args.diffrent[0]

    rnd["type"] = args.type
    rnd["type_num"] = args.type_num
    #rnd["f"] = args.f
#    html_base = ''
#    if args.hbase == "report.html":
#        html_base2 = "report_noAS.html"
#    else:
#        html_base2 = "report.html"

#    tpl = jinja2.Template(open(os.path.join(BASE_DIR, "template", args.hbase)).read())
    tpl = jinja2.Template(open(os.path.join(BASE_DIR,  html_base)).read())
#    with open(os.path.join(report_dir, args.hbase), 'w') as out:
    with open(os.path.join(report_dir, html_base), 'w') as out:
        out.write(tpl.render(rnd))
#    os.system("rm -rf %s_Report/%s" % (args.contract,html_base2))
    os.system("zip -r %s/%s_Report.zip %s/%s_Report" % (args.report_dir,args.contract,args.report_dir,args.contract))


if __name__ == "__main__":
    main()


