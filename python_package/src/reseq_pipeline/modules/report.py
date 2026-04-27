import os
import re
import shutil
import zipfile
import glob
import jinja2
from time import localtime
from collections import defaultdict
from importlib.resources import files
from .report_utils import safecopy,readtbl,schfile,loopcopy,readraw,loopcopy20,readtbl2,readtbl3,readtbl_go
BASE_DIR = str(files("reseq_pipeline").joinpath("scripts"))

def copyfigs(project_dir, report_dir):
    template_dir = os.path.join(BASE_DIR, "template")
    if not os.path.exists(template_dir):
        raise FileNotFoundError(f"Report template directory not found: {template_dir}")
    if os.path.exists(report_dir):
        shutil.rmtree(report_dir)
    shutil.copytree(template_dir, report_dir)

    safecopy("%s/report_noAS_v1.2.html" % BASE_DIR, report_dir)
    imgdir = os.path.join(report_dir, "Result")
    loopcopy("%s/01_QC/QC_result/*.base_quality.png" % project_dir, os.path.join(imgdir, "01_QC"))
    loopcopy("%s/01_QC/QC_result/*.base_content.png" % project_dir, os.path.join(imgdir, "01_QC"))
    
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

    loopcopy("%s/04_Down_analysis/11_association/01_RNA_locr/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "11_association", "01_RNA_locr"))
    loopcopy("%s/04_Down_analysis/11_association/02_RNA_RELB/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "11_association", "02_RNA_RELB"))
    loopcopy("%s/04_Down_analysis/11_association/03_RNA_gener/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "11_association", "03_RNA_gener"))
    loopcopy("%s/04_Down_analysis/11_association/04_RNA_REGB/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "11_association", "04_RNA_REGB"))


def readtbls(project_dir,diff):
    d = {}
    d["table_qc"] = readtbl("%s/01_QC/QC_stat.xls" % project_dir)
    
    d["table_REB"] = readtbl3("%s/04_Down_analysis/09_REB/REB.xls" % project_dir)

    d["table_gene_count_header"],d["table_gene_count"] = readtbl2("%s/04_Down_analysis/06_gener/01_gener_stat/gener_stat.xls" % project_dir, only5=True)
    d["table_gene_TPM_header"],d["table_gene_TPM"] = readtbl2("%s/04_Down_analysis/08_REGB/01_REGB_stat/REGB_stat.xls" % project_dir, only5=True)
    
    d["table_loci_count_header"],d["table_loci_count"] = readtbl2("%s/04_Down_analysis/03_locr/01_locr_stat/locr_stat.xls" % project_dir, only5=True)
    d["table_loci_TPM_header"],d["table_loci_TPM"] = readtbl2("%s/04_Down_analysis/05_RELB/01_RELB_stat/RELB_stat.xls" % project_dir, only5=True)

    d["table_loci_DESeq_header"],d["table_loci_DESeq"] = readtbl2("%s/04_Down_analysis/03_locr/02_locr_diff/%s_DESeq2.xls" % (project_dir, diff),  only5=True)
    d["table_RELB_Per_header"],d["table_RELB_Per"] = readtbl2("%s/04_Down_analysis/05_RELB/02_RELB_diff/%s_Permutation.xls" % (project_dir, diff),  only5=True)
    d["table_gene_DESeq_header"],d["table_gene_DESeq"] = readtbl2("%s/04_Down_analysis/06_gener/02_gener_diff/%s_DESeq2.xls" % (project_dir, diff),  only5=True)
    d["table_REGB_Per_header"],d["table_REGB_Per"] = readtbl2("%s/04_Down_analysis/08_REGB/02_REGB_diff/%s_Permutation.xls" % (project_dir, diff),  only5=True)
    d["table_RNA_Per_header"],d["table_RNA_Per"] = readtbl2("%s/04_Down_analysis/10_RNA/02_RNA_diff/%s_Permutation.xls" % (project_dir, diff),  only5=True)

    d["table_loci_go"] = readtbl_go(glob.glob("%s/04_Down_analysis/03_locr/05_Enrich_R/*/R_Enrichment/*_ALL.xls" % project_dir)[0], only5=True)
    d["table_loci_kegg"] = readtbl_go(glob.glob("%s/04_Down_analysis/03_locr/05_Enrich_R/*/R_Enrichment/*_KK.xls" % project_dir)[0], only5=True)
    d["table_loci_gsea"] = readtbl(glob.glob("%s/04_Down_analysis/03_locr/05_Enrich_R/*/GSEA_Enrichment/*_GO.xls" % project_dir)[0], only5=True)
    d["table_gene_go"] = readtbl_go(glob.glob("%s/04_Down_analysis/06_gener/05_Enrich_R/*/R_Enrichment/*_ALL.xls" % project_dir)[0], only5=True)
    d["table_gene_kegg"] = readtbl_go(glob.glob("%s/04_Down_analysis/06_gener/05_Enrich_R/*/R_Enrichment/*_KK.xls" % project_dir)[0], only5=True)
    d["table_gene_gsea"] = readtbl(glob.glob("%s/04_Down_analysis/06_gener/05_Enrich_R/*/GSEA_Enrichment/*_GO.xls" % project_dir)[0], only5=True)
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

    return d

def listfigs(report_dir):
    d = {}
    imgdir = os.path.join(report_dir, "Result")
    d["figure_base_quality"] = schfile(os.path.join(imgdir, "01_QC"), "*.base_quality.png")
    d["figure_base_content"] = schfile(os.path.join(imgdir, "01_QC"), "*.base_content.png")
    d["figure_Annovar"] = schfile(os.path.join(imgdir, "04_Down_analysis", "01_REseq_detection"), "*png")
    d["figure_Edit_dif_stat1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "02_locr_diff"), "*png")
    d["figure_Edit_dif_Volcanogram1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "03_locr_vol"), "*png")
    d["figure_Edit_dif_heatmap1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "04_locr_pheatmap"), "*png")
    d["figure_Edit_dif_stat3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "02_RELB_diff"), "*png")
    d["figure_Edit_dif_Volcanogram3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "03_RELB_vol"), "*png")
    d["figure_Edit_dif_heatmap3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "04_RELB_pheatmap"), "*png")
    d["figure_Edit_dif_stat4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "02_gener_diff"), "*png")
    d["figure_Edit_dif_Volcanogram4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "03_gener_vol"), "*png")
    d["figure_Edit_dif_heatmap4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "04_gener_pheatmap"), "*png")
    d["figure_Edit_dif_stat6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "02_REGB_diff"), "*png")
    d["figure_Edit_dif_Volcanogram6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "03_REGB_vol"), "*png")
    d["figure_Edit_dif_heatmap6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "04_REGB_pheatmap"), "*png")
    d["figure_Edit_dif_stat7"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "02_RNA_diff"), "*png")
    d["figure_Edit_dif_Volcanogram7"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "03_RNA_vol"), "*png")
    d["figure_Edit_dif_heatmap7"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "04_RNA_pheatmap"), "*png")

    d["figure_veen1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "02_locr_diff", "veen"), "*png")
    d["figure_veen3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "02_gener_diff", "veen"), "*png")
    d["figure_veen5"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "02_RNA_diff", "veen"), "*png")

    d["figure_Edit_dif_GO1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "05_Enrich_R", "R_Enrichment"), "*GO1.png")
    d["figure_Edit_dif_KEGG1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "05_Enrich_R", "R_Enrichment"), "*KK1.png")
    d["figure_Edit_dif_KEGG2"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "05_Enrich_R", "R_Enrichment"), "*KK2.png")
    d["figure_Edit_dif_GSEA1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "03_locr", "05_Enrich_R", "GSEA_Enrichment"), "*png")
    d["figure_Edit_dif_GO3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "05_Enrich_R", "R_Enrichment"), "*GO1.png")
    d["figure_Edit_dif_KEGG5"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "05_Enrich_R", "R_Enrichment"), "*KK1.png")
    d["figure_Edit_dif_KEGG6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "05_Enrich_R", "R_Enrichment"), "*KK2.png")
    d["figure_Edit_dif_GSEA3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "05_RELB", "05_Enrich_R", "GSEA_Enrichment"), "*png")
    d["figure_Edit_dif_GO4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "05_Enrich_R", "R_Enrichment"), "*GO1.png")
    d["figure_Edit_dif_KEGG7"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "05_Enrich_R", "R_Enrichment"), "*KK1.png")
    d["figure_Edit_dif_KEGG8"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "05_Enrich_R", "R_Enrichment"), "*KK2.png")
    d["figure_Edit_dif_GSEA4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "06_gener", "05_Enrich_R", "GSEA_Enrichment"), "*png")
    d["figure_Edit_dif_GO6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"), "*GO1.png")
    d["figure_Edit_dif_KEGG11"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"), "*KK1.png")
    d["figure_Edit_dif_KEGG12"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"), "*KK2.png")
    d["figure_Edit_dif_GSEA6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "GSEA_Enrichment"), "*png")
    d["figure_Edit_dif_GO7"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "05_Enrich_R", "R_Enrichment"), "*GO1.png")
    d["figure_Edit_dif_KEGG13"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "05_Enrich_R", "R_Enrichment"), "*KK1.png")
    d["figure_Edit_dif_KEGG14"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "05_Enrich_R", "R_Enrichment"), "*KK2.png")
    d["figure_Edit_dif_GSEA7"] = schfile(os.path.join(imgdir, "04_Down_analysis", "10_RNA", "05_Enrich_R", "GSEA_Enrichment"), "*png")

    d["figure_association1"] = schfile(os.path.join(imgdir, "04_Down_analysis", "11_association", "01_RNA_locr"), "*png")
    d["figure_association3"] = schfile(os.path.join(imgdir, "04_Down_analysis", "11_association", "02_RNA_RELB"), "*png")
    d["figure_association4"] = schfile(os.path.join(imgdir, "04_Down_analysis", "11_association", "03_RNA_gener"), "*png")
    d["figure_association6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "11_association", "04_RNA_REGB"), "*png")

    return d

def run_report(contract, name, ref, work_dir, out_dir, fz_file, group_file,
                type, type_num):

    html_base = "report_noAS_v1.2.html" 

    report_dir = os.path.join(work_dir, "../03_report", "%s_Report" % contract)
    img_dir = os.path.join(report_dir, "Result")
    copyfigs(out_dir, report_dir)
    fz_fa = []
    with open(fz_file, 'r') as fl:
        for line in fl:
            line=line.strip()
            fz_fa.append(line)
    rnd = readtbls(out_dir, fz_fa[0])
    figs = listfigs(report_dir)
    rnd.update(figs)
    rnd["project_name"] = name
    rnd["project_number"] = contract
    rnd["report_time"] = "%s年%s月%s日" % (localtime()[0], localtime()[1], localtime()[2])
    species = []
    dect1 = defaultdict(str)
    with open(group_file) as fl:
        lines = fl.readlines()[1:]
        for line in lines:
            list1 = line.strip().split("\t")
            species.append(list1[0])
            dect1[list1[1]] += list1[0] + "/"
    group = []
    for line in dect1:
        line2 = dect1[line]
        line3 = line2[:-1]
        line4 = line+":"+line3
        group.append(line4)
    g1 = "&nbsp;&nbsp;&nbsp;&nbsp;".join(species)
    g2 = ""
    for i in group:
        new1 = re.sub(r"^", "<b>", i)
        new2 = re.sub(r"/", "&nbsp;", new1)
        new3 = re.sub(r":", ":&nbsp;</b>", new2)
        g2 = g2 + new3 + "&nbsp;&nbsp;&nbsp;&nbsp;"
    n_g2 = re.sub(r"&nbsp;&nbsp;&nbsp;&nbsp;", "&nbsp;&nbsp;", g2)
    g2 = "&nbsp;&nbsp;&nbsp;&nbsp;".join(group)
    g3 = "&nbsp;&nbsp;&nbsp;&nbsp;".join(fz_fa)
    rnd["species_name"] = g1
    rnd["reference"] = ref
    rnd["group"] = n_g2
    rnd["diffrent"] = g3
    rnd["diff_project"] = fz_fa[0]

    rnd["type"] = type
    rnd["type_num"] = type_num

    tpl = jinja2.Template(open(os.path.join(BASE_DIR,  html_base)).read())
    with open(os.path.join(report_dir, html_base), 'w') as out:
        out.write(tpl.render(rnd))
    os.system("zip -r %s.zip %s" % (report_dir,report_dir))
    