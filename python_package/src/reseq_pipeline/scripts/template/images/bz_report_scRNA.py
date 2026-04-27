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
        safecopy("%s/bz_report.html" % BASE_DIR, report_dir)
    else:
        safecopy("%s/bz_report_noAS.html" % BASE_DIR, report_dir)
    imgdir = os.path.join(report_dir, "Result")
    loopcopy("%s/01_QC/QC_result/*.base_quality.png" % project_dir, os.path.join(imgdir, "01_QC"))
    loopcopy("%s/01_QC/QC_result/*.base_content.png" % project_dir, os.path.join(imgdir, "01_QC"))
    #loopcopy("%s/02_Hisat/mapping_rate.txt" % project_dir, os.path.join(imgdir, "02_Hisat"))
    #loopcopy("%s/02_Hisat/*.hisat_content.png" % project_dir, os.path.join(imgdir, "02_Hisat"))
    
    loopcopy("%s/04_Down_analysis/01_REseq_detection/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "01_REseq_detection"))
    loopcopy("%s/04_Down_analysis/08_REGB/04_REGB_pheatmap/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "04_REGB_pheatmap"))

    loopcopy("%s/04_Down_analysis/08_REGB/01_REGB_stat/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "01_REGB_stat"))
    loopcopy("%s/04_Down_analysis/08_REGB/02_REGB_diff/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "02_REGB_diff"))
    loopcopy("%s/04_Down_analysis/08_REGB/02_REGB_diff/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "02_REGB_diff"))
    loopcopy("%s/04_Down_analysis/08_REGB/03_REGB_vol/*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "03_REGB_vol"))
    loopcopy("%s/04_Down_analysis/08_REGB/05_Enrich_R/*/R_Enrichment/*GO*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/08_REGB/05_Enrich_R/*/R_Enrichment/*KK*png" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"))
    loopcopy("%s/04_Down_analysis/08_REGB/05_Enrich_R/*/R_Enrichment/*xls" % project_dir, os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"))

    loopcopy("%s/04_Down_analysis/09_REB/*" % project_dir, os.path.join(imgdir, "04_Down_analysis", "09_REB"))


def readtbls(project_dir,diff):
    d = {}
    d["table_qc"] = readtbl("%s/01_QC/QC_stat.xls" % project_dir)
    #d["table_mapping"] = readtbl("%s/02_Hisat/mapping_rate.txt" % project_dir)
    
    d["table_REB"] = readtbl3("%s/04_Down_analysis/09_REB/REB.xls" % project_dir)

    d["table_gene_TPM_header"],d["table_gene_TPM"] = readtbl2("%s/04_Down_analysis/08_REGB/01_REGB_stat/REGB_stat.xls" % project_dir, only5=True)
    
    d["table_REGB_Per_header"],d["table_REGB_Per"] = readtbl2("%s/04_Down_analysis/08_REGB/02_REGB_diff/%s_Permutation.xls" % (project_dir, diff),  only5=True)

    d["table_REGB_go"] = readtbl_go(glob.glob("%s/04_Down_analysis/08_REGB/05_Enrich_R/*/R_Enrichment/*_ALL.xls" % project_dir)[0], only5=True)
    d["table_REGB_kegg"] = readtbl_go(glob.glob("%s/04_Down_analysis/08_REGB/05_Enrich_R/*/R_Enrichment/*_KK.xls" % project_dir)[0], only5=True)

    #d["table_expression"] = readtbl("%s/04_RNA/01_RNA_TPM/Expression_stat.xls" % project_dir)

    return d

def listfigs(report_dir):
    d = {}
    imgdir = os.path.join(report_dir, "Result")
    d["figure_base_quality"] = schfile(os.path.join(imgdir, "01_QC"), "*.base_quality.png")
    d["figure_base_content"] = schfile(os.path.join(imgdir, "01_QC"), "*.base_content.png")
    #d["figure_hisat_content"] = schfile(os.path.join(imgdir, "02_Hisat"))
    d["figure_Annovar"] = schfile(os.path.join(imgdir, "04_Down_analysis", "01_REseq_detection"), "*png")
    d["figure_Edit_dif_stat6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "02_REGB_diff"), "*png")
    d["figure_Edit_dif_Volcanogram6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "03_REGB_vol"), "*png")
    d["figure_Edit_dif_heatmap6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "04_REGB_pheatmap"), "*png")

    d["figure_Edit_dif_GO6"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"), "*GO1.png")
    d["figure_Edit_dif_KEGG11"] = schfile(os.path.join(imgdir, "04_Down_analysis", "08_REGB", "05_Enrich_R", "R_Enrichment"), "*KK1.png")

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
    html_base = "bz_report_scRNA.html" if args.asan == 'no' else "bz_report.html"

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


