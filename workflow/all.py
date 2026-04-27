#!/usr/bin/evn python
# -*- coding:utf-8 -*-

import sys,os,re
import argparse
import logging

from collections import OrderedDict
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "./"))
from dagflow import DAG, Task, ParallelTask, do_dag
from REseq.common import read_config, mkdir, check_paths, check_path, read_files, read_sample_file
from REseq.config import *
from collections import defaultdict
import importlib
import glob
import shutil
import jinja2
from time import localtime
from docxtpl import DocxTemplate
from REseq.report_utils import safecopy,readtbl,schfile,loopcopy,readraw,loopcopy20,readtbl2,readtbl3,readtbl_go

importlib.reload(sys)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
LOG = logging.getLogger(__name__)

__version__ = "v1.0.0"
__author__ = ("Jiacheng Yi",)
__email__ = "..."
__all__ = []

def run_fastp(read1, read2, prefix, trim, thread, job_type,
                      work_dir, out_dir, qvalue, concurrent, refresh):
    json = mkdir("%s/json" % work_dir)
    qc_dir = mkdir("%s/QC_result" % out_dir)

    dag = DAG("reseq_fastp")
    tasks = ParallelTask(
        id="ngs_qc",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{fastp} -i {{read1}} -I {{read2}} \\
    -o {{prefix}}.clean.r1.fq.gz -O {{prefix}}.clean.r2.fq.gz \\
    -w {thread} -n 0 -f {trim} -F {trim} -t {trim} -T {trim} \\
    -q {qvalue} --json {{prefix}}_fastp.json

{fastqc} {{prefix}}.clean.r1.fq.gz {{prefix}}.clean.r2.fq.gz \\
    -t {thread} --extract -o {work_dir}

cp {{prefix}}*.json {json}/.

#{python} {scripts}/stat_fastp.py {{prefix}}_fastp.json > QC_stat.xls
#cp QC_stat.xls {out_dir}

{python} {scripts}/plot_fastqc.py \\
    -r1 {{prefix}}.clean.r1_fastqc/fastqc_data.txt \\
    -r2 {{prefix}}.clean.r2_fastqc/fastqc_data.txt \\
    --name {{prefix}}
cp {{prefix}}.base_content.* {qc_dir}
cp {{prefix}}.base_quality.* {qc_dir}
""".format(fastp=FASTP_BIN,
            python=PYTHON_BIN,
            fastqc=FASTQC_BIN,
            scripts=SCRIPTS,
            thread=thread,
            trim=trim,
            qvalue=qvalue,
            json=json,
            work_dir=work_dir,
            out_dir=out_dir,
            qc_dir=qc_dir
        ),
        read1=read1,
        read2=read2,
        prefix=prefix
    )
    clean1 = []
    clean2 = []
    for i in prefix:
        clean1.append(os.path.join(work_dir, "%s.clean.r1.fq.gz" % i))
        clean2.append(os.path.join(work_dir, "%s.clean.r2.fq.gz" % i))
    dag.add_task(*tasks)
    do_dag(dag, concurrent, refresh)
    return clean1, clean2

def stat_fastp(input, work_dir, out_dir, job_type, concurrent, refresh):

    json = mkdir("%s/json" % work_dir)

    dag = DAG("reseq_stat_fastp")
    task = Task(
        id="stat_qc",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{python} {scripts}/stat_fastp.py {input} {json} > QC_stat.xls
cp QC_stat.xls {out_dir}

""".format(python=PYTHON_BIN,
            scripts=SCRIPTS,
            input=input,
            out_dir=out_dir,
            json=json
        )
    )
    dag.add_task(task)
    do_dag(dag, concurrent, refresh)
    return 0

def create_hisat2tasks(prefix, read1, read2, ref, thread,
                        job_type, work_dir, out_dir, concurrent, refresh):

    dag = DAG("reseq_hisat")
    task = ParallelTask(
        id="hisat2",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{hisat} --threads {thread} -x {ref} \\
  -1 {{read1}} -2 {{read2}} | \\
  {samtools} sort -O bam -@ {thread} -o {{prefix}}.sort.bam

""".format(hisat=HISAT2_BIN,
            samtools=SAMTOOLS_BIN,
            ref=ref,
            thread=thread
        ),
        read1=read1,
        read2=read2,
        prefix=prefix
    )
    dag.add_task(*task)
    do_dag(dag, concurrent, refresh)
    hisat_bam = []
    for i in prefix:
        hisat_bam.append(os.path.join(work_dir, "%s.sort.bam" % i))
    return hisat_bam

def create_gene_count(gtf, hisat_bam, prefix, thread, job_type,
                      work_dir, concurrent, refresh):

    dag = DAG("reseq_rna_count")
    tasks = ParallelTask(
        id="work_gene_count_htseq",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{htseq} -f bam -s no -m intersection-nonempty \\
    {{hisat_bam}} {gtf} > {{prefix}}.gene.xls

""".format(htseq=HTSEQ,
            gtf=gtf
        ),
        hisat_bam=hisat_bam,
        prefix=prefix
    )
    dag.add_task(*tasks)
    do_dag(dag, concurrent, refresh)
    gene_list = []
    for i in prefix:
        gene_list.append(os.path.join(work_dir, "%s.gene.xls" % i))
    return gene_list

def merge_gene_count(input, work_dir, out_dir, job_type, concurrent, refresh):

    dag = DAG("reseq_merge_rna_count")
    task = Task(
        id="work_merge_gene_count",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{perl} {scripts}/merge_gene.pl {input} > RNA_Count.xls
cp RNA_Count.xls {out_dir}

""".format(perl=PERL_BIN,
            scripts=SCRIPTS,
            input=input,
            out_dir=out_dir
        )
    )
    dag.add_task(task)
    do_dag(dag, concurrent, refresh)
    gene_count_xls = os.path.join(work_dir, "RNA_Count.xls")
    return gene_count_xls

def work_bam_gatk(hisat_bam, prefix, ref, input, thread, job_type, 
                annovar_ref, work_dir, out_dir, concurrent, refresh):

    tmp = mkdir(os.path.join(work_dir, "tmp"))
    tasks = ParallelTask(
        id="work_bam",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
export PATH={java}:$PATH

{GATK} MarkDuplicates -I {{hisat_bam}} \\
    -O {{prefix}}.deduped.bam -M {{prefix}}.marked_dup_metrics.txt \\
    --CREATE_INDEX true --REMOVE_DUPLICATES true \\
    --TAG_DUPLICATE_SET_MEMBERS true --TMP_DIR {tmp}

{GATK} AddOrReplaceReadGroups -I {{prefix}}.deduped.bam \\
    -O {{prefix}}.picard.bam -LB {{prefix}} -PL illumina \\
    -PU {{prefix}} -SM {{prefix}}

{samtools} index {{prefix}}.picard.bam

{GATK} SplitNCigarReads -R {ref} -I {{prefix}}.picard.bam \\
    -O {{prefix}}.dedup_split.bam

{GATK} HaplotypeCaller -R {ref} -ERC GVCF -I {{prefix}}.dedup_split.bam \\
    -O {{prefix}}.g.vcf

{GATK} GenotypeGVCFs -R {ref} -V {{prefix}}.g.vcf -O {{prefix}}.vcf

{GATK} IndexFeatureFile -F {{prefix}}.vcf

{GATK} SelectVariants -V {{prefix}}.vcf -select-type SNP -O {{prefix}}.snp.vcf
{GATK} SelectVariants -V {{prefix}}.vcf -select-type INDEL -O {{prefix}}.indel.vcf

{GATK} IndexFeatureFile -F {{prefix}}.snp.vcf
{GATK} IndexFeatureFile -F {{prefix}}.indel.vcf

{GATK} BaseRecalibrator -R {ref} -I {{prefix}}.dedup_split.bam \\
    --known-sites {{prefix}}.snp.vcf \\
    --known-sites {{prefix}}.indel.vcf \\
    -O {{prefix}}.recal_data.table

{GATK} ApplyBQSR -R {ref} -I {{prefix}}.dedup_split.bam \\
    --bqsr-recal-file {{prefix}}.recal_data.table \\
    -O {{prefix}}.bqsr.bam

{samtools} index {{prefix}}.bqsr.bam

{GATK} SplitNCigarReads -R {ref} -I {{prefix}}.bqsr.bam \\
    -O {{prefix}}.bqsr.bed.bam 

#GIREMI
##giremi prefix.bqsr.bed.bam -f ref -o prefix.res
####

{GATK} HaplotypeCaller -R {ref} -ERC GVCF -I {{prefix}}.bqsr.bed.bam \\
    -O {{prefix}}.bed.g.vcf

{GATK} GenotypeGVCFs -R {ref} -V {{prefix}}.bed.g.vcf \\
    -O {{prefix}}.bed.vcf

rm -rf {{prefix}}.dedup_split.bam {{prefix}}.picard.bam {{prefix}}.deduped.bam
rm -rf {{prefix}}.snp* {{prefix}}.indel*

{bgzip} {{prefix}}.bed.vcf
{tabix} {{prefix}}.bed.vcf.gz

{GATK} SelectVariants -V {{prefix}}.bed.vcf.gz -select-type SNP \\
    -O {{prefix}}.snp.vcf.gz
{GATK} SelectVariants -V {{prefix}}.bed.vcf.gz -select-type INDEL \\
    -O {{prefix}}.indel.vcf.gz

{GATK} VariantFiltration -V {{prefix}}.snp.vcf.gz \\
    -filter "QD < 2.0" --filter-name "QD2" \\
    -filter "QUAL < 30.0" --filter-name "QUAL30" \\
    -filter "FS > 60.0" --filter-name "FS60" \\
    -filter "MQ < 40.0" --filter-name "MQ40" \\
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
    -O {{prefix}}.vcf

{perl} {scripts}/filter_vcf.pl {{prefix}}.vcf > {{prefix}}.hardfiltered.biallelic.vcf

{perl} {annovar}/convert2annovar.pl --format vcf4 \\
    --includeinfo {{prefix}}.hardfiltered.biallelic.vcf \\
    --allsample --withfreq > {{prefix}}.avinput

{perl} {annovar}/table_annovar.pl {{prefix}}.avinput \\
    {annovar}/{annovar_ref}/ -buildver {annovar_ref} -out {{prefix}} \\
    -remove -protocol refGene -operation g -nastring . \\
    -csvout --thread {thread}

{perl} {scripts}/Reads.pl {{prefix}}
{python} {scripts}/Reads.py {{prefix}}.Reads.csv \\
    {{prefix}}.{annovar_ref}_multianno.csv {{prefix}}

""".format(GATK=GATK_BIN,
            samtools=SAMTOOLS_BIN,
            perl=PERL_BIN,
            python=PYTHON_BIN,
            java=JAVA_BIN,
            annovar=ANNOVAR_BIN,
            scripts=SCRIPTS,
            bgzip=BGZIP,
            tabix=TABIX,
            #giremi=GIREMI,
            thread=thread,
            ref=ref,
            annovar_ref=annovar_ref,
            tmp=tmp,
            work_dir=work_dir
        ),
        prefix=prefix,
        hisat_bam=hisat_bam
    )
    GATK_result = []
    bam_list = []
    work_dir = os.path.abspath(work_dir)
    for i in prefix:
        GATK_result.append(os.path.join(work_dir, "%s.starded.csv" % i))
        bam_list.append(os.path.join(work_dir, "%s.bqsr.bed.bam" % i))

    task = Task(
        id="work_merge_GATK",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{python} {scripts}/merge_sample.py {input} 
cp GATK_information.csv {out_dir}
""".format(python=PYTHON_BIN,
            scripts=SCRIPTS,
            out_dir=out_dir,
            input=input
        )
    )
    gatk_file = os.path.join(work_dir, "GATK_information.csv")

    dag = DAG("REseq_gatk")

    dag.add_task(*tasks)
    dag.add_task(task)
    task.set_upstream(*tasks)
    do_dag(dag, concurrent, refresh)

    return gatk_file, bam_list


def run_reseq_identify_REDItools(gatk_bam_list, prefix, ref, thread, work_dir,
                                job_type, concurrent, refresh):

    tasks = ParallelTask(
        id="work_REDItools",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{python} {reditools} -f {{bam}} -r {ref} \\
    -o {{prefix}}_RED.txt
""".format(python=PYTHON_BIN,
            reditools=REDITOOL_BIN,
            ref=ref
        ),
        bam=gatk_bam_list,
        prefix=prefix
    )

    REDItool_result = []
    for i in prefix:
        REDItool_result.append(os.path.join(work_dir, "%s_RED.txt" % i))

    dag = DAG("REseq_REDItools")
    dag.add_task(*tasks)
    do_dag(dag, concurrent, refresh)
    return REDItools_result

def run_reseq_identify_VarScan(gatk_bam_list, bed, ref, prefix, job_type, work_dir,
                                concurrent, refresh):
    tasks = ParallelTask(
        id="work_VarScan",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{samtools} mpileup -q 1 -f {ref} -l {bed} \\
    {{bam}} --output {{prefix}}.normal.pileup

{java} -Xmx10g -jar {scripts}/VarScan.v2.3.9.jar somatic \\
    {{prefix}}.normal.pileup \\
    {{prefix}}.normal.pileup \\
    {{prefix}}.tumor
""".format(samtools=SAMTOOLS_BIN,
            scripts=SCRIPTS,
            java=JAVA_BIN,
            bed=bed,
            ref=ref
        ),
        bam=gatk_bam_list,
        prefix=prefix
    )

    VarScan_result = []
    for i in prefix:
        VarScan_result.append(os.path.join(work_dir, "%s.tumor" % i))

    dag = DAG("REseq_VarScan")
    dag.add_task(*tasks)
    do_dag(dag, concurrent, refresh)
    return VarScan_result

def run_reseq_identify_Sprint(gatk_bam_list, prefix, ref, gtf, ref_repat, thread, job_type,
                concurrent, refresh, work_dir):
    task = Task(
        id="work_sprint_index",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{sprint} prepare -t {gtf} {ref} {bwa} 

""".format(sprint=SPRINT_BIN,
            bwa=BWA_BIN,
            gtf=gtf,
            ref=ref
        )
    )

    tasks = ParallelTask(
        id="work_sprint",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{sprint_from_bam} -rp {ref_repat} -b {{bam}} \\
    -f {ref} -o {{prefix}} {samtools} 
""".format(samtools=SAMTOOLS_BIN,
            sprint_from_bam=SPRINT_BAM,
            ref_repat=ref_repat,
            ref=ref
        ),
        bam=gatk_bam_list,
        prefix=prefix
    )

    sprint_result = []
    for i in prefix:
        sprint_result.append(os.path.join(work_dir, "%s.txt" % i))

    dag = DAG("REseq_Sprint")
    dag.add_task(task)
    dag.add_task(*tasks)
    task.set_upstream(*tasks)

    do_dag(dag, concurrent, refresh)
    return sprint_result

def run_reseq_identify_RED_ML(gatk_bam_list, prefix, ref, dbsnp, ref_repat,
                ref_alu, ref_snp,
                thread, job_type, concurrent, refresh, work_dir):
    tasks = ParallelTask(
        id="work_RED_ML",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{perl} {red_ml} --rnabam {{bam}} \\
    --reference {ref} \\
    --dbsnp {dbsnp} \\
    --simpleRepeat {ref_repat} \\ 
    --alu {ref_alu} \\
    --snplist {ref_snp} \\
    --outdir {work_dir}
""".format(perl=PERL_BIN,
            red_ml=RED_ML,
            ref_repat=ref_repat,
            dbsnp=dbsnp,
            ref_alu=ref_alu,
            ref_snp=ref_snp,
            ref=ref,
            work_dir=work_dir
        ),
        bam=gatk_bam_list
    )
    dag = DAG("REseq_RED_ML")
    dag.add_task(*tasks)
    do_dag(dag, concurrent, refresh)
    return 0

def reseq_work(input, method, reference, prefix, hisat_bam,
                work_dir, out_dir,
                thread, job_type, concurrent, refresh):

    GATK_result, bam_list = work_bam_gatk(
        input=input,
        prefix=prefix,
        hisat_bam=hisat_bam,
        ref=ref_fasta[reference],
        annovar_ref=reference,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=work_dir,
        out_dir=out_dir
    )

    work_list = method.split(";")
    if "REDItools" in work_list:
        REDItools_result = run_reseq_identify_REDItools(
            gatk_bam_list=bam_list,
            prefix=prefix,
            ref=ref_path[reference],
            thread=thread,
            job_type=job_type,
            concurrent=concurrent,
            refresh=refresh,
            work_dir=work_dir
        )

    if "VarScan" in work_list:
        VarScan_result = run_reseq_identify_VarScan(
            gatk_bam_list=bam_list,
            prefix=prefix,
            ref=ref_path[reference],
            bed=ref_bed[reference],
            thread=thread,
            job_type=job_type,
            concurrent=concurrent,
            refresh=refresh,
            work_dir=work_dir
        ) 

    if "Sprint" in work_list:
        sprint_result = run_reseq_identify_Sprint(
            gatk_bam_list=bam_list,
            prefix=prefix,
            ref=ref_path[reference],
            gtf=ref_gtf[reference],
            ref_repat=ref_repeat[reference],
            thread=thread,
            job_type=job_type,
            concurrent=concurrent,
            refresh=refresh,
            work_dir=work_dir
        ) 

    if "RED_ML" in work_list:
        sprint_result = run_reseq_identify_RED_ML(
            gatk_bam_list=bam_list,
            prefix=prefix,
            ref=ref_path[reference],
            dbsnp=ref_dbsnp[reference],
            ref_repat=ref_repeat[reference],
            ref_alu=ref_alu[reference],
            ref_snp=ref_snp_d[reference],
            thread=thread,
            job_type=job_type,
            concurrent=concurrent,
            refresh=refresh,
            work_dir=work_dir
        )

    return GATK_result

def down_analysis_stat(input, gatk_file, tf_file, gene_length, rna_count,
                      job_type, work_dir, out_dir):

    task = Task(
        id="work_down_analysis",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{perl} {scripts}/get_group1.pl {input} > {work_dir}/group.xls

{Rscript} {scripts}/get_fz.R {work_dir}/group.xls {work_dir}/fz.txt
{Rscript} {scripts}/01_Edit_readGATK.R {out_dir}/01_REseq_detection {gatk_file}
#{Rscript} {scripts}/02_Eidt_GATK_cor.R {out_dir}/02_REseq_PCA \\
#    {out_dir}/01_REseq_detection/GATK_information.csv {work_dir}/group.xls

{Rscript} {scripts}/03_locr_stat.R {out_dir}/01_REseq_detection/GATK_information.csv \\
    {out_dir}/03_locr/01_locr_stat
{perl} {scripts}/get_group.pl {work_dir}/group.xls {work_dir}/fz.txt \\
    {out_dir}/03_locr/01_locr_stat/locr_stat.xls {out_dir}/03_locr/01_locr_stat

#{Rscript} {scripts}/04_locr_f_stat.R {tf_file} \\
#    {out_dir}/01_REseq_detection/GATK_information.csv {out_dir}/04_locr_f/01_locr_f_stat 
#{perl} {scripts}/get_group.pl {work_dir}/group.xls {work_dir}/fz.txt \\
#    {out_dir}/04_locr_f/01_locr_f_stat/locr_f_stat.xls {out_dir}/04_locr_f/01_locr_f_stat

{Rscript} {scripts}/05_RELB_stat.R {gene_length} \\
    {out_dir}/01_REseq_detection/GATK_information.csv {out_dir}/05_RELB/01_RELB_stat {rna_count}
{perl} {scripts}/get_group.pl {work_dir}/group.xls {work_dir}/fz.txt \\
    {out_dir}/05_RELB/01_RELB_stat/RELB_stat.xls {out_dir}/05_RELB/01_RELB_stat

{Rscript} {scripts}/06_gener_stat.R {out_dir}/01_REseq_detection/GATK_information.csv \\
    {out_dir}/06_gener/01_gener_stat
{perl} {scripts}/get_group.pl {work_dir}/group.xls {work_dir}/fz.txt \\
    {out_dir}/06_gener/01_gener_stat/gener_stat.xls {out_dir}/06_gener/01_gener_stat

#{Rscript} {scripts}/07_gener_f_stat.R {tf_file} \\
#    {out_dir}/01_REseq_detection/GATK_information.csv {out_dir}/07_gener_f/01_gener_f_stat 
#{perl} {scripts}/get_group.pl {work_dir}/group.xls {work_dir}/fz.txt \\
#    {out_dir}/07_gener_f/01_gener_f_stat/gener_f_stat.xls {out_dir}/07_gener_f/01_gener_f_stat

{Rscript} {scripts}/08_REGB_stat.R {gene_length} \\
    {out_dir}/01_REseq_detection/GATK_information.csv {out_dir}/08_REGB/01_REGB_stat {rna_count}
{perl} {scripts}/get_group.pl {work_dir}/group.xls {work_dir}/fz.txt \\
    {out_dir}/08_REGB/01_REGB_stat/REGB_stat.xls {out_dir}/08_REGB/01_REGB_stat

""".format(perl=PERL_BIN,
            scripts=SCRIPTS,
            Rscript=RSCRIPT,
            out_dir=out_dir,
            work_dir=work_dir,
            gatk_file=gatk_file,
            tf_file=tf_file,
            input=input,
            gene_length=gene_length,
            rna_count=rna_count
        )
    )

    fz_file = os.path.join(work_dir, "fz.txt")
    group_file = os.path.join(work_dir, "group.xls")
    return task, fz_file, group_file

def pcr_analysis(work_dir, out_dir, job_type):
    task = Task(
        id="work_PCR_analysis",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{Rscript} {scripts}/02_Eidt_GATK_cor.R {out_dir}/02_REseq_PCA \\
    {out_dir}/01_REseq_detection/GATK_information.csv {work_dir}/group.xls
""".format(
            scripts=SCRIPTS,
            Rscript=RSCRIPT,
            out_dir=out_dir,
            work_dir=work_dir
        )
    )

    return task
"""
def down_analysis_deseq(tf_file, rna_count, gene_length,
                      job_type, work_dir, out_dir):

    task = Task(
        id="work_down_analysis_deseq",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script=\"""
{Rscript} {scripts}/Limma.R {out_dir}/03_locr/02_locr_diff \\
    {out_dir}/03_locr/01_locr_stat/locr_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/DEseq2.R {out_dir}/03_locr/02_locr_diff \\
    {out_dir}/03_locr/01_locr_stat/locr_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/Limma.R {out_dir}/04_locr_f/02_locr_f_diff \\
    {out_dir}/04_locr_f/01_locr_f_stat/locr_f_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/Permutation.R {out_dir}/04_locr_f/02_locr_f_diff \\
    {out_dir}/04_locr_f/01_locr_f_stat/locr_f_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/Permutation.R {out_dir}/05_RELB/02_RELB_diff \\
    {out_dir}/05_RELB/01_RELB_stat/RELB_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/Limma.R {out_dir}/06_gener/02_gener_diff \\
    {out_dir}/06_gener/01_gener_stat/gener_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/DEseq2.R {out_dir}/06_gener/02_gener_diff \\
    {out_dir}/06_gener/01_gener_stat/gener_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/Limma.R {out_dir}/07_gener_f/02_gener_f_diff \\
    {out_dir}/07_gener_f/01_gener_f_stat/gener_f_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/Permutation.R {out_dir}/07_gener_f/02_gener_f_diff \\
    {out_dir}/07_gener_f/01_gener_f_stat/gener_f_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/Permutation.R {out_dir}/08_REGB/02_REGB_diff \\
    {out_dir}/08_REGB/01_REGB_stat/REGB_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/09_RNA_stat.R {gene_length} {rna_count} {out_dir}/10_RNA/01_RNA_stat

{perl} {scripts}/get_group.pl {work_dir}/group.xls {work_dir}/fz.txt \\
    {out_dir}/10_RNA/01_RNA_stat/RNA_TPM.xls {out_dir}/10_RNA/01_RNA_stat
{perl} {scripts}/get_group.pl {work_dir}/group.xls {work_dir}/fz.txt \\
    {out_dir}/10_RNA/01_RNA_stat/RNA_Count/RNA_Count.xls {out_dir}/10_RNA/01_RNA_stat/RNA_Count

{Rscript} {scripts}/10_RNA_stat_cor.R {out_dir}/10_RNA/01_RNA_stat/RNA_PCA/ \\
    {out_dir}/10_RNA/01_RNA_stat/RNA_TPM.csv {work_dir}/group.xls

{Rscript} {scripts}/Limma.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/01_RNA_stat/RNA_TPM.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/Permutation.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/01_RNA_stat/RNA_TPM.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/DEseq2.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/01_RNA_stat/RNA_Count/RNA_Count.csv {tf_file} {work_dir}/group.xls

\""".format(scripts=SCRIPTS,
            Rscript=RSCRIPT,
            perl=PERL_BIN,
            out_dir=out_dir,
            work_dir=work_dir,
            gene_length=gene_length,
            rna_count=rna_count,
            tf_file=tf_file
        )
    )

    return task
"""
def down_analysis_locr(tf_file, ref,
                      job_type, work_dir, out_dir):

    task = Task(
        id="work_down_analysis_locr",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{Rscript} {scripts}/Limma.R {out_dir}/03_locr/02_locr_diff \\
    {out_dir}/03_locr/01_locr_stat/locr_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/DEseq2.R {out_dir}/03_locr/02_locr_diff \\
    {out_dir}/03_locr/01_locr_stat/locr_stat.csv {tf_file} {work_dir}/group.xls

ls {out_dir}/03_locr/02_locr_diff/*csv |sed 's/.csv//' | while read i;do \\
    {Rscript} {scripts}/stat_diff.R ${{i}}.csv $i;done

mkdir -p {out_dir}/03_locr/02_locr_diff/veen

cat {work_dir}/fz.txt | sed 's/+/_/'| while read i;do \\
    {Rscript} {scripts}/diff_ggvenn.R {out_dir}/03_locr/02_locr_diff/${{i}}_DESeq2.xls \\
    {out_dir}/03_locr/02_locr_diff/${{i}}_Limma.xls DESeq2 Limma {out_dir}/03_locr/02_locr_diff/veen/${{i}}_veen.png;done

{Rscript} {scripts}/Vol.R {out_dir}/03_locr/02_locr_diff \\
    {out_dir}/03_locr/03_locr_vol Limma {work_dir}/group.xls
{Rscript} {scripts}/Vol.R {out_dir}/03_locr/02_locr_diff \\
    {out_dir}/03_locr/03_locr_vol DESeq2 {work_dir}/group.xls

{Rscript} {scripts}/pheatmap.R {out_dir}/03_locr/02_locr_diff/ \\
    {out_dir}/03_locr/01_locr_stat/ {out_dir}/03_locr/04_locr_pheatmap Limma {work_dir}/group.xls
{Rscript} {scripts}/pheatmap.R {out_dir}/03_locr/02_locr_diff/ \\
    {out_dir}/03_locr/01_locr_stat/ {out_dir}/03_locr/04_locr_pheatmap DESeq2 {work_dir}/group.xls

{Rscript} {scripts}/Enrich_R.R {out_dir}/03_locr/02_locr_diff/ \\
    {out_dir}/03_locr/05_Enrich_R Limma {work_dir}/group.xls {ref} y Up
{Rscript} {scripts}/Enrich_R.R {out_dir}/03_locr/02_locr_diff/ \\
    {out_dir}/03_locr/05_Enrich_R Limma {work_dir}/group.xls {ref} y Down
{Rscript} {scripts}/Enrich_R.R {out_dir}/03_locr/02_locr_diff/ \\
    {out_dir}/03_locr/05_Enrich_R DESeq2 {work_dir}/group.xls {ref} y Up
{Rscript} {scripts}/Enrich_R.R {out_dir}/03_locr/02_locr_diff/ \\
    {out_dir}/03_locr/05_Enrich_R DESeq2 {work_dir}/group.xls {ref} y Down

{Rscript} {scripts}/plot_GO_KEGG_all.R {work_dir}/group.xls Limma {out_dir}/03_locr/05_Enrich_R
{Rscript} {scripts}/plot_GO_KEGG_all.R {work_dir}/group.xls DESeq2 {out_dir}/03_locr/05_Enrich_R

{Rscript} {scripts}/Enrich_GASE.R {out_dir}/03_locr/02_locr_diff/ \\
    {out_dir}/03_locr/05_Enrich_R Limma {work_dir}/group.xls {ref} y
{Rscript} {scripts}/Enrich_GASE.R {out_dir}/03_locr/02_locr_diff/ \\
    {out_dir}/03_locr/05_Enrich_R DESeq2 {work_dir}/group.xls {ref} y
""".format(scripts=SCRIPTS,
            Rscript=RSCRIPT,
            out_dir=out_dir,
            work_dir=work_dir,
            ref=ref,
            tf_file=tf_file
        )
    )

    return task

def down_analysis_locr_f(tf_file, ref,
                        job_type, work_dir, out_dir):

    task = Task(
        id="work_down_analysis_locr_f",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{Rscript} {scripts}/Limma.R {out_dir}/04_locr_f/02_locr_f_diff \\
    {out_dir}/04_locr_f/01_locr_f_stat/locr_f_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/Permutation.R {out_dir}/04_locr_f/02_locr_f_diff \\
    {out_dir}/04_locr_f/01_locr_f_stat/locr_f_stat.csv {tf_file} {work_dir}/group.xls

ls {out_dir}/04_locr_f/02_locr_f_diff/*csv |sed 's/.csv//' | while read i;do \\
    {Rscript} {scripts}/stat_diff.R ${{i}}.csv $i;done

mkdir -p {out_dir}/04_locr_f/02_locr_f_diff/veen

cat {work_dir}/fz.txt | sed 's/+/_/'| while read i;do \\
    {Rscript} {scripts}/diff_ggvenn.R {out_dir}/04_locr_f/02_locr_f_diff/${{i}}_Permutation.xls \\
    {out_dir}/04_locr_f/02_locr_f_diff/${{i}}_Limma.xls Permutation Limma {out_dir}/04_locr_f/02_locr_f_diff/veen/${{i}}_veen.png;done

{Rscript} {scripts}/Vol.R {out_dir}/04_locr_f/02_locr_f_diff/ \\
    {out_dir}/04_locr_f/03_locr_f_vol Limma {work_dir}/group.xls
{Rscript} {scripts}/Vol.R {out_dir}/04_locr_f/02_locr_f_diff/ \\
    {out_dir}/04_locr_f/03_locr_f_vol Permutation {work_dir}/group.xls

{Rscript} {scripts}/pheatmap.R {out_dir}/04_locr_f/02_locr_f_diff/ \\
    {out_dir}/04_locr_f/01_locr_f_stat/ {out_dir}/04_locr_f/04_locr_f_pheatmap Limma {work_dir}/group.xls
{Rscript} {scripts}/pheatmap.R {out_dir}/04_locr_f/02_locr_f_diff/ \\
    {out_dir}/04_locr_f/01_locr_f_stat/ {out_dir}/04_locr_f/04_locr_f_pheatmap Permutation {work_dir}/group.xls

{Rscript} {scripts}/Enrich_R.R {out_dir}/04_locr_f/02_locr_f_diff/ \\
    {out_dir}/04_locr_f/05_Enrich_R Limma {work_dir}/group.xls {ref} y Up
{Rscript} {scripts}/Enrich_R.R {out_dir}/04_locr_f/02_locr_f_diff/ \\
    {out_dir}/04_locr_f/05_Enrich_R Limma {work_dir}/group.xls {ref} y Down
{Rscript} {scripts}/Enrich_R.R {out_dir}/04_locr_f/02_locr_f_diff/ \\
    {out_dir}/04_locr_f/05_Enrich_R Permutation {work_dir}/group.xls {ref} y Up
{Rscript} {scripts}/Enrich_R.R {out_dir}/04_locr_f/02_locr_f_diff/ \\
    {out_dir}/04_locr_f/05_Enrich_R Permutation {work_dir}/group.xls {ref} y Down

{Rscript} {scripts}/plot_GO_KEGG_all.R {work_dir}/group.xls Limma {out_dir}/04_locr_f/05_Enrich_R
{Rscript} {scripts}/plot_GO_KEGG_all.R {work_dir}/group.xls Permutation {out_dir}/04_locr_f/05_Enrich_R

{Rscript} {scripts}/Enrich_GASE.R {out_dir}/04_locr_f/02_locr_f_diff/ \\
    {out_dir}/04_locr_f/05_Enrich_R Limma {work_dir}/group.xls {ref} y
{Rscript} {scripts}/Enrich_GASE.R {out_dir}/04_locr_f/02_locr_f_diff/ \\
    {out_dir}/04_locr_f/05_Enrich_R Permutation {work_dir}/group.xls {ref} y
""".format(scripts=SCRIPTS,
            Rscript=RSCRIPT,
            out_dir=out_dir,
            work_dir=work_dir,
            ref=ref,
            tf_file=tf_file
        )
    )

    return task

def down_analysis_RELB(tf_file, ref,
                      job_type, work_dir, out_dir):

    task = Task(
        id="work_down_analysis_RELB",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{Rscript} {scripts}/Permutation.R {out_dir}/05_RELB/02_RELB_diff \\
    {out_dir}/05_RELB/01_RELB_stat/RELB_stat.csv {tf_file} {work_dir}/group.xls

ls {out_dir}/05_RELB/02_RELB_diff/*csv |sed 's/.csv//' | while read i;do \\
    {Rscript} {scripts}/stat_diff.R ${{i}}.csv $i;done

{Rscript} {scripts}/Vol.R {out_dir}/05_RELB/02_RELB_diff \\
    {out_dir}/05_RELB/03_RELB_vol Permutation {work_dir}/group.xls

{Rscript} {scripts}/pheatmap.R {out_dir}/05_RELB/02_RELB_diff \\
    {out_dir}/05_RELB/01_RELB_stat/ {out_dir}/05_RELB/04_RELB_pheatmap Permutation {work_dir}/group.xls

{Rscript} {scripts}/Enrich_R.R {out_dir}/05_RELB/02_RELB_diff \\
    {out_dir}/05_RELB/05_Enrich_R Permutation {work_dir}/group.xls {ref} y Up
{Rscript} {scripts}/Enrich_R.R {out_dir}/05_RELB/02_RELB_diff \\
    {out_dir}/05_RELB/05_Enrich_R Permutation {work_dir}/group.xls {ref} y Down

{Rscript} {scripts}/plot_GO_KEGG_all.R {work_dir}/group.xls Permutation {out_dir}/05_RELB/05_Enrich_R

{Rscript} {scripts}/Enrich_GASE.R {out_dir}/05_RELB/02_RELB_diff \\
    {out_dir}/05_RELB/05_Enrich_R Permutation {work_dir}/group.xls {ref} y
""".format(scripts=SCRIPTS,
            Rscript=RSCRIPT,
            out_dir=out_dir,
            work_dir=work_dir,
            ref=ref,
            tf_file=tf_file
        )
    )

    return task

def down_analysis_gener(tf_file, ref,
                      job_type, work_dir, out_dir):

    task = Task(
        id="work_down_analysis_gener",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{Rscript} {scripts}/Limma.R {out_dir}/06_gener/02_gener_diff \\
    {out_dir}/06_gener/01_gener_stat/gener_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/DEseq2.R {out_dir}/06_gener/02_gener_diff \\
    {out_dir}/06_gener/01_gener_stat/gener_stat.csv {tf_file} {work_dir}/group.xls

ls {out_dir}/06_gener/02_gener_diff/*csv |sed 's/.csv//' | while read i;do \\
    {Rscript} {scripts}/stat_diff.R ${{i}}.csv $i;done

mkdir -p {out_dir}/06_gener/02_gener_diff/veen

cat {work_dir}/fz.txt | sed 's/+/_/'| while read i;do \\
    {Rscript} {scripts}/diff_ggvenn.R {out_dir}/06_gener/02_gener_diff/${{i}}_DESeq2.xls \\
    {out_dir}/06_gener/02_gener_diff/${{i}}_Limma.xls DESeq2 Limma {out_dir}/06_gener/02_gener_diff/veen/${{i}}_veen.png;done

{Rscript} {scripts}/Vol.R {out_dir}/06_gener/02_gener_diff/ \\
    {out_dir}/06_gener/03_gener_vol Limma {work_dir}/group.xls
{Rscript} {scripts}/Vol.R {out_dir}/06_gener/02_gener_diff/ \\
    {out_dir}/06_gener/03_gener_vol DESeq2 {work_dir}/group.xls

{Rscript} {scripts}/pheatmap.R {out_dir}/06_gener/02_gener_diff/ \\
    {out_dir}/06_gener/01_gener_stat/ {out_dir}/06_gener/04_gener_pheatmap Limma {work_dir}/group.xls
{Rscript} {scripts}/pheatmap.R {out_dir}/06_gener/02_gener_diff/ \\
    {out_dir}/06_gener/01_gener_stat/ {out_dir}/06_gener/04_gener_pheatmap DESeq2 {work_dir}/group.xls

{Rscript} {scripts}/Enrich_R.R {out_dir}/06_gener/02_gener_diff/ \\
    {out_dir}/06_gener/05_Enrich_R Limma {work_dir}/group.xls {ref} n Up
{Rscript} {scripts}/Enrich_R.R {out_dir}/06_gener/02_gener_diff/ \\
    {out_dir}/06_gener/05_Enrich_R Limma {work_dir}/group.xls {ref} n Down
{Rscript} {scripts}/Enrich_R.R {out_dir}/06_gener/02_gener_diff/ \\
    {out_dir}/06_gener/05_Enrich_R DESeq2 {work_dir}/group.xls {ref} n Up
{Rscript} {scripts}/Enrich_R.R {out_dir}/06_gener/02_gener_diff/ \\
    {out_dir}/06_gener/05_Enrich_R DESeq2 {work_dir}/group.xls {ref} n Down

{Rscript} {scripts}/plot_GO_KEGG_all.R {work_dir}/group.xls Limma {out_dir}/06_gener/05_Enrich_R
{Rscript} {scripts}/plot_GO_KEGG_all.R {work_dir}/group.xls DESeq2 {out_dir}/06_gener/05_Enrich_R

{Rscript} {scripts}/Enrich_GASE.R {out_dir}/06_gener/02_gener_diff/ \\
    {out_dir}/06_gener/05_Enrich_R Limma {work_dir}/group.xls {ref} n
{Rscript} {scripts}/Enrich_GASE.R {out_dir}/06_gener/02_gener_diff/ \\
    {out_dir}/06_gener/05_Enrich_R DESeq2 {work_dir}/group.xls {ref} n
""".format(scripts=SCRIPTS,
            Rscript=RSCRIPT,
            out_dir=out_dir,
            work_dir=work_dir,
            ref=ref,
            tf_file=tf_file
        )
    )

    return task

def down_analysis_gener_f(tf_file, ref,
                      job_type, work_dir, out_dir):

    task = Task(
        id="work_down_analysis_gener_f",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{Rscript} {scripts}/Limma.R {out_dir}/07_gener_f/02_gener_f_diff \\
    {out_dir}/07_gener_f/01_gener_f_stat/gener_f_stat.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/Permutation.R {out_dir}/07_gener_f/02_gener_f_diff \\
    {out_dir}/07_gener_f/01_gener_f_stat/gener_f_stat.csv {tf_file} {work_dir}/group.xls

ls {out_dir}/07_gener_f/02_gener_f_diff/*csv |sed 's/.csv//' | while read i;do \\
    {Rscript} {scripts}/stat_diff.R ${{i}}.csv $i;done

mkdir -p {out_dir}/07_gener_f/02_gener_f_diff/veen

cat {work_dir}/fz.txt | sed 's/+/_/'| while read i;do \\
    {Rscript} {scripts}/diff_ggvenn.R {out_dir}/07_gener_f/02_gener_f_diff/${{i}}_Permutation.xls \\
    {out_dir}/07_gener_f/02_gener_f_diff/${{i}}_Limma.xls Permutation Limma {out_dir}/07_gener_f/02_gener_f_diff/veen/${{i}}_veen.png;done

{Rscript} {scripts}/Vol.R {out_dir}/07_gener_f/02_gener_f_diff/ \\
    {out_dir}/07_gener_f/03_gener_f_vol Limma {work_dir}/group.xls
{Rscript} {scripts}/Vol.R {out_dir}/07_gener_f/02_gener_f_diff/ \\
    {out_dir}/07_gener_f/03_gener_f_vol Permutation {work_dir}/group.xls

{Rscript} {scripts}/pheatmap.R {out_dir}/07_gener_f/02_gener_f_diff/ \\
    {out_dir}/07_gener_f/01_gener_f_stat/ {out_dir}/07_gener_f/04_gener_f_pheatmap Limma {work_dir}/group.xls
{Rscript} {scripts}/pheatmap.R {out_dir}/07_gener_f/02_gener_f_diff/ \\
    {out_dir}/07_gener_f/01_gener_f_stat/ {out_dir}/07_gener_f/04_gener_f_pheatmap Permutation {work_dir}/group.xls

{Rscript} {scripts}/Enrich_R.R {out_dir}/07_gener_f/02_gener_f_diff/ \\
    {out_dir}/07_gener_f/05_Enrich_R Limma {work_dir}/group.xls {ref} n Up
{Rscript} {scripts}/Enrich_R.R {out_dir}/07_gener_f/02_gener_f_diff/ \\
    {out_dir}/07_gener_f/05_Enrich_R Limma {work_dir}/group.xls {ref} n Down
{Rscript} {scripts}/Enrich_R.R {out_dir}/07_gener_f/02_gener_f_diff/ \\
    {out_dir}/07_gener_f/05_Enrich_R Permutation {work_dir}/group.xls {ref} n Up
{Rscript} {scripts}/Enrich_R.R {out_dir}/07_gener_f/02_gener_f_diff/ \\
    {out_dir}/07_gener_f/05_Enrich_R Permutation {work_dir}/group.xls {ref} n Down

{Rscript} {scripts}/plot_GO_KEGG_all.R {work_dir}/group.xls Limma {out_dir}/07_gener_f/05_Enrich_R
{Rscript} {scripts}/plot_GO_KEGG_all.R {work_dir}/group.xls Permutation {out_dir}/07_gener_f/05_Enrich_R

{Rscript} {scripts}/Enrich_GASE.R {out_dir}/07_gener_f/02_gener_f_diff/ \\
    {out_dir}/07_gener_f/05_Enrich_R Limma {work_dir}/group.xls {ref} n
{Rscript} {scripts}/Enrich_GASE.R {out_dir}/07_gener_f/02_gener_f_diff/ \\
    {out_dir}/07_gener_f/05_Enrich_R Permutation {work_dir}/group.xls {ref} n
""".format(scripts=SCRIPTS,
            Rscript=RSCRIPT,
            out_dir=out_dir,
            work_dir=work_dir,
            ref=ref,
            tf_file=tf_file
        )
    )

    return task

def down_analysis_REGB(tf_file, ref,
                      job_type, work_dir, out_dir):

    task = Task(
        id="work_down_analysis_REGB",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{Rscript} {scripts}/Permutation.R {out_dir}/08_REGB/02_REGB_diff \\
    {out_dir}/08_REGB/01_REGB_stat/REGB_stat.csv {tf_file} {work_dir}/group.xls

ls {out_dir}/08_REGB/02_REGB_diff/*csv |sed 's/.csv//' | while read i;do \\
    {Rscript} {scripts}/stat_diff.R ${{i}}.csv $i;done

{Rscript} {scripts}/Vol.R {out_dir}/08_REGB/02_REGB_diff \\
    {out_dir}/08_REGB/03_REGB_vol Permutation {work_dir}/group.xls

{Rscript} {scripts}/pheatmap.R {out_dir}/08_REGB/02_REGB_diff \\
    {out_dir}/08_REGB/01_REGB_stat/ {out_dir}/08_REGB/04_REGB_pheatmap Permutation {work_dir}/group.xls

{Rscript} {scripts}/Enrich_R.R {out_dir}/08_REGB/02_REGB_diff \\
    {out_dir}/08_REGB/05_Enrich_R Permutation {work_dir}/group.xls {ref} n Up
{Rscript} {scripts}/Enrich_R.R {out_dir}/08_REGB/02_REGB_diff \\
    {out_dir}/08_REGB/05_Enrich_R Permutation {work_dir}/group.xls {ref} n Down

{Rscript} {scripts}/plot_GO_KEGG_all.R {work_dir}/group.xls Permutation {out_dir}/08_REGB/05_Enrich_R

{Rscript} {scripts}/Enrich_GASE.R {out_dir}/08_REGB/02_REGB_diff \\
    {out_dir}/08_REGB/05_Enrich_R Permutation {work_dir}/group.xls {ref} n
""".format(scripts=SCRIPTS,
            Rscript=RSCRIPT,
            out_dir=out_dir,
            work_dir=work_dir,
            ref=ref,
            tf_file=tf_file
        )
    )

    return task

def down_analysis_REB(work_dir, out_dir, job_type, ref):

    task = Task(
        id="work_down_analysis_REB",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{Rscript} {scripts}/REB.R {out_dir}/01_REseq_detection/GATK_information.csv \\
    {work_dir}/group.xls {out_dir}/09_REB {ref}

""".format(Rscript=RSCRIPT,
            scripts=SCRIPTS,
            work_dir=work_dir,
            out_dir=out_dir,
            ref=ref
        )
    )

    return task

def down_analysis_RNA(tf_file, ref, gene_length, rna_count,
                      job_type, work_dir, out_dir):

    task = Task(
        id="work_down_analysis_RNA",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{Rscript} {scripts}/09_RNA_stat.R {gene_length} {rna_count} {out_dir}/10_RNA/01_RNA_stat

{perl} {scripts}/get_group.pl {work_dir}/group.xls {work_dir}/fz.txt \\
    {out_dir}/10_RNA/01_RNA_stat/RNA_TPM.xls {out_dir}/10_RNA/01_RNA_stat
{perl} {scripts}/get_group.pl {work_dir}/group.xls {work_dir}/fz.txt \\
    {out_dir}/10_RNA/01_RNA_stat/RNA_Count/RNA_Count.xls {out_dir}/10_RNA/01_RNA_stat/RNA_Count

{Rscript} {scripts}/10_RNA_stat_cor.R {out_dir}/10_RNA/01_RNA_stat/RNA_PCA/ \\
    {out_dir}/10_RNA/01_RNA_stat/RNA_TPM.csv {work_dir}/group.xls

{Rscript} {scripts}/Limma.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/01_RNA_stat/RNA_TPM.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/Permutation.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/01_RNA_stat/RNA_TPM.csv {tf_file} {work_dir}/group.xls

{Rscript} {scripts}/DEseq2.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/01_RNA_stat/RNA_Count/RNA_Count.csv {tf_file} {work_dir}/group.xls

ls {out_dir}/10_RNA/02_RNA_diff/*csv |sed 's/.csv//' | while read i;do \\
    {Rscript} {scripts}/stat_diff.R ${{i}}.csv $i;done

mkdir -p {out_dir}/10_RNA/02_RNA_diff/veen

cat {work_dir}/fz.txt | sed 's/+/_/'| while read i;do \\
    {Rscript} {scripts}/diff_ggvenn.R {out_dir}/10_RNA/02_RNA_diff/${{i}}_Permutation.xls \\
    {out_dir}/10_RNA/02_RNA_diff/${{i}}_DESeq2.xls Permutation DESeq2 {out_dir}/10_RNA/02_RNA_diff/veen/${{i}}_veen.png;done

{Rscript} {scripts}/Vol.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/03_RNA_vol/ Limma {work_dir}/group.xls
{Rscript} {scripts}/Vol.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/03_RNA_vol/ Permutation {work_dir}/group.xls
{Rscript} {scripts}/Vol.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/03_RNA_vol/ DESeq2 {work_dir}/group.xls

{Rscript} {scripts}/pheatmap.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/01_RNA_stat/ {out_dir}/10_RNA/04_RNA_pheatmap Limma {work_dir}/group.xls
{Rscript} {scripts}/pheatmap.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/01_RNA_stat/ {out_dir}/10_RNA/04_RNA_pheatmap Permutation {work_dir}/group.xls
{Rscript} {scripts}/pheatmap.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/01_RNA_stat/ {out_dir}/10_RNA/04_RNA_pheatmap DESeq2 {work_dir}/group.xls

{Rscript} {scripts}/Enrich_R.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/05_Enrich_R Limma {work_dir}/group.xls {ref} n Up
{Rscript} {scripts}/Enrich_R.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/05_Enrich_R Limma {work_dir}/group.xls {ref} n Down
{Rscript} {scripts}/Enrich_R.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/05_Enrich_R Permutation {work_dir}/group.xls {ref} n Up
{Rscript} {scripts}/Enrich_R.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/05_Enrich_R Permutation {work_dir}/group.xls {ref} n Down
{Rscript} {scripts}/Enrich_R.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/05_Enrich_R DESeq2 {work_dir}/group.xls {ref} n Up
{Rscript} {scripts}/Enrich_R.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/05_Enrich_R DESeq2 {work_dir}/group.xls {ref} n Down

{Rscript} {scripts}/plot_GO_KEGG_all.R {work_dir}/group.xls Limma {out_dir}/10_RNA/05_Enrich_R
{Rscript} {scripts}/plot_GO_KEGG_all.R {work_dir}/group.xls DESeq2 {out_dir}/10_RNA/05_Enrich_R
{Rscript} {scripts}/plot_GO_KEGG_all.R {work_dir}/group.xls Permutation {out_dir}/10_RNA/05_Enrich_R

{Rscript} {scripts}/Enrich_GASE.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/05_Enrich_R Limma {work_dir}/group.xls {ref} n
{Rscript} {scripts}/Enrich_GASE.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/05_Enrich_R Permutation {work_dir}/group.xls {ref} n
{Rscript} {scripts}/Enrich_GASE.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/10_RNA/05_Enrich_R DESeq2 {work_dir}/group.xls {ref} n
""".format(scripts=SCRIPTS,
            Rscript=RSCRIPT,
            perl=PERL_BIN,
            out_dir=out_dir,
            work_dir=work_dir,
            ref=ref,
            gene_length=gene_length,
            rna_count=rna_count,
            tf_file=tf_file
        )
    )

    return task

def down_analysis_association(work_dir, out_dir, job_type):

    task = Task(
        id="work_down_analysis_association",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 4",
        script="""
{Rscript} {scripts}/association.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/03_locr/02_locr_diff/ {out_dir}/11_association/01_RNA_locr DESeq2 {work_dir}/group.xls y

#{Rscript} {scripts}/association.R {out_dir}/10_RNA/02_RNA_diff/ \\
#    {out_dir}/04_locr_f/02_locr_f_diff/ {out_dir}/11_association/02_RNA_locr_f Permutation {work_dir}/group.xls y

{Rscript} {scripts}/association.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/05_RELB/02_RELB_diff/ {out_dir}/11_association/02_RNA_RELB Permutation {work_dir}/group.xls y

{Rscript} {scripts}/association.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/06_gener/02_gener_diff/ {out_dir}/11_association/03_RNA_gener DESeq2 {work_dir}/group.xls n

#{Rscript} {scripts}/association.R {out_dir}/10_RNA/02_RNA_diff/ \\
#    {out_dir}/07_gener_f/02_gener_f_diff/ {out_dir}/11_association/05_RNA_gener_f Permutation {work_dir}/group.xls n

{Rscript} {scripts}/association.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/08_REGB/02_REGB_diff/ {out_dir}/11_association/04_RNA_REGB Permutation {work_dir}/group.xls n
""".format(Rscript=RSCRIPT,
            scripts=SCRIPTS,
            work_dir=work_dir,
            out_dir=out_dir
        )
    )

    return task

def run_down_analysis(input, gatk_file, tf_file, job_type,
                gene_length, rna_count, ref,
                concurrent, refresh, work_dir, out_dir):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    dag = DAG("run_down_analysis")
    input_f = os.path.abspath(input)
    rna_count = os.path.abspath(rna_count)
    gatk_file = os.path.abspath(gatk_file)
    ref2 = ref
    if ref == "rn6":
        ref = "rno"

    stat_task, fz_file, group_file = down_analysis_stat(
        input=input_f, 
        gatk_file=gatk_file, 
        tf_file=tf_file, 
        gene_length=gene_length, 
        rna_count=rna_count,
        job_type=job_type, 
        work_dir=work_dir, 
        out_dir=out_dir
    )
    dag.add_task(stat_task)
    pcr_task = pcr_analysis(
        work_dir = work_dir, 
        out_dir = out_dir, 
        job_type = job_type
    )
    dag.add_task(pcr_task)
    pcr_task.set_upstream(stat_task)
    """
    diff_task = down_analysis_deseq(
        tf_file=tf_file,
        gene_length=gene_length, 
        rna_count=rna_count,
        job_type=job_type, 
        work_dir=work_dir, 
        out_dir=out_dir
    )
    dag.add_task(diff_task)
    diff_task.set_upstream(stat_task)
    """
    locr_task = down_analysis_locr(
        tf_file=tf_file,
        ref=ref,
        job_type=job_type, 
        work_dir=work_dir, 
        out_dir=out_dir
    )
    dag.add_task(locr_task)
    locr_task.set_upstream(stat_task)
    """
    locr_f_task = down_analysis_locr_f(
        tf_file=tf_file,
        ref=ref,
        job_type=job_type, 
        work_dir=work_dir, 
        out_dir=out_dir
    )
    dag.add_task(locr_f_task)
    locr_f_task.set_upstream(stat_task)
    """
    RELB_task = down_analysis_RELB(
        tf_file=tf_file,
        ref=ref,
        job_type=job_type, 
        work_dir=work_dir, 
        out_dir=out_dir
    )
    dag.add_task(RELB_task)
    RELB_task.set_upstream(stat_task)

    gener_task = down_analysis_gener(
        tf_file=tf_file,
        ref=ref,
        job_type=job_type, 
        work_dir=work_dir, 
        out_dir=out_dir
    )
    dag.add_task(gener_task)
    gener_task.set_upstream(stat_task)
    """
    gener_f_task = down_analysis_gener_f(
        tf_file=tf_file,
        ref=ref,
        job_type=job_type, 
        work_dir=work_dir, 
        out_dir=out_dir
    )
    dag.add_task(gener_f_task)
    gener_f_task.set_upstream(stat_task)
    """
    REGB_task = down_analysis_REGB(
        tf_file=tf_file,
        ref=ref,
        job_type=job_type, 
        work_dir=work_dir, 
        out_dir=out_dir
    )
    dag.add_task(REGB_task)
    REGB_task.set_upstream(stat_task)
    REB_task = down_analysis_REB(
        work_dir=work_dir, 
        out_dir=out_dir, 
        job_type=job_type,
        ref=ref2)
    dag.add_task(REB_task)
    REB_task.set_upstream(stat_task)
    
    RNA_task = down_analysis_RNA(
        tf_file=tf_file, 
        ref=ref, 
        gene_length=gene_length, 
        rna_count=rna_count,
        job_type=job_type, 
        work_dir=work_dir, 
        out_dir=out_dir
    )
    dag.add_task(RNA_task)
    RNA_task.set_upstream(stat_task)
    do_dag(dag, concurrent, refresh)

    dag = DAG("run_down_analysis2")

    association_task = down_analysis_association(
        work_dir=work_dir, 
        out_dir=out_dir,  
        job_type=job_type
    )
    dag.add_task(association_task)
    do_dag(dag, concurrent, refresh)
    return fz_file, group_file

def copyfigs(project_dir, report_dir):
    if os.path.exists(report_dir):
        shutil.rmtree(report_dir)
    shutil.copytree(os.path.join(BASE_DIR, "template"), report_dir)
    if ASAN:
        safecopy("%s/report.html" % BASE_DIR, report_dir)
    else:
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
    #d["table_loci_maf_go"] = readtbl_go(glob.glob("%s/04_Down_analysis/04_locr_f/05_Enrich_R/*/R_Enrichment/*_ALL.xls" % project_dir)[0], only5=True)
    #d["table_loci_maf_kegg"] = readtbl_go(glob.glob("%s/04_Down_analysis/04_locr_f/05_Enrich_R/*/R_Enrichment/*_KK.xls" % project_dir)[0], only5=True)
    #d["table_loci_maf_gsea"] = readtbl(glob.glob("%s/04_Down_analysis/04_locr_f/05_Enrich_R/*/GSEA_Enrichment/*_GO.xls" % project_dir)[0], only5=True)
    d["table_gene_go"] = readtbl_go(glob.glob("%s/04_Down_analysis/06_gener/05_Enrich_R/*/R_Enrichment/*_ALL.xls" % project_dir)[0], only5=True)
    d["table_gene_kegg"] = readtbl_go(glob.glob("%s/04_Down_analysis/06_gener/05_Enrich_R/*/R_Enrichment/*_KK.xls" % project_dir)[0], only5=True)
    d["table_gene_gsea"] = readtbl(glob.glob("%s/04_Down_analysis/06_gener/05_Enrich_R/*/GSEA_Enrichment/*_GO.xls" % project_dir)[0], only5=True)
    #d["table_gene_maf_go"] = readtbl_go(glob.glob("%s/04_Down_analysis/07_gener_f/05_Enrich_R/*/R_Enrichment/*_ALL.xls" % project_dir)[0], only5=True)
    #d["table_gene_maf_kegg"] = readtbl_go(glob.glob("%s/04_Down_analysis/07_gener_f/05_Enrich_R/*/R_Enrichment/*_KK.xls" % project_dir)[0], only5=True)
    #d["table_gene_maf_gsea"] = readtbl(glob.glob("%s/04_Down_analysis/07_gener_f/05_Enrich_R/*/GSEA_Enrichment/*_GO.xls" % project_dir)[0], only5=True)
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

def run_report(asan, contract, name, ref, work_dir, out_dir, fz_file, group_file,
                type, type_num):

    global ASAN
    ASAN = False if asan == 'no' else True
    html_base = "report_noAS_v1.2.html" if asan == 'no' else "report.html"

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

def all(input, ref, trim, thread, job_type, method,
                asan, contract, name, type, type_num,
                concurrent, refresh, work_dir, out_dir, qvalue):
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    prefix, read1, read2 = read_sample_file(input=input)
    input_path = os.path.abspath(input)

    clean1, clean2 = run_fastp(
        read1=read1,
        read2=read2,
        prefix=prefix,
        trim=trim,
        thread=thread,
        job_type=job_type,
        work_dir=mkdir(os.path.join(work_dir, "01_QC")),
        out_dir=mkdir(os.path.join(out_dir, "01_QC")),
        concurrent=concurrent,
        refresh=refresh,
        qvalue=qvalue
    )

    stat_fastp(
        input=input_path,
        work_dir=mkdir(os.path.join(work_dir, "01_QC")),
        out_dir=mkdir(os.path.join(out_dir, "01_QC")),
        concurrent=concurrent,
        refresh=refresh,
        job_type=job_type
    )

    hisat_bam = create_hisat2tasks(
        prefix=prefix,
        read1=clean1,
        read2=clean2,
        ref=ref_path[ref],
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir= mkdir(os.path.join(work_dir, "02_Hisat")),
        out_dir=mkdir(os.path.join(out_dir, "01_QC"))
    )

    gene_list = create_gene_count(
        gtf=ref_gtf[ref],
        hisat_bam=hisat_bam,
        prefix=prefix,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=mkdir(os.path.join(work_dir, "04_RNA_Count"))
    )

    gene_count_xls = merge_gene_count(
        input=input_path,
        work_dir=mkdir(os.path.join(work_dir, "04_RNA_Count")),
        out_dir=mkdir(os.path.join(out_dir, "03_RESeq_file")),
        concurrent=concurrent,
        refresh=refresh,
        job_type=job_type
    )

    gatk_result_xls = reseq_work(
        input=input_path, 
        method=method, 
        hisat_bam=hisat_bam,
        reference=ref, 
        prefix=prefix,
        work_dir=mkdir(os.path.join(work_dir, "03_RESeq")),
        out_dir=mkdir(os.path.join(out_dir, "03_RESeq_file")),
        thread=thread, 
        job_type=job_type, 
        concurrent=concurrent, 
        refresh=refresh
    )

    fz_file, group_file = run_down_analysis(
        input=input_path, 
        gatk_file=gatk_result_xls, 
        tf_file=ref_tf[ref], 
        job_type=job_type,
        gene_length=ref_gene_length[ref], 
        rna_count=gene_count_xls, 
        ref=ref,
        concurrent=concurrent, 
        refresh=refresh, 
        work_dir=mkdir(os.path.join(work_dir, "05_Down_analysis")),
        out_dir=mkdir(os.path.join(out_dir, "04_Down_analysis")))

    run_report(
        asan=asan, 
        contract=contract, 
        name=name, 
        ref=ref_name[ref], 
        work_dir=work_dir, 
        out_dir=out_dir, 
        fz_file=fz_file, 
        group_file=group_file,
        type=type, 
        type_num=type_num)

def run_all(args):
    all(
        input=args.input, 
        ref=args.reference, 
        trim=args.trim, 
        thread=args.thread, 
        method=args.method,
        asan=args.asan,
        contract=args.contract,
        name=args.name,
        type=args.type,
        type_num=args.type_num,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh, 
        work_dir=args.work_dir, 
        out_dir=args.out_dir, 
        qvalue=args.qvalue
    )

def add_args(parser):
    parser.add_argument("-i", "--input", metavar="STR", type=str, required=True,
        help="Input the name of the sample and Second generation sequencing path.")
    parser.add_argument("-ref", "--reference", metavar="FILE", type=str, 
        default="hg38",
        help="""Input the host's reference database.\
        Default: hg38.\
        you can choose mm10 or rn6""")
    parser.add_argument("--trim", metavar="INT", type=int, default=5,
        help="Set trim length, default=5")
    parser.add_argument("--qvalue", metavar="INT", type=int, default=20,
        help="The quality value that a base is qualified, default=20")
    parser.add_argument("-m", "--method", metavar="FILE", type=str, 
        default="gatk",
        help="""RE-seq method.\
        Default: gatk.\
        you can choose gatk; ; REDItools ; VarScan ; Sprint ; RED_ML.\
        if you want choose Multiple methods; you can use: gatk;REDItools;VarScan""")
    parser.add_argument("--thread", metavar="INT", type=int, default=4,
        help="Analysis of the number of threads used, default=4")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    parser.add_argument("--work_dir", metavar="DIR", default=".",
        help="Work directory (default: current directory)")
    parser.add_argument("--out_dir", metavar="DIR", default=".",
        help="Output directory (default: current directory)")

    parser.add_argument('-as', '--asan', help="AS analysis or not, yes or no, default no",
                        choices=['yes', 'no'], default='no')
    parser.add_argument('-cn', '--contract', help="contract number", required=True)
    parser.add_argument('-pn', '--name', help="project name", required=True)
    parser.add_argument('-type', '--type', help="padj or pval", default="padj")
    parser.add_argument('-type_num', '--type_num', 
        help="padj or pval threshold value", default="0.1")    

    return parser

def main():
        
    logging.basicConfig(
            stream=sys.stderr,
            level=logging.INFO,
            format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="""
work RE-Seq

version: %s
contact:  %s <%s>\
            """ % (__version__, " ".join(__author__), __email__))

    parser = add_args(parser)
    args = parser.parse_args()
    run_all(args)

if __name__ == '__main__':
    main()
