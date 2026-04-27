#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import logging
import argparse

from collections import OrderedDict
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../"))
from REseq.config import *
from dagflow import DAG, Task, ParallelTask, do_dag
from REseq.common import check_path, mkdir, read_tsv, read_files, read_sample_file


LOG = logging.getLogger(__name__)

__author__ = ("Jiacheng Yi",)
__email__ = "..."
__version__ = "v1.0.0"

def work_bam_gatk(hisat_bam, prefix, ref, thread, job_type, 
                annovar_ref, work_dir):

    #hisat_bam = []
    #for i in prefix:
    #   hisat_bam.append(os.path.join(work_hisat_dir, "%s.sort.bam" % i))

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
    #stat_qc = os.path.join(out_dir, "%s.stat_qc.tsv" % prefix)
    GATK_result = []
    bam_list = []
    work_dir = os.path.abspath(work_dir)
    for i in prefix:
        GATK_result.append(os.path.join(work_dir, "%s.starded.csv" % i))
        bam_list.append(os.path.join(work_dir, "%s.bqsr.bed.bam" % i))

    return tasks, GATK_result, bam_list

def merge_gatk(input, work_dir, out_dir, job_type):

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
    return task


def work_REDItools(bam, prefix, ref, thread, work_dir, job_type):

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
        bam=bam,
        prefix=prefix
    )

    REDItool_result = []
    for i in prefix:
        REDItool_result.append(os.path.join(work_dir, "%s_RED.txt" % i))
    return tasks, REDItool_result

def work_VarScan(bam, bed, ref, prefix, job_type, work_dir):

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
        bam=bam,
        prefix=prefix
    )

    VarScan_result = []
    for i in prefix:
        VarScan_result.append(os.path.join(work_dir, "%s.tumor" % i))

    return tasks, VarScan_result

def work_Sprint(bam, ref, prefix, gtf, ref_repat, job_type, work_dir):

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
        bam=bam,
        prefix=prefix
    )

    sprint_result = []
    for i in prefix:
        sprint_result.append(os.path.join(work_dir, "%s.txt" % i))

    return task, tasks, sprint_result

def work_RED_ML(bam, ref, dbsnp, ref_repat, ref_alu, ref_snp, 
                job_type, work_dir):

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
        bam=bam
    )

    return tasks

def run_reseq_identify_gatk(input, prefix, ref, thread, job_type, annovar_ref,
                concurrent, refresh, work_dir, out_dir):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    prefix, read1, read2 = read_sample_file(input=input)
    hisat_bam=[]
    input_path = os.path.abspath(work_dir)
    for i in prefix:
        hisat_bam.append(os.path.join(input_path, "02_Hisat", "%s.sort.bam" % i))


    #prefix, read1, read2 = read_sample_file(input=input)
    dag = DAG("REseq_gatk")
    gatk_tasks, GATK_result, bam_list = work_bam_gatk(
        hisat_bam=hisat_bam,
        prefix=prefix, 
        ref=ref, 
        thread=thread, 
        job_type=job_type, 
        annovar_ref=annovar_ref, 
        work_dir=mkdir(os.path.join(work_dir, "03_RESEQ", "GATK"))
    )

    merge_task = merge_gatk(
        input=input,
        work_dir=mkdir(os.path.join(work_dir, "03_RESEQ", "GATK")),
        out_dir=mkdir(os.path.join(out_dir, "03_RESeq_file")),
        job_type=job_type
    )
        
    dag.add_task(*gatk_tasks)
    dag.add_task(merge_task)
    merge_task.set_upstream(*gatk_tasks)
    do_dag(dag, concurrent, refresh)

    return GATK_result, bam_list

def run_reseq_identify_REDItools(gatk_bam_list, prefix, ref, thread, job_type,
                concurrent, refresh, work_dir):

    
    dag = DAG("REseq_REDItools")
    REDItools_tasks, REDItools_result = work_REDItools(
        prefix=prefix,
        bam=gatk_bam_list,
        ref=ref,
        thread=thread,
        job_type=job_type,
        work_dir= mkdir(os.path.join(work_dir, "03_RESEQ", "REDItools"))
    )
    dag.add_task(*REDItools_tasks)
    do_dag(dag, concurrent, refresh)
    return REDItools_result

def run_reseq_identify_VarScan(gatk_bam_list, prefix, ref, bed, thread, job_type,
                concurrent, refresh, work_dir):

    
    dag = DAG("REseq_VarScan")
    VarScan_tasks, VarScan_result = work_VarScan(
        prefix=prefix,
        bam=gatk_bam_list,
        ref=ref,
        bed=bed,
        thread=thread,
        job_type=job_type,
        work_dir= mkdir(os.path.join(work_dir, "03_RESEQ", "VarScan"))
    )
    dag.add_task(*VarScan_tasks)
    do_dag(dag, concurrent, refresh)
    return VarScan_result

def run_reseq_identify_Sprint(gatk_bam_list, prefix, ref, gtf, ref_repat, thread, job_type,
                concurrent, refresh, work_dir):

    dag = DAG("REseq_Sprint")
    index_task, Sprint_tasks, sprint_result = work_Sprint(
        prefix=prefix,
        bam=gatk_bam_list,
        ref=ref,
        gtf=gtf,
        ref_repat=ref_repat,
        thread=thread,
        job_type=job_type,
        work_dir= mkdir(os.path.join(work_dir, "03_RESEQ", "Sprint"))
    )
    dag.add_task(index_task)
    dag.add_task(*Sprint_tasks)
    index_task.set_upstream(*Sprint_tasks)

    do_dag(dag, concurrent, refresh)
    return sprint_result

def run_reseq_identify_RED_ML(gatk_bam_list, prefix, ref, dbsnp, ref_repat,
                ref_alu, ref_snp,
                thread, job_type, concurrent, refresh, work_dir):

    dag = DAG("REseq_RED_ML")
    RED_ML_tasks = work_RED_ML(
        prefix=prefix,
        bam=gatk_bam_list,
        ref=ref,
        dbsnp=dbsnp,
        ref_repat=ref_repat,
        ref_alu=ref_alu,
        ref_snp=ref_snp,
        thread=thread,
        job_type=job_type,
        work_dir= mkdir(os.path.join(work_dir, "03_RESEQ", "RED_ML"))
    )
    dag.add_task(*RED_ML_tasks)

    do_dag(dag, concurrent, refresh)
    return 0

def reseq_work(args):

    GATK_result, bam_list = run_reseq_identify_gatk(
        input=args.input,
        prefix=prefix,
        ref=ref_fasta[args.reference],
        annovar_ref=args.reference,
        thread=args.thread,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir
    )

    work_list = args.method.split(";")
    if "REDItools" in work_list:
        REDItools_result = run_reseq_identify_REDItools(
            gatk_bam_list=bam_list,
            prefix=prefix,
            ref=ref_path[args.reference],
            thread=args.thread,
            job_type=args.job_type,
            concurrent=args.concurrent,
            refresh=args.refresh,
            work_dir=args.work_dir
        )

    if "VarScan" in work_list:
        VarScan_result = run_reseq_identify_VarScan(
            gatk_bam_list=bam_list,
            prefix=prefix,
            ref=ref_path[args.reference],
            bed=ref_bed[args.reference],
            thread=args.thread,
            job_type=args.job_type,
            concurrent=args.concurrent,
            refresh=args.refresh,
            work_dir=args.work_dir
        ) 

    if "Sprint" in work_list:
        sprint_result = run_reseq_identify_Sprint(
            gatk_bam_list=bam_list,
            prefix=prefix,
            ref=ref_path[args.reference],
            gtf=ref_gtf[args.reference],
            ref_repat=ref_repeat[args.reference],
            thread=args.thread,
            job_type=args.job_type,
            concurrent=args.concurrent,
            refresh=args.refresh,
            work_dir=args.work_dir
        ) 

    if "RED_ML" in work_list:
        sprint_result = run_reseq_identify_RED_ML(
            gatk_bam_list=bam_list,
            prefix=prefix,
            ref=ref_path[args.reference],
            dbsnp=ref_dbsnp[args.reference],
            ref_repat=ref_repeat[args.reference],
            ref_alu=ref_alu[args.reference],
            ref_snp=ref_snp_d[args.reference],
            thread=args.thread,
            job_type=args.job_type,
            concurrent=args.concurrent,
            refresh=args.refresh,
            work_dir=args.work_dir
        )

    return 0

def add_reseq_args(parser):

    parser.add_argument("-i", "--input", metavar="STR", type=str, required=True,
        help="Input the name of the sample and Second generation sequencing path.")
    parser.add_argument("-ref", "--reference", metavar="FILE", type=str, 
        default="hg38",
        help="""Input the host's reference database.\
        Default: hg38.\
        you can choose mm10 or rn6""")
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

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_reseq_args(parser)
    args = parser.parse_args()
    reseq_work(args)


if __name__ == "__main__":
    main()