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

def create_fastp_task(read1, read2, prefix, trim, thread, job_type,
                      work_dir, out_dir, qvalue=20):

    json = mkdir("%s/json" % work_dir)
    qc_dir = mkdir("%s/QC_result" % out_dir)

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
    #stat_qc = os.path.join(out_dir, "%s.stat_qc.tsv" % prefix)
    clean1 = []
    clean2 = []
    for i in prefix:
        clean1.append(os.path.join(work_dir, "%s.clean.r1.fq.gz" % i))
        clean2.append(os.path.join(work_dir, "%s.clean.r2.fq.gz" % i))

    return tasks, clean1, clean2

def stat_fastp(input, work_dir, out_dir, job_type):

    json = mkdir("%s/json" % work_dir)

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

    return task

def create_hisat2_tasks(prefix, read1, read2, ref, thread,
                        job_type, work_dir, out_dir):

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
    hisat_bam = []
    for i in prefix:
        hisat_bam.append(os.path.join(work_dir, "%s.sort.bam" % i))
    return task, hisat_bam

def run_reseq_qc(input, ref, trim, thread, job_type,
                concurrent, refresh, work_dir, out_dir, qvalue=20):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    dag = DAG("reseq_qc")

    prefix, read1, read2 = read_sample_file(input=input)
    input_path = os.path.abspath(input)

    qc_tasks, clean1, clean2 = create_fastp_task(
        read1=read1,
        read2=read2,
        prefix=prefix,
        trim=trim,
        thread=thread,
        job_type=job_type,
        work_dir=mkdir(os.path.join(work_dir, "01_QC")),
        out_dir=out_dir,
        qvalue=qvalue
    )
    dag.add_task(*qc_tasks)

    stat_task = stat_fastp(
        input=input_path,
        work_dir=mkdir(os.path.join(work_dir, "01_QC")),
        out_dir=out_dir,
        job_type=job_type
    )
    dag.add_task(stat_task)
    stat_task.set_upstream(*qc_tasks)

    hisat_task, hisat_bam = create_hisat2_tasks(
        prefix=prefix,
        read1=clean1,
        read2=clean2,
        ref=ref,
        thread=thread,
        job_type=job_type,
        work_dir= mkdir(os.path.join(work_dir, "02_Hisat")),
        out_dir=out_dir
    )
    dag.add_task(*hisat_task)
    stat_task.set_downstream(*hisat_task)

    do_dag(dag, concurrent, refresh)

    return prefix, read1, read2, hisat_bam

def reseq_qc(args):

    prefix, read1, read2, hisat_bam = run_reseq_qc(
        input=args.input,
        ref=ref_path[args.reference],
        trim=args.trim,
        thread=args.thread,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        qvalue=args.qvalue
    )

    return 0

def add_reseq_qc_args(parser):

    parser.add_argument("-i", "--input", metavar="STR", type=str, required=True,
        help="Input the name of the sample and Second generation sequencing path.")
    parser.add_argument("-ref", "--reference", metavar="FILE", type=str, 
        default="hg38",
        help="""Input the host's reference database.\
        Default: hg38""")
    parser.add_argument("--trim", metavar="INT", type=int, default=5,
        help="Set trim length, default=5")
    parser.add_argument("--qvalue", metavar="INT", type=int, default=20,
        help="The quality value that a base is qualified, default=20")
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

    parser = add_reseq_qc_args(parser)
    args = parser.parse_args()
    reseq_qc(args)


if __name__ == "__main__":
    main()
