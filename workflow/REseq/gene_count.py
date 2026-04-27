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

def create_gene_count(gtf, hisat_bam, prefix, thread, job_type,
                      work_dir):

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
    gene_list = []
    for i in prefix:
        gene_list.append(os.path.join(work_dir, "%s.gene.xls" % i))

    return tasks, gene_list

def merge_gene_count(input, work_dir, out_dir, job_type):

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

    return task

def run_gene_count(input, gtf, thread, job_type,
                concurrent, refresh, work_dir, out_dir):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    input_file = os.path.abspath(input)

    dag = DAG("stat_gene_count")

    prefix, read1, read2 = read_sample_file(input=input_file)
    hisat_bam=[]
    for i in prefix:
        hisat_bam.append(os.path.join(work_dir, "02_Hisat", "%s.sort.bam" % i))


    htseq_tasks, gene_list = create_gene_count(
        gtf=gtf,
        hisat_bam=hisat_bam,
        prefix=prefix,
        thread=thread,
        job_type=job_type,
        work_dir=mkdir(os.path.join(work_dir, "04_RNA_Count"))
    )
    dag.add_task(*htseq_tasks)

    merge_task = merge_gene_count(
        input=input_file,
        work_dir=mkdir(os.path.join(work_dir, "04_RNA_Count")),
        out_dir=mkdir(os.path.join(out_dir, "03_RESeq_file")),
        job_type=job_type
    )
    dag.add_task(merge_task)
    merge_task.set_upstream(*htseq_tasks)

    do_dag(dag, concurrent, refresh)

    return prefix

def gene_count(args):

    prefix = run_gene_count(
        input=args.input,
        gtf=ref_gtf[args.reference],
        thread=args.thread,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir
    )

    return 0

def add_gene_args(parser):

    parser.add_argument("-i", "--input", metavar="STR", type=str, required=True,
        help="Input the name of the sample and Second generation sequencing path.")
    parser.add_argument("-ref", "--reference", metavar="FILE", type=str, 
        default="hg38",
        help="""Input the host's reference database.\
        Default: hg38.\
        you can choose mm10 or rn6""")
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

    parser = add_gene_args(parser)
    args = parser.parse_args()
    gene_count(args)


if __name__ == "__main__":
    main()
