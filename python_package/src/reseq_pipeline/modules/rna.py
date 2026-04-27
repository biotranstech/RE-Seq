import os

from .dagflow import DAG, Task, ParallelTask, do_dag
from ..toolchain import env_or_default, packaged_script

HTSEQ = env_or_default('HTSEQ', 'htseq-count')
PERL_BIN = env_or_default('PERL_BIN', 'perl')


def create_gene_count(gtf, hisat_bam, input, prefix, thread, job_type,
                      work_dir, out_dir, concurrent, refresh):

    input_path = os.path.abspath(input)
    merge_gene = packaged_script('merge_gene.pl')
    dag = DAG('reseq_rna_count')
    tasks = ParallelTask(
        id='work_gene_count_htseq',
        work_dir=work_dir,
        type=job_type,
        option='-pe smp 4',
        script=f'''
{HTSEQ} -f bam -s no -m intersection-nonempty \
    {{hisat_bam}} {gtf} > {{prefix}}.gene.xls
''',
        hisat_bam=hisat_bam,
        prefix=prefix,
    )
    dag.add_task(*tasks)
    do_dag(dag, concurrent, refresh)

    dag = DAG('reseq_merge_rna_count')
    task = Task(
        id='work_merge_gene_count',
        work_dir=work_dir,
        type=job_type,
        option='-pe smp 4',
        script=f'''
{PERL_BIN} {merge_gene} {input_path} > RNA_Count.xls
cp RNA_Count.xls {out_dir}
''',
    )
    dag.add_task(task)
    do_dag(dag, concurrent, refresh)
    return os.path.join(work_dir, 'RNA_Count.xls')
