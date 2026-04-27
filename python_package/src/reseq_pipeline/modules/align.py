import os

from .dagflow import DAG, ParallelTask, do_dag
from ..toolchain import env_or_default

HISAT2_BIN = env_or_default('HISAT2_BIN', 'hisat2')
SAMTOOLS_BIN = env_or_default('SAMTOOLS_BIN', 'samtools')


def create_hisat2tasks(prefix, read1, read2, ref, thread,
                       job_type, work_dir, out_dir, concurrent, refresh):

    dag = DAG('reseq_hisat')
    task = ParallelTask(
        id='hisat2',
        work_dir=work_dir,
        type=job_type,
        option='-pe smp 4',
        script=f'''
{HISAT2_BIN} --threads {thread} -x {ref} \
  -1 {{read1}} -2 {{read2}} | \
  {SAMTOOLS_BIN} sort -O bam -@ {thread} -o {{prefix}}.sort.bam
''',
        read1=read1,
        read2=read2,
        prefix=prefix,
    )
    dag.add_task(*task)
    do_dag(dag, concurrent, refresh)
    return [os.path.join(work_dir, f'{i}.sort.bam') for i in prefix]
