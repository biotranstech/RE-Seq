import os

from .dagflow import DAG, Task, ParallelTask, do_dag
from .utils import mkdir, read_sample_file
from ..toolchain import env_or_default, packaged_script

FASTP_BIN = env_or_default('FASTP_BIN', 'fastp')
FASTQC_BIN = env_or_default('FASTQC_BIN', 'fastqc')
PYTHON_BIN = env_or_default('PYTHON_BIN', 'python')


def run_fastp(input, trim, thread, job_type,
              work_dir, out_dir, qvalue, concurrent, refresh):
    print('Running fastp...')
    prefix, read1, read2 = read_sample_file(input=input)
    input_path = os.path.abspath(input)
    json = mkdir(f'{work_dir}/json')
    qc_dir = mkdir(f'{out_dir}/QC_result')
    plot_fastqc = packaged_script('plot_fastqc.py')
    stat_fastp = packaged_script('stat_fastp.py')

    dag = DAG('reseq_fastp')
    tasks = ParallelTask(
        id='ngs_qc',
        work_dir=work_dir,
        type=job_type,
        option='-pe smp 4',
        script=f'''
{FASTP_BIN} -i {{read1}} -I {{read2}} \
    -o {{prefix}}.clean.r1.fq.gz -O {{prefix}}.clean.r2.fq.gz \
    -w {thread} -n 0 -f {trim} -F {trim} -t {trim} -T {trim} \
    -q {qvalue} --json {{prefix}}_fastp.json

{FASTQC_BIN} {{prefix}}.clean.r1.fq.gz {{prefix}}.clean.r2.fq.gz \
    -t {thread} --extract -o {work_dir}

cp {{prefix}}*.json {json}/.

{PYTHON_BIN} {plot_fastqc} \
    -r1 {{prefix}}.clean.r1_fastqc/fastqc_data.txt \
    -r2 {{prefix}}.clean.r2_fastqc/fastqc_data.txt \
    --name {{prefix}}
cp {{prefix}}.base_content.* {qc_dir}
cp {{prefix}}.base_quality.* {qc_dir}
''',
        read1=read1,
        read2=read2,
        prefix=prefix,
    )
    clean1 = [os.path.join(work_dir, f'{i}.clean.r1.fq.gz') for i in prefix]
    clean2 = [os.path.join(work_dir, f'{i}.clean.r2.fq.gz') for i in prefix]
    dag.add_task(*tasks)
    do_dag(dag, concurrent, refresh)

    dag = DAG('reseq_stat_fastp')
    task = Task(
        id='stat_qc',
        work_dir=work_dir,
        type=job_type,
        option='-pe smp 4',
        script=f'''
{PYTHON_BIN} {stat_fastp} {input_path} {json} > QC_stat.xls
cp QC_stat.xls {out_dir}
''',
    )
    dag.add_task(task)
    do_dag(dag, concurrent, refresh)

    return prefix, clean1, clean2
