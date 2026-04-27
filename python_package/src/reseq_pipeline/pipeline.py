from __future__ import annotations

import os
from pathlib import Path

import yaml

from .modules.align import create_hisat2tasks
from .modules.down import run_down_analysis
from .modules.qc import run_fastp
from .modules.report import run_report
from .modules.rna import create_gene_count
from .modules.utils import mkdir
from .modules.variant import work_bam_gatk
from .toolchain import validate_toolchain


def load_reference_config(ref_yaml: str, genome: str) -> dict:
    yaml_file = Path(ref_yaml).expanduser().resolve()
    with yaml_file.open('r', encoding='utf-8') as f:
        data = yaml.safe_load(f) or {}
    if genome not in data:
        raise KeyError(f"Genome '{genome}' not found in reference yaml: {yaml_file}")
    ref_info = data[genome]
    required = ['path', 'fasta', 'gtf', 'tf', 'length', 'name']
    missing = [k for k in required if k not in ref_info]
    if missing:
        raise KeyError(f"Missing keys in reference config for {genome}: {', '.join(missing)}")
    return ref_info


def run_pipeline(
    input: str,
    ref_yaml: str,
    genome: str = 'hg38',
    out_dir: str = '.',
    trim: int = 5,
    qvalue: int = 20,
    thread: int = 4,
    concurrent: int = 10,
    refresh: int = 30,
    job_type: str = 'local',
    work_dir: str = '.',
    contract: str | None = None,
    name: str | None = None,
    diff_type: str = 'padj',
    diff_type_num: str = '0.1',
):
    validate_toolchain(strict=True)
    ref_info = load_reference_config(ref_yaml=ref_yaml, genome=genome)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    input_path = Path(input).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f'Input manifest not found: {input_path}')

    print('[1/6] QC...')
    prefix, clean1, clean2 = run_fastp(
        input=str(input_path),
        trim=trim,
        qvalue=qvalue,
        thread=thread,
        concurrent=concurrent,
        refresh=refresh,
        job_type=job_type,
        work_dir=mkdir(os.path.join(work_dir, '01_QC')),
        out_dir=mkdir(os.path.join(out_dir, '01_QC')),
    )

    print('[2/6] Alignment...')
    bam = create_hisat2tasks(
        prefix=prefix,
        read1=clean1,
        read2=clean2,
        ref=ref_info['path'],
        thread=thread,
        concurrent=concurrent,
        refresh=refresh,
        job_type=job_type,
        work_dir=mkdir(os.path.join(work_dir, '02_Hisat')),
        out_dir=mkdir(os.path.join(out_dir, '02_Align')),
    )

    print('[3/6] Gene count...')
    gene_count = create_gene_count(
        gtf=ref_info['gtf'],
        hisat_bam=bam,
        input=str(input_path),
        prefix=prefix,
        thread=thread,
        concurrent=concurrent,
        refresh=refresh,
        job_type=job_type,
        work_dir=mkdir(os.path.join(work_dir, '03_RNA_Count')),
        out_dir=mkdir(os.path.join(out_dir, '03_RESeq_file')),
    )

    print('[4/6] Variant calling...')
    result = work_bam_gatk(
        hisat_bam=bam,
        prefix=prefix,
        ref=ref_info['fasta'],
        input=str(input_path),
        annovar_ref=genome,
        thread=thread,
        concurrent=concurrent,
        refresh=refresh,
        job_type=job_type,
        work_dir=mkdir(os.path.join(work_dir, '04_RESeq')),
        out_dir=mkdir(os.path.join(out_dir, '03_RESeq_file')),
    )

    print('[5/6] Downstream analysis...')
    fz_file, group_file = run_down_analysis(
        input=str(input_path),
        gatk_file=result,
        tf_file=ref_info['tf'],
        gene_length=ref_info['length'],
        rna_count=gene_count,
        ref=genome,
        concurrent=concurrent,
        refresh=refresh,
        job_type=job_type,
        work_dir=mkdir(os.path.join(work_dir, '05_Down_analysis')),
        out_dir=mkdir(os.path.join(out_dir, '04_Down_analysis')),
    )

    print('[6/6] Report...')
    if contract and name:
        try:
            run_report(
                contract=contract,
                name=name,
                ref=ref_info['name'],
                fz_file=fz_file,
                group_file=group_file,
                type=diff_type,
                type_num=diff_type_num,
                work_dir=work_dir,
                out_dir=out_dir,
            )
        except FileNotFoundError as exc:
            print(f'Skip report generation: {exc}')
    else:
        print('Skip report generation because contract or name was not provided.')

    print('Pipeline finished!')
    return result
