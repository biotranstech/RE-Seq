import os

from .dagflow import DAG, Task, ParallelTask, do_dag
from .utils import mkdir
from ..toolchain import env_or_default, packaged_script

GATK_BIN = env_or_default('GATK_BIN', 'gatk')
SAMTOOLS_BIN = env_or_default('SAMTOOLS_BIN', 'samtools')
PERL_BIN = env_or_default('PERL_BIN', 'perl')
PYTHON_BIN = env_or_default('PYTHON_BIN', 'python')
JAVA_BIN = env_or_default('JAVA_BIN', 'java')
ANNOVAR_BIN = env_or_default('ANNOVAR_BIN', 'annovar')
BGZIP = env_or_default('BGZIP', 'bgzip')
TABIX = env_or_default('TABIX', 'tabix')


def work_bam_gatk(hisat_bam, prefix, ref, input, thread, job_type,
                  annovar_ref, work_dir, out_dir, concurrent, refresh):
    tmp = mkdir(os.path.join(work_dir, 'tmp'))
    filter_vcf = packaged_script('filter_vcf.pl')
    reads_pl = packaged_script('Reads.pl')
    reads_py = packaged_script('Reads.py')
    merge_sample = packaged_script('merge_sample.py')
    tasks = ParallelTask(
        id='work_bam',
        work_dir=work_dir,
        type=job_type,
        option='-pe smp 4',
        script=f'''
export PATH={JAVA_BIN}:$PATH

{GATK_BIN} MarkDuplicates -I {{hisat_bam}} \
    -O {{prefix}}.deduped.bam -M {{prefix}}.marked_dup_metrics.txt \
    --CREATE_INDEX true --REMOVE_DUPLICATES true \
    --TAG_DUPLICATE_SET_MEMBERS true --TMP_DIR {tmp}

{GATK_BIN} AddOrReplaceReadGroups -I {{prefix}}.deduped.bam \
    -O {{prefix}}.picard.bam -LB {{prefix}} -PL illumina \
    -PU {{prefix}} -SM {{prefix}}

{SAMTOOLS_BIN} index {{prefix}}.picard.bam

{GATK_BIN} SplitNCigarReads -R {ref} -I {{prefix}}.picard.bam \
    -O {{prefix}}.dedup_split.bam

{GATK_BIN} HaplotypeCaller -R {ref} -ERC GVCF -I {{prefix}}.dedup_split.bam \
    -O {{prefix}}.g.vcf

{GATK_BIN} GenotypeGVCFs -R {ref} -V {{prefix}}.g.vcf -O {{prefix}}.vcf

{GATK_BIN} IndexFeatureFile -F {{prefix}}.vcf

{GATK_BIN} SelectVariants -V {{prefix}}.vcf -select-type SNP -O {{prefix}}.snp.vcf
{GATK_BIN} SelectVariants -V {{prefix}}.vcf -select-type INDEL -O {{prefix}}.indel.vcf

{GATK_BIN} IndexFeatureFile -F {{prefix}}.snp.vcf
{GATK_BIN} IndexFeatureFile -F {{prefix}}.indel.vcf

{GATK_BIN} BaseRecalibrator -R {ref} -I {{prefix}}.dedup_split.bam \
    --known-sites {{prefix}}.snp.vcf \
    --known-sites {{prefix}}.indel.vcf \
    -O {{prefix}}.recal_data.table

{GATK_BIN} ApplyBQSR -R {ref} -I {{prefix}}.dedup_split.bam \
    --bqsr-recal-file {{prefix}}.recal_data.table \
    -O {{prefix}}.bqsr.bam

{SAMTOOLS_BIN} index {{prefix}}.bqsr.bam

{GATK_BIN} SplitNCigarReads -R {ref} -I {{prefix}}.bqsr.bam \
    -O {{prefix}}.bqsr.bed.bam

{GATK_BIN} HaplotypeCaller -R {ref} -ERC GVCF -I {{prefix}}.bqsr.bed.bam \
    -O {{prefix}}.bed.g.vcf

{GATK_BIN} GenotypeGVCFs -R {ref} -V {{prefix}}.bed.g.vcf \
    -O {{prefix}}.bed.vcf

rm -rf {{prefix}}.dedup_split.bam {{prefix}}.picard.bam {{prefix}}.deduped.bam
rm -rf {{prefix}}.snp* {{prefix}}.indel*

{BGZIP} {{prefix}}.bed.vcf
{TABIX} {{prefix}}.bed.vcf.gz

{GATK_BIN} SelectVariants -V {{prefix}}.bed.vcf.gz -select-type SNP \
    -O {{prefix}}.snp.vcf.gz
{GATK_BIN} SelectVariants -V {{prefix}}.bed.vcf.gz -select-type INDEL \
    -O {{prefix}}.indel.vcf.gz

{GATK_BIN} VariantFiltration -V {{prefix}}.snp.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O {{prefix}}.vcf

{PERL_BIN} {filter_vcf} {{prefix}}.vcf > {{prefix}}.hardfiltered.biallelic.vcf

{PERL_BIN} {ANNOVAR_BIN}/convert2annovar.pl --format vcf4 \
    --includeinfo {{prefix}}.hardfiltered.biallelic.vcf \
    --allsample --withfreq > {{prefix}}.avinput

{PERL_BIN} {ANNOVAR_BIN}/table_annovar.pl {{prefix}}.avinput \
    {ANNOVAR_BIN}/{annovar_ref}/ -buildver {annovar_ref} -out {{prefix}} \
    -remove -protocol refGene -operation g -nastring . \
    -csvout --thread {thread}

{PERL_BIN} {reads_pl} {{prefix}}
{PYTHON_BIN} {reads_py} {{prefix}}.Reads.csv \
    {{prefix}}.{annovar_ref}_multianno.csv {{prefix}}
''',
        prefix=prefix,
        hisat_bam=hisat_bam,
    )
    work_dir = os.path.abspath(work_dir)
    task = Task(
        id='work_merge_GATK',
        work_dir=work_dir,
        type=job_type,
        option='-pe smp 4',
        script=f'''
{PYTHON_BIN} {merge_sample} {input}
cp GATK_information.csv {out_dir}
''',
    )
    gatk_file = os.path.join(work_dir, 'GATK_information.csv')
    dag = DAG('REseq_gatk')
    dag.add_task(*tasks)
    dag.add_task(task)
    task.set_upstream(*tasks)
    do_dag(dag, concurrent, refresh)
    return gatk_file
