import os
from .dagflow import DAG, Task, ParallelTask, do_dag
from importlib.resources import files
from .utils import mkdir

RSCRIPT = os.environ.get("RSCRIPT", "Rscript")
PERL_BIN = os.environ.get("PERL_BIN", "perl")
SCRIPTS = files("reseq_pipeline").joinpath("scripts")

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

{Rscript} {scripts}/03_locr_stat.R {out_dir}/01_REseq_detection/GATK_information.csv \\
    {out_dir}/03_locr/01_locr_stat
{perl} {scripts}/get_group.pl {work_dir}/group.xls {work_dir}/fz.txt \\
    {out_dir}/03_locr/01_locr_stat/locr_stat.xls {out_dir}/03_locr/01_locr_stat

{Rscript} {scripts}/05_RELB_stat.R {gene_length} \\
    {out_dir}/01_REseq_detection/GATK_information.csv {out_dir}/05_RELB/01_RELB_stat {rna_count}
{perl} {scripts}/get_group.pl {work_dir}/group.xls {work_dir}/fz.txt \\
    {out_dir}/05_RELB/01_RELB_stat/RELB_stat.xls {out_dir}/05_RELB/01_RELB_stat

{Rscript} {scripts}/06_gener_stat.R {out_dir}/01_REseq_detection/GATK_information.csv \\
    {out_dir}/06_gener/01_gener_stat
{perl} {scripts}/get_group.pl {work_dir}/group.xls {work_dir}/fz.txt \\
    {out_dir}/06_gener/01_gener_stat/gener_stat.xls {out_dir}/06_gener/01_gener_stat

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

{Rscript} {scripts}/association.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/05_RELB/02_RELB_diff/ {out_dir}/11_association/02_RNA_RELB Permutation {work_dir}/group.xls y

{Rscript} {scripts}/association.R {out_dir}/10_RNA/02_RNA_diff/ \\
    {out_dir}/06_gener/02_gener_diff/ {out_dir}/11_association/03_RNA_gener DESeq2 {work_dir}/group.xls n

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

    locr_task = down_analysis_locr(
        tf_file=tf_file,
        ref=ref,
        job_type=job_type, 
        work_dir=work_dir, 
        out_dir=out_dir
    )
    dag.add_task(locr_task)
    locr_task.set_upstream(stat_task)

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