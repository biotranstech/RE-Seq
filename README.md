RE-Seq: RNA Editing Analysis Workflow
=====
## RE-Seq Overview
RE-Seq is an automated RNA editing analysis pipeline built with Python DAGflow. It provides a full end-to-end solution from raw RNA-seq data to biological interpretation.
### What is RNA editing
RNA editing is a dynamic post-transcriptional modification with significant implications for gene regulation and disease mechanisms.
### Key Features
1. RNA-seq quality control
2. Alignment to reference genome
3. RNA editing detection
4. Multi-level statistical analysis
5. Functional enrichment analysis
6. Publication-ready visualization
7. Automatic HTML report generation
### Why Use RE-Seq
RNA editing analysis is typically complex and fragmented:
| Problem                 | Traditional Workflow                   |
| ----------------------- | -------------------------------------- |
| Multiple tools chaining | fastp + hisat2 + GATK + custom scripts |
| Poor reproducibility    | Manual operations are error-prone      |
| Complex pipeline        | Multiple steps and format conversions  |
| Difficult visualization | Requires additional coding             |

👉 RE-Seq addresses these challenges by:
1. Automating the entire workflow (from FASTQ to final report)
2. Standardizing analysis steps (minimizing human-induced errors)
3. Providing built-in statistical analysis and visualization
4. Automatically generating deliverable-ready reports
## Workflow Structure
The project structure is as follows (core module description) :
```
├── all.py                    #Workflow Complete process
├── dagflow                   #python DAGflow package
├── report_noAS.html          #Report template
├── report_RNA_edic.py        #Automated report generation scripts
├── report_utils.py           #Report generating python modules
├── REseq                     #RE-Seq python modules
│   ├── common.py             #Public files handle python packages
│   ├── config.py             #Process software configuration file
│   ├── down_analysis.py      #RE-Seq downstream analysis module
│   ├── gene_count.py         #RNA expression analysis module
│   ├── qc_hista.py           #Quality control and comparison module
│   ├── report_utils.py       #Report generating python modules
│   └── work_bam_reseq.py     #Get the RNA editing information module
├── scripts                   #Analysis script path
├── template                  #html configuration file
│   ├── asset                 #html related js ccs file
│   └── images                #html related graphics files
```

Installation
-----

### Install Conda (if not already installed) 
It is recommended to use:
```
Miniconda / Anaconda
```
### Create Python environment

```
conda create -n reseq python=3.9
conda activate RESeq
```
Install dependencies:
```
conda install -y pandas matplotlib jinja2
```

### Set up the R environment (critical step)
Create environment.yml:
```
name: r_environment
channels:
  - conda-forge
  - defaults
dependencies:
  - r-base=4.2.3
  - bioconda::bioconductor-annotationdbi
  - r-corrplot
  - r-coin
  - bioconda::bioconductor-complexheatmap
  - bioconda::bioconductor-clusterprofiler
  - r-dplyr
  - r-data.table
  - r-foreach
  - bioconda::bioconductor-deseq2
  - r-doparallel
  - bioconda::bioconductor-edger
  - bioconda::bioconductor-enrichplot
  - r-extrafont
  - r-factoextra
  - r-foreach
  - r-fs
  - r-glue
  - conda-forge::r-gridextra
  - r-gdata
  - r-ggplot2
  - r-ggcorrplot
  - r-ggextra
  - r-ggrepel
  - r-igraph
  - bioconda::bioconductor-limma
  - bioconductor-org.hs.eg.db
  - bioconductor-org.mm.eg.db
  - bioconductor-org.rn.eg.db
  - r-rmisc
  - r-rcolorbrewer
  - r-pheatmap
  - bioconda::bioconductor-biocparallel
  - r-pbapply
  - bioconda::bioconductor-pathview
  - r-patchwork
  - r-statmod
  - r-stringr
  - r-tibble
  - r-tidyr
  - r-tidyverse
  - r-reshape2
  - fonts-anaconda
  - pip:
    - custom-functions
```
Install:
```
conda env create -f environment.yml
conda activate r_environment
```

### Install bioinformatics software
```
conda install bioconda::fastp
conda install bioconda::fastqc
conda install bioconda::gatk4
conda install conda-forge::perl
conda install bioconda::hisat2
conda install bioconda::samtools
conda install bioconda::tabix
conda install bioconda::htseq
```
### Install Annovar
Official download： 
```
[Annovar](https://annovar.openbioinformatics.org/en/latest/user-guide/download/)
```
Then configure the path.
### Modify the flow configuration file
Edit：
```
REseq/config.py
```
Set the reference genome and the paths of related software:
```
ref_path = {"hg38": ".../Refer/hg38/hg38",
                        "mm10": ".../Refer/mm10/mm10",
                        "rn6":  ".../Refer/rn6/rn6"
}
ref_fasta = {"hg38": ".../Refer/hg38/hg38.fa",
                        "mm10": ".../Refer/mm10/mm10.fa",
                        "rn6":  ".../Refer/rn6/rn6.fa"
}
...
```

### Input Data
A sample.xls file is required：
```
SampleID   R1.fastq.gz   R2.fastq.gz   Group
S1         /path/1.fq.gz /path/2.fq.gz Tumor
S2         ...           ...           Tumor
S3         ...           ...           Normal
```
Usage
-----
### Help information
You can enter a command to get the help message, which contains information about all the optional parameters.

```python
python all.py -h 
```

### Example
You need to edit a "sample.xls" file with four columns of content, the first column is the sample ID, the second column is the absolute path of the sample reads R1 file, the third column is the absolute path of the sample reads R2 file, and the fourth column is the sample grouping information.
"sample.xls" file example:
```
Rop-1   /home/data/jc1/ZQ-176/CG-1_L1_1.fq.gz   /home/data/jc1/ZQ-176/CG-1_L1_2.fq.gz   Rop
Rop-2   /home/data/jc1/ZQ-176/CG-2_L1_1.fq.gz   /home/data/jc1/ZQ-176/CG-2_L1_2.fq.gz   Rop
Rop-3   /home/data/jc1/ZQ-176/CG-3_L1_1.fq.gz   /home/data/jc1/ZQ-176/CG-3_L1_2.fq.gz   Rop
Rop_Dex-1       /home/data/jc1/ZQ-176/EG-1_L1_1.fq.gz   /home/data/jc1/ZQ-176/EG-1_L1_2.fq.gz   Rop_Dex
Rop_Dex-2       /home/data/jc1/ZQ-176/EG-2_L1_1.fq.gz   /home/data/jc1/ZQ-176/EG-2_L1_2.fq.gz   Rop_Dex
Rop_Dex-3       /home/data/jc1/ZQ-176/EG-3_L1_1.fq.gz   /home/data/jc1/ZQ-176/EG-3_L1_2.fq.gz   Rop_Dex
Control-1       /home/data/jc1/ZQ-176/NG-1_L1_1.fq.gz   /home/data/jc1/ZQ-176/NG-1_L1_2.fq.gz   Control
Control-2       /home/data/jc1/ZQ-176/NG-2_L1_1.fq.gz   /home/data/jc1/ZQ-176/NG-2_L1_2.fq.gz   Control
Control-3       /home/data/jc1/ZQ-176/NG-3_L1_1.fq.gz   /home/data/jc1/ZQ-176/NG-3_L1_2.fq.gz   Control
```
You can then perform RE-Seq analysis with the following command.
```
python you_path.all.py --work_dir ./01_work --out_dir ./02_Result -ref rn6 -i sample.xls -cn ZQ-176 -pn "Rattus norvegicus 9 RE-seq"
```

## Status tracking and restart

The system displays the task progress and whether the task is running successfully. If the workflow interruption needs to be restarted, the system will automatically check the file, skip the completed task, and continue to execute from the breakpoint, improving work efficiency.
