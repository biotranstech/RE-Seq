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

### Running Analysis

```
python all.py \
  --work_dir ./work \
  --out_dir ./result \
  -ref hg38 \
  -i sample.xls \
  -cn PROJECT001 \
  -pn "RNA editing project"
```
### Output
After the operation is completed:
```
02_result/
├── 01_QC
├── 03_RESeq_file
└── 04_Down_analysis
03_report/test_Report
├── asset
├── images
├── Result
└── report_noAS_v1.2.html
```
Directly open:
```
report_noAS_v1.2.html
```
You can then view the complete analysis report.

### Resume
If interrupted:
```
Just re-run the same command.
```
The system will:
1. Automatically skip completed steps
2. Resume from the last checkpoint

### Common Issues
1. Software path error
Must modify:
```
REseq/config.py
```
2. Reference not configured
Check:
```
fasta
index
```
3. Permission issues
```
chmod -R 755 your_project
```
4. Environment conflicts
Recommended:
```
conda clean -a
```
## Workflow Logic (Core Principle)
Actual RE-Seq workflow:
```
FASTQ
 ↓
Quality Control (fastp)
 ↓
Alignment (hisat2)
 ↓
BAM processing
 ↓
RNA editing detection
 ↓
Filtering & annotation
 ↓
Statistical analysis
 ↓
Visualization
 ↓
HTML report
```
