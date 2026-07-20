# Tumor Deconvolution Using Methylation Data
### Note: This project was developed and tested on Ubuntu Linux. The shell scripts assume an Ubuntu environment, and some commands or package installation procedures may differ on other operating systems.

## Introduction

## Data Processing
**Scripts using in this part are stored at `01 data_processing`.**

### Data collecting
The raw sequencing and DNA methylation data used in this study consist of four datasets adopted from MethylBERT. They are downloaded from [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) with the given accession number:

- DLBCL cancer: PRJNA565006
- gastrointestinal cancer: PRJNA628686
- pancreatic cancer: PRJNA266735
- B-cell: PRJNA574550

The original sequencing data are provided in FASTQ (.fastq.gz) format.

In addition, a reference genome is required for read alignment. In this study, the human reference genome **GRCh37 (hg19)** was used. The required reference genome files in FASTA (.fa.gz) format can be downloaded from [National Center of Biotechnology Information (NCBI)](https://ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.22/).

### Environment
For data processing, it is highly recommended to use a virtual environment to manage the required software and dependencies. This project uses Miniconda for creating virtual environments and managing package installations. The following steps describe how to set up the environment for data processing:

1. Install [conda](https://www.anaconda.com/docs/getting-started/installation) and activate it as `base`.
2. Install the environment with `environment.yml` by running the following command:
```bash
conda install -f environment.yml
```
3. Activate the environment:
```bash
conda activate dataprocessing
```

### Alignment, Deduplication and Extraction
The input files required by MethylBERT include BAM files, methylation coverage files (COV files), and DMR result files (CSV files). MethyLYZR requires methylation coverage files to generate the input feature matrix. The following steps describe how to generate these required files.

1. Prepare reference genome:
```bash
bismark_genome_preparation --verbose /path/to/genome/folder
```
2. Align FASTQ files to reference genome. For paired-end sequencing, two FASTQ files (e.g., Sample_1.fastq and Sample_2.fastq) should be supplied as input using the -1 and -2 options:
```bash
cd /your/fastq/file/folder
for f1 in *_1.fastq.gz; do
    f2="${f1%_1.fastq.gz}_2.fastq.gz"
    bismark --genome /your/ref_genome/folder -1 $f1 -2 $f2 -o /folder/for/bam/files
done
```
3. After alignment, deduplication is required before methylation extraction to remove PCR duplicate read pairs:
```bash
cd /your/bam/file/folder
for f in *_1_bismark_bt2_pe.bam; do
    nf="${f%_1_bismark_bt2_pe.bam}_namesorted.bam"
    samtools sort -n -o $nf $f
    deduplicate_bismark -p $nf
done
```
4. Extrace methylation information:
```bash
cd /your/deduplicated/bam/folder
for f in *_namesorted.deduplicated.bam; do
    rn="${f%_namesorted.deduplicated.bam}"
    if [ ! -e extraction ]; do
        mkdir extraction
    fi
    mkdir /extraction/$raw_name
    bismark_methylation_extractor -p --gzip --bedGraph -o extraction/$rn $f
done
```
5. Clean up intermediate files. After methylation extraction, a subdirectory is created for each sample under `extraction/`. Several output files are generated for each sample, including:
- CpG context (TXT file)
- CpG report (TXT file)
- CHG context (TXT file)
- CHG report (TXT file)
- CHH context (TXT file)
- CHH report (TXT file)
- M-bias report (TXT file)
- Splitting report (TXT file)
- Methylation rate report (BEDGRAPH file)
- Methylation coverage report (COV file)

Among these files, the COV file is required for downstream DMR analysis (using DSS) and can also be converted into the feature matrix required by MethyLYZR. The remaining files are not used in the subsequent steps of this pipeline and can be safely removed to save disk space.

To reduce memory consumption and avoid out-of-memory (OOM) errors during DMR analysis, each COV file should split into separate files by chromosome. Chromosomes X, Y, and M are excluded from the analysis. Copy `split_chr.py` to `extraction` folder and run the following shell:
```bash
cd extraction
rm -f */*.txt*
mkdir cov
mkdir bedGraph
mv */*.cov* cov
mv */*.bedGraph* bedGraph
ls cov/*.cov.gz > files.txt

conda activate dataprocessing
python split_chr.py
conda deactivate
```