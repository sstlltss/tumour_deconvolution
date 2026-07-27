# MethylBERT

## Environment
To set up and activate the environment for data preparation and MethylBERT, you might run the following command in `bash`:
```bash
conda install -f methylbert_environment.yml -y
```
DMR analysis is implemented using DSS, an R package for differential methylation analysis. Before running the analysis, create an R environment and install the required dependencies. The following Bash commands create the environment and install the necessary packages:
```bash
conda install -f dss_environment.yml -y
```

## DMR analysis
To reduce memory consumption and avoid out-of-memory (OOM) errors during DMR analysis, each COV file should split into separate files by chromosome. Chromosomes X, Y, and M are excluded from the analysis. Copy `split_chr.py` to `extraction` folder and run the following shell:
```bash
cd extraction
mkdir ../cov
mkdir ../bedGraph
mv */*.cov.gz ../cov
mv */*.bedGraph.gz ../bedGraph
cd ../cov
for f in *_namesorted.deduplicated.bismark.cov.gz; do
  rn="${f%_namesorted.deduplicated.bismark.cov.gz}.cov.gz"
  mv $f $rn
done
ls *.cov* > cov_files.txt

conda activate dataprocessing
python split_chr.py
conda deactivate
```
The splited files are named as `CANCER_chrN.csv.gz`.

List all the samples and their cell type you want to use for training in `train_labels.txt`. You can also use the file mentioned in `tutorial_MethyLYZR.md` since they are in the same format:
```
SRR10099840,DLBCL
SRR10099825,DLBCL
SRR11615800,GC
SRR11615880,GC
SRR3427332,PC
SRR3427335,PC
SRR10165874,normal
SRR10165875,normal
```
The R-script `dss_dmr.R` is for DMR analysis. The line `ctypes <- list("DLBCL","GC","PC")` specified the cancer type you want to analysis. The script will detect the normal cells automatical. Start DMR analysis with the following command:
```
conda activate DSS
Rscript --no-save dss_dmr.R
conda deactivate
```
The output files are named with the cancer name, each content the comparasion result between the cancer and normal cells. Run `merge.py` to merge them up and add the `ctype` column as the requirment of fine-tuning:
```
python merge.py
```
The `ctypes` list in `merge.py` is also to specified to the cancer types your implement.

The output files are the merged DMR report with the name `CANCER.csv`.

## Fine-tuning data preparation
Run `prefinetuning.py` to prepare the dataset for fine-tuning:
```python
from methylbert.data import finetune_data_generate as fdg

bam_folder_path = "path/to/bam/files"
f_bam_file_list = "train_labels_bam.txt"
dmr_folder_path = "path/to/dmr/files"
f_dmr = {
    "DLBCL": dmr_folder_path+"DLBCL.csv",
    "CC": dmr_folder_path+"CC.csv",
    "GC": dmr_folder_path+"GC.csv"
    }
f_ref = "path/to/reference/genome/hg19.fa"
out_dir = "path/to/output/folder"

fdg.finetune_data_generate(
    sc_dataset = bam_folder_path+f_bam_file_list,
    f_dmr = f_dmr,
    f_ref = f_ref,
    output_dir=out_dir,
    split_ratio=1, # Split ratio to make training and validation data
    n_mers=3, # 3-mer DNA sequences 
    n_cores=20,
    use_file_name = True,
    use_existed_files = False
)
```
The variable `f_bam_file_list` point to the file that list the BAM files using in the fine-tuning, while `bam_folder_path` and `dmr_folder_path` are the path to BAM and DMR folders. `f_dmr` is a Dictionary, with keys the cancer type of the DMR file, and values the path to the files. `f_ref` specifies the path to reference genome file and `out_dir` the path to a folder to store the output `data.csv`.

The output file `data.csv` contents the information extracted from BAM files and from the regions specified by DMR reports.

##