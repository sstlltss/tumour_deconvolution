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
    "PC": dmr_folder_path+"PC.csv",
    "GC": dmr_folder_path+"GC.csv"
    }
f_ref = "path/to/reference/genome/hg19.fa"
out_dir = "path/to/output/folder"
suffix_bam = "_1_deduplicated.bam"

fdg.finetune_data_generate(
    sc_dataset = bam_folder_path+f_bam_file_list,
    f_dmr = f_dmr,
    f_ref = f_ref,
    output_dir=out_dir,
    split_ratio=1, # Split ratio to make training and validation data
    n_mers=3, # 3-mer DNA sequences 
    n_cores=20,
    use_file_name = True,
    use_existed_files = False,
    keep_rate = {"DLBCL": 0.1, "GC": 0.2, "PC": 0.5}
)
```
The `keep_rate` parameter requires a dictionary where keys represent cell types and values ​​are floating-point numbers indicating the probability of retaining reads for that cell type. This parameter helps balance the number of extracted reads when there are significant differences in sample number and/or sample size across cell types.

The `suffix_bam` Specifies the common suffix of the input BAM files. All input BAM files must use the same suffix. For example, if the input files are named:
```
SRR001_deduplicated.bam
SRR002_deduplicated.bam
SRR003_deduplicated.bam
```
the suffix_bam parameter should be set to:
```
_deduplicated.bam
```
The suffix is removed from each BAM filename to identify the corresponding sample name.

The variable `f_bam_file_list` point to the file that list the BAM files using in the fine-tuning, while `bam_folder_path` and `dmr_folder_path` are the path to BAM and DMR folders. `f_dmr` is a Dictionary, with keys the cancer type of the DMR file, and values the path to the files. `f_ref` specifies the path to reference genome file and `out_dir` the path to a folder to store the output file(s).

When `split_ratio` and `train_valid_test_ratio` both are not given, or `split_ratio = 1`, or `train_valid_test_ratio = [1,0,0]`, the output file will be a single file `data.csv`, otherwise the output files will be three files named with `train_seq.csv`, `test_seq.csv` and `val_seq.csv`. The file(s) content(s) the information extracted from BAM files and from the regions specified by DMR reports.

## Fine-tuning and evaluation
Once the fine-tuning data are prepared, you might run `finetuning.py` to start model fine-tuning phase. You might run the script with arguments, for example:
```
python finetuning.py --eval_freq 300
```
For more detail of arguments, please run `python finetuning.py --help`.

## Output
During the fine-tuning phase, the model will be trained and evaluated. Several files will be saved to the output directory you specified. They are:
- **acc.jpg**: A plot contents train and evaluate loss as well as evaluation accuracy during the training. The x-axis is training steps.
- **config.json**: The config file of the model.
- **confusion_matrix_step_*step*.jpg**: The classification confusion matrix in **read level** over evaluation steps.
- **sample_confusion_matrix_step_*step*.jpg**: The classification confusion matrix in **sample level** over evaluation steps.
- **eval.csv & train.csv**: A CSV table contents evalutaion & train loss and accuracy over evaluation steps.
- **predictions_step_*step*.csv**: The raw prediction results over evaluation steps.
- **report_step_*step*.csv**: The classification reports over evaluation steps.
- **roc_step_*step*.csv**: The ROC curve over evaluation steps.
- **dmr_encoder.pickle, model.safetensors** and **read_classification_model.pickle**: Saved model information.