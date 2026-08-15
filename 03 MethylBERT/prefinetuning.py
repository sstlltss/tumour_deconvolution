import pandas as pd
from methylbert.data import finetune_data_generate_random_keep as fdg

bam_folder_path = "/home/wuyuwei/nobackup/other/mydata/cancer/bam/"
f_bam_file_list = "labels.txt"
dmr_folder_path = "/home/wuyuwei/nobackup/other/mydata/dmrs/"
f_dmr = {
    "DLBCL": dmr_folder_path+"intersect_dmrs.csv",
    "PC": dmr_folder_path+"intersect_dmrs.csv",
    "GC": dmr_folder_path+"intersect_dmrs.csv",
    "normal": dmr_folder_path+"intersect_dmrs.csv"
    }
f_ref = "/home/wuyuwei/nobackup/other/hg19/hg19.fa"
out_dir = "/home/wuyuwei/nobackup/other/methylbert/tmp/"

fdg.finetune_data_generate(
    sc_dataset = bam_folder_path+f_bam_file_list,
    f_dmr = f_dmr,
    f_ref = f_ref,
    output_dir=out_dir,
    split_ratio=1, # Split ratio to make training and validation data
    n_mers=3, # 3-mer DNA sequences 
    n_cores=1,
    n_dmrs=10000,
    use_file_name = True,
    use_existed_files = True,
    keep_rate = {"DLBCL": 0.12, "GC": 0.015, "PC": 1, "normal": 1}
)

all_data_df = pd.read_csv(out_dir+"data.csv",sep="\t")
n_folds = 3
for i in range(1,n_folds+1):
    for t in ["test","train"]:
        is_first_write = True
        labels = f"tmp/fold{i}_{t}.csv"
        labels_df = pd.read_csv(labels,sep="\t")
        for label in labels_df.itertuples(index=False):
            print(f"Now working on: {label.sample}_1_deduplicated.bam")
            current_reads = all_data_df[all_data_df["filename"]==label.sample+"_1_deduplicated.bam"]
            current_reads.to_csv(f"data{i}_{t}.csv",sep="\t",header=is_first_write,index=False,mode="a")
            is_first_write = False