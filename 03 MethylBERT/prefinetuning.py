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
