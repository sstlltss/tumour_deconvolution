import pandas as pd
import os

SITE_THERSHOLD_RATE = 0.5
SAMPLE_THERSHOLD = 20000
INPUT_PATH = "/home/wuyuwei/nobackup/other/mydata/cancer/cov/"
OUTPUT_PATH = "/home/wuyuwei/nobackup/MethyLYZR/training/prepared_data/"
LABEL_FILE = "train_labels.txt"
TEST_MODE = False

ann = pd.read_csv(
    "/home/wuyuwei/nobackup/other/mydata/450k_annotation_hg19.csv",
    sep=",",
    header=0,
    usecols=["chr","pos","Name"]
)
ann["key"] = (ann["chr"].astype(str) + "_" + ann["pos"].astype(str))
probe_map = dict(zip(ann["key"], ann["Name"]))

files = pd.read_csv(INPUT_PATH+LABEL_FILE, sep=",", header=None)
files.columns = ["sample","cancer"]
cancers = files.iloc[:,1].unique()
sample_df = pd.DataFrame()
counts_sum = pd.DataFrame()

for i,cancer in enumerate(cancers):
    samples = files[files["cancer"]==cancer]["sample"].tolist()
    for sample in samples:
        chunks = []
        counts = []
        sample_path = os.path.join(INPUT_PATH, sample+".cov.gz")
        print(f"Now working on: {sample_path}")
        for chunk in pd.read_csv(sample_path, sep="\t", header=None,chunksize=500000):
            chunk.columns = ["chr","start","end","methylated_rate","methylated_count","unmethylated_count"]
            chunk["pos"] = chunk["start"] - 1
            chunk["key"] = chunk["chr"].astype(str) + "_" + chunk["pos"].astype(str)
            chunk = chunk[chunk["key"].isin(probe_map)]
            chunk["cpg"] = chunk["key"].map(probe_map)
            chunk["count"] = chunk["methylated_count"] + chunk["unmethylated_count"]
            chunk["methylated_rate"] = chunk["methylated_rate"] / 100.0
            chunks.append(chunk.loc[:,["cpg","methylated_rate"]])
            counts.append(chunk.loc[:,["cpg","count"]])
        chunks = pd.concat(chunks)
        counts = pd.concat(counts)
        
        chunks = chunks.set_index("cpg")["methylated_rate"]
        chunks = chunks.rename(sample)
        sample_df = sample_df.join(chunks, how="outer")
        print(f"Sample {sample} is filtered. Matched probes: {chunks.shape[0]}")

if TEST_MODE:
    print("--------Test mode--------")
    for col in sample_df.columns:
        print(f"Now working on sample: {col}")
        single_sample = sample_df[col]
        single_sample = single_sample.dropna()
        single_sample = single_sample.to_frame(name="methylation")
        single_sample = single_sample.reset_index(names=["epic_id"])
        single_sample.to_feather(os.path.join(OUTPUT_PATH, f"{col}_betas.feather"))
    files.to_csv(os.path.join(OUTPUT_PATH, "ClassID.csv"), sep=',', index=False, header=['sample','methylation_class'])
    print("Completed.")
    exit()

sample_cov = sample_df.notna().sum(axis=0)
keep_sample = sample_cov >= SAMPLE_THERSHOLD
sample_df = sample_df.loc[:,keep_sample]    # samples match too few probes will be discarded
filtered_sample = sample_df.columns
filtered_files = files.loc[files["sample"].isin(filtered_sample),:]
filtered_files.to_csv(os.path.join(OUTPUT_PATH,"ClassID.csv"), sep=',', index=False, header=['sample','methylation_class'])
print(f"Samples after filtering: {filtered_sample}")

site_cov = sample_df.notna().sum(axis=1)
keep = site_cov >= SITE_THERSHOLD_RATE * files.shape[0]
sample_df = sample_df.loc[keep,:]   # probes appear in too few samples will be discarded
sample_df.to_csv(os.path.join(OUTPUT_PATH,"merged_betas_with_NA.csv.gz"), sep=",",compression="gzip",index = False)
print(f"Probes number after filtering: {sample_df.index.nunique()}")

# The presence of CpG sites is highly cancer-specific. Therefore, missing values are imputed using the mean value of other probes within the same sample rather than the mean of the same probe across samples of the same cancer type.
sample_df = sample_df.fillna(sample_df.mean())
sample_df.to_feather(os.path.join(OUTPUT_PATH,"merged_betas.feather"))
print(f"Mats:{sample_df}")
#sample_df = sample_df.reset_index()
#sample_df.to_csv(os.path.join(OUTPUT_PATH,"merged_betas.csv.gz"), sep=",",compression="gzip",index = False)

print("Completed.")
