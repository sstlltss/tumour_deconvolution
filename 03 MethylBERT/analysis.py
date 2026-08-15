import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report
import json

seq_len=150
n_mers=3
batch_size=16
num_workers=0
enable_dmr = True
lr=2e-5
warmup_step=200
steps=3000
seed = 88

output_path = "/home/wuyuwei/nobackup/other/methylbert/test11/"
"""with open(output_path+"log.json", "r") as f:
    log = json.load(f)"""

df_train  = pd.read_csv(output_path+"train.csv", sep="\t")
df_eval  = pd.read_csv(output_path+"eval.csv", sep="\t")
sns.lineplot(data=df_train, x="step", y="loss", label="train loss")
sns.lineplot(data=df_eval, x="step", y="loss", label="eval loss")
#plt.title(f'warm up: {log["warmup_step"]}, DMR: {log["enable_dmr"]}, seed: {log["seed"]}')
plt.title(f'warm up: {warmup_step}, DMR: {enable_dmr}, seed: {seed}')
sns.lineplot(data=df_eval, x="step", y="ctype_acc", label="Accuracy")
plt.savefig(output_path+"acc.jpg")
