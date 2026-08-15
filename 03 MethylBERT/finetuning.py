from methylbert.utils import set_seed
from torch.utils.data import DataLoader
from methylbert.data.vocab import MethylVocab
from methylbert.data.dataset import MethylBertFinetuneDataset
from methylbert.trainer import MethylBertFinetuneTrainerWithClassifier
import os
import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import torch
import sys
import argparse
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, classification_report
warnings.filterwarnings("ignore") # Ignore warnings for a clear notebook
torch.multiprocessing.set_sharing_strategy('file_system')

train_data_loader = None
test_data_loader = None

parser = argparse.ArgumentParser(
    usage="%(prog)s [options]",
    description=(
        "Fine-tune and evaluate the MethylBERT model for DNA methylation classification."
    ),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("--seq_len", type=int, default=150, help="Length of the DNA sequence used as model input.")
parser.add_argument("--n_mers", type=int, default=3, help="Size of k-mers used for DNA sequence tokenization.")
parser.add_argument("--batch_size", type=int, default=16, help="Number of samples processed in each training batch.")
parser.add_argument("--num_workers", type=int, default=0, help="Number of worker processes used by the DataLoader.")
parser.add_argument("--warmup_step", type=int, default=50, help="Number of warm-up training steps for the learning rate scheduler.")
parser.add_argument("--steps", type=int, default=6000, help="Total number of training steps.")
parser.add_argument("--lr", type=float, default=2e-5, help="Learning rate used for model training.")
parser.add_argument("--fold", type=int, default=1, help="Cross-validation fold to use for training and evaluation.")
parser.add_argument("--dmr", action=argparse.BooleanOptionalAction, default=True, help="Enable the use of DMR information during model training.")
parser.add_argument("--output_path", type=str, required=True, help="Directory where the model checkpoints and evaluation results will be saved.")
parser.add_argument("--seed", type=int, default=42, help="Random seed used to ensure reproducibility.")
parser.add_argument("--eval_freq", type=int, default=100, help="Evaluate the model every N training steps.")
args = parser.parse_args()

seq_len = args.seq_len
n_mers = args.n_mers
batch_size = args.batch_size
num_workers = args.num_workers
warmup_step = args.warmup_step
steps = args.steps
lr = args.lr
fold = args.fold
dmr = args.dmr
output_path = args.output_path
seed = args.seed
eval_freq = args.eval_freq

set_seed(seed)
log = {"seed": seed,
       "seq_len": seq_len,
       "n_mers": n_mers,
       "batch_size": batch_size,
       "num_workers": num_workers,
       "output_path": output_path,
       "enable_dmr": dmr,
       "lr": lr,
       "warmup_step": warmup_step,
       "steps": steps,
       "fold": fold,
       "eval_freq": eval_freq
       }
with open("log.json", "w") as f:
    f.write(json.dumps(log))

# Creat a look-up table
tokenizer = MethylVocab(n_mers)

# Read number of classes(labels) from labels.csv, and generate "label2id" & "id2label" for config
df_dmr = pd.read_csv("tmp/dmrs.txt", header=None)
if df_dmr.shape[0] > 0:    
    id2label = df_dmr.to_dict()[0]
    label2id = dict([(value, key) for key, value in id2label.items()])
    print(f"DMRs are: {id2label}")
else:
    raise ValueError('Can\'t find any DMRs. Please check "tmp/dmrs.txt"!')

# Load the data files int a data set object
train_data_path = f"data{fold}_train.csv"
test_data_path = f"data{fold}_test.csv"
train_dataset = MethylBertFinetuneDataset(train_data_path, 
                                          tokenizer, 
                                          seq_len=seq_len,
                                          id2label=id2label,
                                          label2id=label2id)
test_dataset = MethylBertFinetuneDataset(test_data_path, 
                                         tokenizer,
                                         seq_len=seq_len,
                                         id2label=id2label,
                                         label2id=label2id) 
# Load the data into a data loader
train_data_loader = DataLoader(train_dataset, batch_size= batch_size, 
                               num_workers= num_workers, pin_memory=False,  
                               shuffle=True)
test_data_loader = DataLoader(test_dataset, batch_size= batch_size, 
                              num_workers= num_workers, pin_memory=False,  
                              shuffle=False)

trainer = MethylBertFinetuneTrainerWithClassifier(
                      len(tokenizer), 
                      save_path=output_path, 
                      train_dataloader=train_data_loader, 
                      test_dataloader=test_data_loader,
                      id2label=id2label,
                      label2id=label2id,
                      enable_dmr=dmr,
                      lr=lr,
                      with_cuda=True, 
                      log_freq=1,
                      #eval_freq=10, #activate this only when you want to evaluate the model with test_data_loader
                      warmup_step=warmup_step,
                      loss="cross_entropy",
                      ignore_mismatched_sizes=True,
                      eval_freq=eval_freq)

# model fine-tuning
trainer.load("hanyangii/methylbert_hg19_4l")
trainer.train(steps=steps)

# result analysis
df_train  = pd.read_csv(output_path+"train.csv", sep="\t")
df_train.head()
df_eval  = pd.read_csv(output_path+"eval.csv", sep="\t")
df_eval.head()
sns.lineplot(data=df_train, x="step", y="loss", label="train loss")
sns.lineplot(data=df_eval, x="step", y="loss", label="eval loss")
plt.title(f"warm up: {warmup_step}, DMR: {dmr}, seed: {seed}")
sns.lineplot(data=df_eval, x="step", y="ctype_acc", label="Accuracy")
plt.savefig(output_path+"acc.jpg")
plt.close()

for i in range(0,steps, eval_freq):
    step = i if i==0 else i-1
    df_pred = pd.read_csv(output_path+f"predictions_step_{step}.csv")
    classes = list(label2id.keys())
    sample_prob = (df_pred.groupby("sample")[classes].mean())
    sample_true = (df_pred.groupby("sample")["true_label"].first())
    sample_pred = sample_prob.idxmax(axis=1)
    report = classification_report(y_true=sample_true, y_pred=sample_pred, labels=classes)
    with open(output_path+f"report_sample_{step}.csv", "w") as f:
        f.write(report)

    cm_labels = sorted(set(sample_true))
    cm = confusion_matrix(y_pred=sample_pred, y_true=sample_true)
    disp = ConfusionMatrixDisplay(cm, display_labels=cm_labels)
    disp.plot(cmap="Blues")
    plt.tight_layout()
    plt.savefig(os.path.join(output_path+f"sample_confusion_matrix_step_{i-1}.png"))
    plt.close()

