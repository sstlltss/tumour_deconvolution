import pandas as pd
from sklearn.metrics import classification_report

step = 199
df_pred = pd.read_csv(f"predictions_step_{step}.csv")
classes = ["DLBCL","GC","PC","normal"]
sample_prob = (df_pred.groupby("sample")[classes].mean())
sample_true = (df_pred.groupby("sample")["true_label"].first())
sample_pred = sample_prob.idxmax(axis=1)
report = classification_report(y_true=sample_true, y_pred=sample_pred, labels=classes)
with open(f"report_sample_{step}.txt", "w") as f:
    f.write(report)