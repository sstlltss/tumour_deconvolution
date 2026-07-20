# coding: utf-8
import pandas as pd
from pathlib import Path

files = []
with open("files.txt", "r") as f:
	for line in f:
		files.append(line.rstrip())

for cov in files:
	rn = cov.split("_1_bismark_bt2_pe.bismark.cov.gz")[0]
	Path(rn).mkdir(exist_ok=True)
	df = pd.read_csv(cov, compression="gzip", header=None, sep="\t")
	uchr = df.iloc[:,0].unique()
	for c in uchr:
		if c[-1] is "MXY":
			continue
		ndf = df[df[0] == c]
		ndf.to_csv(f"{rn}/{rn}_{c}.csv.gz",sep='\t',header=False,index=False,compression='gzip')
