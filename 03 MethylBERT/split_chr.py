# coding: utf-8
import pandas as pd
from pathlib import Path

files = []
with open("files.txt", "r") as f:
	for line in f:
		files.append(line.rstrip())

chrs = [f"chr{i}" for i in range(1,22)]
for cov in files:
	print(f"========{cov}========")
	rn = cov.split(".cov.gz")[0]
	df = pd.read_csv(cov, compression="gzip", header=None, sep="\t")
	for c in chrs:
		print(f"Now working on: {c}.")
		ndf = df[df[0] == c]
		ndf.to_csv(f"{rn}_{c}.csv.gz",sep='\t',header=False,index=False,compression='gzip')
