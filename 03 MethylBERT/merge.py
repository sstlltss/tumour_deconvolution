import pandas as pd
ctypes = ["DLBCL", "GC", "PC"]

for c in ctypes:
    print(f"===={c}====")
    for n in range(1,22):
        file = f"{c}_chr{n}.csv"
        dmr_single = pd.read_csv(file)
        dmr_single["ctype"] = c
        print(f"Adding chr{n} to {file}...")
        dmr_single.to_csv(f"{c}.csv", header=bool(n==1), index=False, mode="a", sep="\t")

print("Completed.")
