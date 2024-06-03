#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd


from hplus_utils_6j3b import edges6j3b
from hplus_utils_6j3b import getHplusDf
from utils_hplus      import masspoints



parser = argparse.ArgumentParser()
parser.add_argument("txtfile")
#parser.add_argument("-m", "--mass", type=int, required=True)
parser.add_argument("-o", "--outfile", required=True)


args = parser.parse_args()


with open(args.txtfile) as f:
	trees = f.read().split(',')


ntuplename = 'tree'

dfs = []
for tree in trees:
	dfs.append(getHplusDf(ntuplename,tree))
df = pd.concat(dfs)
del dfs


print(df)

# Yield tables
print("Yield table")
print(
    df.groupby("sample")["weight"].agg(
        Entries="count",
        Integral="sum")
)


# Get bin edges
# Discretize MVA scores according to binning
for mass in masspoints:
	edges = None
	edges = edges6j3b[mass]
	idx = np.digitize(df[f"scores_{mass}"], bins=edges)

	# Fits uint8
	assert np.all((idx >= 0) & (idx <= 255))
	#None in underflow
	assert np.all(idx != 0)
	#None in overflow
	assert np.all(idx != len(edges))

	df[f"scores_{mass}Bin"] = idx.astype(np.uint8)

df.to_hdf(args.outfile, "df_hplus_6j3b", complevel=9, format="table")

print(df)


