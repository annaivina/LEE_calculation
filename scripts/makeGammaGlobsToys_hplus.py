#!/usr/bin/env python3
from scipy import stats
from tqdm import tqdm
import argparse
import numpy as np
import os
import pandas as pd

from utils_hplus import masspoints

import ROOT as R
R.gROOT.SetBatch(True)


parser = argparse.ArgumentParser()
parser.add_argument("dataframe")
parser.add_argument("asimov")
parser.add_argument("-c", "--channel", choices=["4j3b", "4j4b", "5j3b", "5j4b", "6j3b", "6j4b"], required=True)
parser.add_argument("-o", "--outdir", default="")
args = parser.parse_args()




# Use different seed for different channels for reproducibility and
# independent RVS
if args.channel == "4j3b":
    seed = 9124568754
elif args.channel == "4j4b":
    seed = 4324357532
elif args.channel == "5j3b":
    seed = 2947673012
elif args.channel == "5j4b":
    seed = 8322305006
elif args.channel == "6j3b":
    seed = 1893475993
elif args.channel == "6j4b":
    seed = 3439078643

rng = np.random.default_rng(seed)


df = pd.read_hdf(args.dataframe)
df["weight"] = df["weight"].astype(np.float64)
df["weightSquared"] = df["weight"]**2


df["sample"] = df["sample"].astype('category')
#print(df["sample"].cat.categories)

#df = df.loc[df["sample"] != "Data"]
#df["sample"] = df["sample"].cat.remove_categories(["Data"])

df = df.loc[df["sample"] != "Data"]
if "Data" in df["sample"].cat.categories:
    df["sample"] = df["sample"].cat.remove_categories(["Data"])

# No scaling of Z+HF and ttbar intentional to be consistent with
# treatment in workspaces

# Sum of weights / weights^2 from ntuple
sumw = {}
sumw2 = {}
for mass in masspoints:
    mva = df[f"scores_{mass}Bin"]

    # MVA bin starts counting at 1
    edges = [0] + sorted(mva.unique())
    hsumw, _ = np.histogram(mva - 1, bins=edges, weights=df["weight"])
    hsumw2, _ = np.histogram(mva - 1, bins=edges, weights=df["weightSquared"])

    sumw[mass] = hsumw
    sumw2[mass] = hsumw2

# Calculate tau using ntuple
tau = {}
for mass in masspoints:
    tau[mass] = sumw[mass]**2 / sumw2[mass]

# Compare with tau from workspace
channel_str = None
if args.channel == "4j3b":
    channel_str = "4j3b"
elif args.channel == "4j4b":
    channel_str = "4j4b"
elif args.channel == "5j3b":
    channel_str = "5j3b"
elif args.channel == "5j4b":
    channel_str = "5j4b"
elif args.channel == "6j3b":
    channel_str = "6j3b"
elif args.channel == "6j4b":
    channel_str = "6j4b"
else:
    raise RuntimeError("Unknown channel")


tau_ws = {}
fin = R.TFile.Open(args.asimov)
for mass in masspoints:
	h = fin.Get(f"tau_{channel_str}_m{mass}")
	nbins = h.GetNbinsX()

	assert len(tau[mass]) == nbins

	tau_ws[mass] = np.array([h.GetBinContent(i) for i in range(1, nbins + 1)])
	assert len(tau_ws[mass]) == nbins

	if np.any(np.abs(tau_ws[mass] / tau[mass] - 1) > 5e-2):
		print(f'Possibly problematic histograms: mH+ with {mass} GeV')
		print(f'Relative deviation of WS tau and tau from ntuple:\n{tau_ws[mass] / tau[mass] - 1}\n')

fin.Close()

# Scaling factor
sf = {}
for mass in masspoints:
    sf[mass] = tau_ws[mass] / sumw[mass]

# Calculate random global observables
# This is relatively slow but should be good enough for now
globs = {}
for i in tqdm(range(20000)):
    pois_weight = stats.poisson.rvs(mu=1, size=len(df), random_state=rng)
    df["toy_weight"] = pois_weight * df["weight"]

    for mass in masspoints:
        hist = df.groupby(f"scores_{mass}Bin")["toy_weight"].sum()

        # Make sure that the bins are sorted properly
        hist.sort_index(inplace=True)

        toy_globs = hist.values * sf[mass]
        globs.setdefault(mass, []).append((i, toy_globs))

# Sanity checks
for mass in globs:
    arr = np.array([x for i, x in globs[mass]])
    mean = arr.mean(axis=0)
    var = arr.var(axis=0, ddof=1)

    if np.any(np.abs(tau_ws[mass] / mean) - 1 > 1e-2):
        print("Possibly problematic toy:")
        print(f"Mass: {mass}")
        print(f"{tau_ws[mass] / mean}")

# Write trees
for mass in globs:
    fn_out = os.path.join(args.outdir, f"toy_globs_{args.channel.lower()}_{mass}.root")
    fout = R.TFile.Open(fn_out, "RECREATE")

    tree = R.TTree(f"globs_{args.channel.lower()}", f"globs_{args.channel.lower()}")
    tree.SetDirectory(fout)

    nbins = len(tau_ws[mass])

    tree_index = np.zeros(1, dtype=np.int32)
    tree_globs = np.zeros(nbins, dtype=np.float32)

    tree.Branch("index", tree_index, "index/I")
    tree.Branch("globs", tree_globs, f"globs[{nbins}]/F")

    # Iterate over index, global observables and fill tree
    for i, g in sorted(globs[mass], key=lambda x: x[0]):
        tree_index[0] = i
        np.copyto(tree_globs, g)
        tree.Fill()

    tree.Write()
    fout.Close()
