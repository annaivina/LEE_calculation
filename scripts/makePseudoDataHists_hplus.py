#!/usr/bin/env python
from math import sqrt
from tqdm import tqdm
import argparse
import h5py
import numpy as np
import os

from utils_hplus import masspoints
from hplus_utils_4j3b import edgesHplusPrebin, edges4j3b
from hplus_utils_4j4b import edges4j4b
from hplus_utils_5j3b import edges5j3b
from hplus_utils_5j4b import edges5j4b
from hplus_utils_6j3b import edges6j3b
from hplus_utils_6j4b import edges6j4b


import ROOT as R
R.gROOT.SetBatch(True)


parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("-o", "--outdir", required=True)
parser.add_argument("-c", "--channel", choices=["4j3b", "4j4b", "5j3b", "5j4b", "6j3b", "6j4b"], required=True)
parser.add_argument("--nToys", default=20000, type=int)
args = parser.parse_args()


# Figure out what binning to use
binning = None
binningPostRebin = None

if args.channel == "4j3b":
    binning = edgesHplusPrebin
    binningPostRebin = edges4j3b

elif args.channel == "4j4b":
    binning = edgesHplusPrebin
    binningPostRebin = edges4j4b

elif args.channel == "5j3b":
    binning = edgesHplusPrebin
    binningPostRebin = edges5j3b

elif args.channel == "5j4b":
    binning = edgesHplusPrebin
    binningPostRebin = edges5j4b

elif args.channel == "6j3b":
    binning = edgesHplusPrebin
    binningPostRebin = edges6j3b

elif args.channel == "6j4b":
    binning = edgesHplusPrebin
    binningPostRebin = edges6j4b

else:
    raise RuntimeError("Unknown channel")


# Histogram prototype for cloning
proto_hist = R.TH1F("proto", "", len(binning) - 1, binning)

# Center of the bins after rebinning
centers = {}
for mass in masspoints:
    centers[mass] = 0.5 * (binningPostRebin[mass][:-1] + binningPostRebin[mass][1:])

# Mapping between post-rebin bin index to pre-rebin bin center
#
# Can just put data at center of bin since the rebinning will again be
# performed when building the workspace.
bin_idx_map = {}
for mass in masspoints:
    # Using `side="right"` to mimic ROOT convention of 0 being underflow bin
    bin_idx_map[mass] = np.searchsorted(binning, centers[mass], side="right")


# Poisson random variables to use for WS inputs
with h5py.File(args.infile, "r") as fin:
    # Only use the first couple of toys
    rvs = np.array(fin.get("poisson_rvs")[:args.nToys])
    bin_labels = np.array(fin.get("bin_labels"))


# Making empty histograms to fill
hists = {mass: [] for mass in masspoints}
for mass in hists:
    for itoy in range(args.nToys):
        h = proto_hist.Clone(f"PseudoData{itoy}")
        h.SetTitle(f"PseudoData{itoy}")
        h.SetDirectory(0)
        hists[mass].append(h)


print("Filling histograms...")
for itoy in tqdm(range(args.nToys)):
    row = rvs[itoy, :]
    for entries, (mass, ibin) in zip(row, bin_labels):
        target_index = int(bin_idx_map[mass][ibin - 1])
        hists[mass][itoy].SetBinContent(target_index, entries)
        hists[mass][itoy].SetBinError(target_index, sqrt(entries))


print("Writing histograms...")
for mass in tqdm(hists):
    fn_out = os.path.join(args.outdir, f"pseudodata_{args.channel.lower()}_{mass}.root")
    fout = R.TFile.Open(fn_out, "RECREATE")

    for hist in hists[mass]:
        hist.Write()

    fout.Close()
