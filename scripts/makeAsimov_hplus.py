#!/usr/bin/env python3
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("workspace")
parser.add_argument("-m", "--mass", type=int, required=True)
parser.add_argument("-o", "--outfile", required=True)
args = parser.parse_args()


import ROOT as R
R.gROOT.SetBatch(True)


def fixDataset(dataset):
    fixed_dataset = R.RooDataSet("ds", "", dataset.get(0), "weightVar")
    for i in range(dataset.numEntries()):
        argset = dataset.get(i)
        fixed_dataset.add(argset, dataset.weight())

    return fixed_dataset


def getHistFromAsimov(dataset, observable, cut):
    dataset_reduced = dataset.reduce(R.RooArgSet(observable), cut)
    return dataset_reduced.createHistogram("h", observable)


f = R.TFile.Open(args.workspace)
w = f.Get("combined")
model = w.obj("ModelConfig")
print("Printing MODEL CONFIG: \n")
print(model)
# Observables
obs_4j3b = w.obj("obs_x_c1l4jex3bex")
obs_4j4b = w.obj("obs_x_c1l4jex4bin")
obs_5j3b = w.obj("obs_x_c1l5jex3bex")
obs_5j4b = w.obj("obs_x_c1l5jex4bin")
obs_6j3b = w.obj("obs_x_c1l6jex3bex")
obs_6j4b = w.obj("obs_x_c1l6jex4bin")



# POI
mu = model.GetParametersOfInterest().first()
mu.setVal(0.0)
mu.setConstant()


asimov = R.RooStats.AsymptoticCalculator.MakeAsimovData(
    model,
    R.RooArgSet(mu),
    model.GetGlobalObservables())

# Makes the weighting work
asimov = fixDataset(asimov)


# Histograms of Asimov observables
h_obs_4j3b = getHistFromAsimov(asimov, obs_4j3b,"channelCat == channelCat::c1l4jex3bex")
h_obs_4j4b = getHistFromAsimov(asimov, obs_4j4b,"channelCat == channelCat::c1l4jex4bin")
h_obs_5j3b = getHistFromAsimov(asimov, obs_5j3b,"channelCat == channelCat::c1l5jex3bex")
h_obs_5j4b = getHistFromAsimov(asimov, obs_5j4b,"channelCat == channelCat::c1l5jex4bin")
h_obs_6j3b = getHistFromAsimov(asimov, obs_6j3b,"channelCat == channelCat::c1l6jex3bex")
h_obs_6j4b = getHistFromAsimov(asimov, obs_6j4b,"channelCat == channelCat::c1l6jex4bin")



# Asimov global observables (construct the whole string for the gamma factors)
pattern = re.compile(
    r"^nom_gamma_stat_.*"
    r"(c1l4jex3bex|c1l4jex4bin|c1l5jex3bex|c1l5jex4bin|c1l6jex3bex|c1l6jex4bin)"
    r".*_bin_(\d+)$"
)

print(pattern)
print(model.GetGlobalObservables())
gamma_globs = {}
for param in model.GetGlobalObservables():
    name = param.GetName()
    m = pattern.match(name)
    if not m:
        continue

    region, ibin = m.groups()
    ibin = int(ibin)

    if region == "c1l4jex3bex":
        region = "4j3b"
    elif region == "c1l4jex4bin":
        region = "4j4b"
    elif region == "c1l5jex3bex":
        region = "5j3b"
    elif region == "c1l5jex4bin":
    	region = "5j4b"
    elif region == "c1l6jex3bex":
        region = "6j3b"
    elif region == "c1l6jex4bin":
        region = "6j4b"
    else:
        raise RuntimeError("Unknown value encountered for region")

    gamma_globs.setdefault(region, []).append((ibin, param.getVal()))

print("Printing gamma globs")
print(gamma_globs)

# Make histograms from gamma observables
glob_hists = {}
for key in gamma_globs:
    hbins = [glob for ibin, glob in sorted(gamma_globs[key], key=lambda x: x[0])]
    nbins = len(hbins)

    h = R.TH1F(f"tau_{key}", "tau_{key}", nbins, 0, 1)

    for i, content in enumerate(hbins, start=1):
        h.SetBinContent(i, content)
        h.SetBinError(i, R.TMath.Sqrt(content))

    glob_hists[key] = h


def writeRootObject(obj, name, directory):
    obj.SetName(name)
    obj.SetTitle(name)

    directory.cd()
    obj.Write()


fout = R.TFile.Open(args.outfile, "RECREATE")

# Observables
writeRootObject(h_obs_4j3b, f"h_obs_4j3b_m{args.mass}", fout)
writeRootObject(h_obs_4j4b, f"h_obs_4j4b_m{args.mass}", fout)
writeRootObject(h_obs_5j3b, f"h_obs_5j3b_m{args.mass}", fout)
writeRootObject(h_obs_5j4b, f"h_obs_5j4b_m{args.mass}", fout)
writeRootObject(h_obs_6j3b, f"h_obs_6j3b_m{args.mass}", fout)
writeRootObject(h_obs_6j4b, f"h_obs_6j4b_m{args.mass}", fout)


# Global Observables
writeRootObject(glob_hists["4j3b"], f"tau_4j3b_m{args.mass}", fout)
writeRootObject(glob_hists["4j4b"], f"tau_4j4b_m{args.mass}", fout)
writeRootObject(glob_hists["5j3b"], f"tau_5j3b_m{args.mass}", fout)
writeRootObject(glob_hists["5j4b"], f"tau_5j4b_m{args.mass}", fout)
writeRootObject(glob_hists["6j3b"], f"tau_6j3b_m{args.mass}", fout)
writeRootObject(glob_hists["6j4b"], f"tau_6j4b_m{args.mass}", fout)

fout.Close()

