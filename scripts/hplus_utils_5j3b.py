import numpy as np
import pandas as pd
import pickle
import uproot

from utils_hplus import masspoints


# Bin edges before rebinning (i.e. original histograms)
edgesHplusPrebin = np.linspace(0, 1, 1001, dtype=np.float64)


# Bin edges after rebinning
with open("edges_5jex3bex.pkl", "rb") as f1:
	edges5j3b = pickle.load(f1)


def getHplusDf(treename, filename):
	print(filename)
	flag = 0
	variables = ['nomWeight_weight_btag','nomWeight_weight_jvt','nomWeight_weight_leptonSF','nomWeight_weight_mc','nomWeight_weight_norm','nomWeight_weight_pu']
	file = filename.split("/")[-1].split(".")[0]
	if( file == 'ttbarbb' or file == 'ttbarcc' or file == 'ttbarlight' or file == 'Wtocb'):
		variables += ['branch_rew_func_5j']
		flag = 1
	elif(file=="Data"):
		flag = 2
	elif(file == 'Hpluscb_60' or file == 'Hpluscb_70' or file == 'Hpluscb_80' or file == 'Hpluscb_90' or file == 'Hpluscb_100' or file == 'Hpluscb_110' or file == 'Hpluscb_120' or file == 'Hpluscb_130' or file == 'Hpluscb_140' or file == 'Hpluscb_150'):
		variables += ['branch_rew_func_5j']
		flag = 3
	elif(file == 'Hpluscb_160'):
		variables += ['branch_rew_func_5j']
		flag = 4

	variables += ['jets_n','bjets_n']
	variables += [f"scores_{mass}" for mass in masspoints]



	with uproot.open(filename) as f:
		t = f[treename]
		num_events = t.num_entries-1
		#num_events = 100
		df = t.arrays(variables, entry_stop=num_events, library="pd")


	if(flag==1):
		df["weight"] = df["nomWeight_weight_btag"]*df["nomWeight_weight_jvt"]*df["nomWeight_weight_leptonSF"]*df["nomWeight_weight_mc"]*df["nomWeight_weight_norm"]*df["nomWeight_weight_pu"]*df['branch_rew_func_5j']
		df.drop(columns=['nomWeight_weight_btag','nomWeight_weight_jvt','nomWeight_weight_leptonSF','nomWeight_weight_mc','nomWeight_weight_norm','nomWeight_weight_pu','branch_rew_func_5j'], inplace=True)
	elif(flag==2):
		df["weight"] = 1
		df.drop(columns=['nomWeight_weight_btag','nomWeight_weight_jvt','nomWeight_weight_leptonSF','nomWeight_weight_mc','nomWeight_weight_norm','nomWeight_weight_pu'], inplace=True)
	elif(flag==3):
		df["weight"] = df["nomWeight_weight_btag"]*df["nomWeight_weight_jvt"]*df["nomWeight_weight_leptonSF"]*df["nomWeight_weight_mc"]*df["nomWeight_weight_norm"]*df["nomWeight_weight_pu"]*df['branch_rew_func_5j']
		df.drop(columns=['nomWeight_weight_btag','nomWeight_weight_jvt','nomWeight_weight_leptonSF','nomWeight_weight_mc','nomWeight_weight_norm','nomWeight_weight_pu','branch_rew_func_5j'], inplace=True)
	elif(flag==4):
		sels = (df["nomWeight_weight_mc"]>-100) & (df["nomWeight_weight_mc"]<100)
		df = df.loc[sels]
		df["weight"] = df["nomWeight_weight_btag"]*df["nomWeight_weight_jvt"]*df["nomWeight_weight_leptonSF"]*df["nomWeight_weight_mc"]*df["nomWeight_weight_norm"]*df["nomWeight_weight_pu"]*df['branch_rew_func_5j']*22059247.969726562/185760.6083984375
		df.drop(columns=['nomWeight_weight_btag','nomWeight_weight_jvt','nomWeight_weight_leptonSF','nomWeight_weight_mc','nomWeight_weight_norm','nomWeight_weight_pu','branch_rew_func_5j'], inplace=True)
	else:
		df["weight"] = df["nomWeight_weight_btag"]*df["nomWeight_weight_jvt"]*df["nomWeight_weight_leptonSF"]*df["nomWeight_weight_mc"]*df["nomWeight_weight_norm"]*df["nomWeight_weight_pu"]
		df.drop(columns=['nomWeight_weight_btag','nomWeight_weight_jvt','nomWeight_weight_leptonSF','nomWeight_weight_mc','nomWeight_weight_norm','nomWeight_weight_pu'], inplace=True)


	# Add sample columns
	df["sample"] = filename.split("/")[-1].split(".")[0]

	# Apply selection
	sel = (df["jets_n"]==5) & (df["bjets_n"] == 3)
	df = df.loc[sel]

	#df[f"scores_4j3b_{mass}"] = df[f"scores_{mass}"]

	df.drop(columns=["jets_n", "bjets_n"], inplace=True)

	return df.copy()
