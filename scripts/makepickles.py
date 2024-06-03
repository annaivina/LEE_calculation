#!/usr/bin/env python3
import argparse
import numpy as np
import pickle
import re

parser = argparse.ArgumentParser()
#parser.add_argument("-m", "--mass", type=int, required=True)
parser.add_argument("-o", "--outdir", default="")
args = parser.parse_args()


masses = np.array(['60','70','80','90','100','110','120','130','140','150','160'])
# Output dictionaries
edges_4jex3bex_dict = {}
edges_4jex4bin_dict = {}
edges_5jex3bex_dict = {}
edges_5jex4bin_dict = {}
edges_6jex3bex_dict = {}
edges_6jex4bin_dict = {}



for i in masses:
	edges_4jex3bex_dict[int(i)] = np.array([0.,0.2717,0.3385,0.4184,0.4936,0.5624,0.632,0.711,0.8,0.8909,1.0], dtype=np.float64)
	edges_4jex4bin_dict[int(i)] = np.array([0.,1.0], dtype=np.float64)
	edges_5jex3bex_dict[int(i)] = np.array([0.,0.2178,0.2904,0.3641,0.4365,0.5021,0.5647,0.6331,0.7223,0.8376,1.], dtype=np.float64)
	edges_5jex4bin_dict[int(i)] = np.array([0.,1.0], dtype=np.float64)
	edges_6jex3bex_dict[int(i)] = np.array([0.,0.1479,0.2041,0.2672,0.3355,0.3992,0.4595,0.5186,0.5887,0.702,1.], dtype=np.float64)
	edges_6jex4bin_dict[int(i)] = np.array([0.0,1.], dtype=np.float64)


	with open("edges_4jex3bex.pkl", "wb") as fout:
		pickle.dump(edges_4jex3bex_dict, fout)

	with open("edges_4jex4bin.pkl", "wb") as fout:
		pickle.dump(edges_4jex4bin_dict, fout)

	with open("edges_5jex3bex.pkl", "wb") as fout:
		pickle.dump(edges_5jex3bex_dict, fout)

	with open("edges_5jex4bin.pkl", "wb") as fout:
		pickle.dump(edges_5jex4bin_dict, fout)

	with open("edges_6jex3bex.pkl", "wb") as fout:
		pickle.dump(edges_6jex3bex_dict, fout)

	with open("edges_6jex4bin.pkl", "wb") as fout:
		pickle.dump(edges_6jex4bin_dict, fout)



