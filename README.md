# LEE_calculation
The code is adapted to be able to use for the HPlus analysis. 


# Toy Generation for Hplus ->cb analysis 

## Requirements

- ROOT 6.22.06
- Python 3 (see `requirements.txt` for additional packages)


## Step 0: Preparation

1. Create a working directory
2. Get the workspaces and workspace building log files
3. Get ntuples for hadhad, lephad SLT + LTT


```

## Step 1: Discriminant Binning

Parse the WSMaker workspace building log files to get the bin edges
used for the final discriminant.

```bash
makepickles.py (for our studies)
```

Six pickle files should be produced that store the bin edges indices
used for rebinning -> 6 files -> 6 categories


## Step 2: Create Asimov Datasets

Produce histograms of the Asimov dataset providing:

1. The background expectation
2. The effective number of MC events (tau) for the Barlow-Beeston method


```bash
mkdir asimov

for mass in 60 70 80 90 100 110 120 130 140 150 160; do
    python3 scripts/makeAsimov_hplus.py data/workspaces/${mass}.root -m ${mass} -o asimov/asimov_${mass}.root
done

(
    cd asimov
    [[ -e asimov_merged.root ]] && rm asimov_merged.root
    hadd asimov_merged.root asimov_*.root
)
```

The histograms will store the value of the observables `obs_*`,
effective MC statistics `tau_*` for all channels:
`{4jex3bex,4jex4bex,5jex3bex,5jex4bex,6jex3bex,6jex4bex,}` and for all discriminants `_mass*`.


## Step 3: Make Dataframes

Ntuples are converted to pandas dataframes for ease of use during toy
generation. Additionally, the MVA scores are already discretized
according to the binning used in the final fit.


```bash
python3 scripts/makeHPlusDataframes_4j3b.py data/datainputs.txt -o dataframes/dataframe_4j3b.h5
python3 scripts/makeHPlusDataframes_4j4b.py data/datainputs.txt -o dataframes/dataframe_4j4b.h5
python3 scripts/makeHPlusDataframes_5j3b.py data/datainputs.txt -o dataframes/dataframe_5j3b.h5
python3 scripts/makeHPlusDataframes_5j4b.py data/datainputs.txt -o dataframes/dataframe_5j4b.h5
python3 scripts/makeHPlusDataframes_6j3b.py data/datainputs.txt -o dataframes/dataframe_6j3b.h5
python3 scripts/makeHPlusDataframes_6j4b.py data/datainputs.txt -o dataframes/dataframe_6j4b.h5
```


## Step 4: Expected Correlation Between Bins in Data

Estimates the expected correlation matrix of data yields per bin using
a bivariate Poisson model.

```bash
mkdir correlation_matrices

python3 scripts/makeCorr_hplus.py dataframes/dataframe_4j3b.h5 -o correlation_matrices/corr_4j3b.h5
python3 scripts/makeCorr_hplus.py dataframes/dataframe_4j4b.h5 -o correlation_matrices/corr_4j4b.h5
python3 scripts/makeCorr_hplus.py dataframes/dataframe_5j3b.h5 -o correlation_matrices/corr_5j3b.h5
python3 scripts/makeCorr_hplus.py dataframes/dataframe_5j4b.h5 -o correlation_matrices/corr_5j4b.h5
python3 scripts/makeCorr_hplus.py dataframes/dataframe_6j3b.h5 -o correlation_matrices/corr_6j3b.h5
python3 scripts/makeCorr_hplus.py dataframes/dataframe_6j4b.h5 -o correlation_matrices/corr_6j4b.h5
```


## Step 5: Generate Poisson RVS

This generates correlated Poisson random variables for the observables
in the 6 category regions channels:

```bash
mkdir poisson_rvs

python3 scripts/generateFromCorr_hplus.py correlation_matrices/corr_4j3b.h5 asimov/asimov_merged.root \
    -c 4j3b -o poisson_rvs/rvs_4j3b.h5

python3 scripts/generateFromCorr_hplus.py correlation_matrices/corr_4j4b.h5 asimov/asimov_merged.root \
    -c 4j4b -o poisson_rvs/rvs_4j4b.h5

python3 scripts/generateFromCorr_hplus.py correlation_matrices/corr_5j3b.h5 asimov/asimov_merged.root \
    -c 5j3b -o poisson_rvs/rvs_5j3b.h5

    ... 
```


## Step 6: Generate Global Observables (Barlow-Beeston)

This will create multiple root-files containing trees with the values
of the global observables.


```bash
mkdir gamma_globs

makeGammaGlobsToys_hplus.py dataframes/dataframe_4j3b.h5 asimov/asimov_merged.root \
    -o gamma_globs -c 4j3b

... ( all regions)
```


## Step 7: Generate Global Observables (Others)

All other (Gaussian-constrained) global observables are treated as
fully correlated. The script takes *all* workspaces so that it can
figure what the names of all NPs are.

```bash
mkdir other_globs

python3 scripts/makeAlphaGlobsToys_hplus.py \
 data/workspaces/{60,70,80,90,100,110,120,130,140,150,160}.root \
 -o other_globs/alphas.root
```



## Step 9: Build Workspace Inputs

This (mostly technical) step translates the pseudo-data that was
generated in the final fit binning to the initial histogram binning so
that it can be used in the workspace building as if it were real data.

**Pseudo-data:**

```bash
mkdir ws_inputs

makePseudoDataHists_hplus.py poisson_rvs/rvs_4j3b.h5 -c 4j3b -o ws_inputs/
...
```

**Global observables:**

The global observables are stored as trees where the index of the
entry corresponds to the toy experiment. We can just merge the ROOT
files for all global observables:

```bash
for mass in 60 70 80 90  \
            100 110 120 130 140 150 160; do
         
    hadd ws_inputs/toy_globs_${mass}.root \
        gamma_globs/toy_globs_{4j3b,4j4b,5j3b,5j4b,6j3b,6j4b}_${mass}.root \
        other_globs/alphas.root
done
```

