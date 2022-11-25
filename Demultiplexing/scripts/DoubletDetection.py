#!/usr/bin/env python
import numpy as np
import doubletdetection
import tarfile
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import os
import argparse
import sys
import pandas as pd

# Load read10x function from mods directory

mods_path = "/opt/WG1-pipeline-QC_GoogleCloud/Demultiplexing/mods"
sys.path.append(mods_path)
import read10x

parser = argparse.ArgumentParser( 
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("-m", "--counts_matrix", required = True, help = "cell ranger counts matrix directory containing matrix files or full path to matrix.mtx. Can also also provide the 10x h5.")
parser.add_argument("-b", "--barcodes", required = False, default = None, help = "File containing droplet barcodes. Use barcodes from provided 10x dir by default.")
parser.add_argument("-f", "--filtered_barcodes", required = False, default = None, help = "File containing a filtered list of droplet barcodes. This may be used if you want to use a filtered list of barcodes for doublet detection (ie need to remove droplets that are empty or high in ambient RNA).")
parser.add_argument("-o", "--outdir", required = False, default = os.getcwd(), help = "The output directory; default is current working directory")
parser.add_argument("-p", "--boost_rate", required = False, default = 0.25, type = float, help = "Proportion of cells used to generate synthetic doublets; default is 0.25.")
parser.add_argument("-c", "--n_components", required = False, default = 30, type = int, help = "Number of principal components to use; default is 30.")
parser.add_argument("-g", "--n_top_var_genes", required = False, default = 1000, type = int, help = "Number of top variable genes to use; default is 1000.")
parser.add_argument("-r", "--replace", required = False, default = False, help = "Whether to replace cells when generating synthetic doublets; default is False.")
parser.add_argument("-a", "--clustering_algorithm", required = False, default = 'phenograph', help = "Which clustering algorithm to use; default is 'phenograph'")
parser.add_argument("-k", "--clustering_kwargs", required = False, default = None, help = "Keyword arguments to pass to clustering algorithm; default is None.")
parser.add_argument("-i", "--n_iterations", required = False, default = 50, type = int, help = "Number of iterations to use; default is 50")
parser.add_argument("-e", "--pseudocount", required = False, default = 0.1, type = float, help = "Pseudocount used to normalize counts; default is 0.1.")
parser.add_argument("-n", "--normalizer", required = False, default = None, help = "Method for raw counts normalization; default is None.")
parser.add_argument("-d", "--random_state", required = False, default = 0, type = int, help = "Number to use to seed random state for PCA; default is 0.")
parser.add_argument("-s", "--standard_scaling", required = False, default = True, help = "Whether to use standard scaling of normalized count matrix prior to PCA (True) or not (False); default is True")
parser.add_argument("-j", "--n_jobs", required = False, default = 1, type = int, help = "Number of jobs to to use; default is 1")
parser.add_argument("-t", "--p_thresh", required = False, default = 1e-16, type = float, help = "P-value threshold for doublet calling; default is 1e-16")
parser.add_argument("-v", "--voter_thresh", required = False, default = 0.5, type = float, help = "Voter threshold for doublet calling; default is 0.5")
args = parser.parse_args()


print(args.p_thresh)
print(args.voter_thresh)



if args.clustering_kwargs == 'True':
    pheno = True
elif args.clustering_kwargs == 'False':
    pheno = False
else:
    pheno = args.clustering_kwargs
print(pheno)

if args.standard_scaling == 'True':
    standard_scaling = True
elif args.standard_scaling == 'False':
    standard_scaling = False
else:
    standard_scaling = args.standard_scaling
print(standard_scaling)


### Read in data ###
raw_counts = read10x.import_cellranger_mtx(args.counts_matrix)
barcodes_df = read10x.read_barcodes(args.barcodes)

print('Counts matrix shape: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

# Remove columns with all 0s
zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
raw_counts = raw_counts[:, ~zero_genes]
print('Counts matrix shape after removing unexpressed genes: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

clf = doubletdetection.BoostClassifier(boost_rate=args.boost_rate, n_components=args.n_components, n_top_var_genes=args.n_top_var_genes, replace=args.replace, clustering_algorithm=args.clustering_algorithm, clustering_kwargs=args.clustering_kwargs, n_iters=args.n_iterations, normalizer=args.normalizer, pseudocount=args.pseudocount, random_state=args.random_state, standard_scaling=args.standard_scaling, n_jobs=args.n_jobs, verbose = True)
doublets = clf.fit(raw_counts).predict(p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)

results = pd.Series(doublets, name="DoubletDetection_DropletType")
dataframe = pd.concat([barcodes_df, results], axis=1)
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(1.0, "doublet")
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(0.0, "singlet")

dataframe.to_csv(os.path.join(args.outdir,'DoubletDetection_results.txt'), sep = "\t", index = False)


### Figures ###
doubletdetection.plot.convergence(clf, save=os.path.join(args.outdir,'convergence_test.pdf'), show=False, p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)

f3 = doubletdetection.plot.threshold(clf, save=os.path.join(args.outdir,'threshold_test.pdf'), show=False, p_step=6)