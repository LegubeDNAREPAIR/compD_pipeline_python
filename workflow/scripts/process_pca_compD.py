import cooltools
import numpy as np
import scipy
import scipy.stats
import cooler
import pandas as pd
import functools
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from cooltools import numutils
from scipy.sparse import csr_matrix

def make_pca_diff(clr,my_chrs = "chr17"):
    partition = clr.extent(my_chrs)
    lo = partition[0]
    hi = partition[-1]
    bins = clr.bins()[lo:hi]
    C = clr.matrix(balance=False)[lo:hi, lo:hi]
    C = np.nan_to_num(C)
    # Finally, use PCA to reduce the dimensionality of the data
    pca = PCA(n_components=1)
    pca.fit(np.asarray(C))

    # You can access the principal components of the data using the `components_` attribute of the PCA object
    pca_components = pca.components_
    bins["PC1"] = pca_components[0,:].tolist()
    return(bins)


def make_pca_diff_trans(clr):
    bins = clr.bins()[0:len(clr.bins())]
    C = clr.matrix(balance=False)[0:len(clr.bins()), 0:len(clr.bins())]
    C = np.nan_to_num(C)
    # Finally, use PCA to reduce the dimensionality of the data
    pca = PCA(n_components=1)
    pca.fit(np.asarray(C))

    # You can access the principal components of the data using the `components_` attribute of the PCA object
    pca_components = pca.components_
    bins["PC1"] = pca_components[0,:].tolist()
    return(bins)


# Read the input file from the first command-line argument
input_filename = snakemake.input[0]
clr = cooler.Cooler(input_filename)
my_chrs = clr.chromnames


results = map(functools.partial(make_pca_diff,clr), my_chrs)
results = list(results)
results = pd.concat(results)

# Write the output to the file specified by the second command-line argument
results.to_csv(snakemake.output[0],index=False,header=False,sep='\t')


results_trans = make_pca_diff_trans(clr)

# Write the output to the file specified by the third command-line argument
results_trans.to_csv(snakemake.output[1],index=False,header=False,sep='\t')
