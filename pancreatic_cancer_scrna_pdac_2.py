import os
import scanpy as sc
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import doubletdetection
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import legend
from networkx.algorithms.bipartite import color
from phenograph import cluster
from scipy.stats import median_abs_deviation as mad
import numpy as np
import celltypist
from celltypist import models
import scvi
from scvi.autotune import ModelTuner
from ray import tune
import warnings
from pandas.errors import PerformanceWarning
from seaborn.matrix import dendrogram
from sympy import prime
from sympy.abc import alpha
from scipy import stats
import gseapy as gp
import diffxpy as de

warnings.simplefilter("ignore", PerformanceWarning)

pd.set_option('display.max_columns', None)  # to visualize all the columns in the terminal
pd.set_option('display.width', 1500)

# path to local directories
mtx_pdac = '/home/prince/PycharmProjects/pancreatic_scrna_seq/data_PDAC/'
h5ad_pdac = '/home/prince/PycharmProjects/pancreatic_scrna_seq/h5ad_data_PDAC/'
pp_pdac = '/home/prince/PycharmProjects/pancreatic_scrna_seq/pp_data_PDAC/'
models_dir = '/home/prince/PycharmProjects/pancreatic_scrna_seq/models/'
pred_pdac = '/home/prince/PycharmProjects/pancreatic_scrna_seq/predictions_PDAC/'
extra_pdac = '/home/prince/PycharmProjects/pancreatic_scrna_seq/extra_results_PDAC/'
scvi_pdac_model = '/home/prince/PycharmProjects/pancreatic_scrna_seq/models/scvi_PDAC'
cluster_results = '/home/prince/PycharmProjects/pancreatic_scrna_seq/cluster_results/'
analysis_pdac = '/home/prince/PycharmProjects/pancreatic_scrna_seq/analysis_results_PDAC/'




'''
# ++++++++++++++++++++++++++++++++ Extra works to do (depending on data sets present) ++++++++++++++++++++++++++++++++++#


# 1. (if required) changing the file names of _genes.tsv extension to _features.tsv extension which is preferred by scanpy

os.chdir(mtx_pdac)
for filename in os.listdir(mtx_pdac):
    if filename.endswith('_genes.tsv.gz'):
        new_filename = filename.replace('_genes', '_features')
        os.rename(filename, new_filename)
        print(f'Renamed {filename} to {new_filename}')

# 2. printing sample for verification of mtx data

prefix = 'GSM6567157_PDAC1_'  # change it accordingly
adata = sc.read_10x_mtx(mtx_pdac, prefix=prefix, var_names='gene_symbols', cache=True)

print(adata)
print("Do you see gene names? then ok!\n", adata.var_names[:10])  # Print the first 10 gene names
print("Do you see barcode names? then ok!\n", adata.obs_names[:10])  # print the first 10 cells

# 3. Converting the mtx data into hda5 for later analysis

# Creating the output directory if it does not exist
if not os.path.exists(h5ad_pdac):
    os.makedirs(h5ad_pdac)

# Looping through each dataset file
for file in os.listdir(mtx_pdac):
    if file.endswith('_matrix.mtx.gz'):
        prefix = file.split('matrix.mtx.gz')[0]
        try:
            # Load the data
            adata = sc.read_10x_mtx(
                mtx_pdac,
                var_names='gene_symbols',  # Use 'gene_symbols' or omit if you prefer Ensembl IDs
                prefix=prefix,
                cache=True
            )
            # prefix = prefix.rstrip('_')
            # Saving the processed AnnData object
            output_file = os.path.join(h5ad_pdac, f'{prefix}processed.h5ad')
            adata.write(output_file)
            print(f'Processed and saved {output_file}')

        except Exception as e:
            print(f'Error processing {file}: {e}')

# 4. Cellbender to remove ambient RNA (not needed, if the files are pre-processed)

# List all .h5ad files
files = [f for f in os.listdir(h5ad_pdac) if f.endswith('.h5ad')]

for file in files:
    input_file = os.path.join(h5ad_pdac, file)
    output_file = os.path.join(clean_data, f"{os.path.splitext(file)[0]}_denoised.h5ad")

    # Run the CellBender command
    command = (
        f"conda run -n cellbender cellbender remove-background "
        f"--input {input_file} "
        f"--output {output_file} "
        f"--total-droplets-included 20000 "
        f"--cuda"
    )

    subprocess.run(command, shell=True, executable="/bin/bash")
'''




'''
# +++++++++++++++++++++++++++++++++++ QC ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# list all the files in the directory
adatas = [x for x in os.listdir(h5ad_pdac)]


def load_it(adata):
    samp = adata.split('_')[1]  # Extract the sample identifier
    adata = sc.read_h5ad(h5ad_pdac + adata)  # Load the .h5ad file using Scanpy
    adata.obs['Sample'] = samp  # Add a new column 'Patient' to the cell metadata (obs) with the sample identifier
    adata.obs.index = adata.obs.index + '-' + samp  # Modify the cell barcodes to append the sample identifier
    return adata


# applying load_it function
adatas = [load_it(ad) for ad in adatas]
print(adatas)


# QC function to filter, add new variables and calculate qc metrics
def qc(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

    remove = ['total_counts_mt', 'log1p_total_counts_mt', 'total_counts_ribo',
              'log1p_total_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb']

    adata.obs = adata.obs[[x for x in adata.obs.columns if x not in remove]]

    return adata


# applying qc function
adatas = [qc(ad) for ad in adatas]
# print(adatas)

# combining results (from all .obs dataframes in loaded .h5ad files into one dataframe)
df = pd.concat([x.obs for x in adatas])
df = df.sort_values('Sample')
print(df)
'''


'''
# Visualization using seaborn based on plots of different QC metrics (comment out "value" accordingly)

# Defining the list of metrics to plot
metrics = ["pct_counts_mt", "n_genes", "pct_counts_in_top_20_genes", "log1p_total_counts"]

# Looping over the metrics and plot each one
for value in metrics:
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    g = sns.FacetGrid(df, row="Sample", hue="Sample", aspect=15, height=0.5, palette="tab20")

    # Map the kde plot for the current value
    g.map(sns.kdeplot, value, clip_on=False, fill=True, alpha=1, linewidth=1.5)
    g.map(sns.kdeplot, value, clip_on=False, color="w", lw=2)

    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)

    g.map(label, value)

    g.figure.subplots_adjust(hspace=-.6)

    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)

    for ax in g.axes.flat:
        ax.axvline(x=df[value].median(), color='r', linestyle='-')

    plt.show()
'''

'''
#+++++++++++++++++++++++++++++++++++ Pre-processing ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# 1. checking data to know possibly pre-processed or to what extent (change sample accordingly)
a = df[df.Sample == 'PDAC5'].log1p_total_counts
np.median(a)
np.median(a) - 5 * mad(a)   # mad = median absolute deviation
np.median(a) + 5 * mad(a)

ax = sns.displot(a)

plt.axvline(np.median(a) - 5 * mad(a))
plt.axvline(np.median(a) + 5 * mad(a))

#plt.show()

# 2. Main pre-processing steps

# calling a helper function that identifies outliers in given qc metrics
#   adata = input anndata containing sc data
#   metric = the column in adata.obs to check for outliers
#   nmads = number of MAD (median absolute deviation) to use for determining outliers; will be set while calling function
#   upper_only (if TRUE) = only detects outliers above the upper threshold
#   upper_only (if FALSE) = checks for outliers both above and below of a threshold
#   M is the series (column) from adata.obs corresponding to the metric
def mad_outlier(adata, metric, nmads, upper_only=False):
    M = adata.obs[metric]

    if not upper_only:
        return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))

    return (M > np.median(M) + nmads * mad(M))  # this part for pct_counts_mt only

# doublet detection classifier
clf = doubletdetection.BoostClassifier(
    n_iters=10,
    clustering_algorithm="louvain",
    standard_scaling=True,          # standardize input data
    pseudocount=0.1,                # adds a small value to avoid division by zero
    n_jobs=-1)  #  Number of CPU cores to use for parallel processing (-1 uses all available cores)

# pre-processing function
def pp(adata):
    adata = adata[adata.obs.pct_counts_mt < 25]     # high mt == dying cells

#   applying mad_outlier function (prints TRUE if only one of the metrics is TRUE)
    bool_vector = mad_outlier(adata, 'log1p_total_counts', 5) + \
                  mad_outlier(adata, 'log1p_n_genes_by_counts', 5) + \
                  mad_outlier(adata, 'pct_counts_in_top_20_genes', 5) + \
                  mad_outlier(adata, 'pct_counts_mt', 3, upper_only=True)   # lower values aren't a concern
    adata = adata[~bool_vector]     # takes only FALSE of the bool_vector (i.e., removes cells flagged as outliers)

#   storing number of cells removed as outliers
    adata.uns['cells_removed'] = sum(bool_vector)

#   applying doublet detection classifier
#   lower the p_thresh, the more the stringency of filtering (fewer doublets detected)
#   voter_thresh (0 to 1); higher = more conservative (fewer doublets detected)
    doublets = clf.fit(adata.X).predict(p_thresh=1e-3, voter_thresh=0.5)    # tweak threshold after observing results
    doublet_score = clf.doublet_score()

    adata.obs["doublet"] = doublets             # 0 = singlet, 1 = doublet
    adata.obs["doublet_score"] = doublet_score

#   removing doublets
    adata.uns['doublets_removed'] = adata.obs.doublet.sum()
    adata = adata[adata.obs.doublet == 0]

    return adata

# applying the pp function
adatas = [pp(ad) for ad in adatas]

# printing to observe the degree of doublet removal
for adata in adatas:
    print(len(adata), adata.uns['cells_removed'], adata.uns['doublets_removed'])

# Saving each preprocessed AnnData object
if not os.path.exists(pp_pdac):
    os.makedirs(pp_pdac)

for adata in adatas:
    sample_name = adata.obs['Sample'].unique()[0]  # Assuming 'Sample' column is the same for all cells in the adata
    file_name = f"{sample_name}_pp.h5ad"  # Define file name
    file_path = os.path.join(pp_pdac, file_name)  # Full path to save the file
    adata.write(file_path)  # Save the AnnData object
    print(f"Saved preprocessed data for {sample_name} to {file_path}")
'''

'''
#+++++++++++++++++++++++++++++++++++ QC after pre-processing ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

df2 = pd.concat([x.obs for x in adatas])
df2 = df2.sort_values('Sample')

# Defining the list of metrics to plot
metrics = ["pct_counts_mt", "n_genes", "pct_counts_in_top_20_genes", "log1p_total_counts"]

# Looping over the metrics and plot each one
for value in metrics:
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    g = sns.FacetGrid(df2, row="Sample", hue="Sample", aspect=15, height=0.5, palette="tab20")

    # Map the kde plot for the current value
    g.map(sns.kdeplot, value, clip_on=False, fill=True, alpha=1, linewidth=1.5)
    g.map(sns.kdeplot, value, clip_on=False, color="w", lw=2)

    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)

    g.map(label, value)

    g.figure.subplots_adjust(hspace=-.6)

    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)

    for ax in g.axes.flat:
        ax.axvline(x=df2[value].median(), color='r', linestyle='-')

    plt.show()
'''

'''
# observing the result so far
for x in os.listdir(pp_pdac):
    full_path = os.path.join(pp_pdac, x)  # Create the full path to the file
    print(x)
    a = sc.read_h5ad(full_path)  # Use the full path to read the file
    print(a.obs)
'''

'''
#+++++++++++++++++++++++++++++++++++++++++++++ Label transfer ++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# can use either celltypist or SCVI trained models
# celltypist models for download: https://www.celltypist.org/organs

# showing the inbuilt models of celltypist
#all = (models.get_all_models())
#for x in all:
#    print(x)

# path to models
model_1 = models.Model.load(model=os.path.join(models_dir, "Adult_Pancreas_Fasolinoetal.2022.pkl"))
model_2 = models.Model.load(model=os.path.join(models_dir, "Adult_Pancreas_Tostietal.2021.pkl"))
model_3 = models.Model.load(model="Immune_All_Low.pkl")     # this one is downloaded in-built

# functions to predict cell types
def predict_cells(adata):
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)  # not recommended for typical pp (but for celltypist)
    sc.pp.log1p(adata)

    adata.X = adata.X.toarray()

    predictions = celltypist.annotate(adata, model=model_1, majority_voting=False)
    predictions_adata = predictions.to_adata()
    adata.obs["mod1_label"] = predictions_adata.obs.loc[adata.obs.index, "predicted_labels"]
    adata.obs["mod1_score"] = predictions_adata.obs.loc[adata.obs.index, "conf_score"]

    predictions = celltypist.annotate(adata, model=model_2, majority_voting=False)
    predictions_adata = predictions.to_adata()
    adata.obs["mod2_label"] = predictions_adata.obs.loc[adata.obs.index, "predicted_labels"]
    adata.obs["mod2_score"] = predictions_adata.obs.loc[adata.obs.index, "conf_score"]

    predictions = celltypist.annotate(adata, model=model_3, majority_voting=False)
    predictions_adata = predictions.to_adata()
    adata.obs["mod3_label"] = predictions_adata.obs.loc[adata.obs.index, "predicted_labels"]
    adata.obs["mod3_score"] = predictions_adata.obs.loc[adata.obs.index, "conf_score"]

    return adata.obs


adatas = [sc.read_h5ad(pp_pdac + x) for x in os.listdir(pp_pdac)]
predictions = [predict_cells(ad.copy()) for ad in adatas]
predictions = pd.concat(predictions)[['mod1_label', 'mod1_score', 'mod2_label', 'mod2_score', 'mod3_label', 'mod3_score']]

print(predictions)

if not os.path.exists(pred_pdac):
    os.makedirs(pred_pdac)
predictions.to_csv(os.path.join(pred_pdac, 'PREDICTIONS.csv'))

'''

'''
#+++++++++++++++++++++++++++++++++++++++++++++++++ Concatenation (merging) ++++++++++++++++++++++++++++++++++++++++++++#
# Reading in all h5ad files
adatas = [sc.read_h5ad(pp_pdac + x) for x in os.listdir(pp_pdac)]

# Concatenating them into a single AnnData object
adata = sc.concat(adatas)
print(adata)

# Reading the predictions CSV file
predictions = pd.read_csv(pred_pdac + "PREDICTIONS.csv", index_col=0)
print(predictions)

# Merging predictions into the AnnData object based on the index
adata.obs = adata.obs.merge(right=predictions, left_index=True, right_index=True)
#print(adata.obs)

# Saving the updated AnnData object
#adata.write_h5ad(extra_pdac + 'unintegrated_PDAC.h5ad')

sc.pp.filter_genes(adata, min_cells = 100)  # set min_cells accordingly (observe cells to genes ratio)
print(adata)
'''

# +++++++++++++++++++++++++++++++++++++++++++++++++ Integration +++++++++++++++++++++++++++++++++++++++++++++++++++++++#


# Step 1: Tuning the hyperparameters
'''
#================ comment out this part after hyperparameters have been set ===========================================#
# setting up scVI model
model_cls = scvi.model.SCVI
model_cls.setup_anndata(adata,
                        categorical_covariate_keys = ['Sample'],
                        continuous_covariate_keys=['pct_counts_mt', 'pct_counts_ribo'])

# initializing model tuner
tuner = ModelTuner(model_cls)
print(tuner.info())

# defining hyperparameter search space (specifies a list from which to choose from)
search_space = {
    "n_hidden": tune.choice([92, 128, 192, 256]),   # more hidden units == more capacity (but risk of overfitting)
    "n_latent": tune.choice([10, 20, 30, 40, 50, 60]),  # larger n_latent == more complex structures captured (risk ,,)
    "n_layers": tune.choice([1, 2, 3]), # defines depth of network
    "lr": tune.loguniform(1e-4, 1e-2), # smaller lr (learning rate) == gradual updates per epoch (fine tuning)
    "gene_likelihood": tune.choice(["nb", "zinb"])} # nb = negative binomial likelihood; zinb = zero-inflated nb

# fitting the tuner (performs the hyperparameter tuning)
results = tuner.fit(
                    adata,
                    metric="validation_loss",   # specifies the evaluation metric to be minimized
                    resources = {'gpu': 1},
                    search_space = search_space,
                    num_samples = 100,          # number of different hyperparameter combinations to try
                    max_epochs = 20)

# evaluating results (to find the best model based on validation loss)
best_vl = 10000
best_i = 0
for i, res in enumerate(results.results):
    vl = res.metrics['validation_loss']

    if vl < best_vl:
        best_vl = vl
        best_i = i
print(results.results[best_i])
#======================================================================================================================#
'''

'''
Note: results were printed by running (50 minutes for 100 num_samples) up to this part of the script;
      based on the results tuned, the model for adata was set up later in:
      model = scvi.model.SCVI(adata, n_hidden = ?, n_latent = ?, n_layers = ?, gene_likelihood = 'zinb')
      also kwargs was set up observing the result

Result(
  metrics={'validation_loss': 5341.3056640625},
  path='/home/prince/PycharmProjects/pancreatic_scrna_seq/code/scvi_log/autotune/2024-09-18_12-31-29_scvi/_trainable_2536f1ac_30_gene_likelihood=zinb,
  lr=0.0014,n_hidden=192,n_latent=10,n_layers=2_2024-09-18_12-53-39',
  filesystem='local',
  checkpoint=None
)

'''

'''
# Step 2: Setup (based on results of tuning) and Training

scvi.model.SCVI.setup_anndata(adata,
                                categorical_covariate_keys = ['Sample'],
                                continuous_covariate_keys=['pct_counts_mt', 'pct_counts_ribo'])
model = scvi.model.SCVI(adata, n_hidden = 192, n_latent = 10, n_layers = 2, gene_likelihood = 'zinb')   # set up accordingly
kwargs = {'lr': 0.0014}

model.train(max_epochs = 200, early_stopping = True, plan_kwargs = kwargs)

model.save(scvi_pdac_model)

y = model.history['reconstruction_loss_validation']['reconstruction_loss_validation'].min()
plt.plot(model.history['reconstruction_loss_train']['reconstruction_loss_train'], label='train')
plt.plot(model.history['reconstruction_loss_validation']['reconstruction_loss_validation'], label='validation')

plt.axhline(y, c = 'k')

plt.legend()
plt.show()

adata.write_h5ad(extra_pdac + 'integrated_pdac.h5ad')
'''

'''
#+++++++++++++++++++++++++++++++++++++++++++++++++ Clustering  +++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# loading preprocessed anndata object
print("reading adata......")
adata = sc.read_h5ad(extra_pdac + 'integrated_pdac.h5ad')

# loading the pretrained scvi model (plus linking to the same anndata object)
print("reading model......")
model = scvi.model.SCVI.load(scvi_pdac_model, adata)

# latent representation is a compressed, low-dimensional space that captures key biological variation in the data
print("obtaining latent representation......")
adata.obsm['X_scVI'] = model.get_latent_representation()
print(adata.obsm['X_scVI'].shape)

#adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size=1e4) 

# clustering and UMAP
print("computing neighbors.....")
sc.pp.neighbors(adata, use_rep = 'X_scVI')  # computes neighborhood graph based on latent representation
print("using leiden for clustering.....")
sc.tl.leiden(adata, resolution = 3, key_added = 'overcluster')  # uses leiden alg to cluster cells based on neighborhood
sc.tl.umap(adata)

print("normalizing and log conversion....")
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

adata.write_h5ad(extra_pdac + 'normalized_clsutered_pdac.h5ad')

#print(adata.obs)
'''

'''
adata.obs['model_1'] = adata.obs.groupby('overcluster')['mod1_label'].transform(lambda x: x.mode()[0])
sc.pl.umap(adata, color = ['model_1'], s = 5, legend_loc='on data')

adata.obs['model_2'] = adata.obs.groupby('overcluster')['mod2_label'].transform(lambda x: x.mode()[0])
sc.pl.umap(adata, color = ['model_2'], s = 5, legend_loc='on data')

adata.obs['model_3'] = adata.obs.groupby('overcluster')['mod3_label'].transform(lambda x: x.mode()[0])
sc.pl.umap(adata, color = ['model_3'], s = 5, legend_loc='on data')

sc.pl.umap(adata, color = ['mod1_score', 'mod2_score', 'mod3_score'])

plt.show()
'''

# +++++++++++++++++++++++++++++++++++++++++++++++++ Annotation  +++++++++++++++++++++++++++++++++++++++++++++++++++++++#
'''
adata = sc.read_h5ad(extra_pdac + 'normalized_clsutered_pdac.h5ad')

# Converting adata.obs to a DataFrame and add the 'status' column
def assign_status(sample_name):
    if 'ADJ' in sample_name:
        return 'Adjunct'
    elif 'PDAC' in sample_name:
        return 'PDAC'
    else:
        return 'Unknown'

# Applying the function to assign status based on the 'sample' column
adata.obs['status'] = adata.obs['Sample'].apply(assign_status)

#print(adata)
#print(adata.obs)

#sc.pl.umap(adata, color = ['overcluster'], legend_loc = 'on data', s = 5)

np.random.seed(1)   # ensure consistent shuffling of cells
ri = np.random.permutation(list(range(adata.shape[0])))

#sc.pl.umap(adata[ri,:], color = ['status'], vmin = .5, size = 2)

labels = adata.obs[['mod1_label', 'mod2_label', 'mod3_label', 'overcluster']].groupby('overcluster').agg(lambda x: x.mode())
scores = adata.obs[['mod1_score', 'mod2_score', 'mod3_score', 'overcluster']].groupby('overcluster').agg(lambda x: x.mean())

mapping_res = labels.merge(right = scores, left_index=True, right_index=True)
#print(mapping_res)
#print(mapping_res[0:10]) # to check on the terminal (change range accordingly)
#mapping_res.to_csv(f'{cluster_results}/default_label.csv', index=True)  # Saves with the index column


# finding marker genes for each cluster
sc.tl.rank_genes_groups(adata, groupby = 'overcluster')
marks = sc.get.rank_genes_groups_df(adata, group = None)
#print("Genes with lfc ranked:\n", marks)

# get an empty dictionary to fill-up cell types later
#for x in range(len(adata.obs.overcluster.unique())):
#    print(f'"{x}":"",')
'''
# individual clustering and top gene observation
'''
# Ensure output directory exists
if not os.path.exists(cluster_results):
    os.makedirs(cluster_results)

# Loop over each unique cluster in 'overcluster'
unique_clusters = adata.obs['overcluster'].unique()

for cluster in unique_clusters:
    # Filter rows for the current cluster
    marks_group = marks[marks['group'] == str(cluster)]

    # Extract top 25 genes
    gene_list = marks_group['names'].head(25).tolist()

    # Create a text file for the current cluster
    with open(f'{cluster_results}cluster_{cluster}_results.txt', 'w') as f:
        f.write(f"Genes with lfc ranked in target cluster {cluster}:\n")
        f.write(str(marks_group[0:25]) + '\n\n')

        # Write the gene list
        f.write(f"Gene list for DAVID analysis (cluster {cluster}):\n")
        for name in gene_list:
            f.write(name + '\n')
        f.write('\n')

    print(f"Gene list for cluster {cluster}:\n", gene_list)

    # Plot and save the UMAP for the current cluster
    ax = sc.pl.umap(adata, palette='blue', show=False)  # Plot full UMAP
    sc.pl.umap(adata[adata.obs['overcluster'] == str(cluster)],
               color='overcluster',
               ax=ax,
               legend_loc=None,
               palette='k',
               show=False)

    # Save the plot as a PNG file
    ax.figure.savefig(f'{cluster_results}cluster_{cluster}.png')
    ax.clear()
'''

# =========================== trials (edit and run for each cluster separately) ==+=====================================#

'''
#marks_d = marks[marks.names == 'IFITM2'].sort_values('logfoldchanges', ascending = False)
#print(marks_d[0:25])  # to observe nearby cluster's expression of the same marker

# markers searching (using the marker genes from https://panglaodb.se/markers.html?cell_type=%27choose%27)
# also use marker genes from http://117.50.127.228/CellMarker/
# https://ngdc.cncb.ac.cn/databasecommons/database/id/406
# also look at DAVID (https://david.ncifcrf.gov/tools.jsp) to get idea on possible gene ontology

fibroblast_markers = ['VIM', 'COL1A2', 'LUM']
stellate_markers = ['COL6A1', 'RGS5', 'THY1', 'MMP11']
cytT_markers = ['TRAC', 'CD8A', 'GZMB'] # GZMB is more of a NK marker
helperT_markers = ['CD4', 'CCR4']
regT_markers = ['IKZF2', 'FOXP3', 'CCR4', 'IL2RA']
gammadeltaT_markers = ['TUBB', 'STMN1']
monocyte_markers = ['CD14', 'LYZ', 'CFP', 'APOBEC3A']
B_markers = ['PXK', 'MS4A1', 'NPIPB15']
tumor_markers = ['KRAS', 'MUC1', 'TP53', 'SMAD4']   # note MUC1 is also an epithelial marker
CAF_markers = ['ACTA2', 'PDGFRA', 'PDGFRB', 'FAP']
macrophage_markers = ['CD68', 'NAAA', 'MARCH1', 'JAML'] # JAML seems more to be with DCs
dendritic_markers = ['ITGAX', 'ZBTB46', 'CD86', 'LAMP3']
mast_markers = ['KIT', 'TPSAB1', 'HDC', 'IL1RL1']
NK_markers = ['TRDC', 'NKG7', 'KLRF1']
NK_T_markers = ['NCAM1', 'GATA3']
neutrophil_markers = ['CSF3R', 'S100A8', 'CTSG']
plasma_markers = ['MZB1', 'IGHG1', 'JCHAIN', 'SPAG4']
epithelial_markers = ['KRT19', 'EPCAM', 'CDH1']
basal_markers = ['KRT5', 'KRT14', 'ITGA3']  # KRT14 might be in both epithelial and basal cells
endothelial_markers = ['CD93', 'VWF', 'EMCN', 'EGFL7']
ductal_markers = ['CFTR', 'SERPINA5', 'SLPI', 'TFF1']
acinar_markers = ['PRSS1', 'KLK1', 'CTRC']
peri_schwann_markers = ['MPZ', 'OLFML2A', 'GULP1', 'GFRA3']
schwann_markers = ['SOX10', 'S100B', 'CRYAB']
alpha_markers = ['GCG', 'TTR', 'PCSK2']

sc.pl.umap(adata, color = alpha_markers, legend_loc ='on data', s = 5)
#sc.pl.umap(adata, color = gammadeltaT_markers, legend_loc ='on data', s = 5)
#sc.pl.umap(adata, color = ['Sample'], vmin = .5, size = 2)

marks_a = marks[marks['names'].isin(alpha_markers)].sort_values('logfoldchanges', ascending=False)
print("With target gene markers over all clusters:\n", marks_a[0:50])

marks_a = marks[marks['names'].isin(epithelial_markers)].sort_values('logfoldchanges', ascending=False)
print("With target gene markers over all clusters:\n", marks_a[0:50])

'''

# =====================================================================================================================#

'''
# Define the marker sets and names
marker_sets = {
    'fibroblast_markers': ['VIM', 'COL1A2', 'LUM'],
    'stellate_markers': ['COL6A1', 'RGS5', 'THY1', 'MMP11'],
    'cytT_markers': ['TRAC', 'CD8A', 'GZMB'],
    'helperT_markers': ['CD4', 'CCR4', 'IL13'],
    'regT_markers': ['IKZF2', 'FOXP3', 'CCR4', 'IL2RA'],
    'gammadeltaT_markers': ['TUBB', 'STMN1'],
    'monocyte_markers': ['CD14', 'LYZ', 'CFP', 'APOBEC3A'],
    'B_markers': ['PXK', 'MS4A1', 'NPIPB15'],
    'tumor_markers': ['KRAS', 'MUC1', 'TP53', 'SMAD4'],
    'CAF_markers': ['ACTA2', 'PDGFRA', 'PDGFRB', 'FAP'],
    'macrophage_markers': ['CD68', 'NAAA', 'MARCH1', 'JAML'],
    'dendritic_markers': ['ITGAX', 'ZBTB46', 'CD86', 'LAMP3'],
    'mast_markers': ['KIT', 'TPSAB1', 'HDC', 'IL1RL1'],
    'NK_markers': ['TRDC', 'NKG7', 'KLRF1'],
    'NK_T_markers': ['NCAM1', 'GATA3'],
    'neutrophil_markers': ['CSF3R', 'S100A8', 'CTSG'],
    'plasma_markers': ['MZB1', 'IGHG1', 'JCHAIN', 'SPAG4'],
    'epithelial_markers': ['KRT19', 'EPCAM', 'CDH1'],
    'basal_markers': ['KRT5', 'KRT14', 'ITGA3'],
    'endothelial_markers': ['CD93', 'VWF', 'EMCN', 'EGFL7'],
    'ductal_markers': ['CFTR', 'SERPINA5', 'SLPI', 'TFF1'],
    'acinar_markers': ['PRSS1', 'KLK1', 'CTRC'],
    'peri_schwann_markers': ['MPZ', 'OLFML2A', 'GULP1', 'GFRA3'],
    'schwann_markers': ['SOX10', 'S100B', 'CRYAB'],
    'alpha_markers': ['GCG', 'TTR', 'PCSK2']
}

# Set output directories
output_file = os.path.join(cluster_results, "marker_analysis.txt")

# Open the output file in append mode
with open(output_file, "a") as f_output:
    # Loop over the marker sets
    for marker_name, marker_list in marker_sets.items():
        # Create and save UMAP plot
        sc.pl.umap(adata, color=marker_list, legend_loc='on data', s=5, show=False)
        plt.savefig(os.path.join(cluster_results, f"z_{marker_name}.png"))
        plt.close()

        # Filter marker genes from the marker list
        marks_a = marks[marks['names'].isin(marker_list)].sort_values('logfoldchanges', ascending=False)

        # Write the top 50 markers to the output file
        f_output.write(f"Top 50 results for {marker_name}:\n\n")
        f_output.write(marks_a[0:50].to_string())
        f_output.write("\n\n\n")  # Three new lines after each marker block

print("Process completed successfully.")

'''

# ======================================================================================================================#

#print("reading adata......")
#adata = sc.read_h5ad(extra_pdac + 'normalized_clsutered_pdac.h5ad')

# after checking individually, the cell types were filled

cell_type = {"0":"Stellates",
"1":"CD8 T",
"2":"B pro",
"3":"CD4 T pro",
"4":"B pro",
"5":"Stellates",
"6":"Macrophages",
"7":"CD4 T pro",
"8":"CD8 T",
"9":"Ductal",
"10":"Endothelial",
"11":"CD8 T",
"12":"Ductal",
"13":"Reg T",
"14":"CD8 T",
"15":"CD8 T",
"16":"CD4 T pro",
"17":"CD4 T",
"18":"Neutrophils",
"19":"Macrophages",
"20":"Reg T",
"21":"Ductal",
"22":"Endothelial",
"23":"Ductal",
"24":"NK",
"25":"CD8 T",
"26":"Stellates",
"27":"Mast",
"28":"CD8 T",
"29":"Macrophages",
"30":"Stellates",
"31":"Macrophages",
"32":"CD8 T",
"33":"CD8 T",
"34":"CD8 T pro",
"35":"Acinar",
"36":"Neutrophils",
"37":"Acinar",
"38":"B pro",
"39":"Stellates",
"40":"Stellates",
"41":"Acinar",
"42":"Fibroblasts",
"43":"Macrophages",
"44":"Schwann",
"45":"Endothelial",
"46":"Ductal",
"47":"Fibroblasts",
"48":"Macrophages",
"49":"Fibroblasts",
"50":"Ductal",
"51":"Gamma-delta T",
"52":"Fibroblasts",
"53":"Acinar",
"54":"Endothelial",
"55":"Acinar",
"56":"Stellates",
"57":"NK T",
"58":"Alpha",
"59":"Ductal",
"60":"B",
"61":"Endothelial",
"62":"Dendritic",
"63":"Macrophages",
"64":"Acinar",
"65":"Dendritic",
"66":"Basal epithelial",
"67":"Plasma",
"68":"Macrophages",
"69":"Fibroblasts"}

'''
# Define a color palette
color_palette = sns.color_palette("tab20", n_colors=len(cell_type))

#annotated umap generation
adata.obs['cell_type'] = adata.obs.overcluster.map(cell_type)

# Generating the UMAP plot
sc.pl.umap(adata, color=['cell_type'], palette=color_palette, size=20,
           alpha = 0.8, legend_loc = 'on data')
'''

# broad resolution annotation

cell_type_broad = {"0":"Stellates",
"1":"T cells",
"2":"B cells",
"3":"T cells",
"4":"B cells",
"5":"Stellates",
"6":"Macrophages",
"7":"T cells",
"8":"T cells",
"9":"Ductal",
"10":"Endothelial",
"11":"T cells",
"12":"Ductal",
"13":"T cells",
"14":"T cells",
"15":"T cells",
"16":"T cells",
"17":"T cells",
"18":"Neutrophils",
"19":"Macrophages",
"20":"T cells",
"21":"Ductal",
"22":"Endothelial",
"23":"Ductal",
"24":"NK cells",
"25":"T cells",
"26":"Stellates",
"27":"Mast cells",
"28":"T cells",
"29":"Macrophages",
"30":"Stellates",
"31":"Macrophages",
"32":"T cells",
"33":"T cells",
"34":"T cells",
"35":"Acinar",
"36":"Neutrophils",
"37":"Acinar",
"38":"B cells",
"39":"Stellates",
"40":"Stellates",
"41":"Acinar",
"42":"Fibroblasts",
"43":"Macrophages",
"44":"Schwann cells",
"45":"Endothelial",
"46":"Ductal",
"47":"Fibroblasts",
"48":"Macrophages",
"49":"Fibroblasts",
"50":"Ductal",
"51":"T cells",
"52":"Fibroblasts",
"53":"Acinar",
"54":"Endothelial",
"55":"Acinar",
"56":"Stellates",
"57":"T cells",
"58":"Alpha cells",
"59":"Ductal",
"60":"B cells",
"61":"Endothelial",
"62":"Dendritic",
"63":"Macrophages",
"64":"Acinar",
"65":"Dendritic",
"66":"Basal epithelial",
"67":"Plasma cells",
"68":"Macrophages",
"69":"Fibroblasts"}

'''
# Define a color palette
color_palette = sns.color_palette("tab20", n_colors=len(cell_type))

# Annotate the cell types
adata.obs['cell_type_broad'] = adata.obs.overcluster.map(cell_type_broad)

# Generating the UMAP plot without the legend
sc.pl.umap(adata, color=['cell_type_broad'], palette = color_palette,
           size=20, legend_loc='on data')  # Disable show to customize the legend
'''


'''
# Creating the legend below the plot
fig, ax = plt.gcf(), plt.gca()
handles, labels = ax.get_legend_handles_labels()

# Removing the default legend
ax.legend_.remove()

# Adding a new legend below the plot
fig.legend(handles, labels, loc='lower center', ncol=4, fontsize='small')

# Adjusting the space for the legend to fit below the plot
plt.subplots_adjust(bottom=0.2)
plt.show()
'''


'''
# saving the gene marks in unstructured (so that if needed they can be accessed later easily)
adata.uns['rank_genes_groups'] = marks

adata.write_h5ad(extra_pdac + 'annotated_pdac.h5ad')

print(adata)

'''

# +++++++++++++++++++++++++++++++++++++++++++++++++ Analysis  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# loading preprocessed anndata object
print("reading adata......")
adata = sc.read_h5ad(extra_pdac + 'annotated_pdac.h5ad')

# loading the pretrained scvi model (plus linking to the same anndata object) (?? not needed for analysis ig)
print("reading model......")
model = scvi.model.SCVI.load(scvi_pdac_model, adata)

print(adata)
#print(dir(model))

#GOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOD
# 1. finding frequency (need total cell counts and count by cell type):-
'''
# counting total cells per sample
num_tot_cells = adata.obs.groupby(['Sample']).count()   # finds the number of cells in a given sample
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.doublet))
print(num_tot_cells)

# counting cell types per sample and condition (options: cell_type_broad or cell_type)
cell_type_counts = adata.obs.groupby(['Sample', 'status', 'cell_type']).count()  # finds number of cells with this unique combination
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()   # keeps only the rows with non-zero count and converts the grouped dataframe back to regular
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:4]]  # selects the first 4 columns
cell_type_counts = cell_type_counts.rename(columns={'n_genes': 'cell_count_per_cell_type'})   # renaming column for user viewers convenience
#print(cell_type_counts)

# calculating frequency
cell_type_counts['total_cell_count'] = cell_type_counts.Sample.map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.cell_count_per_cell_type / cell_type_counts.total_cell_count

print(cell_type_counts)

# saving the dataframe of frequency (highres or broad)
output_file = analysis_pdac + 'frequency_data.csv'
cell_type_counts.to_csv(output_file, index=False)
print(f"Frequency data saved to {output_file}")

# plot forming for frequency (options: cell_type_broad or cell_type; depending on per sample change above)
plt.figure(figsize = (7, 7))

custom_palette = {
    'Adjunct': '#7fffd4',  # aquamarine
    'PDAC': '#b94e48',  # deep chestnut
    }
ax = sns.boxplot(data = cell_type_counts, x = 'cell_type', y = 'frequency',
                 hue = 'status', palette=custom_palette)

unique_cell_types = adata.obs['cell_type'].unique()

for cell_type in unique_cell_types:
    # Filter the data for the specific cell type
    temp = cell_type_counts[cell_type_counts['cell_type'] == cell_type]

    # Extract the data for each condition
    a = temp[temp['status'] == 'Adjunct']['frequency']
    b = temp[temp['status'] == 'PDAC']['frequency']

    # Performing Mann-Whitney U test (non-parametric test)
    if len(a) > 0 and len(b) > 0:  # Check if both groups are non-empty
        mannwhitney_result = stats.mannwhitneyu(a, b)
        p_value = mannwhitney_result.pvalue
        statistic = mannwhitney_result.statistic

        # Annotating the p-value on the plot
        # Adjusting the coordinates (x, y) based on your specific plot layout
        plt.text(x=cell_type, y=max(temp['frequency']) + 0.005,  # Adjust y-position as needed
                 s=f'p={p_value:.3f}', ha='center', va='bottom', rotation=90)

plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')   # rotates the x-axis labels to avoid overlapping
plt.tight_layout()
plt.show()
'''


# 2. top differentially expressed genes in PDAC vs Adjunct

#print("Obtaining normalized data.....")
#adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size=1e4)    # no need to do this if already done while clustering
#adata.write_h5ad(extra_pdac + 'annotated_pdac.h5ad')
'''
np.random.seed(42)

# subsetting Adjunct and PDAC cell types as a copy from adata
subset = adata[adata.obs['status'].isin(['Adjunct', 'PDAC'])].copy()

# converting to dense matrix
subset.X = subset.X.toarray()

# filtering based on cells
sc.pp.filter_genes(subset, min_cells=100)

# scvi training for DE
scvi_de = model.differential_expression(
    idx1 = [adata.obs['status'] == 'PDAC'], # test here
    idx2 = [adata.obs['status'] == 'Adjunct']   # put reference here
    )

# filtering of SCVI trained DE data
scvi_de = scvi_de[(scvi_de['is_de_fdr_0.05']) & (abs(scvi_de.lfc_mean) > .5)]   # false discovery rate and log fold changes considered for de analysis
scvi_de = scvi_de.sort_values('lfc_mean')
scvi_de = scvi_de[(scvi_de.raw_normalized_mean1 > .5) | (scvi_de.raw_normalized_mean2 > .5)]
#print(scvi_de)

# sorting top DEGs and plotting in heatmap
diff_genes_overall = scvi_de[-100:].index.tolist() + scvi_de[:100].index.tolist() #top 25 and bottom 25 from sorted df
sc.pl.heatmap(subset, diff_genes_overall, groupby='status', swap_axes=True, show_gene_labels=True,
              layer ='scvi_normalized', log = True)  # tried with both layer='counts' and layer='scvi_normalized'



# 3. Functional Enrichment [GO enrichment and Pathway enrichment]

# Upregulated genes (lfc_mean > 0)
upregulated_genes = scvi_de[scvi_de['lfc_mean'] > 0]
upregulated_gene_names = upregulated_genes.index.tolist()

# Downregulated genes (lfc_mean < 0)
downregulated_genes = scvi_de[scvi_de['lfc_mean'] < 0]
downregulated_gene_names = downregulated_genes.index.tolist()

# Perform enrichment for upregulated genes
enr_up = gp.enrichr(gene_list=upregulated_gene_names,
                    gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2023', 'Cancer_Cell_Line_Encyclopedia', 'DisGeNET'],
                    organism='human',
                    outdir=None,
                    background=subset.var_names.tolist())

# Perform enrichment for downregulated genes
enr_down = gp.enrichr(gene_list=downregulated_gene_names,
                      gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2023', 'Cancer_Cell_Line_Encyclopedia', 'DisGeNET'],
                      organism='human',
                      outdir=None,
                      background=subset.var_names.tolist())

# Print and save results for upregulated genes
print("Upregulated Genes Enrichment Results")
print(enr_up.results)
enr_up.results.to_csv(extra_pdac + '/upregulated_enrichment.csv', index=False)

# Print and save results for downregulated genes
print("Downregulated Genes Enrichment Results")
print(enr_down.results)
enr_down.results.to_csv(extra_pdac + '/downregulated_enrichment.csv', index=False)
'''



''' 
# 2. top differentially expressed genes in T cells of PDAC vs Adjunct (???????????????????????????????????????)

#print("Obtaining normalized data.....")
#adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size=1e4)    # no need to do this if already done while clustering
#adata.write_h5ad(extra_pdac + 'annotated_pdac.h5ad')

np.random.seed(42)

# subsetting ajunct vs padc cell types as a copy from adata
subset = adata[adata.obs['status'].isin(['Adjunct', 'PDAC'])].copy()

# converting to dense matrix
subset.X = subset.X.toarray()

# filtering based on cells
sc.pp.filter_genes(subset, min_cells=100)

# scvi training for DE
scvi_de = model.differential_expression(
     idx1 = [(subset.obs.status == 'PDAC') & (adata.obs.cell_type_broad == 'T cells')],
     idx2 = [(subset.obs.status == 'Adjunct') & (adata.obs.cell_type_broad == 'T cells')]
     )

# filtering of SCVI trained DE data
scvi_de = scvi_de[(scvi_de['is_de_fdr_0.05']) & (abs(scvi_de.lfc_mean) > .5)]   # false discovery rate and log fold changes considered for de analysis
scvi_de = scvi_de.sort_values('lfc_mean')
scvi_de = scvi_de[(scvi_de.raw_normalized_mean1 > .5) | (scvi_de.raw_normalized_mean2 > .5)]
#print(scvi_de)

# sorting top DEGs and plotting in heatmap
diff_genes_overall = scvi_de[-30:].index.tolist() + scvi_de[:30].index.tolist() #top 25 and bottom 25 from sorted df
sc.pl.heatmap(subset, diff_genes_overall, groupby='status', swap_axes=True, show_gene_labels=True,
              layer ='scvi_normalized', log = True)  # tried with both layer='counts' and layer='scvi_normalized'

'''


#unique_cell_types = adata.obs['cell_type'].unique()
#print(unique_cell_types)


 # GOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOD
'''
# Dictionary to store top genes
top_genes_per_cell_type = {}

# Loop through each cell type category to find top genes
for cell_type in adata.obs['cell_type'].cat.categories:
    scvi_de = model.differential_expression(
        groupby="cell_type",
        group1=cell_type,
    )
    # Manually select top genes if necessary
    top_genes = scvi_de.sort_values(by='lfc_mean', ascending=False).head(5)
    top_genes_per_cell_type[cell_type] = top_genes.index.tolist()

# Print top genes for each cell type
for cell_type, genes in top_genes_per_cell_type.items():
    print(f"Top genes for {cell_type}: {genes}")

# Define the list of top genes for all cell types (combine top 5 genes per cell type)
top_genes = [gene for genes in top_genes_per_cell_type.values() for gene in genes]

# Plot heatmap for the top genes across all cell types
sc.pl.heatmap(adata, var_names=top_genes, groupby='cell_type', swap_axes=True, show_gene_labels=True, layer='scvi_normalized')

# Create a dot plot
sc.pl.dotplot(adata, var_names=top_genes, groupby='cell_type', layer='scvi_normalized', standard_scale='var')
'''

#++++++++++++++++++++++++++++++++++++++++++++++++++++ Enrichment analysis +++++++++++++++++++++++++++++++++++++++++++++#

'''
# Define immune cell types
immune_cell_types = ['T cells', 'B cells', 'NK cells', 'Macrophages', 'Neutrophils', 'Mast cells']  # Adjust based on your annotations

# Subset the adata object to include only immune cells using the broad cell type
adata_immune = adata[adata.obs['cell_type_broad'].isin(immune_cell_types)]

# Define the statuses to compare
statuses = ['PDAC', 'Adjunct']  # Adjust based on your actual status labels

# Dictionary to store top genes for each status comparison
top_genes_per_status = {}

# Loop through each status to find top genes
for status in statuses:
    scvi_de = model.differential_expression(
        groupby="status",  # Change to your column name for status
        group1=status,
    )
    # Select the top 100 genes based on log fold change
    top_genes = scvi_de.sort_values(by='lfc_mean', ascending=False).head(100)
    top_genes_per_status[status] = top_genes.index.tolist()

# Print top genes for each status
for status, genes in top_genes_per_status.items():
    print(f"Top genes for {status}: {genes}")

# Automate the extraction of upregulated and downregulated genes
upregulated_gene_names = top_genes_per_status['PDAC']  # Genes upregulated in PDAC
downregulated_gene_names = top_genes_per_status['Adjunct']  # Genes upregulated in Adjunct

# Print the gene lists for verification
print(f"Upregulated genes: {upregulated_gene_names}")
print(f"Downregulated genes: {downregulated_gene_names}")

# Define the list of top genes for all statuses (combine top 100 genes)
top_genes = [gene for genes in top_genes_per_status.values() for gene in genes]

# Plot heatmap for the top genes across statuses in immune cells
sc.pl.heatmap(adata_immune, var_names=top_genes, groupby='status',
               swap_axes=True, show_gene_labels=True, layer='scvi_normalized')

# Create a dot plot for the immune cells by status
sc.pl.dotplot(adata_immune, var_names=top_genes, groupby='status', layer='scvi_normalized', standard_scale='var')


# Perform enrichment for upregulated genes
enr_up = gp.enrichr(gene_list=upregulated_gene_names,
                    gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2023', 'BioCarta_2013',
                               'Cancer_Cell_Line_Encyclopedia', 'DisGeNET'],
                    organism='human',
                    outdir=None,
                    background=adata_immune.var_names.tolist())

# Perform enrichment for downregulated genes
enr_down = gp.enrichr(gene_list=downregulated_gene_names,
                      gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2023', 'BioCarta_2013',
                                 'Cancer_Cell_Line_Encyclopedia', 'DisGeNET'],
                      organism='human',
                      outdir=None,
                      background=adata_immune.var_names.tolist())

enr_up.results.to_csv(extra_pdac + '/upregulated_immune_pdac_enrichment.csv', index=False)
enr_down.results.to_csv(extra_pdac + '/downregulated_immune_pdac_enrichment.csv', index=False)
'''

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


#+++++++++++++++++++++++++++++++++++++++++++ PLOTTING OF ENRICHED DATA +++++++++++++++++++++++++++++++++++++++++++++++#
'''
upregulated_enrichment = pd.read_csv(extra_pdac + '/upregulated_immune_pdac_enrichment.csv')
downregulated_enrichment = pd.read_csv(extra_pdac + '/downregulated_immune_pdac_enrichment.csv')
# Define a function to plot enrichment results for a given gene set
def plot_enrichment_results(enrichment_df, title):
    # Set the plot style
    sns.set(style='whitegrid')

    # Filter for significant pathways (optional)
    significant_terms = enrichment_df[enrichment_df['Adjusted P-value'] < 0.05]

    # Plotting each pathway in separate figures
    for gene_set in significant_terms['Gene_set'].unique():
        pathway_data = significant_terms[significant_terms['Gene_set'] == gene_set].sort_values(by='Adjusted P-value')

        plt.figure(figsize=(10, 6))
        sns.barplot(x='Adjusted P-value', y='Term', data=pathway_data, color='skyblue')
        plt.title(f'{title} - {gene_set}')
        plt.xlabel('Adjusted P-value')
        plt.ylabel('Enriched Term')
        plt.tight_layout()
        plt.show()

# Plot for upregulated genes
plot_enrichment_results(upregulated_enrichment, 'Enrichment Results for Upregulated Genes')

# Plot for downregulated genes
plot_enrichment_results(downregulated_enrichment, 'Enrichment Results for Downregulated Genes')
'''


'''
upregulated_enrichment = pd.read_csv(extra_pdac + '/upregulated_enrichment.csv')
downregulated_enrichment = pd.read_csv(extra_pdac + '/downregulated_enrichment.csv')

# Define a function to plot enrichment results for a given gene set
def plot_enrichment_results(enrichment_df, title, pval_threshold=0.05):
    # Set the plot style
    sns.set(style='whitegrid')

    # Filter for significant pathways based on adjusted p-value
    significant_terms = enrichment_df[enrichment_df['Adjusted P-value'] < pval_threshold]

    # Plotting each pathway in separate figures
    for gene_set in significant_terms['Gene_set'].unique():
        pathway_data = significant_terms[significant_terms['Gene_set'] == gene_set].sort_values(by='Adjusted P-value')

        plt.figure(figsize=(10, 6))
        sns.barplot(x='Adjusted P-value', y='Term', data=pathway_data, color='red')
        plt.title(f'{title} - {gene_set} (Significant Results)')
        plt.xlabel('Adjusted P-value')
        plt.ylabel('Enriched Term')
        plt.tight_layout()
        plt.show()

# Plot for upregulated genes
plot_enrichment_results(upregulated_enrichment, 'Enrichment Results for Upregulated Genes')

# Plot for downregulated genes
plot_enrichment_results(downregulated_enrichment, 'Enrichment Results for Downregulated Genes')
'''


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


'''
# Dictionary to store top genes for T cells
top_genes_per_immune_cells = {}

# Define T cell populations you want to analyze
immune_cell_types = ['B cells', 'Dendritic', 'T cells', 'NK cells', 'Mast cells', 'Macrophages', 'Neutrophils', 'Plasma cells']

# Loop through each T cell type to find top genes across PDAC vs. Adjunct
for t_cell in immune_cell_types:
    scvi_de = model.differential_expression(
        groupby="cell_type_broad",
        group1="cell_type_broad"
#        group1="PDAC",
#        group2="Adjunct"
    )
    # Manually select top genes based on log fold change
    top_genes = scvi_de.sort_values(by='lfc_mean', ascending=False).head(10)
    top_genes_per_immune_cells[t_cell] = top_genes.index.tolist()

# Print top genes for each T cell type
for t_cell, genes in top_genes_per_immune_cells.items():
    print(f"Top genes for {t_cell}: {genes}")

# Combine top genes for plotting
top_genes_t_cells = [gene for genes in top_genes_per_immune_cells.values() for gene in genes]

# Plot heatmap for the top genes across T cell types
sc.pl.heatmap(adata, var_names=top_genes_t_cells, groupby='status', swap_axes=True, show_gene_labels=True, layer='scvi_normalized')

# Create a dot plot for the top genes
sc.pl.dotplot(adata, var_names=top_genes_t_cells, groupby='status', layer='scvi_normalized', standard_scale='var')
'''



'''
# Dictionary to store top genes for immune cells
top_genes_per_immune_cells = {}

# Define immune cell populations you want to analyze
immune_cell_types = ['B cells', 'Dendritic', 'T cells', 'NK cells', 'Mast cells', 'Macrophages', 'Neutrophils', 'Plasma cells']

# Loop through each immune cell type to find top genes
for immune_cell in immune_cell_types:
    scvi_de = model.differential_expression(
        groupby="cell_type_broad",  # Now grouping by broad immune cell types
        group1=immune_cell,  # Compare each immune cell type
    )
    # Select top genes based on log fold change
    top_genes = scvi_de.sort_values(by='lfc_mean', ascending=False).head(10)
    top_genes_per_immune_cells[immune_cell] = top_genes.index.tolist()

# Print top genes for each immune cell type
for immune_cell, genes in top_genes_per_immune_cells.items():
    print(f"Top genes for {immune_cell}: {genes}")

# Combine top genes for plotting
top_genes_immune_cells = [gene for genes in top_genes_per_immune_cells.values() for gene in genes]

# Plot heatmap for the top genes across immune cell types
sc.pl.heatmap(adata, var_names=top_genes_immune_cells, groupby='cell_type_broad', swap_axes=True, show_gene_labels=True, layer='scvi_normalized')

# Create a dot plot for the top genes
sc.pl.dotplot(adata, var_names=top_genes_immune_cells, groupby='cell_type_broad', layer='scvi_normalized', standard_scale='var')

'''



'''
# subsetting Adjunct and PDAC cell types as a copy from adata
subset = adata[adata.obs['cell_type_broad'].isin(['Adjunct', 'PDAC'])].copy()

# converting to dense matrix
subset.X = subset.X.toarray()

# filtering based on cells
sc.pp.filter_genes(subset, min_cells=100)

# scvi training for DE
scvi_de = model.differential_expression(
     idx1 = [(subset.obs.status == 'PDAC') & (adata.obs.cell_type_broad == 'T cells')],
     idx2 = [(subset.obs.status == 'Adjunct') & (adata.obs.cell_type_broad == 'T cells')]
     )

# filtering of SCVI trained DE data
scvi_de = scvi_de[(scvi_de['is_de_fdr_0.05']) & (abs(scvi_de.lfc_mean) > .5)]   # false discovery rate and log fold changes considered for de analysis
scvi_de = scvi_de.sort_values('lfc_mean')
scvi_de = scvi_de[(scvi_de.raw_normalized_mean1 > .5) | (scvi_de.raw_normalized_mean2 > .5)]
#print(scvi_de)

# sorting top DEGs and plotting in heatmap
diff_genes = scvi_de[-30:].index.tolist() + scvi_de[:30].index.tolist() #top 25 and bottom 25 from sorted df
sc.pl.heatmap(subset, diff_genes, groupby='status', swap_axes=True, show_gene_labels=True,
              layer ='scvi_normalized', log = True)  # tried with both layer='counts' and layer='scvi_normalized'

# Create a dot plot for the top genes
sc.pl.dotplot(adata, var_names=diff_genes, groupby='cell_type_broad', layer='scvi_normalized', standard_scale='var')


'''


'''
# Filter to keep only immune cells in the original AnnData
print("Filtering for immune cells...")
immune_cells = adata[adata.obs['cell_type_broad'].isin(['B cells', 'T cells', 'NK cells', 'Dendritic', 'Macrophages', 'Plasma cells', 'Mast cells'])].copy()

# Check the size of the filtered immune cells
print(f"Subset size: {immune_cells.n_obs} observations and {immune_cells.n_vars} variables")

# Create boolean masks for PDAC and Adjunct using the original AnnData object
pdac_mask = adata.obs['status'] == 'PDAC'
adjunct_mask = adata.obs['status'] == 'Adjunct'

# Ensure the masks are correct
print(f'PDAC cells: {pdac_mask.sum()}')
print(f'Adjunct cells: {adjunct_mask.sum()}')

# Perform differential expression analysis using original indices
print("Performing differential expression analysis...")
scvi_de = model.differential_expression(
    idx1=pdac_mask,
    idx2=adjunct_mask
)

# Filter the DE results
scvi_de = scvi_de[(scvi_de['is_de_fdr_0.05']) & (abs(scvi_de.lfc_mean) > 0.5)]
scvi_de = scvi_de.sort_values('lfc_mean')
scvi_de = scvi_de[(scvi_de.raw_normalized_mean1 > 0.5) | (scvi_de.raw_normalized_mean2 > 0.5)]

# Sorting top DEGs and plotting in heatmap
diff_genes = scvi_de[-50:].index.tolist() + scvi_de[:50].index.tolist()  # Top 30 and bottom 30 from sorted df
sc.pl.heatmap(immune_cells, diff_genes, groupby='status', swap_axes=True, show_gene_labels=True,
               layer='scvi_normalized', log=True)

# Create a dot plot for the top genes
sc.pl.dotplot(immune_cells, var_names=diff_genes, groupby='cell_type_broad', layer='scvi_normalized', standard_scale='var')

print("Analysis complete.")
'''


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


immune_marker_sets = {
    'cytT_markers': ['TRAC', 'CD8A', 'GZMB'],
    'helperT_markers': ['CD4', 'CCR4'],
    'regT_markers': ['IKZF2', 'FOXP3', 'CCR4', 'IL2RA'],
    'B_markers': ['PXK', 'MS4A1', 'NPIPB15'],
    'macrophage_markers': ['CD68', 'NAAA', 'MARCH1', 'JAML'],
    'dendritic_markers': ['ITGAX', 'ZBTB46', 'CD86', 'LAMP3'],
    'mast_markers': ['KIT', 'TPSAB1', 'HDC', 'IL1RL1'],
    'NK_markers': ['TRDC', 'NKG7', 'KLRF1'],
    'NK_T_markers': ['NCAM1', 'GATA3'],
    'neutrophil_markers': ['CSF3R', 'S100A8', 'CTSG'],
}




# PDAC vs Adj (not so significant)
'''
# Open a PDF file to save the plots
with PdfPages('immune_marker_plots.pdf') as pdf:
    # Loop through each marker set and create a separate violin plot
    for marker_set, markers in immune_marker_sets.items():
        # Extract the expression data for the markers
        if adata.raw is not None:
            data = adata.raw[:, markers].X
        else:
            data = adata[:, markers].X

        # Convert to a DataFrame
        if isinstance(data, np.ndarray):
            data = pd.DataFrame(data, columns=markers)
        else:
            data = pd.DataFrame(data.toarray(), columns=markers)  # In case of sparse matrix

        data['status'] = adata.obs['status'].values

        # Melt the DataFrame for seaborn
        melted_data = data.melt(id_vars='status', var_name='marker', value_name='expression')

        # Create a new figure for each marker set
        plt.figure(figsize=(10, 6))  # Adjust size as needed

        # Create the violin plot
        sns.violinplot(x='marker', y='expression', hue='status', data=melted_data,
                       split=True, inner=None)  # Excluded jitter and mean overlay

        plt.title(marker_set)
        plt.xlabel('')
        plt.ylabel('Expression')

        # Rotate x-axis labels if necessary
        plt.xticks(rotation=45)

        # Save the current figure into the PDF
        pdf.savefig()  # Save the current figure into the PDF
        plt.close()  # Close the figure to avoid display

print("Plots saved to immune_marker_plots.pdf")
'''



''' # GOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOD
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Open a PDF file to save the heatmaps
with PdfPages('immune_marker_heatmaps_by_cell_type.pdf') as pdf:
    # Loop through each marker set and create a heatmap
    for marker_set, markers in immune_marker_sets.items():
        # Extract the expression data for the markers
        if adata.raw is not None:
            data = adata.raw[:, markers].X
        else:
            data = adata[:, markers].X

        # Convert to a DataFrame
        if isinstance(data, np.ndarray):
            data = pd.DataFrame(data, columns=markers)
        else:
            data = pd.DataFrame(data.toarray(), columns=markers)  # In case of sparse matrix

        data['cell_type'] = adata.obs['cell_type'].values  # Use appropriate column for cell types

        # Melt the DataFrame for seaborn
        melted_data = data.melt(id_vars='cell_type', var_name='marker', value_name='expression')

        # Aggregate the data to ensure no duplicates
        aggregated_data = melted_data.groupby(['cell_type', 'marker']).mean().reset_index()

        # Pivot the data for heatmap
        heatmap_data = aggregated_data.pivot('marker', 'cell_type', 'expression')

        # Create a new figure for each marker set
        plt.figure(figsize=(12, 8))  # Adjust size as needed

        # Create the heatmap
        sns.heatmap(heatmap_data, cmap='viridis', annot=True, fmt='.1f', cbar=True, linewidths=0.5)

        plt.title(f'Heatmap of Marker Expression by Cell Type - {marker_set}')
        plt.xlabel('Cell Types')
        plt.ylabel('Markers')

        # Save the current figure into the PDF
        pdf.savefig()  # Save the current figure into the PDF
        plt.close()  # Close the figure to avoid display

print("Heatmaps saved to immune_marker_heatmaps_by_cell_type.pdf")
'''





#GOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOD



'''
# Creating dotplot for anna
dotplot = sc.pl.dotplot(adata, var_names=immune_marker_sets, groupby='cell_type',
                        var_group_rotation= 80,  # Rotate the variable (gene) names by 45 degrees
                        figsize=(14, 10))
'''

# Define the list of cell types to exclude
exclude_list = ['Plasma cells', 'Fibroblasts', 'Basal epithelial', 'Schwann cells', 'Alpha cells']  # add more cell types as needed

# Filter the adata by excluding cells with the specified types in `exclude_list`
filtered_adata = adata[~adata.obs['cell_type_broad'].isin(exclude_list)].copy()

marker_sets_filtered = {
    'acinar_markers': ['PRSS1', 'KLK1', 'CTRC'],
    'B_markers': ['PXK', 'MS4A1', 'NPIPB15'],
    'dendritic_markers': ['ITGAX', 'ZBTB46', 'CD86', 'LAMP3'],
    'ductal_markers': ['CFTR', 'SERPINA5', 'SLPI', 'TFF1'],
    'endothelial_markers': ['CD93', 'VWF', 'EMCN', 'EGFL7'],
    'macrophage_markers': ['CD68', 'NAAA', 'MARCH1', 'JAML'],
    'mast_markers': ['KIT', 'TPSAB1', 'HDC', 'IL1RL1'],
    'NK_markers': ['TRDC', 'NKG7', 'KLRF1'],
    'neutrophil_markers': ['CSF3R', 'S100A8', 'CTSG'],
    'stellate_markers': ['COL6A1', 'RGS5', 'THY1', 'MMP11'],
    'T_markers': ['TRAC', 'CD8A', 'CD3D', 'TRBC2'],
}


# Creating dotplot for the filtered adata
dotplot = sc.pl.dotplot(
    filtered_adata,
    var_names=marker_sets_filtered,
    groupby='cell_type_broad',
    var_group_rotation=70,  # Rotate the variable (gene) names by 80 degrees
    figsize=(30, 12),
#    dendrogram=True
)



'''
# ++++++++++++++++++++++ DEG between immune cell types of PDAC vs Adj +++++++++++++++++++++++++++++++++++++++++++++++#

## T cells

# Subset T cells from the entire dataset based on the cell type annotation
t_cells = adata[adata.obs['cell_type_broad'] == 'T cells'].copy()

# Ensure that both Adjunct and PDAC statuses are present
print(t_cells.obs['status'].value_counts())

# Subset T cells into Adjunct and PDAC samples
t_cells_adjunct = t_cells[t_cells.obs['status'] == 'Adjunct'].copy()
t_cells_pdac = t_cells[t_cells.obs['status'] == 'PDAC'].copy()

# Compare expression of key T cell markers (e.g., CD3D, CD4, CD8A) between Adjunct and PDAC
markers = ['CD3D', 'CD4', 'CD8A', 'TRAC', 'FOXP3', 'CCR4']

# Create subplots for each marker
n = len(markers)
fig, axes = plt.subplots(n, 1, figsize=(8, 4 * n), sharex=True)

# Iterate through markers to plot and calculate p-values
for i, marker in enumerate(markers):
    # Plotting the violin plot
    sc.pl.violin(t_cells, marker, groupby='status', ax=axes[i], show=False)

    # Finding corresponding p-value for significance test
    a = t_cells[t_cells.obs.status == 'PDAC'].X[:, t_cells.var_names.get_loc(marker)].toarray().flatten()  # Convert sparse matrix to dense array
    b = t_cells[t_cells.obs.status == 'Adjunct'].X[:, t_cells.var_names.get_loc(marker)].toarray().flatten()

    # Perform Mann-Whitney U test
    mannwhitney_result = stats.mannwhitneyu(a, b)

    # Add p-value and statistic to the plot
    p_value = mannwhitney_result.pvalue
    statistic = mannwhitney_result.statistic
    axes[i].text(0.5, 0.95, f'Mann-Whitney U test\nStatistic: {statistic:.2f}\nP-value: {p_value:.2e}',
                  ha='center', va='top', fontsize=12, transform=axes[i].transAxes)

# Adjust layout
plt.tight_layout()
plt.show()

# Perform differential gene expression analysis between Adjunct and PDAC T cells
sc.tl.rank_genes_groups(t_cells, groupby='status', reference='Adjunct', method='wilcoxon')

# View the top differentially expressed genes
sc.pl.rank_genes_groups(t_cells, n_genes=50, sharey=False)

# Flatten the structured array to get a 1D list of gene names
top_genes = [gene[0] for gene in t_cells.uns['rank_genes_groups']['names'][:50]]
print(top_genes)  # Check the extracted top genes

# Generate dot plot for the top 50 differentially expressed genes
sc.pl.dotplot(t_cells, var_names=top_genes, groupby='status', dendrogram=False)
plt.show()

# Optionally, create a heatmap of the top differentially expressed genes
sc.pl.heatmap(t_cells, var_names=top_genes, groupby='status', cmap='viridis', standard_scale='var')
plt.show()
'''


'''
# Subset B cells from the entire dataset based on the cell type annotation
b_cells = adata[adata.obs['cell_type_broad'] == 'B cells'].copy()

# Ensure that both Adjunct and PDAC statuses are present
print(b_cells.obs['status'].value_counts())

# Subset T cells into Adjunct and PDAC samples
b_cells_adjunct = b_cells[b_cells.obs['status'] == 'Adjunct'].copy()
b_cells_pdac = b_cells[b_cells.obs['status'] == 'PDAC'].copy()

# Compare expression of 
markers = ['PXK', 'MS4A1', 'NPIPB15']

# Create subplots for each marker
n = len(markers)
fig, axes = plt.subplots(n, 1, figsize=(8, 4 * n), sharex=True)

# Iterate through markers to plot and calculate p-values
for i, marker in enumerate(markers):
    # Plotting the violin plot
    sc.pl.violin(b_cells, marker, groupby='status', ax=axes[i], show=False)

    # Finding corresponding p-value for significance test
    a = b_cells[b_cells.obs.status == 'PDAC'].X[:, b_cells.var_names.get_loc(marker)].toarray().flatten()  # Convert sparse matrix to dense array
    b = b_cells[b_cells.obs.status == 'Adjunct'].X[:, b_cells.var_names.get_loc(marker)].toarray().flatten()

    # Perform Mann-Whitney U test
    mannwhitney_result = stats.mannwhitneyu(a, b)

    # Add p-value and statistic to the plot
    p_value = mannwhitney_result.pvalue
    statistic = mannwhitney_result.statistic
    axes[i].text(0.5, 0.95, f'Mann-Whitney U test\nStatistic: {statistic:.2f}\nP-value: {p_value:.2e}',
                  ha='center', va='top', fontsize=12, transform=axes[i].transAxes)

# Adjust layout
plt.tight_layout()
plt.show()

# Perform differential gene expression analysis between Adjunct and PDAC b cells
sc.tl.rank_genes_groups(b_cells, groupby='status', reference='Adjunct', method='wilcoxon')

# View the top differentially expressed genes
sc.pl.rank_genes_groups(b_cells, n_genes=50, sharey=False)

# Flatten the structured array to get a 1D list of gene names
top_genes = [gene[0] for gene in b_cells.uns['rank_genes_groups']['names'][:50]]
print(top_genes)  # Check the extracted top genes

# Generate dot plot for the top 50 differentially expressed genes
sc.pl.dotplot(b_cells, var_names=top_genes, groupby='status', dendrogram=False)
plt.show()

# Optionally, create a heatmap of the top differentially expressed genes
sc.pl.heatmap(b_cells, var_names=top_genes, groupby='status', cmap='viridis', standard_scale='var')
plt.show()
'''



'''
# Subset macrophage from the entire dataset based on the cell type annotation
macrophage = adata[adata.obs['cell_type_broad'] == 'Macrophages'].copy()

# Ensure that both Adjunct and PDAC statuses are present
print(macrophage.obs['status'].value_counts())

# Subset given cells into Adjunct and PDAC samples
macrophage_adjunct = macrophage[macrophage.obs['status'] == 'Adjunct'].copy()
macrophage_pdac = macrophage[macrophage.obs['status'] == 'PDAC'].copy()

# Compare expression of key markers between Adjunct and PDAC
markers = ['CD68', 'MARCH1', 'LYZ']

# Create subplots for each marker
n = len(markers)
fig, axes = plt.subplots(n, 1, figsize=(8, 4 * n), sharex=True)

# Iterate through markers to plot and calculate p-values
for i, marker in enumerate(markers):
    # Plotting the violin plot
    sc.pl.violin(macrophage, marker, groupby='status', ax=axes[i], show=False)

    # Finding corresponding p-value for significance test
    a = macrophage[macrophage.obs.status == 'PDAC'].X[:, macrophage.var_names.get_loc(marker)].toarray().flatten()  # Convert sparse matrix to dense array
    b = macrophage[macrophage.obs.status == 'Adjunct'].X[:, macrophage.var_names.get_loc(marker)].toarray().flatten()

    # Perform Mann-Whitney U test
    mannwhitney_result = stats.mannwhitneyu(a, b)

    # Add p-value and statistic to the plot
    p_value = mannwhitney_result.pvalue
    statistic = mannwhitney_result.statistic
    axes[i].text(0.5, 0.95, f'Mann-Whitney U test\nStatistic: {statistic:.2f}\nP-value: {p_value:.2e}',
                  ha='center', va='top', fontsize=12, transform=axes[i].transAxes)

# Adjust layout
plt.tight_layout()
plt.show()

# Perform differential gene expression analysis between Adjunct and PDAC b cells
sc.tl.rank_genes_groups(macrophage, groupby='status', reference='Adjunct', method='wilcoxon')

# View the top differentially expressed genes
sc.pl.rank_genes_groups(macrophage, n_genes=50, sharey=False)

# Flatten the structured array to get a 1D list of gene names
top_genes = [gene[0] for gene in macrophage.uns['rank_genes_groups']['names'][:50]]
print(top_genes)  # Check the extracted top genes

# Generate dot plot for the top 50 differentially expressed genes
sc.pl.dotplot(macrophage, var_names=top_genes, groupby='status', dendrogram=False)
plt.show()

# Optionally, create a heatmap of the top differentially expressed genes
sc.pl.heatmap(macrophage, var_names=top_genes, groupby='status', cmap='viridis', standard_scale='var')
plt.show()
'''




'''
# Subset NK cells from the entire dataset based on the cell type annotation
nk = adata[adata.obs['cell_type_broad'] == 'NK cells'].copy()

# Ensure that both Adjunct and PDAC statuses are present
print(nk.obs['status'].value_counts())

# Subset given cells into Adjunct and PDAC samples
nk_adjunct = nk[nk.obs['status'] == 'Adjunct'].copy()
nk_pdac = nk[nk.obs['status'] == 'PDAC'].copy()

# Compare expression of key markers between Adjunct and PDAC
markers = ['TRDC', 'NKG7', 'KLRF1']

# Create subplots for each marker
n = len(markers)
fig, axes = plt.subplots(n, 1, figsize=(8, 4 * n), sharex=True)

# Iterate through markers to plot and calculate p-values
for i, marker in enumerate(markers):
    # Plotting the violin plot
    sc.pl.violin(nk, marker, groupby='status', ax=axes[i], show=False)

    # Finding corresponding p-value for significance test
    a = nk[nk.obs.status == 'PDAC'].X[:, nk.var_names.get_loc(marker)].toarray().flatten()  # Convert sparse matrix to dense array
    b = nk[nk.obs.status == 'Adjunct'].X[:, nk.var_names.get_loc(marker)].toarray().flatten()

    # Perform Mann-Whitney U test
    mannwhitney_result = stats.mannwhitneyu(a, b)

    # Add p-value and statistic to the plot
    p_value = mannwhitney_result.pvalue
    statistic = mannwhitney_result.statistic
    axes[i].text(0.5, 0.95, f'Mann-Whitney U test\nStatistic: {statistic:.2f}\nP-value: {p_value:.2e}',
                  ha='center', va='top', fontsize=12, transform=axes[i].transAxes)

# Adjust layout
plt.tight_layout()
plt.show()

# Perform differential gene expression analysis between Adjunct and PDAC b cells
sc.tl.rank_genes_groups(nk, groupby='status', reference='Adjunct', method='wilcoxon')

# View the top differentially expressed genes
sc.pl.rank_genes_groups(nk, n_genes=50, sharey=False)

# Flatten the structured array to get a 1D list of gene names
top_genes = [gene[0] for gene in nk.uns['rank_genes_groups']['names'][:50]]
print(top_genes)  # Check the extracted top genes

# Generate dot plot for the top 50 differentially expressed genes
sc.pl.dotplot(nk, var_names=top_genes, groupby='status', dendrogram=False)
plt.show()

# Optionally, create a heatmap of the top differentially expressed genes
sc.pl.heatmap(nk, var_names=top_genes, groupby='status', cmap='viridis', standard_scale='var')
plt.show()
'''



'''
# Subset given cells from the entire dataset based on the cell type annotation
neutrophil = adata[adata.obs['cell_type_broad'] == 'Neutrophils'].copy()

# Ensure that both Adjunct and PDAC statuses are present
print(neutrophil.obs['status'].value_counts())

# Subset given cells into Adjunct and PDAC samples
neutrophil_adjunct = neutrophil[neutrophil.obs['status'] == 'Adjunct'].copy()
neutrophil_pdac = neutrophil[neutrophil.obs['status'] == 'PDAC'].copy()

# Compare expression of key markers between Adjunct and PDAC
markers = ['CSF3R', 'S100A8', 'CTSG']

# Create subplots for each marker
n = len(markers)
fig, axes = plt.subplots(n, 1, figsize=(8, 4 * n), sharex=True)

# Iterate through markers to plot and calculate p-values
for i, marker in enumerate(markers):
    # Plotting the violin plot
    sc.pl.violin(neutrophil, marker, groupby='status', ax=axes[i], show=False)

    # Finding corresponding p-value for significance test
    a = neutrophil[neutrophil.obs.status == 'PDAC'].X[:, neutrophil.var_names.get_loc(marker)].toarray().flatten()  # Convert sparse matrix to dense array
    b = neutrophil[neutrophil.obs.status == 'Adjunct'].X[:, neutrophil.var_names.get_loc(marker)].toarray().flatten()

    # Perform Mann-Whitney U test
    mannwhitney_result = stats.mannwhitneyu(a, b)

    # Add p-value and statistic to the plot
    p_value = mannwhitney_result.pvalue
    statistic = mannwhitney_result.statistic
    axes[i].text(0.5, 0.95, f'Mann-Whitney U test\nStatistic: {statistic:.2f}\nP-value: {p_value:.2e}',
                  ha='center', va='top', fontsize=12, transform=axes[i].transAxes)

# Adjust layout
plt.tight_layout()
plt.show()

# Perform differential gene expression analysis between Adjunct and PDAC b cells
sc.tl.rank_genes_groups(neutrophil, groupby='status', reference='Adjunct', method='wilcoxon')

# View the top differentially expressed genes
sc.pl.rank_genes_groups(neutrophil, n_genes=50, sharey=False)

# Flatten the structured array to get a 1D list of gene names
top_genes = [gene[0] for gene in neutrophil.uns['rank_genes_groups']['names'][:50]]
print(top_genes)  # Check the extracted top genes

# Generate dot plot for the top 50 differentially expressed genes
sc.pl.dotplot(neutrophil, var_names=top_genes, groupby='status', dendrogram=False)
plt.show()

# Optionally, create a heatmap of the top differentially expressed genes
sc.pl.heatmap(neutrophil, var_names=top_genes, groupby='status', cmap='viridis', standard_scale='var')
plt.show()
'''



'''
# Subset given cells from the entire dataset based on the cell type annotation
mast = adata[adata.obs['cell_type_broad'] == 'Mast cells'].copy()

# Ensure that both Adjunct and PDAC statuses are present
print(mast.obs['status'].value_counts())

# Subset given cells into Adjunct and PDAC samples
mast_adjunct = mast[mast.obs['status'] == 'Adjunct'].copy()
mast_pdac = mast[mast.obs['status'] == 'PDAC'].copy()

# Compare expression of key markers between Adjunct and PDAC
markers = ['KIT', 'TPSAB1', 'HDC', 'IL1RL1']

# Create subplots for each marker
n = len(markers)
fig, axes = plt.subplots(n, 1, figsize=(8, 4 * n), sharex=True)

# Iterate through markers to plot and calculate p-values
for i, marker in enumerate(markers):
    # Plotting the violin plot
    sc.pl.violin(mast, marker, groupby='status', ax=axes[i], show=False)

    # Finding corresponding p-value for significance test
    a = mast[mast.obs.status == 'PDAC'].X[:, mast.var_names.get_loc(marker)].toarray().flatten()  # Convert sparse matrix to dense array
    b = mast[mast.obs.status == 'Adjunct'].X[:, mast.var_names.get_loc(marker)].toarray().flatten()

    # Perform Mann-Whitney U test
    mannwhitney_result = stats.mannwhitneyu(a, b)

    # Add p-value and statistic to the plot
    p_value = mannwhitney_result.pvalue
    statistic = mannwhitney_result.statistic
    axes[i].text(0.5, 0.95, f'Mann-Whitney U test\nStatistic: {statistic:.2f}\nP-value: {p_value:.2e}',
                  ha='center', va='top', fontsize=12, transform=axes[i].transAxes)

# Adjust layout
plt.tight_layout()
plt.show()

# Perform differential gene expression analysis between Adjunct and PDAC b cells
sc.tl.rank_genes_groups(mast, groupby='status', reference='Adjunct', method='wilcoxon')

# View the top differentially expressed genes
sc.pl.rank_genes_groups(mast, n_genes=50, sharey=False)

# Flatten the structured array to get a 1D list of gene names
top_genes = [gene[0] for gene in mast.uns['rank_genes_groups']['names'][:50]]
print(top_genes)  # Check the extracted top genes

# Generate dot plot for the top 50 differentially expressed genes
sc.pl.dotplot(mast, var_names=top_genes, groupby='status', dendrogram=False)
plt.show()

# Optionally, create a heatmap of the top differentially expressed genes
sc.pl.heatmap(mast, var_names=top_genes, groupby='status', cmap='viridis', standard_scale='var')
plt.show()
'''




# +++++++++++++++++++++++++++++++++++++++++++++++++ Immune checkpoint analysis +++++++++++++++++++++++++++++++++++++++#


# ICPs in t cells (significant) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'''
# Define the ICI markers you want to analyze
checkpoint_markers = ['CTLA4', 'PDCD1', 'CD274', 'LAG3', 'HAVCR2', 'TIGIT', 'SIRPA']

# Subset T cells from the entire dataset based on the cell type annotation
t_cells = adata[adata.obs['cell_type_broad'] == 'T cells'].copy()

# Create a DataFrame to hold mean expression data
mean_expression = pd.DataFrame(index=checkpoint_markers)

# Calculate mean expression for each checkpoint marker based on status
for marker in checkpoint_markers:
    mean_expression.loc[marker, 'Adjunct'] = t_cells[t_cells.obs['status'] == 'Adjunct'].X[:, t_cells.var_names.get_loc(marker)].mean()
    mean_expression.loc[marker, 'PDAC'] = t_cells[t_cells.obs['status'] == 'PDAC'].X[:, t_cells.var_names.get_loc(marker)].mean()

# List to hold p-values
p_values = []

# Perform statistical tests for each checkpoint marker
for marker in checkpoint_markers:
    adj_expr = t_cells[t_cells.obs['status'] == 'Adjunct'].X[:, t_cells.var_names.get_loc(marker)].toarray().flatten()
    pdac_expr = t_cells[t_cells.obs['status'] == 'PDAC'].X[:, t_cells.var_names.get_loc(marker)].toarray().flatten()
    _, p_val = stats.mannwhitneyu(adj_expr, pdac_expr)
    p_values.append(p_val)

# Add p-values to the DataFrame
mean_expression['p-value'] = p_values

# Create bar plot for mean expression levels
ax = mean_expression.drop(columns='p-value').plot(kind='bar', figsize=(12, 6), colormap='viridis')
plt.title('Mean Expression of Immune Checkpoint Markers in T Cells (PDAC vs Adjunct)')
plt.ylabel('Mean Expression Level')
plt.axhline(0, color='gray', linewidth=0.8, linestyle='--')
plt.xticks(rotation=45)

# Annotate with p-values
for index, p_value in enumerate(p_values):
    ax.text(index, ax.patches[index].get_height() + 0.05, f'p={p_value:.2e}', ha='center')

plt.tight_layout()
plt.show()
'''


# immune checkpoint markers along the entire dataset

'''
# Define the ICI markers you want to analyze
checkpoint_markers = ['CTLA4', 'PDCD1', 'CD274', 'LAG3', 'HAVCR2', 'TIGIT', 'SIRPA']

# Create a DataFrame to hold mean expression data
mean_expression = pd.DataFrame(index=checkpoint_markers)

# Calculate mean expression for each checkpoint marker based on status for the entire dataset
for marker in checkpoint_markers:
    mean_expression.loc[marker, 'Adjunct'] = adata[adata.obs['status'] == 'Adjunct'].X[:, adata.var_names.get_loc(marker)].mean()
    mean_expression.loc[marker, 'PDAC'] = adata[adata.obs['status'] == 'PDAC'].X[:, adata.var_names.get_loc(marker)].mean()

# List to hold p-values
p_values = []

# Perform statistical tests for each checkpoint marker across the entire dataset
for marker in checkpoint_markers:
    adj_expr = adata[adata.obs['status'] == 'Adjunct'].X[:, adata.var_names.get_loc(marker)].toarray().flatten()
    pdac_expr = adata[adata.obs['status'] == 'PDAC'].X[:, adata.var_names.get_loc(marker)].toarray().flatten()
    _, p_val = stats.mannwhitneyu(adj_expr, pdac_expr)
    p_values.append(p_val)

# Add p-values to the DataFrame
mean_expression['p-value'] = p_values

# Create a bar plot for mean expression levels
ax = mean_expression.drop(columns='p-value').plot(kind='bar', figsize=(12, 6), colormap='viridis')
plt.title('Mean Expression of Immune Checkpoint Markers (PDAC vs Adjunct)')
plt.ylabel('Mean Expression Level')
plt.axhline(0, color='gray', linewidth=0.8, linestyle='--')
plt.xticks(rotation=45)

# Annotate with p-values
for index, p_value in enumerate(p_values):
    ax.text(index, ax.patches[index].get_height() + 0.05, f'p={p_value:.2e}', ha='center')

plt.tight_layout()
plt.show()
'''



#  immune checkpoint analysis along the cell_types
'''
checkpoint_markers = ['CTLA4', 'PDCD1', 'CD274', 'LAG3', 'HAVCR2', 'TIGIT', 'SIRPA']
dotplot = sc.pl.dotplot(adata, var_names=checkpoint_markers, groupby='status',
                        var_group_rotation= 80,
                        figsize=(14, 10))
'''



'''

# optional violin plotting of icp markers

# Define the ICI markers you want to analyze
checkpoint_markers = ['CTLA4', 'PDCD1', 'CD274', 'LAG3', 'HAVCR2', 'TIGIT']

# Subset T cells from the entire dataset based on the cell type annotation
t_cells = adata[adata.obs['cell_type_broad'] == 'T cells'].copy()

# Check the expression of ICI markers
sc.pl.violin(t_cells, keys=checkpoint_markers, groupby='status', jitter=0.4, multi_panel=True)
plt.title('Expression of Immune Checkpoint Inhibitors in T Cells')
plt.show()

# Perform differential gene expression analysis between Adjunct and PDAC T cells for checkpoint markers
sc.tl.rank_genes_groups(t_cells, groupby='status', reference='Adjunct', method='wilcoxon', key='checkpoint_markers')

# View the results of differential expression analysis
sc.pl.rank_genes_groups(t_cells, n_genes=50, sharey=False, key='checkpoint_markers')
plt.show()

# Extracting the top differentially expressed checkpoint genes
top_checkpoint_genes = [gene[0] for gene in t_cells.uns['checkpoint_markers']['names'][:50]]
print("Top differentially expressed checkpoint genes:", top_checkpoint_genes)
'''



#+++++++++++++++++++++++++++++++++++++++++++++ innate subtyping in PDAC vs ADJ ++++++++++++++++++++++++++++++++++++++++#



''' # done on entire cellsssss
# Define the M1 and M2 markers you want to analyze
m1_markers = ['CD80', 'CD86', 'TNF', 'IL1B', 'CXCL10', 'IL12A']
m2_markers = ['CD163', 'MRC1', 'IL10', 'TGFBI', 'CCL22', 'ARG1', 'VEGFA']

# Create a DataFrame to hold mean expression data for M1 and M2 markers
mean_expression_m1 = pd.DataFrame(index=m1_markers)
mean_expression_m2 = pd.DataFrame(index=m2_markers)

# Calculate mean expression for each M1 marker based on status
for marker in m1_markers:
    mean_expression_m1.loc[marker, 'Adjunct'] = adata[adata.obs['status'] == 'Adjunct'].X[:, adata.var_names.get_loc(marker)].mean()
    mean_expression_m1.loc[marker, 'PDAC'] = adata[adata.obs['status'] == 'PDAC'].X[:, adata.var_names.get_loc(marker)].mean()

# Calculate mean expression for each M2 marker based on status
for marker in m2_markers:
    mean_expression_m2.loc[marker, 'Adjunct'] = adata[adata.obs['status'] == 'Adjunct'].X[:, adata.var_names.get_loc(marker)].mean()
    mean_expression_m2.loc[marker, 'PDAC'] = adata[adata.obs['status'] == 'PDAC'].X[:, adata.var_names.get_loc(marker)].mean()

# List to hold p-values for M1 markers
p_values_m1 = []

# Perform statistical tests for each M1 marker
for marker in m1_markers:
    adj_expr = adata[adata.obs['status'] == 'Adjunct'].X[:, adata.var_names.get_loc(marker)].toarray().flatten()
    pdac_expr = adata[adata.obs['status'] == 'PDAC'].X[:, adata.var_names.get_loc(marker)].toarray().flatten()
    _, p_val = stats.mannwhitneyu(adj_expr, pdac_expr)
    p_values_m1.append(p_val)

# Add p-values to the DataFrame for M1
mean_expression_m1['p-value'] = p_values_m1

# Create bar plot for mean expression levels of M1 markers
ax_m1 = mean_expression_m1.drop(columns='p-value').plot(kind='bar', figsize=(12, 6), colormap='viridis')
plt.title('Mean Expression of M1 Markers (PDAC vs Adjunct)')
plt.ylabel('Mean Expression Level')
plt.axhline(0, color='gray', linewidth=0.8, linestyle='--')
plt.xticks(rotation=45)

# Annotate with p-values for M1
for index, p_value in enumerate(p_values_m1):
    ax_m1.text(index, ax_m1.patches[index].get_height() + 0.05, f'p={p_value:.2e}', ha='center')

plt.tight_layout()
plt.show()

# List to hold p-values for M2 markers
p_values_m2 = []

# Perform statistical tests for each M2 marker
for marker in m2_markers:
    adj_expr = adata[adata.obs['status'] == 'Adjunct'].X[:, adata.var_names.get_loc(marker)].toarray().flatten()
    pdac_expr = adata[adata.obs['status'] == 'PDAC'].X[:, adata.var_names.get_loc(marker)].toarray().flatten()
    _, p_val = stats.mannwhitneyu(adj_expr, pdac_expr)
    p_values_m2.append(p_val)

# Add p-values to the DataFrame for M2
mean_expression_m2['p-value'] = p_values_m2

# Create bar plot for mean expression levels of M2 markers
ax_m2 = mean_expression_m2.drop(columns='p-value').plot(kind='bar', figsize=(12, 6), colormap='viridis')
plt.title('Mean Expression of M2 Markers (PDAC vs Adjunct)')
plt.ylabel('Mean Expression Level')
plt.axhline(0, color='gray', linewidth=0.8, linestyle='--')
plt.xticks(rotation=45)

# Annotate with p-values for M2
for index, p_value in enumerate(p_values_m2):
    ax_m2.text(index, ax_m2.patches[index].get_height() + 0.05, f'p={p_value:.2e}', ha='center')

plt.tight_layout()
plt.show()
'''



# lets focus on the macrophage population only
'''
# Define the M1 and M2 markers you want to analyze
m1_markers = ['CD80', 'CD86', 'TNF', 'IL1B', 'CXCL10', 'IL12A']
m2_markers = ['CD163', 'MRC1', 'IL10', 'TGFBI', 'CCL22', 'ARG1', 'VEGFA']

# Subset macrophages from the entire dataset based on the cell type annotation
macrophages = adata[adata.obs['cell_type_broad'] == 'Macrophages'].copy()

# Create a DataFrame to hold mean expression data for M1 and M2 markers
mean_expression_m1 = pd.DataFrame(index=m1_markers)
mean_expression_m2 = pd.DataFrame(index=m2_markers)

# Calculate mean expression for each M1 marker based on status
for marker in m1_markers:
    mean_expression_m1.loc[marker, 'Adjunct'] = macrophages[macrophages.obs['status'] == 'Adjunct'].X[:, macrophages.var_names.get_loc(marker)].mean()
    mean_expression_m1.loc[marker, 'PDAC'] = macrophages[macrophages.obs['status'] == 'PDAC'].X[:, macrophages.var_names.get_loc(marker)].mean()

# Calculate mean expression for each M2 marker based on status
for marker in m2_markers:
    mean_expression_m2.loc[marker, 'Adjunct'] = macrophages[macrophages.obs['status'] == 'Adjunct'].X[:, macrophages.var_names.get_loc(marker)].mean()
    mean_expression_m2.loc[marker, 'PDAC'] = macrophages[macrophages.obs['status'] == 'PDAC'].X[:, macrophages.var_names.get_loc(marker)].mean()

# List to hold p-values for M1 markers
p_values_m1 = []

# Perform statistical tests for each M1 marker
for marker in m1_markers:
    adj_expr = macrophages[macrophages.obs['status'] == 'Adjunct'].X[:, macrophages.var_names.get_loc(marker)].toarray().flatten()
    pdac_expr = macrophages[macrophages.obs['status'] == 'PDAC'].X[:, macrophages.var_names.get_loc(marker)].toarray().flatten()
    _, p_val = stats.mannwhitneyu(adj_expr, pdac_expr)
    p_values_m1.append(p_val)

# Add p-values to the DataFrame for M1
mean_expression_m1['p-value'] = p_values_m1

# Create bar plot for mean expression levels of M1 markers
ax_m1 = mean_expression_m1.drop(columns='p-value').plot(kind='bar', figsize=(12, 6), colormap='viridis')
plt.title('Mean Expression of M1 Markers in Macrophages (PDAC vs Adjunct)')
plt.ylabel('Mean Expression Level')
plt.axhline(0, color='gray', linewidth=0.8, linestyle='--')
plt.xticks(rotation=45)

# Annotate with p-values for M1
for index, p_value in enumerate(p_values_m1):
    ax_m1.text(index, ax_m1.patches[index].get_height() + 0.05, f'p={p_value:.2e}', ha='center')

plt.tight_layout()
plt.show()

# List to hold p-values for M2 markers
p_values_m2 = []

# Perform statistical tests for each M2 marker
for marker in m2_markers:
    adj_expr = macrophages[macrophages.obs['status'] == 'Adjunct'].X[:, macrophages.var_names.get_loc(marker)].toarray().flatten()
    pdac_expr = macrophages[macrophages.obs['status'] == 'PDAC'].X[:, macrophages.var_names.get_loc(marker)].toarray().flatten()
    _, p_val = stats.mannwhitneyu(adj_expr, pdac_expr)
    p_values_m2.append(p_val)

# Add p-values to the DataFrame for M2
mean_expression_m2['p-value'] = p_values_m2

# Create bar plot for mean expression levels of M2 markers
ax_m2 = mean_expression_m2.drop(columns='p-value').plot(kind='bar', figsize=(12, 6), colormap='viridis')
plt.title('Mean Expression of M2 Markers in Macrophages (PDAC vs Adjunct)')
plt.ylabel('Mean Expression Level')
plt.axhline(0, color='gray', linewidth=0.8, linestyle='--')
plt.xticks(rotation=45)

# Annotate with p-values for M2
for index, p_value in enumerate(p_values_m2):
    ax_m2.text(index, ax_m2.patches[index].get_height() + 0.05, f'p={p_value:.2e}', ha='center')

plt.tight_layout()
plt.show()
'''

'''
m1_markers = ['CD80', 'CD86', 'TNF', 'IL1B', 'CXCL10', 'IL12A']
m2_markers = ['CD163', 'MRC1', 'IL10', 'TGFBI', 'CCL22', 'ARG1', 'VEGFA']
dotplot_1 = sc.pl.dotplot(adata, var_names=m1_markers, groupby='cell_type',
                        var_group_rotation= 80,
                        figsize=(14, 10))
dotplot_2 = sc.pl.dotplot(adata, var_names=m2_markers, groupby='cell_type',
                        var_group_rotation= 80,
                        figsize=(14, 10))
'''

'''
# NK plotting


# Define the inflammatory and cytotoxic markers for NK cells
inflammatory_markers = ['FCGR3A', 'NCAM1', 'KLRC1', 'CCR7', 'CXCR4', 'CXCR3']
cytotoxic_markers = ['FCGR3A', 'PRF1']

# Subset NK cells from the entire dataset based on the cell type annotation
nk_cells = adata[adata.obs['cell_type_broad'] == 'NK cells'].copy()

# Create DataFrames to hold mean expression data for inflammatory and suppressive markers
mean_expression_inflammatory = pd.DataFrame(index=inflammatory_markers)
mean_expression_cytotoxic = pd.DataFrame(index=cytotoxic_markers)

# Calculate mean expression for each inflammatory marker based on status
for marker in inflammatory_markers:
    mean_expression_inflammatory.loc[marker, 'Adjunct'] = nk_cells[nk_cells.obs['status'] == 'Adjunct'].X[:, nk_cells.var_names.get_loc(marker)].mean()
    mean_expression_inflammatory.loc[marker, 'PDAC'] = nk_cells[nk_cells.obs['status'] == 'PDAC'].X[:, nk_cells.var_names.get_loc(marker)].mean()

# Calculate mean expression for each suppressive marker based on status
for marker in cytotoxic_markers:
    mean_expression_cytotoxic.loc[marker, 'Adjunct'] = nk_cells[nk_cells.obs['status'] == 'Adjunct'].X[:, nk_cells.var_names.get_loc(marker)].mean()
    mean_expression_cytotoxic.loc[marker, 'PDAC'] = nk_cells[nk_cells.obs['status'] == 'PDAC'].X[:, nk_cells.var_names.get_loc(marker)].mean()

# List to hold p-values for inflammatory markers
p_values_inflammatory = []

# Perform statistical tests for each inflammatory marker
for marker in inflammatory_markers:
    adj_expr = nk_cells[nk_cells.obs['status'] == 'Adjunct'].X[:, nk_cells.var_names.get_loc(marker)].toarray().flatten()
    pdac_expr = nk_cells[nk_cells.obs['status'] == 'PDAC'].X[:, nk_cells.var_names.get_loc(marker)].toarray().flatten()
    _, p_val = stats.mannwhitneyu(adj_expr, pdac_expr)
    p_values_inflammatory.append(p_val)

# Add p-values to the DataFrame for inflammatory markers
mean_expression_inflammatory['p-value'] = p_values_inflammatory

# Create bar plot for mean expression levels of inflammatory markers
ax_inflammatory = mean_expression_inflammatory.drop(columns='p-value').plot(kind='bar', figsize=(12, 6), colormap='viridis')
plt.title('Mean Expression of Inflammatory Markers in NK cells (PDAC vs Adjunct)')
plt.ylabel('Mean Expression Level')
plt.axhline(0, color='gray', linewidth=0.8, linestyle='--')
plt.xticks(rotation=45)

# Annotate with p-values for inflammatory markers
for index, p_value in enumerate(p_values_inflammatory):
    ax_inflammatory.text(index, ax_inflammatory.patches[index].get_height() + 0.05, f'p={p_value:.2e}', ha='center')

plt.tight_layout()
plt.show()

# List to hold p-values for suppressive markers
p_values_cytotoxic = []

# Perform statistical tests for each suppressive marker
for marker in cytotoxic_markers:
    adj_expr = nk_cells[nk_cells.obs['status'] == 'Adjunct'].X[:, nk_cells.var_names.get_loc(marker)].toarray().flatten()
    pdac_expr = nk_cells[nk_cells.obs['status'] == 'PDAC'].X[:, nk_cells.var_names.get_loc(marker)].toarray().flatten()
    _, p_val = stats.mannwhitneyu(adj_expr, pdac_expr)
    p_values_cytotoxic.append(p_val)

# Add p-values to the DataFrame for suppressive markers
mean_expression_cytotoxic['p-value'] = p_values_cytotoxic

# Create bar plot for mean expression levels of suppressive markers
ax_cytotoxic = mean_expression_cytotoxic.drop(columns='p-value').plot(kind='bar', figsize=(12, 6), colormap='viridis')
plt.title('Mean Expression of Suppressive Markers in NK cells (PDAC vs Adjunct)')
plt.ylabel('Mean Expression Level')
plt.axhline(0, color='gray', linewidth=0.8, linestyle='--')
plt.xticks(rotation=45)

# Annotate with p-values for suppressive markers
for index, p_value in enumerate(p_values_cytotoxic):
    ax_cytotoxic.text(index, ax_cytotoxic.patches[index].get_height() + 0.05, f'p={p_value:.2e}', ha='center')

plt.tight_layout()
plt.show()
'''


'''
# Define N1 and N2 markers for neutrophils (with official gene symbols)
n1_markers = ['TNF', 'ICAM1', 'FAS', 'SOD1']
n2_markers = ['ARG1', 'CCL2', 'CCL5']

# Subset neutrophils from the dataset based on cell type annotation
neutrophils = adata[adata.obs['cell_type_broad'] == 'Neutrophils'].copy()

# Create a DataFrame to hold mean expression data for N1 and N2 markers
mean_expression_n1 = pd.DataFrame(index=n1_markers)
mean_expression_n2 = pd.DataFrame(index=n2_markers)

# Calculate mean expression for each N1 marker based on status
for marker in n1_markers:
    mean_expression_n1.loc[marker, 'Adjunct'] = neutrophils[neutrophils.obs['status'] == 'Adjunct'].X[:, neutrophils.var_names.get_loc(marker)].mean()
    mean_expression_n1.loc[marker, 'PDAC'] = neutrophils[neutrophils.obs['status'] == 'PDAC'].X[:, neutrophils.var_names.get_loc(marker)].mean()

# Calculate mean expression for each N2 marker based on status
for marker in n2_markers:
    mean_expression_n2.loc[marker, 'Adjunct'] = neutrophils[neutrophils.obs['status'] == 'Adjunct'].X[:, neutrophils.var_names.get_loc(marker)].mean()
    mean_expression_n2.loc[marker, 'PDAC'] = neutrophils[neutrophils.obs['status'] == 'PDAC'].X[:, neutrophils.var_names.get_loc(marker)].mean()

# List to hold p-values for N1 markers
p_values_n1 = []

# Perform statistical tests for each N1 marker
for marker in n1_markers:
    adj_expr = neutrophils[neutrophils.obs['status'] == 'Adjunct'].X[:, neutrophils.var_names.get_loc(marker)].toarray().flatten()
    pdac_expr = neutrophils[neutrophils.obs['status'] == 'PDAC'].X[:, neutrophils.var_names.get_loc(marker)].toarray().flatten()
    _, p_val = stats.mannwhitneyu(adj_expr, pdac_expr)
    p_values_n1.append(p_val)

# Add p-values to the DataFrame for N1
mean_expression_n1['p-value'] = p_values_n1

# Create bar plot for mean expression levels of N1 markers
ax_n1 = mean_expression_n1.drop(columns='p-value').plot(kind='bar', figsize=(12, 6), colormap='viridis')
plt.title('Mean Expression of N1 Markers in Neutrophils (PDAC vs Adjunct)')
plt.ylabel('Mean Expression Level')
plt.axhline(0, color='gray', linewidth=0.8, linestyle='--')
plt.xticks(rotation=45)

# Annotate with p-values for N1
for index, p_value in enumerate(p_values_n1):
    ax_n1.text(index, ax_n1.patches[index].get_height() + 0.05, f'p={p_value:.2e}', ha='center')

plt.tight_layout()
plt.show()

# List to hold p-values for N2 markers
p_values_n2 = []

# Perform statistical tests for each N2 marker
for marker in n2_markers:
    adj_expr = neutrophils[neutrophils.obs['status'] == 'Adjunct'].X[:, neutrophils.var_names.get_loc(marker)].toarray().flatten()
    pdac_expr = neutrophils[neutrophils.obs['status'] == 'PDAC'].X[:, neutrophils.var_names.get_loc(marker)].toarray().flatten()
    _, p_val = stats.mannwhitneyu(adj_expr, pdac_expr)
    p_values_n2.append(p_val)

# Add p-values to the DataFrame for N2
mean_expression_n2['p-value'] = p_values_n2

# Create bar plot for mean expression levels of N2 markers
ax_n2 = mean_expression_n2.drop(columns='p-value').plot(kind='bar', figsize=(12, 6), colormap='viridis')
plt.title('Mean Expression of N2 Markers in Neutrophils (PDAC vs Adjunct)')
plt.ylabel('Mean Expression Level')
plt.axhline(0, color='gray', linewidth=0.8, linestyle='--')
plt.xticks(rotation=45)

# Annotate with p-values for N2
for index, p_value in enumerate(p_values_n2):
    ax_n2.text(index, ax_n2.patches[index].get_height() + 0.05, f'p={p_value:.2e}', ha='center')

plt.tight_layout()
plt.show()
'''

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#



'''
# manual enrichment of the top 50 DEGs x 5 immune cell types (in PDAC vs ADjunct)

# Define immune cell types
immune_cell_types = ['T cells', 'B cells', 'NK cells', 'Macrophages', 'Neutrophils', 'Mast cells']  # Adjust based on your annotations

# Subset the adata object to include only immune cells using the broad cell type
adata_immune = adata[adata.obs['cell_type_broad'].isin(immune_cell_types)]

de = ['HSP90AA1', 'HSP90AB1', 'CREM', 'DNAJA1', 'FTH1', 'RGCC', 'H3F3B', 'TUBB4B', 'HSPA8', 'JUND', 'EZR', 'CREM', 'CD55', 'JUND', 'CALM1', 'TUBB4B', 'PTMA', 'SNHG5', 'TMSB10', 'CRIP1', 'NUPR1', 'TMEM176B', 'FOSB', 'PMP22', 'TMEM176A', 'LGALS3BP', 'FCGR3A', 'FN1', 'SLC16A10', 'FTL', 'CXCR4', 'CREM', 'HSP90AA1', 'HSP90AB1', 'FTH1', 'RPL13A', 'JUND', 'PTMA', 'TNFAIP3', 'LDHA', 'CXCL8', 'FTH1', 'PLAUR', 'IER3', 'TPT1', 'G0S2', 'NFKBIA', 'PPIF', 'C15orf48', 'TNFAIP3']

enr_de = gp.enrichr(gene_list=de,
                      gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2023', 'BioCarta_2013',
                                 'Cancer_Cell_Line_Encyclopedia', 'DisGeNET'],
                      organism='human',
                      outdir=None,
                      background=adata_immune.var_names.tolist())

enr_de.results.to_csv(extra_pdac + '/de_immune_pdac_enrichment.csv', index=False)


de_enrichment = pd.read_csv(extra_pdac + '/de_immune_pdac_enrichment.csv')

# Define a function to plot enrichment results for a given gene set
def plot_enrichment_results(enrichment_df, title):
    # Set the plot style
    sns.set(style='whitegrid')

    # Filter for significant pathways (optional)
    significant_terms = enrichment_df[enrichment_df['Adjusted P-value'] < 0.02]

    # Plotting each pathway in separate figures
    for gene_set in significant_terms['Gene_set'].unique():
        pathway_data = significant_terms[significant_terms['Gene_set'] == gene_set].sort_values(by='Adjusted P-value')

        plt.figure(figsize=(10, 6))
        sns.barplot(x='Adjusted P-value', y='Term', data=pathway_data, color='red')
        plt.title(f'{title} - {gene_set}')
        plt.xlabel('Adjusted P-value')
        plt.ylabel('Enriched Term')
        plt.tight_layout()
        plt.show()

# Plot for upregulated genes
plot_enrichment_results(de_enrichment, 'Enrichment Results for Upregulated Genes')



de = ['HSP90AA1', 'HSP90AB1', 'CREM', 'DNAJA1', 'FTH1', 'RGCC', 'H3F3B', 'TUBB4B', 'HSPA8', 'JUND', 'EZR', 'CREM', 'CD55', 'JUND', 'CALM1', 'TUBB4B', 'PTMA', 'SNHG5', 'TMSB10', 'CRIP1', 'NUPR1', 'TMEM176B', 'FOSB', 'PMP22', 'TMEM176A', 'LGALS3BP', 'FCGR3A', 'FN1', 'SLC16A10', 'FTL', 'CXCR4', 'CREM', 'HSP90AA1', 'HSP90AB1', 'FTH1', 'RPL13A', 'JUND', 'PTMA', 'TNFAIP3', 'LDHA', 'CXCL8', 'FTH1', 'PLAUR', 'IER3', 'TPT1', 'G0S2', 'NFKBIA', 'PPIF', 'C15orf48', 'TNFAIP3']
# Creating dotplot for anna
dotplot = sc.pl.dotplot(adata, var_names=de, groupby='status',
                        var_group_rotation= 80,  # Rotate the variable (gene) names by 45 degrees
                        figsize=(14, 10))
'''

'''
# Define immune cell types based on your annotations (adjust as necessary)
immune_cell_types = ['T cells', 'B cells', 'NK cells', 'Macrophages', 'Neutrophils', 'Mast cells']

# Subset the adata object to include only immune cells
adata_immune = adata[adata.obs['cell_type_broad'].isin(immune_cell_types)].copy()

# Ensure the subsetted data is prepared for scVI (if not already done)
scvi.model.SCVI.setup_anndata(adata_immune, layer='counts', categorical_covariate_keys=['status'])

# Train scVI model for the immune cell subset
scvi_immune_model = scvi.model
scvi_immune_model.train()

# Perform differential expression between immune cells in PDAC vs. Adjunct
de_results_immune = scvi_immune_model.differential_expression(
    groupby="status",
    group1="PDAC",
    group2="Adjunct",
    n_top_genes=50
)

# View the top 50 differentially expressed genes for immune cells
print(de_results_immune.head(50))

# Save the results to a CSV file
de_results_immune.to_csv(extra_pdac + '/immune_de_results_scvi.csv', index=False)

# Optional: Plot top differentially expressed genes
top_genes_immune = de_results_immune['gene'][:50].tolist()

# Dot plot for top 50 genes
sc.pl.dotplot(adata_immune, var_names=top_genes_immune, groupby='status', dendrogram=False)

# Heatmap for top 50 genes
sc.pl.heatmap(adata_immune, var_names=top_genes_immune, groupby='status', cmap='viridis', standard_scale='var')
'''