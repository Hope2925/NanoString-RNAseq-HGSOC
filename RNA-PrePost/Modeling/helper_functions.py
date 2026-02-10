import numpy as np
import pandas as pd
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, roc_auc_score, confusion_matrix
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import LeaveOneOut
from statsmodels.stats.multitest import multipletests
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import RobustScaler, StandardScaler

from collections import defaultdict, Counter
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

color_dict = {
        "Nano_M": "#4876FF", 
        "Nano_J": "#27408B", 
        "Nano_B": "#009ACD",
        "RNAseq_A": "#C03830", 
        "RNAseq_J": "#860E07", 
        "RNAseq_EGA": "#F85D48", 
        "Nano": "#0000EE", 
        "RNAseq": "#EE4000", 
        "All": "#7D26CD", 
				"Tr/Val": "#7D26CD", 
        "All_Assay": "#9966ff",
        "Hold_Out": "#CDCD00",
        
		# NanoString and RNAseq paired data
        "Manso_Adzib": "#000066", # blue
        "Manso_Jav": "#0000ff", 
        "Manso_EGA": "#8080ff",
        
        "James_Adzib": "#800080",  # pink
        "James_Jav": "#cc00cc", 
        "James_EGA": "#ff66ff",

        "Bitler_Adzib": "#800000", # red
        "Bitler_Jav": "#ff3333", 
        "Bitler_EGA": "#ff9999",
				
		"PFS>12": "#008B45", 
		"PFS<=12": "grey"
        
    }

color_dict = {
        "Nano_M": "#FDAE6B", 
        "Nano_J": "#D94801", 
        "Nano_B": "#7F2704",
        "RNAseq_A": "#3f007d", 
        "RNAseq_J": "#bcbddc", 
        "RNAseq_EGA": "#807dba", 
        "Nano": "#FD8D3C", 
        "RNAseq": "#6A51A3", 
        "All": "#ae017e", 
				"Tr/Val": "#ae017e", 
        "All_Assay": "#ae017e",
        "Hold_Out": "#008B45",
        
		# NanoString and RNAseq paired data
        "Manso_Adzib": "#000066", # blue
        "Manso_Jav": "#0000ff", 
        "Manso_EGA": "#8080ff",
        
        "James_Adzib": "#800080",  # pink
        "James_Jav": "#cc00cc", 
        "James_EGA": "#ff66ff",

        "Bitler_Adzib": "#000000", # red
        "Bitler_Jav": "#525252", 
        "Bitler_EGA": "#969696",
				
		"PFS>12": "#008B45", 
		"PFS<=12": "grey"
        
    }

# EGA Patients split by SC, Bulk, or Both
sc_only = ["EOC1005", "EOC153", "EOC227", "EOC349", "EOC540", "EOC733"] # 6
bulk_only = ["EOC1129", "EOC183", "EOC218", "EOC26", "EOC376", "EOC423", "EOC568", "EOC587", "EOC649", "EOC677", "EOC883", "EOC933", "EOC868", "EOC1133", "EOC167", "EOC740", "EOC551", "EOC60"]
sc_bulk = ["EOC136", "EOC3", "EOC372", "EOC443", "EOC87"] # all in training dataset

### TEST DATA
# third lowest and third highest AND middle
NanovRNAseq_pat_ho_dict = {"Nano_M":["P4M", "P7M", "P6M"], "Nano_J":["P22J", "P7J", "P21J"], "Nano_B":["P8B", "P17B", "P18B"], 
                           "RNAseq_A":["USF3", "USF16", "WSU_F"], "RNAseq_J":["PJ_25", "PJ_17", "PJ_21"], "RNAseq_EGA":["EOC26", "EOC1129", "EOC1133"]}
# third lowest and third highest
PFS_pat_ho_dict = {"Nano_M":["P4M", "P6M"], "Nano_J":["P22J", "P21J"], "Nano_B":["P8B", "P18B"], 
					"RNAseq_A":["USF3", "WSU_F"], "RNAseq_J":["PJ_25", "PJ_21"], "RNAseq_EGA":["EOC26", "EOC1133"]}
# third from lowest/highest AND middle
PFS_pat_ho_dict_large = {"Nano_M":["P4M", "P6M"], "Nano_J":["P22J", "P21J", "P1J", "P15J"], "Nano_B":["P8B", "P18B", "P30B", "P12B"], 
					"RNAseq_A":["USF3", "WSU_F", "USF14", "USF34"], "RNAseq_J":["PJ_25", "PJ_21", "PJ_5", "PJ_29"], "RNAseq_EGA":["EOC26", "EOC1133", "EOC740", "EOC568"]}

#####################################
##### DISPLAY #######################
def display_summary_ns(summary_df):
	# Flatten all values to compute global min/max
	all_vals = summary_df.values.flatten()
	global_min, global_max = all_vals.min(), all_vals.max()

	# Normalize values across the whole table
	cell_colors = np.zeros((summary_df.shape[0], summary_df.shape[1], 4))  # RGBA
	cmap = plt.cm.Blues

	for i in range(summary_df.shape[0]):
		for j in range(summary_df.shape[1]):
			val = summary_df.iloc[i, j]
			norm = (val - global_min) / (global_max - global_min + 1e-6)  # avoid div0
			cell_colors[i, j, :] = cmap(norm)

	# Create figure
	fig, ax = plt.subplots(figsize=(10, len(summary_df)*0.5 + 1))
	ax.axis('off')

	# Build table
	the_table = ax.table(
		cellText=summary_df.values,
		rowLabels=summary_df.index,
		colLabels=summary_df.columns,
		cellColours=cell_colors,
		cellLoc='center',
		rowLoc='center',
		colLoc='center',
		loc='center'
	)

	# Adjust text color based on brightness
	for (row_idx, col_idx), cell in the_table.get_celld().items():
		if row_idx == 0:  # header row
			cell.get_text().set_color("black")
			cell.set_fontsize(12)
		else:
			r, g, b, _ = cell.get_facecolor()
			brightness = 0.299*r + 0.587*g + 0.114*b
			cell.get_text().set_color("white" if brightness < 0.5 else "black")
			cell.set_fontsize(12)

	the_table.auto_set_font_size(False)
	the_table.scale(1, 1.5)
	plt.title(f"Dataset Sample Sizes and Class Counts", fontsize=14, pad=20)
	plt.show()


#####################################################
#### DATA PREP ###########
#####################################################
# Get the pairwise datasets with appropriate labels
def get_final_nr_dataset(nr_data_set_dict, NanovRNAseq_pat_ho_dict, nano_genes):
	"""
	Takes in a dictionary of just dfs and validation patients to get the 
	dictionary  of datasets for NanoString vs RNAseq comparison

	Inputs:
	* nr_data_set_dict: dictionary with Nano_M, Nano_J, Nano_B, RNAseq_A, RNAseq_J, RNAseq_EGA keys
    and values are the full dataframes (includes columns Patient and genes)
	* NanovRNAseq_pat_ho_dict: dictionary with the same keys as nr_data_set_dict that determines which 
    of the patients are reserved for the hold out dataset
	* nano_genes: list of genes to be kept
	
	Returns:
	* The new dictionary ready for training/analysis (nr_new_data_set_dict)
	* the df with all the values for possible future usage (nr_all_df) (no hold out)
	"""
	## STEP 1: Get the datasets with the hold out removed
	X_final_val_nr = []
	for dataset, df in nr_data_set_dict.items():
		# split up according to patients
		X_final_val_nr.append(df[df['Patient'].isin(NanovRNAseq_pat_ho_dict[dataset])])
		nr_data_set_dict[dataset] = df[~df['Patient'].isin(NanovRNAseq_pat_ho_dict[dataset])]
	
	## STEP 2: GET THE PAIRWISE datasets
	man_adzib_df = pd.concat([nr_data_set_dict["Nano_M"],nr_data_set_dict["RNAseq_A"]], axis=0, ignore_index=True)
	man_jav_df = pd.concat([nr_data_set_dict["Nano_M"],nr_data_set_dict["RNAseq_J"]], axis=0, ignore_index=True)
	man_ega_df = pd.concat([nr_data_set_dict["Nano_M"],nr_data_set_dict["RNAseq_EGA"]], axis=0, ignore_index=True)

	bitler_adzib_df = pd.concat([nr_data_set_dict["Nano_B"],nr_data_set_dict["RNAseq_A"]], axis=0, ignore_index=True)
	bitler_jav_df = pd.concat([nr_data_set_dict["Nano_B"],nr_data_set_dict["RNAseq_J"]], axis=0, ignore_index=True)
	bitler_ega_df = pd.concat([nr_data_set_dict["Nano_B"],nr_data_set_dict["RNAseq_EGA"]], axis=0, ignore_index=True)

	james_adzib_df = pd.concat([nr_data_set_dict["Nano_J"],nr_data_set_dict["RNAseq_A"]], axis=0, ignore_index=True)
	james_jav_df = pd.concat([nr_data_set_dict["Nano_J"],nr_data_set_dict["RNAseq_J"]], axis=0, ignore_index=True)
	james_ega_df = pd.concat([nr_data_set_dict["Nano_J"],nr_data_set_dict["RNAseq_EGA"]], axis=0, ignore_index=True)

	# make all df
	nr_all_df = pd.concat([nr_data_set_dict["Nano_M"], nr_data_set_dict["Nano_B"], nr_data_set_dict["Nano_J"], 
						nr_data_set_dict["RNAseq_A"], nr_data_set_dict["RNAseq_J"], nr_data_set_dict["RNAseq_EGA"]], axis=0, ignore_index=True)

	# save data in dict for easy access
	nr_new_data_set_dict = {"Manso_Adzib":man_adzib_df, "Manso_Jav":man_jav_df,  "Manso_EGA":man_ega_df, 
					"Bitler_Adzib":bitler_adzib_df, "Bitler_Jav":bitler_jav_df, "Bitler_EGA":bitler_ega_df,
					"James_Adzib":james_adzib_df, "James_Jav":james_jav_df, "James_EGA":james_ega_df, 
					"Tr/Val":nr_all_df, "Hold_Out":pd.concat(X_final_val_nr, axis=0)}
      
	## STEP 3: GET THE y (RNASEQ VS NANOSTRING)
	for dataset, df in nr_new_data_set_dict.items():
		# Also get the info from the original FULL data for comparison
		# Numeric binary labels
		y = (df["Experiment"].isin(['Nano Manso', 'Nano Bitler', 'Nano James'])).astype(int)
		# Labeled version
		y_clean = y.map({1: 'Nanostring', 0: 'RNAseq'})
		X =  df.drop(columns=["Patient", "Experiment", "PFS_mths", "PFS_bool"]) 
		# only include genes considered in Nanostring
		cols_to_keep_existing = [c for c in nano_genes if c in X.columns]
		X = X[cols_to_keep_existing]
		# Remove columns containing any NaN values
		X = X.dropna(axis=1)
		# get the metadata too
		metadata_keep = df.loc[:, ~df.columns.isin(cols_to_keep_existing)]
		# remove any genes that are constant for all patients (otherwise will get warning later0
		X = X.loc[:, X.nunique() > 1]
		print(f"Number of genes left for {dataset}:", X.shape[1])
		nr_new_data_set_dict[dataset] = [X, y, y_clean, metadata_keep]
	return(nr_new_data_set_dict, nr_all_df)

def get_final_pfs_dataset(data_set_dict, PFS_pat_ho_dict, nano_genes, prediction_mth=12, RNAseqonly=False):
	"""
	Takes in a dictionary of dfs and validation patients to get the 
	dictionary  of datasets for PFS comparison

	Inputs:
	* data_set_dict: dictionary with Nano_M, Nano_J, Nano_B, RNAseq_A, RNAseq_J, RNAseq_EGA keys
    and values are the full dataframes (includes columns Patient and genes)
	* PFS_pat_ho_dict: dictionary with the same keys as data_set_dict that determines which 
    of the patients are reserved for the hold out dataset
	* nano_genes: list of genes to be kept
	
	Returns:
	* The new dictionary ready for training/analysis (new_data_set_dict)
	* the df with all the values for possible future usage (all_df) (no hold out)
	"""
	
	new_data_set_dict = data_set_dict.copy()
	# if the PFS_pat_ho_dict is "Adzib", that means just use Adzib as the hold out group
	if PFS_pat_ho_dict == "Adzib":
		new_data_set_dict["Hold_Out"] = new_data_set_dict["RNAseq_A"]
		new_data_set_dict["Tr/Val"] = pd.concat([new_data_set_dict["RNAseq_J"], new_data_set_dict["RNAseq_EGA"]], axis=0)
		del new_data_set_dict['RNAseq_A']
	else:
		## STEP 1: Get the datasets with the hold out removed
		X_final_val = []
		for dataset, df in data_set_dict.items():
			# split up according to patients
			X_final_val.append(df[df['Patient'].isin(PFS_pat_ho_dict[dataset])])
			new_data_set_dict[dataset] = df[~df['Patient'].isin(PFS_pat_ho_dict[dataset])]
		# save the final hold out set
		new_data_set_dict["Hold_Out"] = pd.concat(X_final_val, axis=0)
		if RNAseqonly is False:
			new_data_set_dict["Nano"] = pd.concat([new_data_set_dict["Nano_M"], new_data_set_dict["Nano_J"], new_data_set_dict["Nano_B"]], axis=0)
			new_data_set_dict["RNAseq"] = pd.concat([new_data_set_dict["RNAseq_A"], new_data_set_dict["RNAseq_J"], new_data_set_dict["RNAseq_EGA"]], axis=0)
			new_data_set_dict["Tr/Val"] = pd.concat([new_data_set_dict["Nano"], new_data_set_dict["RNAseq"]], axis=0)
			RNA_copy = new_data_set_dict["RNAseq"].copy()
			Nano_copy = new_data_set_dict["Nano"].copy()
			Nano_copy["Nano"] = 1
			RNA_copy["Nano"] = 0
			#new_data_set_dict["All_Assay"] = pd.concat([Nano_copy, RNA_copy], axis=0)
		else:
			new_data_set_dict["Tr/Val"] = pd.concat([new_data_set_dict["RNAseq_A"], new_data_set_dict["RNAseq_J"], new_data_set_dict["RNAseq_EGA"]], axis=0)

 
	

	## STEP 2: Get the dictionary set up for training/analysis
	for dataset, df in new_data_set_dict.items():
		# Numeric or labeled binary labels
		y = (df["PFS_mths"].astype(float) <= prediction_mth).astype(int)
		y_clean = y.map({1: f'PFS<={prediction_mth}', 0: f'PFS>{prediction_mth}'})
		# only keep gene information (and whether NanoString or not)
		cols_to_keep_existing = [c for c in nano_genes if c in df.columns]
		if dataset == "All_Assay":
			cols_to_keep_existing.append("Nano")
		X = df[cols_to_keep_existing]
		# Remove columns containing any NaN values
		X = X.dropna(axis=1)
		metadata_keep = df.loc[:, ~df.columns.isin(cols_to_keep_existing)]
		# remove any genes that are constant for all patients (otherwise will get warning later0
		X = X.loc[:, (X.nunique() > 1) | (X.columns == "Nano")]
		print(f"Number of genes left for {dataset}:", X.shape[1])
		new_data_set_dict[dataset] = [X, y, y_clean, metadata_keep]
	return(new_data_set_dict)
	


#####################################################
#### FEATURE SELECTION ###########
#####################################################

####### BOOTSTRAPPING APPROACH #################
def graph_bootstrap_samples(freq_df, sample_size, n_bootstraps):
	"""
	Graph the sampling of the boostrapped samples with "expectation" based on completely random values
	"""
	n_samples_total = len(freq_df)
        # total number of unique samples
	expected_per_sample = (n_bootstraps * sample_size) / n_samples_total
	# --- Plot ---
	plt.figure(figsize=(7, 5))
	sns.histplot(freq_df['Count'], bins=20, kde=True, color="skyblue", label='Observed')
	plt.axvline(expected_per_sample, color='red', linestyle='--', lw=2, label='Expected mean frequency')
	# Add ±1 SD region (optional, helps visualize variance)
	std_count = np.std(freq_df['Count'])
	plt.axvspan(expected_per_sample - std_count, expected_per_sample + std_count,
				color='red', alpha=0.15, label='±1 SD region')
	# Style
	plt.title(f"Bootstrap Sampling Frequency (SampleSize={sample_size}, #Btstps={n_bootstraps})", fontsize=14)
	plt.xlabel("Number of times sample was selected", fontsize=12)
	plt.ylabel("Number of samples", fontsize=12)
	plt.legend()
	plt.tight_layout()
	plt.show()

def plot_convergence(top_lists, n_bootstraps_to_plot=None, top_display=20):
	"""
	Plots the convergence towards top genes for bootstrapping
	top_display: Number of genes consider
	"""
	if n_bootstraps_to_plot is None:
		n_bootstraps_to_plot = len(top_lists)
	all_genes = sorted({g for lst in top_lists for g in lst})
	counts = Counter()
	cum_prop = {g: [] for g in all_genes}
	for i, lst in enumerate(top_lists[:n_bootstraps_to_plot], start=1):
		counts.update(lst)
		for g in all_genes:
			cum_prop[g].append(counts[g] / i)
	# pick top by final proportion
	final_prop = {g: cum_prop[g][-1] for g in all_genes}
	top_genes = sorted(final_prop.keys(), key=lambda x: final_prop[x], reverse=True)[:top_display]

	plt.figure(figsize=(8,4))
	num = 0
	for g in top_genes:
		if num < 11:
 			plt.plot(range(1, n_bootstraps_to_plot+1), cum_prop[g], label=g)
		else:
			plt.plot(range(1, n_bootstraps_to_plot+1), cum_prop[g])
		num = num + 1
	plt.xlabel('Number of bootstraps')
	plt.ylabel(f'Cumulative proportion in top-{top_display}')
	plt.legend(bbox_to_anchor=(1.05,1), loc='upper left', fontsize='small')
	plt.title(f"Convergence of top-{top_display} selection frequencies")
	plt.tight_layout()
	plt.show()
	return cum_prop

def btsp_selection_replacement_weighted(X_use, y_use, meta_use,
n_bootstraps=200, top_k=50, random_state=30):
	"""
    Stability selection via repeated sampling WITH replacement.
    - X_use: pandas DataFrame (samples x features).
    - y_use: pandas Series or array-like of binary labels (length = n_samples).
	- meta_use: pandas DataFrame with columns including Patient, and Experiment
    - n_bootstraps: number of subsamples to draw
    - top_k: track top_k features by F-score each bootstrap
    Returns:
      - stability_df: DataFrame with Gene, Frequency, Mean_F, Median_padj, Freq_pct
      - top_lists: list of top_k lists for each bootstrap
      - cum_prop: dict -- convergence tracking output
    """
	# get a random seed fixed for reproducibility
	np.random.seed(random_state)

	# Median imputation in case of NAs
	X_use = X_use.copy()
	for column in X_use.columns:
		if X_use[column].dtype in ['int64', 'float64']:
			X_use.loc[:, column] = X_use[column].fillna(X_use[column].median())

	# Remove any 0 variance cases
	std_devs = X_use.std()
	cols_to_drop = std_devs[std_devs == 0].index.tolist()
	X_use_cleaned = X_use.drop(columns=cols_to_drop)


	# Align indices
	y_use = y_use.reset_index(drop=True)
	X_use_cleaned = X_use_cleaned.reset_index(drop=True)
	meta_use = meta_use.reset_index(drop=True)

	# Get experiment weighting so that each experiment is likely equally represented
	exp_counts = meta_use['Experiment'].value_counts()
	exp_weights = 1 / exp_counts
	exp_weights = exp_weights / exp_weights.sum() # normalize to sum to one
	# Map weight to each sample
	sample_weights = meta_use['Experiment'].map(exp_weights)

	# Get class indices
	idx_0 = y_use[y_use == 0].index
	idx_1 = y_use[y_use == 1].index
	print("Samples per class:", len(idx_0), len(idx_1))

	# Store top genes across bootstraps
	top_lists = []; top_genes_all = []
	# Store the frequency of samples across bootstraps to ensure one isn't called more than the rest
	sample_freq = []
	fscore_dict = defaultdict(list)
	padj_dict = defaultdict(list)

	for i in range(n_bootstraps):
		# Stratified bootstrap sample (with replacement)
		btsp_idx_0 = np.random.choice(idx_0, size=len(idx_0), replace=True,
									p=sample_weights[idx_0] / sample_weights[idx_0].sum())
		btsp_idx_1 = np.random.choice(idx_1, size=len(idx_1), replace=True, 
									p=sample_weights[idx_1] / sample_weights[idx_1].sum())
		btsp_idx = np.concatenate([btsp_idx_0, btsp_idx_1])
		X_btsp = X_use_cleaned.loc[btsp_idx]
		y_btsp = y_use.loc[btsp_idx]
		meta_btsp = meta_use.loc[btsp_idx]

		# Run ANOVA F-test
		f_vals, p_vals = f_classif(X_btsp, y_btsp)
        # Adjust p-values
		_, p_adj, _, _ = multipletests(p_vals, method='fdr_bh')
		# Collect results
		# Store results per gene
		for gene, f_val, padj in zip(X_use_cleaned.columns, f_vals, p_adj):
			fscore_dict[gene].append(f_val)
			padj_dict[gene].append(padj)
						
		# Get top_k genes this round
		topk_genes = pd.Series(X_use_cleaned.columns)[np.argsort(f_vals)[::-1][:top_k]]
		# Record top_k genes for stability tracking
		top_lists.append(topk_genes)
		top_genes_all.extend(topk_genes)
		# Record indices for tracking
		sample_freq.extend(btsp_idx.tolist())

	# Compute stability frequency
	gene_freq = Counter(top_genes_all)
	sample_freq = Counter(sample_freq)
	all_genes = list(X_use_cleaned.columns)
	stability_df = pd.DataFrame({
        "Gene": all_genes,
        "Frequency": [gene_freq.get(g, 0) for g in all_genes],
        "Mean_F": [np.nanmean(fscore_dict[g]) for g in all_genes],
        "Median_padj": [np.nanmedian(padj_dict[g]) for g in all_genes]
    })

	# Get sample frequency in bootstrap with metadata
	freq_df = pd.DataFrame.from_dict(sample_freq, orient='index', columns=['Count'])
	freq_df = pd.concat([freq_df, meta_use.iloc[freq_df.index]], axis=1)
	freq_df = freq_df.sort_values('Count', ascending=False).reset_index().rename(columns={'index': 'Sample'})
	graph_bootstrap_samples(freq_df, len(btsp_idx), n_bootstraps)
	cum_prop = plot_convergence(top_lists, top_display=100)
	

	# Get the Fraction frequency 
	stability_df["Freq_Percent"] = stability_df["Frequency"] / n_bootstraps
	# Sort by both Frequency and Mean F-score
	stability_df = stability_df.sort_values(by=["Frequency", "Mean_F"], ascending=[False, False]).reset_index(drop=True)
	
	return stability_df, freq_df, cum_prop

def get_stable_genes(cum_prop, num_btsp, prop_limit):
	# cum_prop: dictionrary where keys are genes and values are lists of bootstraps
	# num_btsp: Bootstrap number at which to use limit
	# prop_limit: lower limit of genes that will be considered
	# returns a list of genes that hit the limit
	kept_genes = []
	prop_list = []
	all_genes = []
	all_prop_list = []
	for gene, btsps in cum_prop.items():
		all_genes.append(gene)
		all_prop_list.append(btsps[num_btsp])
		if btsps[num_btsp] > prop_limit:
			kept_genes.append(gene)
			prop_list.append(btsps[num_btsp])
	print(len(kept_genes), "kept genes")
	# sort by the proportion
	kept_genes_sorted = [gene for _, gene in sorted(zip(prop_list, kept_genes), reverse=True)]
	# also return a dictionary with the cumulative proportions at this point
	cum_prop_df = pd.DataFrame({"Gene":all_genes, "Cum_Prop":all_prop_list})
	return kept_genes_sorted, cum_prop_df

## Plotting the F score percentages between the different considerations
def plot_Fscore_comps(merged_df, type_):
	freq_cols = [c for c in merged_df.columns if "Freq_Percent" in c]
	for ycol in freq_cols:
		if ycol == 'Freq_Percent_NR_All':
			continue
		# Get the labeling info for genes above 0.2 in both
		label_df = merged_df[(merged_df['Freq_Percent_NR_All'] > 0.1) & (merged_df[ycol] > 0.2)]
		plt.figure(figsize=(6,6))
		sns.scatterplot(data=merged_df, x='Freq_Percent_NR_All', y=ycol)
		for _, row in label_df.iterrows():
			plt.text(row['Freq_Percent_NR_All'], row[ycol], row["Gene"], fontsize=9,
				ha="left", va="bottom", color="darkblue")
		plt.title(f"{type_}: Freq_Percent_NR_All vs {ycol}")
		plt.xlabel('Freq_Percent_NR_All')
		plt.ylabel(ycol)
		plt.plot([0, merged_df[['Freq_Percent_NR_All',ycol]].max().max()],
				[0, merged_df[['Freq_Percent_NR_All',ycol]].max().max()],
				'r--')
		plt.tight_layout()
		plt.show()

	label_df = merged_df[(merged_df['Freq_Percent_PFS_RNA'] > 0.2) & (merged_df['Freq_Percent_PFS_Nano'] > 0.2)]
	plt.figure(figsize=(6,6))
	sns.scatterplot(data=merged_df, x='Freq_Percent_PFS_RNA', y='Freq_Percent_PFS_Nano')
	for _, row in label_df.iterrows():
		plt.text(row['Freq_Percent_PFS_RNA'], row['Freq_Percent_PFS_Nano'], row["Gene"], fontsize=9,
				ha="left", va="bottom", color="darkblue")
	plt.title(type_)
	plt.xlabel('Freq_Percent_PFS_RNA')
	plt.ylabel('Freq_Percent_PFS_Nano')
	plt.plot([0, merged_df[['Freq_Percent_PFS_RNA','Freq_Percent_PFS_Nano']].max().max()],
				[0, merged_df[['Freq_Percent_PFS_RNA','Freq_Percent_PFS_Nano']].max().max()],
				'r--')
	plt.tight_layout()
	plt.show()


####### CONSISTENT ACROSS EXPERIMENTS APPROACH #################
# From 11.3.25.ipynb
def get_gene_Fscores(data_set_dict_use):
	"""
	Get the F-scores for all genes within each experimental group (in data_set_dict_use) in predicting
	the y variable.
	"""
	gene_score_dict = {}
	for dataset, (X_use, y_use, y_clean_use, meta_use) in data_set_dict_use.items():   
		if dataset in ["Hold_Out", "Tr/Val", "All_Filt", "RNAseq", "Nano", "All_Assay"]:
			continue
		# Replace any NaN with per-gene median
		imputer = SimpleImputer(strategy="median")
		X_filled = pd.DataFrame(imputer.fit_transform(X_use), 
								columns=X_use.columns, index=X_use.index)
		
		# try to capture the top 20 genes with highest F classification scores for Nano vs RNA-seq
		selector = SelectKBest(f_classif, k=20)
		X_sel = selector.fit_transform(X_filled, y_use)
		# make a dataframe with the scores and genes
		gene_score_dict[dataset] = pd.DataFrame({
					"Gene": X_use.columns,
					"Score": selector.scores_
				}).sort_values(by="Score", ascending=False).reset_index(drop=True)
		combined_scores = pd.concat([
    	df.assign(Dataset=name) for name, df in gene_score_dict.items()
		])
		score_matrix = combined_scores.pivot(index="Gene", columns="Dataset", values="Score")
		# Median score per gene across datasets 
		median_scores = score_matrix.median(axis=1).sort_values(ascending=False)
	return(gene_score_dict, median_scores)

def get_top_sets(gene_score_dict_use, N=30):
	# Get the top 30 genes per dataset
	top_sets = {
		ds: set(df["Gene"].head(N))
		for ds, df in gene_score_dict_use.items()
	}
	return(top_sets)

#####################################################
#### TRANSFORMATION AND SCALING ###########
#####################################################
def get_scaler(scale):
    """
	Get the scaler to be used
    """
    if scale is None:
        return None
    if scale.lower() == "standard":
        return StandardScaler()
    if scale.lower() == "robust":
        return RobustScaler(
            with_centering=True,
            with_scaling=True,
            quantile_range=(25, 75)
        )
    raise ValueError(f"Unknown scale option: {scale}")

from sklearn.preprocessing import QuantileTransformer

def fit_cross_platform_transform(X_train, n_quantiles=50):
    """
    Learns a monotonic transform that aligns feature distributions
    without destroying relative ordering or sign.
    """
    qt = QuantileTransformer(
        n_quantiles=min(n_quantiles, X_train.shape[0]),
        output_distribution="normal",
        random_state=0
    )
    X_qt = qt.fit_transform(X_train)
    return qt, X_qt

def apply_model_with_transform(model, scaler, qt, X):
    """
    Applies the quantile transform -> optional scaler -> model
    Inputs: model (model used to predict), scaler (scaler item or None), qt (quantile from transformer object)
    """
    X_t = qt.transform(X)
    if scaler is not None:
        X_t = scaler.transform(X_t)
    probs = model.predict_proba(X_t)[:, 1]
    preds = (probs > 0.5).astype(int)
    return preds, probs

# WARNS IF THE LOGITS ARE HIGHLY COMPRESSED
def logit_qc(model, X, label):
    logits = model.decision_function(X)
    print(f"[{label}] Logit min/max:", logits.min(), logits.max())
    print(f"[{label}] Logit std:", logits.std())

    if logits.std() < 0.2:
        print(f"⚠️ WARNING: {label} logits are highly compressed")

#####################################################
#### TRAINING, VALIDATION, TESTING ###########
#####################################################
def train_and_eval_loocv(model, X, y, k, features_list=None, exhaustive=False, scale=None):
	"""
	Runs LOOCV on a model of interest to figure out the best set of genes to use for a given set of ks. 
	It considers all possible combinations of the features for a given K if exhaustive = True.
	It also saves the gene scores used based on SelectKBest (default is f_classif) but this will change with each LOOCV.
	Input:
	* model: Type of model being used (iwll already be a LogR preprared for example)
	* features_list: the list of features to consider (pre-filtering in the function)
	* all: True indicates that the full list of features in features_list are considered with no additional filtering
    Output:
    * results: pandas dataframe with K, gene subset, and auc
    """
	loo = LeaveOneOut()
	results = []
	
	if exhaustive:
		gene_sets = list(combinations(features_list, k))
	else:
		gene_sets = [features_list[0:k]]

	for gene_subset in gene_sets:
		y_true, y_pred = [], []
		for train_idx, test_idx in loo.split(X):
			X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
			y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
			# feature selection
			X_train = X_train[gene_subset]
			X_test = X_test[gene_subset]
            # get the scaler and scale data
			scaler = get_scaler(scale)
			if scaler is not None:
				X_train = scaler.fit_transform(X_train)
				X_test = scaler.transform(X_test)
            # train and predict
			model.fit(X_train, y_train)
			pred = model.predict_proba(X_test)[:, 1]
			y_true.append(y_test.values[0])
			y_pred.append(pred[0])
		# compute overall AUC
		auc = roc_auc_score(y_true, y_pred)
		results.append({
                "K": k,
                "Genes": gene_subset if gene_subset is not None else "TopK_from_train",
                "AUC": auc
            })
	
	# Get the best performing subset
	results_df = pd.DataFrame(results)

	best_row = results_df.loc[results_df["AUC"].idxmax()]
	best_auc = best_row["AUC"]
	best_genes = best_row["Genes"]

	# Compute summary statistics
	auc_min = results_df["AUC"].min()
	auc_max = results_df["AUC"].max()
	auc_range = auc_max - auc_min
	auc_std = results_df["AUC"].std()

	if exhaustive:
		print("LOOCV Feature Selection Summary")
		print(f"Best K = {k}")
		print(f"Best AUC = {best_auc:.4f}")
		print(f"Best Gene Subset = {best_genes}")
		print(f"\nAUC range: {auc_min:.4f} - {auc_max:.4f} (Δ = {auc_range:.4f})")
		print(f"AUC std: {auc_std:.4f}")
	return best_auc, best_genes

def train_eval_main(X, y, data_use="Nano", plot_data_set_dict_use=None, gene_set_sizes=[5, 10, 20, 25, 30, 35, 40, 45, 50], feature_list=None, exhaustive=False, scale="Robust", file_name=None):
    """
    Runs LOOCV to select best models (LogR L2 and SVM) and top genes.
    Plots ROC curves for all datasets and returns a summary of results for downstream analysis.
    """
    # ---------------------------
    # Run LOOCV for different gene set sizes for both Logistic Regression (L2) and SVM
    # ---------------------------
    results = []
    
    for k in gene_set_sizes:
        for model_name, model in [
            ("LogR-L1", LogisticRegression(max_iter=1000, penalty="l1", solver="liblinear")),
            ("LogR-L2", LogisticRegression(max_iter=1000, penalty="l2", solver="lbfgs")),
            ("SVM", SVC(probability=True, kernel="linear"))
        ]:
            best_auc, best_genes = train_and_eval_loocv(model, X, y, k, feature_list, exhaustive, scale)
            results.append((model_name, k, best_auc, best_genes))

    # Convert results to DataFrame
    res_df = pd.DataFrame(results, columns=["Model", "NumGenes", "AUC", "SelectedGenes"])
    display(res_df)

    # ---------------------------
    # Pick best model for SVM and L2 (top AUC but smallest # genes)
    # ---------------------------
    best_models = (
        res_df.sort_values(by=["AUC", "NumGenes"], ascending=[False, True])
        .groupby("Model")
        .first()
        .reset_index()
    )
    
    print("Best models per type:")
    print(best_models)

      # ---------------------------
    # Retrain best models on full initial dataset to get final coefficients (since LOOCV)
    # ---------------------------
    trained_models = {}
    selected_gene_sets = {}
    for _, row in best_models.iterrows():
        model_type = row["Model"]
        top_genes = row["SelectedGenes"]
        selected_gene_sets[model_type] = top_genes

        if model_type == "LogR-L1":
            model = LogisticRegression(max_iter=1000, penalty="l1", solver="liblinear")
        elif model_type == "LogR_L2":
            model = LogisticRegression(max_iter=1000, penalty="l2", solver="lbfgs")    
        else:
            model = SVC(probability=True, kernel="linear")
        X_sel = X[top_genes]
        scaler = get_scaler(scale)
        if scaler is not None:
             X_sel = scaler.fit_transform(X_sel)
        model.fit(X_sel, y)
        trained_models[model_type] = { "model": model, "scaler": scaler}
    
    # ---------------------------
    # Calculate and Plot AUCs for ALL datasets
    # ---------------------------
    auc_summary = {model_type: {} for model_type in ["LogR-L1", "LogR-L2", "SVM"]}
    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    list_ds_y, list_probs, list_patient, list_PFS, type_list = [], [], [], [], []
    for model_type, ax in zip(["LogR-L1", "LogR-L2", "SVM"], axes):
        model = trained_models[model_type]["model"]
        scaler = trained_models[model_type]["scaler"]
        top_genes = selected_gene_sets[model_type]
        # get nicer model_type name
        if model_type == "LogR-L1":
            model_name = "LASSO"
        elif model_type == "LogR-L2":
            model_name = "Ridge"
        else:
            model_name = "SVM"

        for dataset_name, (X_ds, y_ds, y_clean_ds, meta_ds) in plot_data_set_dict_use.items():
            if dataset_name == "All_Filt":
                continue
            # Only use the selected genes
            X_ds_sel = X_ds[top_genes]
            if scaler is not None:
                 X_ds_sel = scaler.transform(X_ds_sel)
            # Predict probabilities
            probs_test = model.predict_proba(X_ds_sel)[:,1]
						
						# print the differences and metadata
            if dataset_name in ["Hold_Out", "Tr/Val"]:
                list_ds_y.extend(y_ds)
                list_probs.extend(probs_test)
                list_patient.extend(meta_ds.Patient)
                list_PFS.extend(meta_ds.PFS_mths)
                model_name_list = [model_name]*len(y_ds)
                type_list.extend(model_name_list)

            # Compute ROC
            fpr, tpr, _ = roc_curve(y_ds, probs_test)
            auc_test = roc_auc_score(y_ds, probs_test)
            auc_summary[model_type][dataset_name] = auc_test

            # Plot
            if dataset_name in ["Nano", "RNAseq", "Tr/Val", "Hold_Out"]:
                ax.plot(fpr, tpr, label=f"{dataset_name} (AUC={auc_test:.2f})", color=color_dict[dataset_name])
            else:
                ax.plot(fpr, tpr, '--', label=f"{dataset_name} (AUC={auc_test:.2f})", color=color_dict[dataset_name], alpha=0.5)
				
        ax.plot([0,1],[0,1],'k--', linewidth=1, alpha=0.5)
        
        ax.set_xlabel("False Positive Rate", fontsize=17)
        ax.set_ylabel("True Positive Rate", fontsize=17)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_title(f"{model_name} ROC")
        ax.legend(fontsize=13, title_fontsize=14)

    plt.suptitle(f"ROC Curves for {data_use}-based Models on All Datasets")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    if file_name:
      plt.savefig(f"../plots/{file_name}_AUCs.png", transparent=True, bbox_inches='tight')
    plt.show()
		
	# Get the predictions for the samples to see if there is clear bias on where performing well
    res_full_df = pd.DataFrame({"Truth":list_ds_y, "Predicted":list_probs, "Patient":list_patient, "PFS_mths":list_PFS, "log_PFS_mths":np.log(np.array(list_PFS)+0.1), "Model":type_list})
    res_full_df["Predicted_Bool"] = res_full_df["Predicted"].apply(lambda x: 1 if x > 0.5 else 0)
    res_full_df["Match"] = (res_full_df["Predicted_Bool"] == res_full_df["Truth"]).astype(int)

	# 3. Assign confidence levels
    def classify_confidence(p):
        if p < 0.4:
            return "High Conf 0"
        elif p < 0.5:
            return "Low Conf 0"
        elif p < 0.6:
            return "Low Conf 1"
        else:
            return "High Conf 1"

    res_full_df["Confidence_Case"] = res_full_df["Predicted"].apply(classify_confidence)
		
	# ---------------------------
    # Compute sensitivity and specificity
    # ---------------------------
    def compute_sens_spec(df):
        #tn, fp, fn, tp = confusion_matrix(df["Truth"], df["Predicted_Bool"]).ravel()
        tn, fp, fn, tp = confusion_matrix(
						df["Truth"],
						df["Predicted_Bool"],
						labels=[0, 1]
				).ravel()
        sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
        spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
        return sens, spec

    metrics_summary = []
    for model in res_full_df["Model"].unique():
        df_m = res_full_df[res_full_df["Model"] == model]
        sens, spec = compute_sens_spec(df_m)
        metrics_summary.append([model, "Overall", round(sens, 2), round(spec, 2), df_m.shape[0]])
		# High Conf
        df_sub = df_m[(df_m["Predicted"] > 0.6) | (df_m["Predicted"] < 0.4)]
        sens, spec = compute_sens_spec(df_sub)
        metrics_summary.append([model, "High Conf", sens, spec,  df_sub.shape[0]])
        # Low Conf
        df_sub = df_m[(df_m["Predicted"] <= 0.6) & (df_m["Predicted"] >= 0.4)]
        sens, spec = compute_sens_spec(df_sub)
        metrics_summary.append([model, "Low Conf", sens, spec, df_sub.shape[0]])
    metrics_df = pd.DataFrame(metrics_summary, columns=["Model", "Confidence_Case", "Sensitivity", "Specificity", "N"])
		
    # ---------------------------
    # Return summary for downstream analysis
    # ---------------------------
    return {
        "best_models": best_models,               # Best model info per type
        "trained_models": trained_models,         # Fitted models retrained on full train set
        "selected_gene_sets": selected_gene_sets, # Top genes per best model
        "auc_summary": auc_summary,               # AUCs for each model/dataset combination
        "res_df": res_df,                         # All LOOCV results for comparison
		"res_full_df": res_full_df, 
		"metrics_df": metrics_df
    }

from sklearn.metrics import roc_curve, roc_auc_score, confusion_matrix

# Get the coefficient dataframes
def get_coef_dfs(te_results_):
	lasso_coef_df = pd.DataFrame({
		'Feature':     te_results_["selected_gene_sets"]["LogR-L1"],
		'Coefficient': te_results_["trained_models"]["LogR-L1"]["model"].coef_[0]
	})
	lasso_coef_df = lasso_coef_df.loc[lasso_coef_df.Coefficient != 0,]
	lasso_coef_df = lasso_coef_df.sort_values(by='Coefficient', ascending=True)
	lasso_coef_df.iloc[0:2,]

	ridge_coef_df = pd.DataFrame({
		'Feature':     te_results_["selected_gene_sets"]["LogR-L2"],
		'Coefficient': te_results_["trained_models"]["LogR-L2"]["model"].coef_[0]
	})
	ridge_coef_df = ridge_coef_df.loc[ridge_coef_df.Coefficient != 0,]
	ridge_coef_df = ridge_coef_df.sort_values(by='Coefficient', ascending=True)

	svm_coef_df = pd.DataFrame({
		'Feature':     te_results_["selected_gene_sets"]["SVM"],
		'Coefficient': te_results_["trained_models"]["SVM"]["model"].coef_[0]
	})
	svm_coef_df = svm_coef_df.loc[svm_coef_df.Coefficient != 0,]
	svm_coef_df = svm_coef_df.sort_values(by='Coefficient', ascending=True)
	return {"LASSO":lasso_coef_df, "Ridge":ridge_coef_df, "SVM":svm_coef_df}

#####################################################
#### PLOTTING COEFFICIENTS & GENES ###########
#####################################################


def graph_NR_coefs(coef_df, all_full_df, all_txn_df, corr_summary_df, title_use, scale=True, fig_width=15):
	experiment_order = ['RNAseq Adzib', 'RNAseq Jav', 'RNAseq EGA', 'Nano Manso', 'Nano Bitler', 'Nano James']  # user-defined order
	# Select only the features in the coefficient plot
	genes_to_use = coef_df['Feature'].tolist()
	median_df = all_full_df.groupby('Experiment')[genes_to_use].median()
	# get the second highest and lowest to avoid any extreme outliers
	second_highest = all_full_df[genes_to_use].apply(lambda x: x.dropna().nlargest(2).iloc[-1])
	second_lowest  = all_full_df[genes_to_use].apply(lambda x: x.dropna().nsmallest(2).iloc[-1])
	# combine
	maxmin_df = pd.concat([second_lowest, second_highest], axis=1)
	maxmin_df.columns = ["Min", "Max"]
	median_df = median_df.loc[experiment_order]  # reorder experiments
	median_txn_df = all_txn_df.groupby('Experiment')[genes_to_use].median()
	median_log_txn_df = median_txn_df.apply(lambda x: np.log10(x+0.01) if np.issubdtype(x.dtype, np.number) else x)
	corr_summary_filt_df = corr_summary_df.loc[genes_to_use]
	if scale:
		median_scaled = median_df.copy()
		median_scaled = (median_scaled - median_scaled.min()) / (median_scaled.max() - median_scaled.min())
		# Now transpose if needed for heatmap (genes as rows)
		median_plot = median_scaled.T  # genes x experiments

	
	fig, axes = plt.subplots(ncols=5, figsize=(fig_width, 7), gridspec_kw={'width_ratios': [1, 0.7, 0.4, 0.7, 0.4]})
   
	# --- Odds Ratio bar plot ---
	sns.barplot(
		x='Coefficient', 
		y='Feature', 
		data=coef_df, 
		palette='vlag',
		ax=axes[0]
	)
	axes[0].axvline(0, color='black', linestyle='--')  # 1 = no effect
	axes[0].set_xlabel('Coefficient (Effect on Probability of calling Nanostring)')
	axes[0].set_title(title_use)

	# --- Heatmap of median log2fc values ---
	# get the colors
	# Transpose if needed so genes are rows to match coefficient plot
	# Create a diverging color map
	if scale:
		cmap = sns.diverging_palette(240, 10, n=256, as_cmap=True)
		norm = TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=1)  # 0 = min, 0.5 = midpoint, 1 = max
	else:
		cmap = sns.diverging_palette(240, 10, n=256, as_cmap=True)  # blue-white-red
		median_plot = median_df.T  # genes x experiments
		vmin = median_plot.values.min()
		vmax = median_plot.values.max()
		if vmin >= 0:
			vmin = -0.05
		norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
	sns.heatmap(
		median_plot,  # transpose so genes are rows to match coef plot
		ax=axes[1],
		cmap=cmap,
		norm=norm,
		cbar_kws={'label': 'Median Expression'},
		yticklabels=True
	)
	axes[1].axvline(3, color='black', linestyle='--')  # 1 = no effect
	axes[1].set_xlabel('Experiment')
	axes[1].set_title('Median Log2FC by Experiment')

	# Max and Min
	norm = TwoSlopeNorm(vmin=maxmin_df.values.min(), vcenter=0, vmax=maxmin_df.values.max())
	sns.heatmap(
		maxmin_df,  # transpose so genes are rows to match coef plot
		ax=axes[2],
		cmap=cmap,
		norm=norm,
		cbar_kws={'label': 'log2FC'},
		yticklabels=True
	)
	axes[2].axvline(2, color='black', linestyle='--')  # 1 = no effect
	axes[2].set_xlabel('')
	axes[2].set_title('Min & Max L2FC')
	for i, gene in enumerate(maxmin_df.index):
		for j, l2fc in enumerate(maxmin_df.columns):
			if maxmin_df.loc[gene, l2fc]:
				axes[2].text(j + 0.5, i + 0.5, round(maxmin_df.loc[gene, l2fc], 1), color="black", ha="center", va="center", fontsize=9)


	# --- Heatmap of median txn values ---
	# get the colors
	# Transpose if needed so genes are rows to match coefficient plot
	# Create a diverging color map
	cmap_txn = LinearSegmentedColormap.from_list("white_red", ["white", "red", "darkred"])
	median_txn_plot = median_log_txn_df.T  # genes x experiments
	#median_txn_plot = median_txn_df.T  # genes x experiments
	vmin = median_txn_plot.values.min()
	vmax = median_txn_plot.values.max()
	print(vmin, vmax)
	if vmin >= 0:
			vmin = -0.05
	norm_txn = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
	#norm_txn = OneSlopeNorm(vmin=vmin, vmax=vmax)
	# Get Xs to add if they are below the cutoff
	# Create an empty mask of False values (same shape)
	mask = pd.DataFrame(False, index=median_txn_df.index, columns=median_txn_df.columns)
	mask.loc[mask.index.str.contains("RNAseq", case=False),:] = (
    	median_txn_df.loc[mask.index.str.contains("RNAseq", case=False),:] < 0.4
	)
	mask.loc[mask.index.str.contains("Nano", case=False),:] = (
    	median_txn_df.loc[mask.index.str.contains("Nano", case=False),:] < 23
	)
	mask_T = mask.T # transpose to have rows = genes

	sns.heatmap(
		median_txn_plot,  # transpose so genes are rows to match coef plot
		ax=axes[3],
		cmap=cmap_txn,
		#norm=norm_txn,
		cbar_kws={'label': 'log10(Median Expression)'},
		yticklabels=True
	)
	axes[3].axvline(3, color='black', linestyle='--')  # 1 = no effect
	axes[3].set_xlabel('Experiment')
	axes[3].set_title('Median Exp by Experiment')
	# Overlay X marks wherever mask is True
	# have the mask_T be in the same order
	for i, gene in enumerate(mask_T.index):
		for j, exp in enumerate(mask_T.columns):
			if mask_T.loc[gene, exp]:
				axes[3].text(j + 0.5, i + 0.5, "X", color="black", ha="center", va="center", fontsize=9, fontweight="bold")

	# Heatmap for Exon and Isoform Correlation
	corr_cmap = sns.diverging_palette(240, 10, n=256, as_cmap=True)  # blue-white-red
	vmin = corr_summary_filt_df.values.min()
	vmax = corr_summary_filt_df.values.max()
	if vmin < 0:
		corr_norm = TwoSlopeNorm(vmin=min(-0.2, vmin), vcenter=0.5, vmax=vmax)
	else:
		corr_norm = TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=vmax)
	
	sns.heatmap(
		corr_summary_filt_df,  # transpose so genes are rows to match coef plot
		ax=axes[4],
		cmap=corr_cmap,
		norm=corr_norm,
		cbar_kws={'label': 'Spearman Correlation Coefficient'},
		yticklabels=True
	)
	for i, gene in enumerate(corr_summary_filt_df.index):
		for j, feat in enumerate(corr_summary_filt_df.columns):
			if corr_summary_filt_df.loc[gene, feat]:
				axes[4].text(j + 0.5, i + 0.5, round(corr_summary_filt_df.loc[gene, feat], 1), color="black", ha="center", va="center", fontsize=9)

	axes[4].axvline(3, color='black', linestyle='--')  # 1 = no effect
	axes[4].set_xlabel('Count Type')
	axes[4].set_title('Correlation by Count Type')
	plt.tight_layout()
	plt.savefig(f"../plots/Btsp_Features/{title_use}_NR_Feature.pdf", transparent=True, bbox_inches='tight')

	plt.show()

def graph_coefs_PFS(coef_df, all_full_df, all_txn_df, corr_summary_df, title_use, scale=True, fig_width=15, RNAseqonly=False, save_fig=False):
	if RNAseqonly:
		experiment_order = ['RNAseq Adzib', 'RNAseq Jav', 'RNAseq EGA']  # user-defined order
	else:
		experiment_order = ['RNAseq Adzib', 'RNAseq Jav', 'RNAseq EGA', 'Nano Manso', 'Nano Bitler', 'Nano James']  # user-defined order

	# Select only the features in the coefficient plot
	genes_to_use = coef_df['Feature'].tolist()
	median_df = all_full_df.groupby(['Experiment', 'PFS_bool'])[genes_to_use].median()
	# Move PFS_bool into columns for subtraction
	median_unstacked = median_df.unstack('PFS_bool')
	pfs_log2fc_diff_df = median_unstacked.xs(">12mths", level=1, axis=1) - median_unstacked.xs("<=12mths", level=1, axis=1)
	# Transpose so genes are rows, experiments are columns
	pfs_log2fc_diff_df = pfs_log2fc_diff_df.T
	print(pfs_log2fc_diff_df.iloc[0:2,])
	# get the second highest and lowest to avoid any extreme outliers
	second_highest = all_full_df[genes_to_use].apply(lambda x: x.dropna().nlargest(2).iloc[-1])
	second_lowest  = all_full_df[genes_to_use].apply(lambda x: x.dropna().nsmallest(2).iloc[-1])
	# combine
	maxmin_df = pd.concat([second_lowest, second_highest], axis=1)
	maxmin_df.columns = ["Min", "Max"]
	pfs_log2fc_diff_df = pfs_log2fc_diff_df.loc[:,experiment_order]  # reorder experiments
	median_txn_df = all_txn_df.groupby('Experiment')[genes_to_use].median()
	median_log_txn_df = median_txn_df.apply(lambda x: np.log10(x+0.01) if np.issubdtype(x.dtype, np.number) else x)
	if RNAseqonly is False:
		corr_summary_filt_df = corr_summary_df.loc[genes_to_use]
	if scale:
		pfs_log2fc_diff_df_scaled = pfs_log2fc_diff_df.copy()
		pfs_log2fc_diff_df_scaled = (pfs_log2fc_diff_df_scaled - pfs_log2fc_diff_df_scaled.min()) / (pfs_log2fc_diff_df_scaled.max() - pfs_log2fc_diff_df_scaled.min())
		# Now transpose if needed for heatmap (genes as rows)
		pfs_log2fc_diff_df = pfs_log2fc_diff_df_scaled

	
	if RNAseqonly:
		fig, axes = plt.subplots(ncols=4, figsize=(fig_width, 7), gridspec_kw={'width_ratios': [1, 0.7, 0.4, 0.7]})
	else:
		fig, axes = plt.subplots(ncols=5, figsize=(fig_width, 7), gridspec_kw={'width_ratios': [1, 0.7, 0.4, 0.7, 0.4]})
   
	# --- Odds Ratio bar plot ---
	sns.barplot(
		x='Coefficient', 
		y='Feature', 
		data=coef_df, 
		palette='vlag',
		ax=axes[0]
	)
	axes[0].axvline(0, color='black', linestyle='--')  # 1 = no effect
	axes[0].set_xlabel('Coefficient (Effect on Probability of calling Nanostring)')
	axes[0].set_title(title_use)

	# --- Heatmap of >12mth - <=12mth median log2fc values ---
	# get the colors
	# Transpose if needed so genes are rows to match coefficient plot
	# Create a diverging color map
	if scale:
		cmap = sns.diverging_palette(240, 10, n=256, as_cmap=True)
		norm = TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=1)  # 0 = min, 0.5 = midpoint, 1 = max
	else:
		cmap = sns.diverging_palette(240, 10, n=256, as_cmap=True)  # blue-white-red
		vmin = pfs_log2fc_diff_df.values.min()
		vmax = pfs_log2fc_diff_df.values.max()
		if vmin >= 0:
			vmin = -0.05
		norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
	sns.heatmap(
		pfs_log2fc_diff_df,  # transpose so genes are rows to match coef plot
		ax=axes[1],
		cmap=cmap,
		norm=norm,
		cbar_kws={'label': 'PFS >12mth log2FC - PFS <=12mth log2FC'},
		yticklabels=True
	)
	if RNAseqonly is False:
		axes[1].axvline(3, color='black', linestyle='--')  # 1 = no effect
	axes[1].set_xlabel('Experiment')
	axes[1].set_title('Median Difference in Log2FC by Experiment')

	# Max and Min
	norm = TwoSlopeNorm(vmin=maxmin_df.values.min(), vcenter=0, vmax=maxmin_df.values.max())
	sns.heatmap(
		maxmin_df,  # transpose so genes are rows to match coef plot
		ax=axes[2],
		cmap=cmap,
		norm=norm,
		cbar_kws={'label': 'log2FC'},
		yticklabels=True
	)
	axes[2].axvline(2, color='black', linestyle='--')  # 1 = no effect
	axes[2].set_xlabel('')
	axes[2].set_title('Min & Max L2FC')
	for i, gene in enumerate(maxmin_df.index):
		for j, l2fc in enumerate(maxmin_df.columns):
			if maxmin_df.loc[gene, l2fc]:
				axes[2].text(j + 0.5, i + 0.5, round(maxmin_df.loc[gene, l2fc], 1), color="black", ha="center", va="center", fontsize=9)


	# --- Heatmap of median txn values ---
	# get the colors
	# Transpose if needed so genes are rows to match coefficient plot
	# Create a diverging color map
	cmap_txn = LinearSegmentedColormap.from_list("white_red", ["white", "red", "darkred"])
	median_txn_plot = median_log_txn_df.T  # genes x experiments
	#median_txn_plot = median_txn_df.T  # genes x experiments
	vmin = median_txn_plot.values.min()
	if vmin >= 0:
		vmin = -0.05
	vmax = median_txn_plot.values.max()
	print(vmin, vmax)
	norm_txn = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
	#norm_txn = OneSlopeNorm(vmin=vmin, vmax=vmax)
	# Get Xs to add if they are below the cutoff
	# Create an empty mask of False values (same shape)
	mask = pd.DataFrame(False, index=median_txn_df.index, columns=median_txn_df.columns)
	mask.loc[mask.index.str.contains("RNAseq", case=False),:] = (
    	median_txn_df.loc[mask.index.str.contains("RNAseq", case=False),:] < 0.4
	)
	if RNAseqonly is False:
		mask.loc[mask.index.str.contains("Nano", case=False),:] = (
				median_txn_df.loc[mask.index.str.contains("Nano", case=False),:] < 23
		)
	mask_T = mask.T # transpose to have rows = genes

	sns.heatmap(
		median_txn_plot,  # transpose so genes are rows to match coef plot
		ax=axes[3],
		cmap=cmap_txn,
		#norm=norm_txn,
		cbar_kws={'label': 'log10(Median Expression)'},
		yticklabels=True
	)
	if RNAseqonly is False:
		axes[3].axvline(3, color='black', linestyle='--')  # 1 = no effect
	axes[3].set_xlabel('Experiment')
	axes[3].set_title('Median Exp by Experiment')
	# Overlay X marks wherever mask is True
	# have the mask_T be in the same order
	for i, gene in enumerate(mask_T.index):
		for j, exp in enumerate(mask_T.columns):
			if mask_T.loc[gene, exp]:
				axes[3].text(j + 0.5, i + 0.5, "X", color="black", ha="center", va="center", fontsize=9, fontweight="bold")

	# Heatmap for Exon and Isoform Correlation
	if RNAseqonly is False:
		corr_cmap = sns.diverging_palette(240, 10, n=256, as_cmap=True)  # blue-white-red
		vmin = corr_summary_filt_df.values.min()
		vmax = corr_summary_filt_df.values.max()
		if vmin < 0:
			corr_norm = TwoSlopeNorm(vmin=min(-0.2, vmin), vcenter=0.5, vmax=vmax)
		else:
			corr_norm = TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=vmax)
	
		sns.heatmap(
			corr_summary_filt_df,  # transpose so genes are rows to match coef plot
			ax=axes[4],
			cmap=corr_cmap,
			norm=corr_norm,
			cbar_kws={'label': 'Spearman Correlation Coefficient'},
			yticklabels=True
		)
		for i, gene in enumerate(corr_summary_filt_df.index):
			for j, feat in enumerate(corr_summary_filt_df.columns):
				if corr_summary_filt_df.loc[gene, feat]:
					axes[4].text(j + 0.5, i + 0.5, round(corr_summary_filt_df.loc[gene, feat], 1), color="black", ha="center", va="center", fontsize=9)

		axes[4].axvline(3, color='black', linestyle='--')  # 1 = no effect
		axes[4].set_xlabel('Count Type')
		axes[4].set_title('Correlation by Count Type')
	plt.tight_layout()
	if save_fig:
		plt.savefig(f"../plots/Btsp_Features/{title_use}_PFS_Feature.pdf", transparent=True, bbox_inches='tight')
	plt.show()


def plot_conf_PFS(res_full_df, file_name=None):
    """
    Plot PFS distributions faceted by model and confidence category.
    Annotates each boxplot with number of Matches and Mismatches.
    """
    # Compute counts per Confidence_Case and Model
    counts_df = (
        res_full_df
        .groupby(["Model", "Confidence_Case", "Match"])
        .size()
        .reset_index(name="Count")
    )
    print(counts_df)

    # Pivot so we can easily see Matches vs Mismatches
    counts_pivot = (
        counts_df
        .pivot_table(index=["Model", "Confidence_Case"], 
                     columns="Match", values="Count", fill_value=0)
        .reset_index()
    )
    counts_pivot = counts_pivot.rename(columns={0: "Mismatches", 1: "Matches"})

    # Create the plot
    g = sns.catplot(
        data=res_full_df,
        x="Confidence_Case",
        y="log_PFS_mths",
        col="Model",
        kind="box",
        palette={
            "High Conf 1": "grey", 
            "Low Conf 1": "grey", 
            "High Conf 0": "#008B45", 
            "Low Conf 0": "#008B45"
        },
        height=5,
        aspect=0.9,
        order=["High Conf 0", "Low Conf 0", "Low Conf 1", "High Conf 1"]
    )

    # Add horizontal reference line and annotations
    for ax in g.axes.flatten():
        model = ax.get_title().replace("Model = ", "")
        ax.axhline(np.log(12), color="black", linestyle="--", linewidth=1)

        # Subset counts for this model
        model_counts = counts_pivot[counts_pivot["Model"] == model]

        # Annotate each category
        for i, cat in enumerate(["High Conf 0", "Low Conf 0", "Low Conf 1", "High Conf 1"]):
            if cat in model_counts["Confidence_Case"].values:
                row = model_counts[model_counts["Confidence_Case"] == cat].iloc[0]
                text = f'M: {round(row["Matches"])}\nMM: {round(row["Mismatches"])}'
                ax.text(
                    i, 
                    ax.get_ylim()[1] * 1.03,  # place near top
                    text, 
                    ha="center", 
                    va="top", 
                    fontsize=9, 
                    color="black"
                )

        ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.05)  # Adjust y-limit for annotations

    g.set_axis_labels("Prediction Confidence Category", "PFS (months, log scale)")
    plt.tight_layout()
    if file_name:
      plt.savefig(f"../plots/{file_name}.pdf", format="pdf")
    plt.show()
		
def graph_sens_spec(metrics_df, file_name=None):
	# Melt Sensitivity & Specificity into a long format for easier plotting
	metrics_long = metrics_df.melt(
		id_vars=["Model", "Confidence_Case", "N"],
		value_vars=["Sensitivity", "Specificity"],
		var_name="Metric",
		value_name="Value"
	)

	# Sort confidence cases if you want a specific order
	conf_order = ["Overall", "Low Conf", "High Conf"]
	metrics_long["Confidence_Case"] = pd.Categorical(metrics_long["Confidence_Case"], categories=conf_order, ordered=True)

	# Create the plot
	g = sns.catplot(
		data=metrics_long,
		x="Confidence_Case",
		y="Value",
		hue="Metric",
		col="Model",
		kind="bar",
		order=conf_order,
		height=3,
		aspect=0.8,
		palette=["#707071", "#BBBBBC"]
	)
	# Add N annotations on each bar
	for ax in g.axes.flatten():
		sub_df = metrics_long[metrics_long["Model"] == ax.get_title().split(' = ')[-1]]
		num = 0
		for container in ax.containers:
			if num == 1:
				num = 0
				continue
			else:
				num = 1
				for bar, conf_name in zip(container, conf_order):
					x = bar.get_x() + bar.get_width()
					# Match N for annotation
					n_val = sub_df[(sub_df["Confidence_Case"] == conf_name)]["N"].iloc[0]
					ax.text(x, 0.1, f"N={n_val}", ha="center", va="bottom", fontsize=8)
	# aesthetics
	plt.yticks(fontsize=10)
	plt.xticks(fontsize=10)
	if file_name:
		plt.savefig(f"../plots/{file_name}", format="pdf")
		plt.show()

#################################
####### PLOT GENE TRENDS ########
#################################
# both from Cons_Iter_Feature_11.3.25.
def plot_row_violin(matrices, row_name, matrix_names, show_points=True, log_transform=False, FC=False):
    """
    Plot box plots for a given row (by name) across multiple matrices.
    Optional log-transform of counts.
    
    Parameters
    ----------
    matrices : list of pandas.DataFrame
        List of matrices with the same row index labels.
    row_name : str
        The index label of the row to plot.
    matrix_names : list of str
        Names to use for each matrix in the plot.
    show_points : bool
        Whether to overlay raw points on the plot.
    log_transform : bool
        Whether to apply log(value + 0.0001) to counts before plotting.
    """
    if len(matrices) != len(matrix_names):
        raise ValueError("Number of matrices and matrix_names must match")
    
    dfs = []
    for mat, name in zip(matrices, matrix_names):
        if row_name not in mat.index:
            raise ValueError(f"Row '{row_name}' not found in matrix {name}")
        row_values = mat.loc[row_name].values.astype(float)
        if log_transform:
            row_values = np.log(row_values + 0.1)
        temp_df = pd.DataFrame({
            "Counts": row_values,
            "Matrix": [name] * len(row_values)
        })
        dfs.append(temp_df)
    
    df_all = pd.concat(dfs, ignore_index=True)
    
    plt.figure(figsize=(5,4))
    
    sns.boxplot(data=df_all, x="Matrix", y="Counts", palette=color_dict)
    
    if show_points:
        sns.stripplot(data=df_all, x="Matrix", y="Counts", color="black", jitter=True)
    
    title = f"Distribution of Counts for row '{row_name}'"
    if log_transform:
        title += " (log-transformed)"
        if FC:
            plt.ylabel("log(Log2FC+0.1)")
            plt.axhline(1, color="black")
        else:
            plt.ylabel("log(Exp+0.1)")
            plt.axhline(np.log(0 + 0.1), color="black")
            plt.axhline(np.log(23 + 0.1), color="blue")
            plt.axhline(np.log(0.4 + 0.1), color="red")
    else:
        if FC:
            plt.ylabel("Log2FC Post/Pre")
            plt.axhline(0, color="black")
        else:
            plt.ylabel("Exp")
            plt.axhline(0, color="black")
            plt.axhline(23, color="blue")
            plt.axhline(0.4, color="red")
		
    plt.title(title)
    plt.show()



def plot_model_gene_L2FC_facet(te_results, data_set_dict, max_genes=None, genes_include=[]):
    """
    Plot boxplots of Log2FC (Expression) across Classes for each gene,
    faceted per gene, showing all datasets for each gene.
    """
    for dataset_model, res in te_results.items():
        # Get genes from both models
        svm_genes = res["selected_gene_sets"].get("SVM", [])
        lasso_genes = res["selected_gene_sets"].get("LogR-L1", [])
        ridge_genes = res["selected_gene_sets"].get("LogR-L2", [])
        if max_genes is not None:
            svm_genes = svm_genes[:max_genes]
            lasso_genes = lasso_genes[:max_genes]
            ridge_genes = ridge_genes[:max_genes]

        # Build a dict mapping gene -> model source
        gene_model_map = {}
        all_genes = set(svm_genes) | set(lasso_genes) | set(ridge_genes)
        if genes_include != []:
            all_genes = genes_include
        for gene in all_genes:
            in_svm = gene in svm_genes
            in_lasso= gene in lasso_genes
            in_ridge= gene in ridge_genes
            if in_svm and in_lasso and in_ridge:
                gene_model_map[gene] = "All"
            elif in_svm and in_lasso:
                gene_model_map[gene] = "SVM_LASSO"
            elif in_svm and in_ridge:
                gene_model_map[gene] = "SVM_Ridge"
            elif in_lasso and in_ridge:
                gene_model_map[gene] = "LASSO_Ridge"
            elif in_lasso:
                gene_model_map[gene] = "LASSO"
            elif in_svm:
                gene_model_map[gene] = "SVM"
            else:
                gene_model_map[gene] = "Ridge"

        # Collect data across all datasets
        plot_data = []
        for dataset_name, (X, y, y_clean, _) in data_set_dict.items():
            genes_present = [g for g in all_genes if g in X.columns]
            if not genes_present:
                continue
            X_sub = X[genes_present].copy()
            X_sub["Class"] = y_clean.values
            df_melt = X_sub.melt(id_vars="Class", var_name="Gene", value_name="Log2FC")
            df_melt["Dataset"] = dataset_name
            df_melt["GeneLabel"] = df_melt["Gene"].map(lambda g: f"{g} ({gene_model_map[g]})")
            df_melt["ModelGroup"] = df_melt["Gene"].map(lambda g: gene_model_map[g])
            plot_data.append(df_melt)

        if not plot_data:
            continue

                

        df_plot = pd.concat(plot_data, ignore_index=True)
        # Order GeneLabel by ModelGroup: Both -> SVM -> LogR
        model_order = ["All", "LogR-L1", "LogR-L2", "SVM"]
        gene_order = df_plot.groupby("GeneLabel")["ModelGroup"].first().sort_values(key=lambda x: x.map({k:i for i,k in enumerate(model_order)})).index
        df_plot["Gene"] = pd.Categorical(df_plot["GeneLabel"], categories=gene_order, ordered=True)

        # FacetGrid by GeneLabel
        g = sns.catplot(
            data=df_plot,
            x="Dataset",
            y="Log2FC",
            hue="Class",
            col="Gene",
            kind="box",
            palette={"PFS>12":"#008B45", "PFS<=12":"grey"},
            col_wrap=3,
            sharey=False,
            height=3,
            aspect=1
        )

        # Rotate x-axis labels and only show for bottom row
        n_cols = g._col_wrap
        n_rows = (len(df_plot["Gene"].unique()) + n_cols - 1) // n_cols
        for i, ax in enumerate(g.axes.flat):
            ax.set_xticks(range(len(df_plot["Dataset"].unique())))
            ax.set_xticklabels(df_plot["Dataset"].unique(), rotation=45, ha="right")
            row = i // n_cols
            if row < n_rows - 1:
                ax.set_xlabel("")
        # Add horizontal line at y=0 to all facets
        for ax in g.axes.flat:
            ax.axhline(y=0, color='grey', linestyle='--')
        g.fig.suptitle(f"{dataset_model} Genes Selected by SVM/LogR", y=1)
        g.add_legend(title="Class")
        plt.tight_layout()
        #plt.savefig(f"../plots/SVM_RidgeR/{dataset_model}_RNAseqNano_SVMRidge_L2FCs.png", transparent=True, bbox_inches='tight')
        plt.show()

