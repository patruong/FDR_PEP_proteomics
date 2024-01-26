import pandas as pd
import numpy as np 
import os

os.chdir("/hd2/patrick/FDR_PEP_proteomics/reanalysis/percolator")

target = pd.read_csv("percolator.target.peptides.txt", sep = "\t")
decoy = pd.read_csv("percolator.decoy.peptides.txt", sep = "\t")
df = pd.concat([target, decoy], axis = 0)
file_idx_map = {
    0: "reanalysis/mzML/NEG1.mzML",
    1: "reanalysis/mzML/NEG1rep.mzML",
    2: "reanalysis/mzML/NEG2.mzML",
    3: "reanalysis/mzML/NEG2rep.mzML",
    4: "reanalysis/mzML/NEG3.mzML",
    5: "reanalysis/mzML/NEG3rep.mzML",
    6: "reanalysis/mzML/NEG4.mzML",
    7: "reanalysis/mzML/NEG4rep.mzML",
    8: "reanalysis/mzML/NEG5.mzML",
    9: "reanalysis/mzML/NEG5rep.mzML",
    10: "reanalysis/mzML/POS1.mzML",
    11: "reanalysis/mzML/POS1rep.mzML",
    12: "reanalysis/mzML/POS2.mzML",
    13: "reanalysis/mzML/POS2rep.mzML",
    14: "reanalysis/mzML/POS3.mzML",
    15: "reanalysis/mzML/POS3rep.mzML",
    16: "reanalysis/mzML/POS4.mzML",
    17: "reanalysis/mzML/POS4rep.mzML",
    18: "reanalysis/mzML/POS5.mzML"
}
df["run"] = df["file_idx"].map(file_idx_map)
df["target_decoy_FDR"] = calculate_targetdecoy_fdr(df)

#len(target)/len(df)
#len(decoy)/len(df)

def calculate_targetdecoy_fdr(df):
    df.sort_values("percolator score", inplace = True, ascending = False)
    df["labels"] =  df["protein id"].apply(lambda x: 0 if "decoy" in x else 1)
    # Sort scores and labels
    sorted_labels = df["labels"]
    # Calculate FDR
    cumulative_decoy_hits = np.cumsum(sorted_labels == 0)
    cumulative_target_hits = np.cumsum(sorted_labels == 1)
    fdr = cumulative_decoy_hits / cumulative_target_hits
    return fdr

def triqler_qVal_searchScore(df):
    triqler_input = pd.DataFrame()
    triqler_input["run"] = df["run"]
    triqler_input["condition"] = df["run"].apply(lambda x: "positive" if "POS" in x else "negative")
    triqler_input["charge"] = df["charge"]
    triqler_input["searchScore"] = -np.log10(df["percolator q-value"]) # this should be change to different search scores.
    triqler_input["intensity"] = df["spectrum precursor m/z"]
    triqler_input["peptide"] = df["sequence"]
    triqler_input["protein"] = df["protein id"]
    return triqler_input

def triqler_targetDecoy_FDR_searchScore(df):
    triqler_input = pd.DataFrame()
    triqler_input["run"] = df["run"]
    triqler_input["condition"] = df["run"].apply(lambda x: "positive" if "POS" in x else "negative")
    triqler_input["charge"] = df["charge"]
    triqler_input["searchScore"] = -np.log10(df["target_decoy_FDR"]) # this should be change to different search scores.
    triqler_input["intensity"] = df["spectrum precursor m/z"]
    triqler_input["peptide"] = df["sequence"]
    triqler_input["protein"] = df["protein id"]
    return triqler_input

def triqler_PEP_searchScore(df):
    triqler_input = pd.DataFrame()
    triqler_input["run"] = df["run"]
    triqler_input["condition"] = df["run"].apply(lambda x: "positive" if "POS" in x else "negative")
    triqler_input["charge"] = df["charge"]
    triqler_input["searchScore"] = -np.log10(df["percolator PEP"]) # this should be change to different search scores.
    triqler_input["intensity"] = df["spectrum precursor m/z"]
    triqler_input["peptide"] = df["sequence"]
    triqler_input["protein"] = df["protein id"]
    return triqler_input

def triqler_percolatorScore_searchScore(df):
    triqler_input = pd.DataFrame()
    triqler_input["run"] = df["run"]
    triqler_input["condition"] = df["run"].apply(lambda x: "positive" if "POS" in x else "negative")
    triqler_input["charge"] = df["charge"]
    triqler_input["searchScore"] = df["percolator score"] # this should be change to different search scores.
    triqler_input["intensity"] = df["spectrum precursor m/z"]
    triqler_input["peptide"] = df["sequence"]
    triqler_input["protein"] = df["protein id"]
    return triqler_input





np.corrcoef(df["percolator PEP"], df["percolator q-value"])
np.corrcoef(df["percolator PEP"], df["percolator score"])
np.corrcoef(df["percolator q-value"], df["percolator score"])
np.corrcoef(df["percolator q-value"], df["target_decoy_FDR"])




import numpy as np

# Hypothetical target and decoy scores
target_scores = np.array([5, 4, 3, 2, 1])
decoy_scores = np.array([10, 9, 8, 7, 6])

# Concatenate scores and labels (1 for target, 0 for decoy)
all_scores = np.concatenate([target_scores, decoy_scores])
labels = np.concatenate([np.ones_like(target_scores), np.zeros_like(decoy_scores)])

# Sort scores and labels
sorted_indices = np.argsort(all_scores)[::-1]
sorted_labels = labels[sorted_indices]

all_scores
sorted_indices
sorted_labels
# Calculate FDR
cumulative_decoy_hits = np.cumsum(sorted_labels == 0)
cumulative_target_hits = np.cumsum(sorted_labels == 1)
fdr = cumulative_decoy_hits / cumulative_target_hits

print("FDR at each score:", fdr)