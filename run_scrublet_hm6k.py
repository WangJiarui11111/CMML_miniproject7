import scrublet as scr
import scipy.io
import pandas as pd
import numpy as np

# Read matrix
# R exported genes x cells, Scrublet needs cells x genes
counts_matrix = scipy.io.mmread("hm_counts_filtered.mtx").T.tocsc()

genes = pd.read_csv("hm_genes_filtered.csv").iloc[:, 0].values
barcodes = pd.read_csv("hm_barcodes_filtered.csv").iloc[:, 0].values

dbr_df = pd.read_csv("hm_true_doublet_rate.csv")
expected_doublet_rate = float(dbr_df["true_dbr"].iloc[0])

print("Counts matrix shape:", counts_matrix.shape)
print("Number of barcodes:", len(barcodes))
print("Expected doublet rate:", expected_doublet_rate)

# Run Scrublet
scrub = scr.Scrublet(
    counts_matrix,
    expected_doublet_rate=expected_doublet_rate
)

doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=2,
    min_cells=3,
    min_gene_variability_pctl=85,
    n_prin_comps=30
)

# Save result
scrublet_result = pd.DataFrame({
    "barcode": barcodes,
    "Scrublet_score": doublet_scores,
    "Scrublet_pred": predicted_doublets
})

scrublet_result.to_csv("scrublet_result.csv", index=False)

print(scrublet_result.head())
print(scrublet_result["Scrublet_pred"].value_counts())