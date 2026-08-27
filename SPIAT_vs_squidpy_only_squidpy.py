import numpy as np
import pandas as pd
import anndata as ad
import squidpy as sq

df = pd.read_csv("/scratch/pawsey1172/dle/rstudio-dir/R/non_interactive_slurm/SD1_dfs/spe1.csv")

adata = ad.AnnData(
    X=np.zeros((len(df), 1)),  # placeholder expression matrix
    obs=pd.DataFrame(
        {"Cell.Type": df["Cell.Type"].values},
        index=df.index.astype(str)
    )
)

# Store spatial coordinates
adata.obsm["spatial"] = df[["Cell.X.Position", "Cell.Y.Position"]].values

# Optional: make cell type categorical
adata.obs["Cell.Type"] = adata.obs["Cell.Type"].astype("category")

sq.gr.co_occurrence(
    adata,
    cluster_key="Cell.Type"
)

occ_dict = adata.uns["Cell.Type_co_occurrence"]
distances = adata.uns["Cell.Type_co_occurrence"]["interval"][1:]

cell_types = adata.obs["Cell.Type"].cat.categories.tolist()

occ = occ_dict["occ"]

co_df = pd.DataFrame(
    [
        (reference, target, distance, value)
        for i, reference in enumerate(cell_types)
        for j, target in enumerate(cell_types)
        for k, distance in enumerate(distances)
        for value in [occ[i, j, k]]
    ],
    columns=["reference", "target", "distance", "value"]
)

co_df.to_csv('/scratch/pawsey1172/dle/rstudio-dir/R/non_interactive_slurm/co_stat.csv', index=False)

# Calculate L
ripley_L = sq.gr.ripley(
    adata,
    cluster_key="Cell.Type",
    mode="L",
    copy=True
)

# Calculate G
ripley_G = sq.gr.ripley(
    adata,
    cluster_key="Cell.Type",
    mode="G",
    copy=True
)

def extract_ripley(ripley_result, mode):
    df = ripley_result[f"{mode}_stat"].copy()

    # Rename distance column
    df = df.rename(columns={"bins": "distance", "Cell.Type": "reference", "stats": "value"})

    return df

L_df = extract_ripley(ripley_L, "L")
G_df = extract_ripley(ripley_G, "G")

L_df.to_csv('/scratch/pawsey1172/dle/rstudio-dir/R/non_interactive_slurm/L_stat.csv', index=False)

G_df.to_csv('/scratch/pawsey1172/dle/rstudio-dir/R/non_interactive_slurm/G_stat.csv', index=False)
