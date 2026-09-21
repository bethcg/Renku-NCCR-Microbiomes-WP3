"""
Data loading + diversity/ordination for the WP3 wheat-rhizosphere dashboard.

Kept separate from the Streamlit UI (microbiome_app.py) so the analysis logic
can be imported and tested without launching a server. This mirrors the R
phyloseq/vegan workflow in analysis/ using pandas/scipy.

Data: Garrido-Sanz & Keel, "Sequential propagation of a reproducible wheat
rhizosphere microbiome" (Keel lab, UNIL - NCCR Microbiomes WP3).
Zenodo 10.5281/zenodo.14514438 (CC-BY-4.0); mirror github.com/dgarrs/RhizCom.
"""
from __future__ import annotations
import os
import urllib.request
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

# repo-root/data, regardless of where the app is launched from
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
_MIRROR = "https://raw.githubusercontent.com/dgarrs/RhizCom/main/Amplicon_sequence_analyses"
_FILES = {  # local name -> remote name (same mapping as analysis/00_get_data.R)
    "asv_counts.txt": "ASV_sequences.txt",
    "taxonomy.txt": "Taxtable_dada2.txt",
    "metadata.txt": "metadata.txt",
}


def ensure_data(data_dir: str = DATA_DIR) -> None:
    """Download the data files from the archived mirror if they aren't present.

    The repo ships the analysis but not the (already-public) data, so a fresh
    session fetches it once - the same pattern as analysis/00_get_data.R. On
    RenkuLab you can instead mount the Zenodo record as a data connector.
    """
    os.makedirs(data_dir, exist_ok=True)
    for local, remote in _FILES.items():
        dest = os.path.join(data_dir, local)
        if not os.path.exists(dest):
            urllib.request.urlretrieve(f"{_MIRROR}/{remote}", dest)


def load(data_dir: str = DATA_DIR):
    """Return (counts [samples x ASVs], metadata, taxonomy), sample-aligned."""
    ensure_data(data_dir)
    counts = pd.read_csv(os.path.join(data_dir, "asv_counts.txt"), sep=";", index_col=0)
    counts.index = counts.index.str.replace(r"_$", "", regex=True)

    meta = pd.read_csv(os.path.join(data_dir, "metadata.txt"), sep="\t").set_index("id_samples")
    meta["RGBcol"] = meta["RGBcol"].astype(str).str.replace('"', "", regex=False)

    tax = pd.read_csv(os.path.join(data_dir, "taxonomy.txt"), sep=";", index_col=0)

    common = counts.index.intersection(meta.index)
    return counts.loc[common], meta.loc[common], tax


def filter_counts(counts: pd.DataFrame, min_prevalence: int = 2, min_total: int = 10) -> pd.DataFrame:
    """Drop noise ASVs: present in < min_prevalence samples or total < min_total."""
    keep = ((counts > 0).sum(axis=0) >= min_prevalence) & (counts.sum(axis=0) >= min_total)
    return counts.loc[:, keep]


def _shannon(x: np.ndarray) -> float:
    p = x[x > 0] / x.sum()
    return float(-(p * np.log(p)).sum())


def alpha_table(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Per-sample Observed richness and Shannon diversity, joined to metadata."""
    df = pd.DataFrame({
        "Observed": (counts > 0).sum(axis=1),
        "Shannon": counts.apply(_shannon, axis=1),
    })
    df = df.join(meta[["Name", "Order2", "RGBcol"]])
    return df.reset_index().rename(columns={"index": "id_samples"}).sort_values("Order2")


def braycurtis_pcoa(counts: pd.DataFrame, meta: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    """Relative-abundance Bray-Curtis PCoA (classical MDS). Returns (coords_df, %var)."""
    rel = counts.div(counts.sum(axis=1), axis=0)
    D = squareform(pdist(rel.values, metric="braycurtis"))
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J.dot(D ** 2).dot(J)
    vals, vecs = np.linalg.eigh(B)
    idx = np.argsort(vals)[::-1]
    vals, vecs = vals[idx], vecs[:, idx]
    pos = vals > 1e-9
    coords = vecs[:, pos] * np.sqrt(vals[pos])
    var = vals[pos] / vals[pos].sum() * 100
    out = pd.DataFrame({"PCo1": coords[:, 0], "PCo2": coords[:, 1]}, index=counts.index)
    out = out.join(meta[["Name", "Order2", "RGBcol"]])
    return out.reset_index().rename(columns={"index": "id_samples"}), var[:2]


def group_order(meta: pd.DataFrame) -> list[str]:
    return meta.drop_duplicates("Name").sort_values("Order2")["Name"].tolist()


def palette(meta: pd.DataFrame) -> dict[str, str]:
    d = meta.drop_duplicates("Name").set_index("Name")["RGBcol"]
    return {k: (v if str(v).startswith("#") else f"#{v}") for k, v in d.items()}
