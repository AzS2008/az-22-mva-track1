"""
Step 1: Download the proband VCF, run quality checks, and filter to the
three genes hypothesised to cause the proband's phenotype (Mosaic
Variegated Aneuploidy syndrome): BUB1B, CEP57, TRIP13.

Requires: a Hugging Face account with access to SageBio/mva-hackathon-2026-data
(gated dataset -- run `huggingface-cli login` first, or set HF_TOKEN env var).

Output: outputs/candidates.csv (NOT committed to this repo -- see README
for why real patient-derived variant data is excluded from version control).
"""

import gzip
import pandas as pd
from huggingface_hub import hf_hub_download

# --- GRCh38 coordinates for the three MVA-causing genes (UCSC Genome Browser) ---
GENE_COORDS = {
    "BUB1B": {"chrom": "15", "start": 40161023, "end": 40221123},
    "CEP57": {"chrom": "11", "start": 95789965, "end": 95832693},
    "TRIP13": {"chrom": "5", "start": 892884, "end": 919357},
}

VCF_FILENAME = "WGS_EX2312012_HGWCNDSX7.vcf.gz"
DTYPES = {
    "#CHROM": str, "POS": "int64", "ID": str, "REF": str, "ALT": str,
    "QUAL": "float64", "FILTER": str, "INFO": str, "FORMAT": str,
}


def download_vcf():
    return hf_hub_download(
        repo_id="SageBio/mva-hackathon-2026-data",
        filename=VCF_FILENAME,
        repo_type="dataset",
        local_dir="./data",
    )


def get_header_and_build(file_path):
    """Extract column names and confirm the reference genome build."""
    columns = None
    with gzip.open(file_path, "rt") as f:
        for line in f:
            if line.startswith("##reference"):
                print("Reference genome line:", line.strip())
            if line.startswith("#CHROM"):
                columns = line.strip().split("\t")
                break
    return columns


def load_vcf(file_path, columns):
    dtypes = dict(DTYPES)
    dtypes[columns[-1]] = str  # sample column, whatever it's named
    return pd.read_csv(file_path, sep="\t", comment="#", header=None,
                        names=columns, dtype=dtypes)


def filter_to_candidate_genes(df):
    """Build one boolean mask per gene, then combine with OR."""
    mask = pd.Series(False, index=df.index)
    for gene, coords in GENE_COORDS.items():
        gene_mask = (
            (df["#CHROM"].astype(str) == coords["chrom"])
            & (df["POS"] >= coords["start"])
            & (df["POS"] <= coords["end"])
        )
        print(f"{gene}: {gene_mask.sum()} variant(s)")
        mask |= gene_mask
    return df[mask].copy()


if __name__ == "__main__":
    vcf_path = download_vcf()
    columns = get_header_and_build(vcf_path)
    df = load_vcf(vcf_path, columns)

    print(f"\nTotal variants: {df.shape[0]}")
    print("\nFILTER breakdown:")
    print(df["FILTER"].value_counts())

    candidates = filter_to_candidate_genes(df)
    print(f"\nCandidates after gene-coordinate filtering: {candidates.shape[0]}")

    candidates.to_csv("outputs/candidates.csv", index=False)
    print("Saved outputs/candidates.csv")
