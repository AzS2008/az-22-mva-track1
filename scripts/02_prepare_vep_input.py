"""
Step 2: Split gene-filtered candidates by genotype, and prepare a
VCF-format text block for submission to the Ensembl VEP web tool
(https://www.ensembl.org/Tools/VEP).

Rationale: since Mosaic Variegated Aneuploidy is autosomal recessive,
homozygous (1/1) candidates were prioritised first and checked manually
against gnomAD (see README -- this step is NOT automated here, since it
was done interactively via the gnomAD website using exact
chr-pos-ref-alt coordinates, not rsID, after finding rsID lookups could
occasionally redirect to a different nearby variant).

All homozygous candidates in this case turned out to be common
population variants (22-96% allele frequency) upon manual gnomAD lookup
and were excluded. The remaining heterozygous candidates were batch
submitted to VEP instead -- this script prepares that input.
"""

import pandas as pd


def split_by_genotype(candidates_path="outputs/candidates.csv", sample_col=None):
    df = pd.read_csv(candidates_path, dtype={"#CHROM": str})
    if sample_col is None:
        # sample column is whatever comes after FORMAT -- last column
        sample_col = df.columns[-1]

    genotype = df[sample_col].str.split(":", n=1).str[0]
    homozygous = df[genotype.isin(["1/1", "1|1"])]
    heterozygous = df[genotype.isin(["0/1", "0|1", "1|0"])]

    print(f"Homozygous candidates: {len(homozygous)} "
          "(checked manually against gnomAD by exact coordinates -- see README)")
    print(f"Heterozygous candidates: {len(heterozygous)} (submitted to VEP below)")

    return homozygous, heterozygous


def write_vep_input(heterozygous_df, out_path="outputs/vep_input.vcf"):
    """Write a minimal VCF-format block VEP's web tool accepts directly."""
    with open(out_path, "w") as f:
        for _, row in heterozygous_df.iterrows():
            f.write(f"{row['#CHROM']}\t{row['POS']}\t{row['ID']}\t"
                     f"{row['REF']}\t{row['ALT']}\t{row['QUAL']}\t"
                     f"{row['FILTER']}\t.\n")
    print(f"Saved {out_path} -- paste this into the VEP web tool "
          "(Species: Human, Assembly: GRCh38).")


if __name__ == "__main__":
    homozygous, heterozygous = split_by_genotype()
    write_vep_input(heterozygous)
