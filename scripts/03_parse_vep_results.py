"""
Step 3: Parse the annotated VCF returned by Ensembl VEP and extract the
key fields used to rank candidates: consequence, in-silico damage
scores, population frequency, and ClinVar clinical significance.

Usage:
    python 03_parse_vep_results.py path/to/vep_output.vcf

The CSQ field order below must match your own VEP job's header
(`##INFO=<ID=CSQ,...`) -- check this first if VEP's output format
differs (e.g. if you enabled different annotation options).
"""

import sys

TARGET_GENES = {"BUB1B", "CEP57", "TRIP13"}

# Update this to match the ##INFO=<ID=CSQ,...> line in your own VEP output
CSQ_FIELDS = (
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|"
    "INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|"
    "Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|"
    "SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|MANE_SELECT|MANE_PLUS_CLINICAL|"
    "GENCODE_PRIMARY|SOURCE|HGVS_OFFSET|MaveDB_doi|MaveDB_score|MaveDB_urn|"
    "OpenTargets_beta|OpenTargets_gwasDiseases|OpenTargets_gwasGeneId|"
    "OpenTargets_gwasLocusToGeneScore|OpenTargets_pValueExponent|"
    "OpenTargets_pValueMantissa|OpenTargets_qtlBiosampleName|"
    "OpenTargets_qtlGeneId|PHENOTYPES|am_class|am_pathogenicity|CADD_PHRED|"
    "CADD_RAW|REVEL|LOEUF|Geno2MP_HPO_count|SpliceAI_pred_DP_AG|"
    "SpliceAI_pred_DP_AL|SpliceAI_pred_DP_DG|SpliceAI_pred_DP_DL|"
    "SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|"
    "SpliceAI_pred_DS_DL|SpliceAI_pred_SYMBOL|NMD|Conservation|gnomAD_exomes|"
    "gnomAD_exomes_AF|gnomAD_genomes|gnomAD_genomes_AF|gnomAD_genomes_AF_grpmax|"
    "AoU|AoU_gvs_all_af|AoU_gvs_max_af|AoU_gvs_max_subpop|ClinVar|"
    "ClinVar_CLNDN|ClinVar_CLNDNINCL|ClinVar_CLNDISDB|ClinVar_CLNDISDBINCL|"
    "ClinVar_CLNREVSTAT|ClinVar_CLNSIG|ClinVar_CLNSIGCONF|ClinVar_CLNSIGINCL|"
    "ClinVar_ONCDN|ClinVar_ONCDNINCL|ClinVar_ONCDISDB|ClinVar_ONCDISDBINCL|"
    "ClinVar_ONC|ClinVar_ONCINCL|ClinVar_ONCREVSTAT|ClinVar_ONCSCV|"
    "ClinVar_ONCCONF|ClinVar_CLNSUBA|ClinVar_CLNPMID|ClinVar_CLNSUBN|"
    "ClinVar_CLNRCV|ClinVar_SCI|ClinVar_SCIREVSTAT|ClinVar_SCIDN|"
    "ClinVar_SCIDISDB|ClinVar_GENEINFO|SOURCE2"
).split("|")


def best_transcript_for_gene(transcripts, gene_symbols):
    """Prefer the canonical / MANE Select transcript for a target gene."""
    best = None
    for t in transcripts:
        d = dict(zip(CSQ_FIELDS, t.split("|")))
        if d.get("SYMBOL") not in gene_symbols:
            continue
        score = (2 if d.get("CANONICAL") == "YES" else 0) \
            + (2 if d.get("MANE_SELECT") else 0) \
            + (1 if d.get("HGVSc") else 0)
        if best is None or score > best[0]:
            best = (score, d)
    return best[1] if best else None


def parse_vep_vcf(path, gene_symbols=TARGET_GENES):
    results = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt, _, _, info = parts[:8]
            csq = next(x[4:] for x in info.split(";") if x.startswith("CSQ="))
            d = best_transcript_for_gene(csq.split(","), gene_symbols)
            if d is None:
                continue
            results.append({
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "gene": d.get("SYMBOL"), "consequence": d.get("Consequence"),
                "hgvsc": d.get("HGVSc"), "hgvsp": d.get("HGVSp"),
                "cadd_phred": d.get("CADD_PHRED"),
                "am_class": d.get("am_class"),
                "am_pathogenicity": d.get("am_pathogenicity"),
                "revel": d.get("REVEL"),
                "gnomad_exomes_af": d.get("gnomAD_exomes_AF"),
                "gnomad_genomes_af": d.get("gnomAD_genomes_AF"),
                "aou_af": d.get("AoU_gvs_all_af"),
                "clinvar_sig": d.get("ClinVar_CLNSIG"),
                "clinvar_disease": d.get("ClinVar_CLNDN"),
            })
    return results


if __name__ == "__main__":
    vep_path = sys.argv[1] if len(sys.argv) > 1 else "outputs/vep_output.vcf"
    for r in parse_vep_vcf(vep_path):
        print("=" * 90)
        for k, v in r.items():
            print(f"  {k}: {v}")
