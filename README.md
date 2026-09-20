# MVA Hackathon 2026 — Track 1

Pipeline for identifying the causal variant(s) in proband WGS_EX2312012.

## What I did

Started from the proband's VCF (5,012,204 variants) and did a quick PASS/FILTER
check. The clinical phenotype (rhabdomyosarcoma, nephrocalcinosis, short
stature, failure to thrive, severe IUGR, and a history of parental recurrent
miscarriage) pointed toward Mosaic Variegated Aneuploidy syndrome, so I
restricted the search to its three known causative genes — BUB1B, CEP57,
TRIP13 (GRCh38, coordinates from UCSC) — instead of scanning the whole genome.
That got it down to 21 variants.

From there: split by genotype, manually checked every homozygous candidate
against gnomAD using exact coordinates (rsID lookups occasionally redirected
to the wrong nearby variant — worth watching for). All of them turned out to
be common (22–96% frequency) and got ruled out. Batch-annotated the remaining
heterozygous candidates through Ensembl VEP, which surfaced a BUB1B
stop-gained variant (chr15:40209701 T>G, p.Leu737Ter) already flagged
Pathogenic in ClinVar for MVA1, plus a second ultra-rare BUB1B missense
variant (chr15:40220612 T>G, p.Asn1002Lys) as a likely compound-het partner.

**This pair was confirmed correct by the Hackathon's own scoring** (100/100
rank points, F-max 1.000, full match at rank 1).

Full writeup is in the methods report.

## Running it

    pip install -r requirements.txt
    huggingface-cli login
    python scripts/01_filter_vcf.py
    python scripts/02_prepare_vep_input.py

Then manually submit `outputs/vep_input.vcf` to Ensembl VEP
(ensembl.org/Tools/VEP) — Human, GRCh38 — with HGVS, AlphaMissense, REVEL,
SpliceAI, CADD, GERP, gnomAD Exomes/Genomes v4.1.1, All of Us, NMD, LOEUF,
MaveDB, Phenotypes, Geno2MP, and OpenTargets enabled. Save the result as
`outputs/vep_output.vcf`, then:

    python scripts/03_parse_vep_results.py outputs/vep_output.vcf

The gnomAD lookups and VEP submission were manual, not scripted — everything
else runs as-is.

## Data

The gated dataset this pipeline runs on can't be redistributed, so `outputs/`
(besides the final `predictions.csv`) is gitignored. You'll need your own
access to `SageBio/mva-hackathon-2026-data` to regenerate it.

## AI use

Used Claude (Anthropic) for debugging and explaining genomics concepts along
the way — see the methods report for the full disclosure.