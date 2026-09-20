# MVA Hackathon 2026 — Track 1: Variant Prediction

Pipeline used to identify candidate causal variant(s) for proband `WGS_EX2312012`,
using their WGS VCF and provided clinical phenotype.

## Approach summary

1. Loaded the proband VCF and reviewed quality (`FILTER`/`PASS`).
2. Used the provided clinical phenotype (rhabdomyosarcoma, nephrocalcinosis,
   short stature, failure to thrive, severe prematurity/IUGR, parental
   recurrent miscarriage) to hypothesise Mosaic Variegated Aneuploidy (MVA)
   syndrome, and restricted the search to its three known causative genes:
   **BUB1B, CEP57, TRIP13** (GRCh38 coordinates, UCSC Genome Browser).
3. Split candidates by genotype; homozygous candidates were checked manually
   against gnomAD by exact coordinates (not rsID — see note below) and all
   turned out to be common population variants (22–96% frequency), ruling
   them out.
4. Remaining heterozygous candidates were batch-annotated via the
   [Ensembl VEP web tool](https://www.ensembl.org/Tools/VEP) (GRCh38),
   with ClinVar, AlphaMissense, REVEL, CADD, SpliceAI, GERP, LOEUF, NMD,
   MaveDB, Geno2MP, OpenTargets, gnomAD Exomes/Genomes v4.1.1, and NIH
   All of Us annotations enabled.
5. Identified a **BUB1B stop-gained variant** (chr15:40209701 T>G,
   p.Leu737Ter) already classified Pathogenic in ClinVar for MVA syndrome 1,
   and a second ultra-rare **BUB1B missense variant** (chr15:40220612 T>G,
   p.Asn1002Lys) as a candidate compound-heterozygous partner.

Full reasoning and limitations are in the accompanying methods report
(`methods_description_form.xlsx` / exported PDF).

## Setup

```bash
pip install -r requirements.txt
huggingface-cli login   # requires access to the gated SageBio/mva-hackathon-2026-data dataset
```

## Running the pipeline

```bash
python scripts/01_filter_vcf.py          # download VCF, quality-check, filter to BUB1B/CEP57/TRIP13
python scripts/02_prepare_vep_input.py   # split by genotype, prepare VEP web-tool input
```

At this point, **manually** submit `outputs/vep_input.vcf` to the
[Ensembl VEP web tool](https://www.ensembl.org/Tools/VEP):
Species = Human, Assembly = GRCh38, with these annotation options enabled:
HGVS, AlphaMissense, REVEL, SpliceAI, CADD, GERP conservation, gnomAD
Exomes v4.1.1, gnomAD Genomes v4.1.1, NIH All of Us, NMD, LOEUF, MaveDB,
Phenotypes, Geno2MP, OpenTargets. Download the result as a VCF and save it
as `outputs/vep_output.vcf`, then:

```bash
python scripts/03_parse_vep_results.py outputs/vep_output.vcf
```

**Note on manual steps:** the homozygous-candidate gnomAD lookups and the
VEP web submission itself were done interactively, not fully scripted —
this is documented in the methods report. `02_prepare_vep_input.py`
automates the genotype-based split and file preparation, but the actual
frequency lookups and VEP submission are manual steps in this pipeline.

## ⚠️ Data access and redistribution

This pipeline operates on data from a **gated dataset**
(`SageBio/mva-hackathon-2026-data`) governed by a data use agreement that
restricts redistribution and requires deletion of local data after the
Hackathon concludes. Accordingly:

- **No raw or intermediate patient-derived data files are committed to this
  repository** (`outputs/` is gitignored). Running the scripts above will
  regenerate them locally, provided you have your own dataset access.
- The final submitted predictions (`predictions.csv`) are included, since
  the Hackathon's own submission terms explicitly cover sharing of
  submitted results.
- If you have access to this dataset and want to reproduce these results,
  request access via the dataset's own gating process on Hugging Face.

## AI assistance disclosure

This analysis was conducted with assistance from Claude (Anthropic). See
the methods report for the specific provider/plan/data-handling disclosure.
