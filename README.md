# CoreGene-Bayes

CoreGene-Bayes is a reproducible real-data research prototype for ranking genes
that satisfy two model-defined conditions:

1. genetically anchored gene-to-CAD Wald-ratio evidence is non-null; and
2. the gene receives probabilistic influences from at least `m` physically
   distinct source loci.

The case study starts from real LDL-associated loci, uses eQTLGen whole-blood
cis/trans summary statistics and a separate CAD GWAS, and performs a GTEx v8
liver audit of the SORT1/CELSR2/PSRC1 region.

**Project timeline:** the project concept and study design were initiated in
**July 2026**. The frozen real-data snapshot and its manifests retain their
actual September 2026 acquisition and analysis timestamps.

> The reported probabilities are exploratory model-based evidence conditional
> on the candidate set, tested pairs, priors, coverage, and unverified MR
> assumptions. They are not probabilities of biological truth, direct
> regulatory mechanisms, drug efficacy, or clinical utility.

## Study at a glance

| Component | Frozen choice |
|---|---|
| Candidate trait | LDL cholesterol |
| Candidate GWAS | GCST90019512 |
| Candidate rule | Genome-wide significant loci, Open Targets L2G ≥ 0.40, plus verified BayesV2G rank-one seeds |
| Disease outcome | CAD, GCST90132314 |
| Expression instruments | eQTLGen Phase I significant cis-eQTLs |
| Directed influence layer | Complete eQTLGen Phase I trans-tested coefficients |
| Discovery tissue | Whole blood |
| Locus audit | GTEx v8 liver, QTD000266 release 7 |
| Proxy annotation audit | PanelApp familial hypercholesterolaemia panel 772 v2.6, green genes |
| Primary definition | `P(D_j = 1 and K_j ≥ 2 | observed evidence)` |
| Sensitivity | `m = 3/5`, inclusion-prior grid, and slab-width grid |
| Genome builds | GRCh38 for CAD/liver; eQTLGen GRCh37 matched by rsID and signed alleles |

All primary analysis inputs are real public summary data. Synthetic effects are
used only in unit tests.

## Scientific design

```text
real LDL loci
  → prespecified candidate genes
  → cis-expression instruments
  → allele-frequency-scaled gene-to-CAD Wald-ratio evidence
  → complete tested trans coefficients
  → exploratory directed influence graph
  → conservative collapse into 1 Mb physical source loci
  → exact Poisson-binomial convergence probability
  → ratio evidence × convergence core-evidence probability
  → external PanelApp annotation and liver SORT1 audits
```

For gene `j`, `D_j` denotes non-null gene-to-CAD ratio evidence and `E_lj`
denotes non-null influence evidence from physical source locus `l`:

```text
K_j = Σ_l E_lj
C_j(m) = D_j × 1(K_j ≥ m)
```

The main output is `P(C_j(m)=1 | observed data, tested pairs, candidate set,
priors, model)`. Within a physical locus, the maximum edge probability is used
so correlated source genes are not treated as independent replication. Missing
pairs are never converted to negative edges.

See [Methods](docs/METHODS.md), [Limitations](docs/LIMITATIONS.md), and the
[reference list](docs/REFERENCES.md) before interpreting rankings.

## Real data sources

| Role | Dataset |
|---|---|
| LDL loci and candidates | Open Targets Platform release 26.09, GCST90019512 |
| Candidate-seed bridge | Verified frozen sibling BayesV2G `results/ldl_demo` artifact |
| CAD outcome | GWAS Catalog harmonized GCST90132314 |
| Whole-blood cis instruments | eQTLGen Phase I significant cis release |
| Whole-blood target coefficients | eQTLGen Phase I complete trans release |
| Liver locus audit | eQTL Catalogue release 7, GTEx v8 liver QTD000266 |
| Proxy annotation | Genomics England PanelApp panel 772 v2.6 |

Exact URLs, acquisition time, raw and compact file sizes, SHA-256 checksums,
and TLS status are frozen in
`data/snapshots/real_ldl_cad/source_manifest.json`. Dataset details are in
[Data provenance](docs/DATA_PROVENANCE.md).

## Install and test

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -e '.[dev,download]'
.venv/bin/pytest
```

For the pinned environment:

```bash
.venv/bin/python -m pip install -r requirements.lock
.venv/bin/python -m pip install -e . --no-deps
```

## Run the included compact real-data snapshot

```bash
.venv/bin/coregene-bayes run
```

The command verifies every required input against the source manifest before
analysis and writes `results/ldl_cad/`. Custom paths are supported:

```bash
.venv/bin/coregene-bayes run \
  --config config/ldl_cad.json \
  --data data/snapshots/real_ldl_cad \
  --output results/ldl_cad
```

## Rebuild the snapshot from public sources

```bash
.venv/bin/coregene-bayes fetch --cache .cache/coregene-bayes
.venv/bin/coregene-bayes run
```

The main compressed downloads are approximately 6.09 GB for eQTLGen trans,
0.32 GB for eQTLGen cis, and 1.30 GB for CAD. Allow at least 10–12 GB of free
space. Download and sequential compaction can take roughly 20–90 minutes on a
typical workstation. Valid raw downloads and completed compaction stages are
checksum-verified and reused after interruption.

Candidate ascertainment verifies the adjacent BayesV2G artifact. The default
layout is:

```text
workspace/
├── bayesv2g/
│   └── results/ldl_demo/
└── coregene-bayes/
```

An existing liver slice can be supplied with `--existing-liver-snapshot`.
`--insecure` disables TLS certificate verification and is recorded in the
manifest; use it only in a controlled environment with a certificate-chain
problem.

## Frozen real-data results

The verified run included with this project produced:

| Result | Value |
|---|---:|
| Candidate genes | 239 |
| Harmonized gene-to-CAD estimates | 173 |
| Network source genes / target genes | 114 / 182 |
| Tested directed source-target pairs | 20,493 |
| Physical source loci | 97 |
| Genes eligible for the primary ranking | 173 |
| Highest-ranked gene at `m=2` | HP (`P=0.769`) |
| Highest-ranked gene at `m=3` and `m=5` | LDLR (`P=0.543`, `P=0.209`) |
| Edges with model probability ≥0.5 / also FDR ≤0.05 | 19 / 8 |

The default prior alone implies `P(K≥2)=0.253` and
`P(C(m=2))=0.025` across 97 source loci. Prior/slab sensitivity is wide—for
example, HP ranges from 0.029 to 0.997 at `m=2`—so point rankings are
hypothesis-generating. LDLR is the only eligible overlap with the PanelApp
green proxy set and ranks third at `m=2`. In the separate liver audit, classic
single-signal coloc gives `PP(H4)=0.999831` for SORT1, with rs12740374 as the
top shared variant; CELSR2 and PSRC1 show similarly high regional sharing, so
this audit does not by itself distinguish mediation among the three genes.

The complete portable report is `results/ldl_cad/report.html`. Scope is narrow:
the ranking concerns core-like behavior among LDL-locus, whole-blood-testable
candidate genes, not all CAD genes.

## Output guide

Main result:

- `report.html` — narrative report, figures, top genes, edges, audits, and
  interpretation guardrails.
- `core_gene_posteriors.csv` — ranked `m=2/3/5` evidence probabilities,
  coverage, sensitivity range, and PanelApp overlap.
- `analysis_attrition.csv` and `run_summary.json` — analysis population and
  machine-readable run totals.

Evidence tables:

- `gene_disease_mr.csv` — strongest harmonizable cis-instrument Wald ratio per
  candidate, using real CAD allele frequencies for scale alignment.
- `gene_gene_edges.csv.gz` — tested directed pairs and exploratory edge
  probabilities.
- `locus_collapsed_edges.csv`, `network_coverage.csv`, and
  `gene_features.csv` — convergence inputs and missingness diagnostics.
- `prior_sensitivity.csv.gz` — gene × threshold × prior/slab grid.
- `harmonization_audit.csv.gz` — retained and rejected CAD allele matches.

Independent audits:

- `sort1_liver_cad_colocalization.csv` and
  `sort1_liver_cad_shared_variants.csv.gz` — liver regional MR/coloc audit.
- `validation_panelapp_genes.csv` and `validation_metrics.csv` — proxy
  annotation overlap/retrieval metrics, not gold-standard validation.

Reproducibility:

- `model.json` — exact event definition, priors, collapse rule, and prior-only
  baselines.
- `run_manifest.json` and `SHA256SUMS` — code/input/output provenance and file
  integrity.

See the [Data dictionary](docs/DATA_DICTIONARY.md) for column definitions.

Verify results on macOS:

```bash
cd results/ldl_cad
shasum -a 256 -c SHA256SUMS
```

## Interpretation guardrails

- The genome-wide screen uses one cis instrument and lacks regional
  colocalization; linkage and pleiotropy can mimic mediation.
- Whole blood is a scalable discovery context, not a complete tissue model of
  lipid biology or CAD.
- Source and target eQTLGen estimates share cohorts; unknown covariance can
  miscalibrate edge uncertainty.
- eQTLGen's trans SNPs were selected for trait associations, making horizontal
  pleiotropy particularly plausible.
- `Z/√N` eQTL effects are approximate; real CAD allele frequencies align
  genotype units.
- A 1 Mb physical cluster and maximum-edge collapse are conservative
  bookkeeping, not LD-aware joint inference.
- The prior alone creates nonzero convergence probability with many source
  loci; the report shows this baseline.
- PanelApp FH genes provide face validity, not labels for true CAD core genes.
- Lifelong genetic effects do not directly predict acute drug interventions.

This is a complete, auditable, hypothesis-generating MVP—not a clinically
validated gene-discovery system.

## Development, citation, and license

```bash
make test
make coverage
make build
```

The tests cover allele harmonization, effect-scale conversion, Wald ratios,
IVW, weighted median, MR-Egger, spike-and-slab probabilities, exact
Poisson-binomial tails, conservative locus collapse, colocalization,
configuration, reporting, and the BayesV2G artifact contract.

Use `CITATION.cff` to cite the software and cite each source dataset listed in
[Data provenance](docs/DATA_PROVENANCE.md). Project code is MIT licensed;
source-dataset terms remain controlling.
