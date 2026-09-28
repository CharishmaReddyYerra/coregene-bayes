# Methods

## Study question and prespecified scope

CoreGene-Bayes asks a deliberately operational question:

> Among genes surfaced at real LDL-associated loci, which genes both have
> genetically anchored expression-to-coronary-artery-disease (CAD) ratio
> evidence and receive expression associations from multiple physically
> distinct, genetically perturbed source loci?

The primary candidate trait is LDL cholesterol; the outcome is CAD. Candidate
selection is completed without reading the disease-outcome estimates used in
the model. Whole-blood eQTLGen data supply a scalable discovery network. A separate
GTEx v8 liver analysis audits the SORT1/CELSR2/PSRC1 locus. RNA is the only
molecular layer in the model; the adjacent BayesV2G project's plasma pQTL data
are not silently relabelled as eQTLs or used as replication.

This is a research prototype, not a clinical predictor or a therapeutic claim.

## Candidate ascertainment

All Open Targets `GCST90019512` LDL credible sets are retrieved from release
26.09. Candidate locus-gene rows must satisfy both:

- lead association p-value ≤ 5 × 10⁻⁸; and
- Open Targets Locus-to-Gene score ≥ 0.40.

The rank-one gene from each of the ten frozen BayesV2G demonstration loci is
also retained so the two projects have a verified artifact-level connection.
BayesV2G probabilities are used only for ascertainment. They are not entered
again as causal evidence because BayesV2G already includes molecular-QTL
colocalization features; reuse would double-count evidence.

When a gene occurs at more than one locus, the row with the largest L2G score
is retained as node metadata. Candidate selection is performed without reading
the CAD estimates or eQTLGen trans associations.

## Data sources

- CAD: GWAS Catalog `GCST90132314`, 181,522 cases and 984,168 controls,
  harmonized GRCh38 summary statistics.
- Whole-blood QTLs: eQTLGen Phase I, up to 31,684 predominantly European
  participants. The significant cis table is used to identify instruments;
  the **complete**, not significance-filtered, trans table provides target
  coefficients for all tested source-SNP/target-gene pairs.
- Liver locus audit: eQTL Catalogue release 7 `QTD000266`, GTEx v8 liver,
  n=208, GRCh38. The release is pinned because the QTD identifier was reused in
  the 2026 release-8 prerelease for GTEx v10.
- Candidate loci/features: Open Targets Platform release 26.09.
- External proxy annotation audit: Genomics England PanelApp 772, signed-off
  version 2.6, confidence-level-3 genes.

Exact URLs, fetch time, TLS status, byte counts, and SHA-256 checksums are in
`data/snapshots/real_ldl_cad/source_manifest.json`.

## Harmonization

MR matching uses rsID and signed alleles because eQTLGen positions are GRCh37,
whereas the harmonized CAD GWAS and liver QTL coordinates are GRCh38. Coordinate
equality is never used across builds.

For every matched record, the CAD or target-expression effect is oriented to
the exposure assessed allele. Direct, swapped, strand-complemented, and
strand-complemented/swapped states are distinguished. A/T and C/G instruments
are excluded because eQTLGen Phase-I does not provide the effect-allele
frequency required to resolve their orientation safely. Incompatible alleles
remain in the harmonization audit rather than disappearing silently.

The eQTLGen public tables report Z statistics and per-pair sample size rather
than beta and standard error. The analysis uses the genotype-standardized
approximation

```text
beta ≈ Z / sqrt(N)
SE   ≈ 1 / sqrt(N).
```

This is an association per approximate standard deviation of genotype, not a
per-allele effect. For gene-to-CAD MR, the real CAD effect-allele frequency
`p` therefore transforms the per-allele CAD association onto the same scale:

```text
beta_Y,standardized = beta_Y,per-allele × sqrt(2p(1-p))
SE_Y,standardized   = SE_Y,per-allele   × sqrt(2p(1-p)).
```

Rows without a valid CAD allele frequency are not estimated. Every affected
output labels this scale. It is a first-order representation, not a recovered
effect on native expression units.

## Gene-to-CAD evidence

For each candidate gene, up to the 20 strongest non-palindromic significant
cis associations are retained during data compaction. After CAD matching and
allele harmonization, the strongest available instrument by F statistic is
used. The real-data MVP therefore uses the single-instrument Wald ratio:

```text
theta_hat = beta_Y / beta_X.
```

The full first-order delta variance propagates uncertainty in both associations:

```text
Var(theta_hat) = SE_Y² / beta_X²
               + beta_Y² SE_X² / beta_X⁴
               - 2 beta_Y Cov(beta_X,beta_Y) / beta_X³.
```

For the separate-study gene-to-CAD screen, covariance is set to zero; exact
participant non-overlap is not asserted. F is
`(beta_X / SE_X)²`; instruments below F=10 fail QC. Benjamini-Hochberg FDR is
reported over estimable genes. Genome-wide regional colocalization is not
available, so these are Wald-ratio signals compatible with an expression
effect, not confirmed expression mediation. Linkage and horizontal pleiotropy
remain live explanations. This one-instrument screen cannot estimate
MR-Egger pleiotropy, heterogeneity, or leave-one-instrument-out diagnostics.

The package implements IVW, weighted median, MR-Egger, heterogeneity, and
signed harmonization with unit tests, but does not manufacture their use where
the real-data instrument design cannot support them.

## Gene-to-gene directed influence graph

For every candidate source gene, the strongest non-palindromic cis association
that belongs to eQTLGen's fixed trans-tested SNP universe is selected. Complete
target coefficients are then extracted for candidate genes. Consequently,
null target associations remain represented; edges are not ascertained from a
significant-trans-only file.

Each source-target pair is estimated with the same Wald ratio using approximate
standardized eQTLGen effects. The source and target estimates come from
overlapping, often identical cohorts. Their unknown covariance is set to zero,
so standard errors and edge probabilities are exploratory and may be
miscalibrated. An edge is called a **genetically anchored directed influence**,
not a direct gene-regulatory mechanism. Horizontal pleiotropy, source-gene
ambiguity, mediation, and feedback remain unresolved.

Nearby source genes are not counted as independent inputs. Candidate genes are
greedily clustered by chromosome with a 1 Mb linkage window. Within one such
physical source locus, edge probabilities are collapsed conservatively as

```text
q_lj = max_i(q_ij).
```

This avoids treating correlated, and sometimes identical, within-locus edge
evidence as independent replication. The maximum is a conservative evidence
rule, not a joint posterior for the locus event. The physical cluster is
bookkeeping; it is not a measured LD block.

## Spike-and-slab evidence probabilities

For either a gene-to-CAD estimate or an edge estimate,

```text
theta_hat | theta ~ Normal(theta, SE²)
```

is combined with a point-normal prior:

```text
H0: theta = 0
H1: theta ~ Normal(0, tau²)
P(H1) = pi.
```

The closed-form normal-normal Bayes factor yields
`P(H1 | theta_hat, SE, pi, tau)`. Defaults are:

- disease: `pi=0.10`, `tau=0.20`;
- gene-to-gene edge: `pi=0.01`, `tau=0.15`.

These are evidence-inclusion probabilities under a specified model. They do
not include the probability that the instrumental-variable assumptions hold.

## Label-free core definition

Let `D_j` indicate non-null gene-to-CAD Wald-ratio evidence under the specified
spike-and-slab model and `E_lj` non-null directed-influence evidence from
physical source locus l to gene j. Define

```text
K_j = sum_l E_lj
C_j(m) = D_j × 1(K_j >= m).
```

The reported quantity is

```text
P(C_j(m)=1 | observed summary statistics, tested pairs, candidate set, priors).
```

Incoming locus indicators are treated as Bernoulli variables, and their
Poisson-binomial tail is evaluated exactly. The result is multiplied by the
disease ratio-evidence probability. The primary threshold is `m=2`; `m=3` and `m=5`
are prespecified sensitivity definitions.

Because there can be many source loci, even a sparse edge prior creates a
nonzero convergence baseline. The run records prior-only `P(K>=m)` and
prior-only `P(C_j(m))` so rankings can be interpreted as evidence updates from
that baseline.

This avoids the unidentified model in which a logistic regression is written
for an entirely unobserved core label without any observation likelihood.
Genes with network coverage below 50% retain an observed-data diagnostic but
do not receive a primary ranking probability. Missing pairs are not converted
to negative edges.

## Sensitivity analysis

The complete core calculation is repeated over:

- disease prior inclusion: 0.05, 0.10, 0.20;
- disease slab SD: 0.10, 0.20, 0.40;
- edge prior inclusion: 0.001, 0.01, 0.05;
- edge slab SD: 0.075, 0.15, 0.30;
- convergence thresholds: 2, 3, 5.

The displayed range is a prior-sensitivity interval, not a posterior credible
interval.

## SORT1 locus audit

GTEx v8 liver associations for SORT1, CELSR2, and PSRC1 are aligned to the CAD
GWAS within a 1 Mb window around rs12740374. The audit reports a strongest-
instrument Wald ratio and classic H0-H4 Wakefield-ABF colocalization.

Classic coloc assumes at most one causal variant per trait in the region and
does not use LD. PP(H4), PP(H3), `PP4/(PP3+PP4)`, effect-prior-variance
sensitivity, and shared-signal-prior (`p12`) sensitivity are reported together.
Colocalization supports a shared signal; it
does not prove that expression mediates disease risk.

## External proxy annotation audit

PanelApp familial-hypercholesterolaemia green genes are held out from all
fitting. The report gives total panel size, candidate overlap, eligible overlap,
annotation-retrieval average precision, and top-k Fisher enrichment. These
calculations treat unannotated candidates as negatives for metric computation
only. Panel genes are diagnostic monogenic evidence, not known “core genes,”
and genes outside the panel must not be interpreted as true negatives.

## Reporting framework

The project follows the intent of the STROBE-MR checklist by recording
hypotheses, sample sources, instrument construction, harmonization, ancestry,
sample overlap, multiplicity, sensitivity analyses, and limitations. The
frozen run manifest and `SHA256SUMS` make every compact input and result
auditable.
