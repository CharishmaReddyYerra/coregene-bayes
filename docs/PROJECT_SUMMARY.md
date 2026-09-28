# Project summary

CoreGene-Bayes connects real molecular-QTL instruments, a real disease GWAS,
and a label-free probability model to answer a narrower and more defensible
version of the proposed project question.

The project concept and study design were initiated in **July 2026**. The
included real-data snapshot was acquired and analyzed in September 2026, so
its source and run manifests preserve those later provenance timestamps.

The pipeline is:

```text
real LDL loci (Open Targets + BayesV2G seeds)
  → prespecified candidate genes
  → eQTLGen cis instruments
  → separate CAD GWAS, allele-frequency-scaled Wald-ratio evidence
  → complete eQTLGen trans target coefficients
  → directed influence graph
  → physical source-locus convergence
  → P(D_j=1 and K_j≥m | observed evidence)
  → held-out PanelApp and liver SORT1 audits
```

The central improvement over a degree-based score is the explicit event being
estimated: a candidate must have disease-linked expression evidence and
probabilistic inputs from at least `m` distinct source loci. Co-located source
genes cannot inflate convergence simply by being numerous.

The central scientific restraint is equally important. Genome-wide ratio
evidence lacks regional colocalization, and a gene-to-gene Wald
ratio from same-cohort blood eQTLs remains pleiotropy-ambiguous and is called a
directed influence, not a mechanistic regulatory edge. The HTML report and all
tables carry that interpretation boundary.

## Frozen real-data run

The final snapshot contains 239 candidate genes, 173 harmonized gene-to-CAD
ratio estimates, 20,493 tested directed pairs from 114 source genes to 182
targets, and 97 physical source loci. All 173 disease-estimable genes meet the
primary network-coverage rule.

HP ranks first at the primary `m=2` definition (`P=0.769`), followed by LAMC1,
LDLR, and SORT1. LDLR ranks first at the stricter `m=3` and `m=5` definitions,
reflecting the strongest convergence signal in the candidate set. The broad
prior/slab sensitivity ranges prevent these point estimates from being treated
as calibrated biological probabilities.

The separate liver audit identifies rs12740374 as the top shared variant and
gives classic-coloc `PP(H4)` near 1 for SORT1, CELSR2, and PSRC1. Because all
three track the same regional signal, this supports shared genetic architecture
but does not uniquely prove which expression trait mediates CAD risk.
