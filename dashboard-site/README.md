# CoreGene Bayes dashboard

Static analytical dashboard for the real-data LDL → coronary artery disease run.

Project concept and study design initiated: **July 2026**. The analytical
snapshot keeps its actual September 2026 provenance timestamps.

The published site is generated from the parent project's versioned outputs. To refresh the embedded snapshot and downloadable files after rerunning the analysis:

```bash
node scripts/build-data.mjs
```

Serve `dist/` with any static HTTP server. The dashboard has no third-party JavaScript, font, analytics, or runtime API dependency.

## Source contract

- `../results/ldl_cad/core_gene_posteriors.csv`
- `../results/ldl_cad/gene_disease_mr.csv`
- `../results/ldl_cad/gene_gene_edges.csv.gz`
- `../results/ldl_cad/sort1_liver_cad_colocalization.csv`
- `../results/ldl_cad/analysis_attrition.csv`
- `../data/snapshots/real_ldl_cad/source_manifest.json`

The dashboard describes probabilities as conditional model outputs and directed gene-to-gene estimates as exploratory. It does not present them as clinical guidance or confirmed causal effects.
