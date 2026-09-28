# Data dictionary

## `candidate_genes.csv`

- `gene_id`, `gene_symbol`, `gene_name`, `gene_biotype`: Open Targets gene
  identity; Ensembl ID is the primary key.
- `gene_chromosome`, `gene_start`, `gene_end`, `gene_strand`: GRCh38 target
  location reported by Open Targets.
- `source_study_locus_id`, `source_lead_variant*`, `source_*`: best-scoring
  LDL locus used for node metadata.
- `l2g_score`: candidate-ascertainment score; not a core-model feature.
- `bayesv2g_seed`: whether the gene is rank one at a frozen BayesV2G demo locus.
- `selection_reason`: threshold and/or seed rule that retained the gene.

## `gene_disease_mr.csv`

- `instrument_snp`: selected eQTLGen cis instrument.
- `beta_x_approx`, `se_x_approx`: standardized Z/N approximation.
- `beta_y_aligned_per_allele`, `se_y_per_allele`: allele-aligned CAD
  log-odds association on the source per-allele scale.
- `outcome_effect_allele_frequency`, `genotype_sd`: real CAD allele frequency
  and `sqrt(2p(1-p))` scale factor.
- `beta_y_per_genotype_sd`, `se_y_per_genotype_sd`: CAD association placed on
  the same approximate genotype-SD scale as the eQTL effect and used in MR.
- `mr_beta`, `mr_se`, `mr_ci_*`, `mr_p_value`, `mr_fdr`: Wald-ratio result.
- `instrument_f`: squared exposure Z score.
- `disease_probability`: exploratory spike-and-slab probability of a non-null
  Wald-ratio signal under configured priors; missing if QC fails. It does not
  include the probability that MR assumptions or colocalization hold.
- `qc_status`: estimability and instrument-strength result.

## `gene_gene_edges.csv.gz`

- `source_gene_*`, `target_gene_*`: directed candidate pair.
- `source_locus_id`: 1 Mb physical cluster used for convergence collapse.
- `beta_*_approx`, `se_*_approx`: standardized eQTLGen effects.
- `mr_*`: single-instrument Wald-ratio result.
- `edge_probability`: model-based probability of a non-null influence under
  the edge spike-and-slab prior.
- `edge_qc`: `supported_pleiotropy_ambiguous`, `model_supported_not_fdr`,
  `indeterminate`, weak-instrument, or harmonization failure.
- `interpretation`: explicit non-mechanistic edge scope.

## `locus_collapsed_edges.csv`

- `target_gene_id`, `source_locus_id`: convergence unit.
- `locus_probability`: maximum edge probability within the source locus, used
  to avoid multiplying correlated within-locus evidence.
- `n_edges_collapsed`: source-gene edges represented by that unit.
- `collapse_rule`: records the conservative maximum rule explicitly.

## `network_coverage.csv`

- `tested_source_loci`: physical source loci with an allele-aligned tested
  coefficient for the target.
- `source_locus_universe`: loci with usable source instruments.
- `network_coverage`: tested/universe ratio. Missing pairs are not negatives.

## `core_gene_posteriors.csv`

- `disease_probability`: posterior non-null expression-to-CAD evidence.
- `expected_incoming_loci`: sum of collapsed locus probabilities.
- `p_core_m2`, `p_core_m3`, `p_core_m5`: main rankable probabilities after the
  network-coverage rule.
- `p_core_m*_observed`: diagnostic probabilities conditional on whatever
  source loci were observed, even below coverage threshold.
- `p_core_m2_sensitivity_min/max`: range across prespecified prior scenarios;
  not a credible interval.
- `core_rank`: descending primary `m=2` rank among eligible genes.
- `validation_positive`: overlap with the held-out PanelApp green annotation;
  this is not a verified core-gene label.

## `sort1_liver_cad_colocalization.csv`

- `pp_h0`–`pp_h4`: classic coloc hypothesis posterior probabilities.
- `pp_h4_given_h3_h4`: conditional support for a shared rather than distinct
  signal when both traits associate.
- `pp_h4_prior_variance_*`: prior sensitivity.
- `pp_h4_shared_prior_*`: sensitivity to the coloc shared-signal prior `p12`.
- `top_shared_variant*`: top conditional shared variant.
- `mr_*`: strongest liver cis-eQTL Wald-ratio audit.

## Audit and provenance files

- `harmonization_audit.csv.gz`: every attempted gene-to-CAD allele match and
  retain/drop status.
- `prior_sensitivity.csv.gz`: gene × prior × threshold results.
- `validation_panelapp_genes.csv`, `validation_metrics.csv`: external proxy
  annotation, overlap counts, retrieval average precision, and top-k
  enrichment metrics.
- `model.json`: exact latent-event definition and priors.
- `source_manifest.json`, `run_manifest.json`, `SHA256SUMS`: sources, versions,
  environment, counts, and file integrity.
