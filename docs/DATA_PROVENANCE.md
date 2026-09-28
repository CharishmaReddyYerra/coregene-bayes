# Data provenance

All primary analysis inputs are real public summary statistics. The main run
does not substitute simulated effects.

| Role | Resource | Frozen identifier/release | Public source |
|---|---|---|---|
| Candidate loci | Open Targets Platform | release 26.09, LDL study GCST90019512 | `api.platform.opentargets.org/api/v4/graphql` |
| Upstream candidate model | BayesV2G sibling artifact | frozen `results/ldl_demo` checksums | local verified contract |
| Disease outcome | GWAS Catalog | CAD GCST90132314 | EBI harmonized summary-statistics FTP |
| Whole-blood cis instruments | eQTLGen Phase I | 2019-12-11 significant cis release | eQTLGen public downloads |
| Whole-blood target effects | eQTLGen Phase I | 2018-09-04 complete trans release | eQTLGen public downloads |
| Liver locus QTL | eQTL Catalogue | release 7 QTD000266, GTEx v8 liver n=208 | EBI QTL FTP |
| Proxy annotation audit | Genomics England PanelApp | panel 772, version 2.6 | PanelApp API |

Important distinctions:

- GCST90019512 supplies candidate LDL loci; GCST90132314 supplies the independent
  CAD outcome.
- The eQTL Catalogue release-7 QTD000266 file is GTEx v8 liver. A 2026
  prerelease reused QTD000266 for GTEx v10; URL, release, sample size, and
  format are therefore pinned, not merely the dataset ID.
- eQTLGen coordinates are GRCh37; GWAS Catalog and eQTL Catalogue coordinates
  used by the analysis are GRCh38. In the harmonized CAD file,
  `base_pair_location` is the lifted coordinate (`hm_coordinate_conversion=lo`)
  while `markername` preserves the submitted coordinate. Cross-source MR
  matching is by rsID plus allele orientation, not by raw position.
- The adjacent QTD000584 snapshot is a Sun et al. plasma pQTL resource and is
  deliberately not called an eQTL or used in the RNA network model.

The full trans table is 6.09 GB (5.67 GiB) compressed. The acquisition script
streams it and retains complete coefficients only for the prespecified
candidate targets and selected source instruments. This compact derivative is
included with checksums so the analysis runs offline after one acquisition.
Raw-source byte counts and SHA-256 digests are recorded in the source manifest.
Validated compaction-state records allow an interrupted fetch to resume without
rescanning the multi-gigabyte trans table.

Source licenses remain controlling. GWAS Catalog summary statistics are
distributed under CC0; eQTL Catalogue data are distributed under CC BY 4.0;
eQTLGen and PanelApp retain their source terms. The MIT license applies only to
this project's code.
