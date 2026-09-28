# Data layout

`snapshots/real_ldl_cad/` is the frozen, compact real-data input used by the
included analysis. It contains prespecified rows extracted from public Open
Targets, eQTLGen, GWAS Catalog, eQTL Catalogue/GTEx, and PanelApp resources.
No simulated effects are used in the main run.

The multi-gigabyte raw files are intentionally not copied into this project.
Their exact URLs, byte sizes, and SHA-256 digests are recorded in
`source_manifest.json`; the compact derivatives are checksum-listed there as
well. `eqtl_compaction_state.json` and `cad_compaction_state.json` let the fetch
workflow safely reuse finished streaming stages.

Rebuild from the public sources with:

```bash
.venv/bin/coregene-bayes fetch --cache .cache/coregene-bayes
```

Then run the offline analysis with:

```bash
.venv/bin/coregene-bayes run
```

Source licenses and dataset-specific caveats are documented in
`../docs/DATA_PROVENANCE.md`. Do not interpret absence from a compact table as
a negative association: the snapshot is conditional on the candidate and
tested-SNP universes described in the methods.
