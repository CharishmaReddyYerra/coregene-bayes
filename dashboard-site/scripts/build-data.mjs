import { readFile, writeFile, mkdir, copyFile } from "node:fs/promises";
import { gunzipSync } from "node:zlib";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const siteRoot = path.resolve(here, "..");
const projectRoot = path.resolve(siteRoot, "..");
const results = path.join(projectRoot, "results", "ldl_cad");
const snapshot = path.join(projectRoot, "data", "snapshots", "real_ldl_cad");
const dist = path.join(siteRoot, "dist");

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let i = 0; i < text.length; i += 1) {
    const char = text[i];
    if (quoted) {
      if (char === '"' && text[i + 1] === '"') { field += '"'; i += 1; }
      else if (char === '"') quoted = false;
      else field += char;
    } else if (char === '"') quoted = true;
    else if (char === ",") { row.push(field); field = ""; }
    else if (char === "\n") { row.push(field.replace(/\r$/, "")); rows.push(row); row = []; field = ""; }
    else field += char;
  }
  if (field || row.length) { row.push(field); rows.push(row); }
  const headers = rows.shift();
  return rows.filter(item => item.some(Boolean)).map(values => Object.fromEntries(headers.map((header, index) => [header, values[index] ?? ""])));
}

const number = value => value === "" || value == null ? null : Number(value);
const bool = value => String(value).toLowerCase() === "true";
const csv = async file => parseCsv(await readFile(file, "utf8"));
const gzCsv = async file => parseCsv(gunzipSync(await readFile(file)).toString("utf8"));

const [summary, model, sourceManifest, genesRaw, candidatesRaw, diseaseRaw, edgesRaw, colocRaw, attritionRaw, validationRaw] = await Promise.all([
  readFile(path.join(results, "run_summary.json"), "utf8").then(JSON.parse),
  readFile(path.join(results, "model.json"), "utf8").then(JSON.parse),
  readFile(path.join(snapshot, "source_manifest.json"), "utf8").then(JSON.parse),
  csv(path.join(results, "core_gene_posteriors.csv")),
  csv(path.join(results, "candidate_genes.csv")),
  csv(path.join(results, "gene_disease_mr.csv")),
  gzCsv(path.join(results, "gene_gene_edges.csv.gz")),
  csv(path.join(results, "sort1_liver_cad_colocalization.csv")),
  csv(path.join(results, "analysis_attrition.csv")),
  csv(path.join(results, "validation_metrics.csv"))
]);

const candidates = new Map(candidatesRaw.map(row => [row.gene_id, row]));
const disease = new Map(diseaseRaw.filter(row => row.analysis_retained === "True").map(row => [row.gene_id, row]));

const genes = genesRaw.map(row => {
  const candidate = candidates.get(row.gene_id) ?? {};
  const mr = disease.get(row.gene_id) ?? {};
  return {
    id: row.gene_id,
    gene: row.gene_symbol,
    name: candidate.gene_name || row.gene_symbol,
    chromosome: candidate.gene_chromosome || null,
    rank: number(row.core_rank),
    m2: number(row.p_core_m2),
    m3: number(row.p_core_m3),
    m5: number(row.p_core_m5),
    sensitivityMin: number(row.p_core_m2_sensitivity_min),
    sensitivityMax: number(row.p_core_m2_sensitivity_max),
    diseaseProbability: number(row.disease_probability),
    expectedIncoming: number(row.expected_incoming_loci),
    sourceLociTested: number(row.tested_source_loci),
    sourceLocusUniverse: number(row.source_locus_universe),
    networkCoverage: number(row.network_coverage),
    l2g: number(row.l2g_score_for_ascertainment),
    bayesv2gSeed: bool(row.bayesv2g_seed),
    panelAppPositive: bool(row.validation_positive),
    instrument: mr.instrument_snp || null,
    mrBeta: number(mr.mr_beta),
    mrSe: number(mr.mr_se),
    mrP: number(mr.mr_p_value),
    mrFdr: number(mr.mr_fdr),
    instrumentF: number(mr.instrument_f)
  };
});

const edges = edgesRaw
  .filter(row => number(row.edge_probability) >= 0.5)
  .map(row => ({
    source: row.source_gene_symbol,
    target: row.target_gene_symbol,
    sourceId: row.source_gene_id,
    targetId: row.target_gene_id,
    locus: row.source_locus_id,
    instrument: row.instrument_snp,
    probability: number(row.edge_probability),
    beta: number(row.mr_beta),
    se: number(row.mr_se),
    fdr: number(row.mr_fdr),
    qc: row.edge_qc
  }))
  .sort((a, b) => b.probability - a.probability);

const coloc = colocRaw.map(row => ({
  gene: row.gene_symbol,
  sharedVariants: number(row.shared_variant_count),
  ppH4: number(row.pp_h4),
  ppH3: number(row.pp_h3),
  priorLow: number(row.pp_h4_shared_prior_1e_6),
  priorHigh: number(row.pp_h4_shared_prior_1e_4),
  topVariant: row.top_shared_variant,
  variantProbability: number(row.top_shared_variant_probability_given_h4),
  beta: number(row.mr_beta),
  se: number(row.mr_se),
  p: number(row.mr_p_value),
  instrumentF: number(row.instrument_f)
}));

const payload = {
  generatedAt: new Date().toISOString(),
  summary,
  model,
  genes,
  edges,
  coloc,
  attrition: attritionRaw.map(row => ({ stage: row.stage, count: number(row.count), note: row.note })),
  validation: validationRaw.map(row => ({ metric: row.metric, topK: number(row.top_k), value: number(row.value), p: number(row.p_value) })),
  sources: {
    analysis: sourceManifest.analysis_sources,
    selection: sourceManifest.candidate_selection,
    counts: sourceManifest.counts,
    raw: sourceManifest.raw_sources,
    urls: sourceManifest.urls,
    createdAt: sourceManifest.created_at_utc,
    tlsVerified: sourceManifest.tls_verified
  }
};

await mkdir(path.join(dist, "downloads"), { recursive: true });
await writeFile(path.join(dist, "data.json"), `${JSON.stringify(payload)}\n`);

const downloads = [
  [path.join(results, "core_gene_posteriors.csv"), "ranked-genes.csv"],
  [path.join(results, "gene_disease_mr.csv"), "gene-disease-estimates.csv"],
  [path.join(results, "gene_gene_edges.csv.gz"), "directed-edges.csv.gz"],
  [path.join(results, "sort1_liver_cad_colocalization.csv"), "sort1-locus-colocalization.csv"],
  [path.join(results, "run_manifest.json"), "run-manifest.json"],
  [path.join(snapshot, "source_manifest.json"), "source-manifest.json"],
  [path.join(projectRoot, "docs", "METHODS.md"), "methods.md"],
  [path.join(projectRoot, "docs", "LIMITATIONS.md"), "limitations.md"]
];
await Promise.all(downloads.map(([source, name]) => copyFile(source, path.join(dist, "downloads", name))));

console.log(JSON.stringify({ genes: genes.length, edges: edges.length, coloc: coloc.length, output: path.join(dist, "data.json") }));
