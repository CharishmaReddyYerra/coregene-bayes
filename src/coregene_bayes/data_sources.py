"""Acquisition and compaction of the real public datasets used by the project."""

from __future__ import annotations

import csv
import gzip
import heapq
import json
import os
import shutil
import ssl
import tempfile
import urllib.error
import urllib.request
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TextIO

from coregene_bayes.bayesv2g_adapter import load_bayesv2g_candidates, top_seed_genes
from coregene_bayes.provenance import sha256_file, write_json


OPEN_TARGETS_QUERY = r"""
query CoreGeneCandidates($index: Int!) {
  meta { dataVersion { year month iteration } }
  credibleSets(studyIds: ["GCST90019512"], page: {index: $index, size: 100}) {
    count
    rows {
      studyLocusId studyId chromosome position pValueMantissa pValueExponent
      beta standardError sampleSize confidence
      variant { id rsIds chromosome position referenceAllele alternateAllele }
      l2GPredictions(page: {index: 0, size: 50}) {
        count
        rows {
          score
          target {
            id approvedSymbol approvedName biotype
            genomicLocation { chromosome start end strand }
          }
        }
      }
    }
  }
}
"""

EQTLGEN_CIS_URL = (
    "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/"
    "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
)
EQTLGEN_TRANS_URL = (
    "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/"
    "2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz"
)
CAD_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
    "GCST90132001-GCST90133000/GCST90132314/harmonised/GCST90132314.h.tsv.gz"
)
CAD_METADATA_URL = CAD_URL + "-meta.yaml"
GTEX_LIVER_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
    "QTS000015/QTD000266/QTD000266.all.tsv.gz"
)
PANELAPP_URL = "https://panelapp.genomicsengland.co.uk/api/v1/panels/772/?version=2.6"

EQTLGEN_COLUMNS = [
    "Pvalue",
    "SNP",
    "SNPChr",
    "SNPPos",
    "AssessedAllele",
    "OtherAllele",
    "Zscore",
    "Gene",
    "GeneSymbol",
    "GeneChr",
    "GenePos",
    "NrCohorts",
    "NrSamples",
    "FDR",
    "BonferroniP",
]
QTL_COLUMNS = [
    "molecular_trait_id",
    "chromosome",
    "position",
    "ref",
    "alt",
    "variant",
    "ma_samples",
    "maf",
    "pvalue",
    "beta",
    "se",
    "type",
    "ac",
    "an",
    "r2",
    "molecular_trait_object_id",
    "gene_id",
    "median_tpm",
    "rsid",
]
SORT1_GENE_IDS = {"ENSG00000134243", "ENSG00000143126", "ENSG00000134222"}
SORT1_START = 108_774_968
SORT1_END = 109_774_968
EQTL_COMPACT_FILES = (
    "eqtlgen_cis_instrument_candidates.tsv.gz",
    "eqtlgen_network_instruments.tsv.gz",
    "eqtlgen_trans_outcomes.tsv.gz",
)
CAD_COMPACT_FILES = ("cad_instrument_outcomes.tsv.gz", "cad_sort1_region.tsv.gz")


def _ssl_context(verify_tls: bool) -> ssl.SSLContext | None:
    return None if verify_tls else ssl._create_unverified_context()  # noqa: SLF001


def _request_bytes(url: str, *, verify_tls: bool, accept: str = "application/json") -> bytes:
    request = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "CoreGene-Bayes/0.1"})
    with urllib.request.urlopen(request, context=_ssl_context(verify_tls), timeout=180) as response:
        return response.read()


def _graphql(endpoint: str, query: str, variables: dict[str, Any], *, verify_tls: bool) -> dict[str, Any]:
    body = json.dumps({"query": query, "variables": variables}).encode("utf-8")
    request = urllib.request.Request(
        endpoint,
        data=body,
        headers={"Content-Type": "application/json", "User-Agent": "CoreGene-Bayes/0.1"},
        method="POST",
    )
    with urllib.request.urlopen(request, context=_ssl_context(verify_tls), timeout=180) as response:
        payload = json.load(response)
    if payload.get("errors"):
        raise RuntimeError(f"Open Targets GraphQL errors: {payload['errors']}")
    return payload


def fetch_ldl_candidates(
    endpoint: str,
    *,
    l2g_threshold: float,
    p_threshold: float,
    bayesv2g_directory: Path,
    verify_tls: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Fetch all significant LDL loci and create a prespecified candidate universe."""

    rows: list[dict[str, Any]] = []
    release: dict[str, Any] | None = None
    total: int | None = None
    page = 0
    while total is None or len(rows) < total:
        payload = _graphql(endpoint, OPEN_TARGETS_QUERY, {"index": page}, verify_tls=verify_tls)
        release = payload["data"]["meta"]["dataVersion"]
        credible_sets = payload["data"]["credibleSets"]
        total = int(credible_sets["count"])
        batch = credible_sets["rows"]
        if not batch:
            break
        rows.extend(batch)
        page += 1
    if total is None or len(rows) != total:
        raise ValueError(f"Open Targets response was truncated: expected {total}, received {len(rows)}")
    for locus in rows:
        page_data = locus["l2GPredictions"]
        if int(page_data["count"]) != len(page_data["rows"]):
            raise ValueError(f"L2G results were truncated for {locus['studyLocusId']}")

    bayesv2g, bayesv2g_manifest = load_bayesv2g_candidates(bayesv2g_directory)
    seeds = top_seed_genes(bayesv2g)
    seed_ids = set(seeds["geneId"].astype(str))
    selected_locus_rows: list[dict[str, Any]] = []
    for locus in rows:
        pvalue = float(locus["pValueMantissa"]) * 10.0 ** int(locus["pValueExponent"])
        if pvalue > p_threshold:
            continue
        lead = locus.get("variant") or {}
        for prediction in locus["l2GPredictions"]["rows"]:
            target = prediction.get("target") or {}
            gene_id = target.get("id")
            score = float(prediction["score"])
            if not gene_id or (score < l2g_threshold and gene_id not in seed_ids):
                continue
            location = target.get("genomicLocation") or {}
            selected_locus_rows.append(
                {
                    "gene_id": gene_id,
                    "gene_symbol": target.get("approvedSymbol"),
                    "gene_name": target.get("approvedName"),
                    "gene_biotype": target.get("biotype"),
                    "gene_chromosome": location.get("chromosome"),
                    "gene_start": location.get("start"),
                    "gene_end": location.get("end"),
                    "gene_strand": location.get("strand"),
                    "source_study_locus_id": locus["studyLocusId"],
                    "source_lead_variant": (lead.get("rsIds") or [lead.get("id")])[0],
                    "source_lead_variant_id": lead.get("id"),
                    "source_chromosome": locus["chromosome"],
                    "source_position": locus["position"],
                    "source_p_value": pvalue,
                    "l2g_score": score,
                    "bayesv2g_seed": gene_id in seed_ids,
                    "selection_reason": (
                        "l2g_threshold_and_bayesv2g_seed"
                        if score >= l2g_threshold and gene_id in seed_ids
                        else "l2g_threshold"
                        if score >= l2g_threshold
                        else "bayesv2g_rank1_seed"
                    ),
                }
            )
    if not selected_locus_rows:
        raise ValueError("Candidate selection returned no genes")

    best: dict[str, dict[str, Any]] = {}
    for row in selected_locus_rows:
        current = best.get(row["gene_id"])
        if current is None or (row["l2g_score"], -row["source_p_value"]) > (
            current["l2g_score"],
            -current["source_p_value"],
        ):
            best[row["gene_id"]] = row.copy()
    candidate_rows = sorted(best.values(), key=lambda item: (-item["l2g_score"], item["gene_id"]))
    metadata = {
        "open_targets_release": release,
        "total_credible_sets": total,
        "significant_locus_count": len(
            {
                row["source_study_locus_id"]
                for row in selected_locus_rows
            }
        ),
        "candidate_locus_gene_rows": len(selected_locus_rows),
        "unique_candidate_genes": len(candidate_rows),
        "bayesv2g_seed_gene_count": len(seed_ids),
        "bayesv2g_seed_genes_retained": sum(row["gene_id"] in seed_ids for row in candidate_rows),
        "bayesv2g_manifest_created_at": bayesv2g_manifest.get("created_at_utc"),
    }
    return candidate_rows, selected_locus_rows, metadata


def download(url: str, destination: Path, *, verify_tls: bool) -> Path:
    """Download a large source atomically, reusing a complete existing file."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    head = urllib.request.Request(url, headers={"User-Agent": "CoreGene-Bayes/0.1"}, method="HEAD")
    with urllib.request.urlopen(head, context=_ssl_context(verify_tls), timeout=180) as response:
        expected_size_text = response.headers.get("Content-Length")
    expected_size = int(expected_size_text) if expected_size_text else None
    if destination.exists() and destination.stat().st_size > 0:
        if expected_size is None or destination.stat().st_size == expected_size:
            return destination
        raise ValueError(
            f"Cached source is incomplete: {destination} has {destination.stat().st_size} bytes; "
            f"expected {expected_size}. Remove or replace that cache file before retrying."
        )
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "CoreGene-Bayes/0.1"})
    with urllib.request.urlopen(request, context=_ssl_context(verify_tls), timeout=600) as response:
        with temporary.open("wb") as handle:
            shutil.copyfileobj(response, handle, length=8 * 1024 * 1024)
    if expected_size is not None and temporary.stat().st_size != expected_size:
        raise IOError(
            f"Incomplete download for {url}: received {temporary.stat().st_size} of {expected_size} bytes"
        )
    temporary.replace(destination)
    return destination


def _write_dict_rows(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> int:
    count = 0
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def _is_palindromic(first: str, second: str) -> bool:
    return {first.upper(), second.upper()} in ({"A", "T"}, {"C", "G"})


def _source_fingerprint(path: Path, url: str) -> dict[str, Any]:
    return {"url": url, "bytes": path.stat().st_size, "sha256": sha256_file(path)}


def _reusable_compaction(
    state_path: Path,
    *,
    expected: dict[str, Any],
    output: Path,
    required_files: tuple[str, ...],
) -> dict[str, int] | None:
    if not state_path.exists() or any(
        not (output / name).exists() or (output / name).stat().st_size == 0
        for name in required_files
    ):
        return None
    try:
        state = json.loads(state_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if state.get("inputs") != expected:
        return None
    counts = state.get("counts")
    if not isinstance(counts, dict) or not all(isinstance(value, int) for value in counts.values()):
        return None
    return counts


def compact_eqtlgen(
    *,
    cis_path: Path,
    trans_path: Path,
    candidates: list[dict[str, Any]],
    output: Path,
    top_disease_instruments: int = 20,
) -> dict[str, int]:
    """Reduce the public 6 GB trans file to complete candidate-pair coefficients."""

    symbol_to_id = {
        str(row["gene_symbol"]): str(row["gene_id"])
        for row in candidates
        if row.get("gene_symbol")
    }
    candidate_symbols = set(symbol_to_id)
    tested_snps: set[str] = set()
    temporary_path = output / ".trans_candidate_rows.tsv.gz"
    target_rows = 0
    with gzip.open(trans_path, "rt", encoding="utf-8", newline="") as source, gzip.open(
        temporary_path, "wt", encoding="utf-8", newline=""
    ) as destination:
        header = source.readline().rstrip("\n\r").split("\t")
        if header != EQTLGEN_COLUMNS:
            raise ValueError("Unexpected eQTLGen trans schema")
        destination.write("\t".join(EQTLGEN_COLUMNS + ["target_gene_id"]) + "\n")
        for line in source:
            stripped = line.rstrip("\n\r")
            fields = stripped.split("\t")
            if len(fields) != len(EQTLGEN_COLUMNS):
                raise ValueError("Malformed row in eQTLGen trans file")
            tested_snps.add(fields[1])
            symbol = fields[8]
            if symbol in candidate_symbols:
                destination.write(f"{stripped}\t{symbol_to_id[symbol]}\n")
                target_rows += 1

    disease_heaps: dict[str, list[tuple[float, str, int, dict[str, Any]]]] = {}
    network_best: dict[str, dict[str, Any]] = {}
    cis_candidate_rows = 0
    with gzip.open(cis_path, "rt", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames != EQTLGEN_COLUMNS:
            raise ValueError("Unexpected eQTLGen cis schema")
        for row in reader:
            symbol = row["GeneSymbol"]
            if symbol not in candidate_symbols:
                continue
            if _is_palindromic(row["AssessedAllele"], row["OtherAllele"]):
                continue
            cis_candidate_rows += 1
            row = dict(row)
            row["source_gene_id"] = symbol_to_id[symbol]
            score = abs(float(row["Zscore"]))
            heap = disease_heaps.setdefault(symbol, [])
            entry = (score, row["SNP"], cis_candidate_rows, row)
            if len(heap) < top_disease_instruments:
                heapq.heappush(heap, entry)
            elif entry[:2] > heap[0][:2]:
                heapq.heapreplace(heap, entry)
            if row["SNP"] in tested_snps:
                current = network_best.get(symbol)
                if current is None or score > abs(float(current["Zscore"])):
                    network_best[symbol] = row

    disease_rows = [
        entry[3]
        for symbol in sorted(disease_heaps)
        for entry in sorted(disease_heaps[symbol], reverse=True)
    ]
    disease_fields = EQTLGEN_COLUMNS + ["source_gene_id"]
    disease_count = _write_dict_rows(
        output / "eqtlgen_cis_instrument_candidates.tsv.gz", disease_rows, disease_fields
    )
    network_rows = [network_best[symbol] for symbol in sorted(network_best)]
    network_count = _write_dict_rows(
        output / "eqtlgen_network_instruments.tsv.gz", network_rows, disease_fields
    )

    instrument_by_snp: dict[str, list[dict[str, Any]]] = {}
    for row in network_rows:
        instrument_by_snp.setdefault(row["SNP"], []).append(row)
    trans_output = output / "eqtlgen_trans_outcomes.tsv.gz"
    edge_rows = 0
    trans_fields = EQTLGEN_COLUMNS + ["target_gene_id"]
    with gzip.open(temporary_path, "rt", encoding="utf-8", newline="") as source, gzip.open(
        trans_output, "wt", encoding="utf-8", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(destination, delimiter="\t", fieldnames=trans_fields)
        writer.writeheader()
        for row in reader:
            if row["SNP"] in instrument_by_snp:
                writer.writerow(row)
                edge_rows += 1
    temporary_path.unlink()
    return {
        "tested_trans_snps": len(tested_snps),
        "candidate_target_rows_scanned": target_rows,
        "candidate_cis_rows_scanned": cis_candidate_rows,
        "disease_instrument_candidate_rows": disease_count,
        "network_source_instruments": network_count,
        "compact_trans_outcomes": edge_rows,
    }


def compact_cad(
    *,
    cad_path: Path,
    cis_instruments_path: Path,
    output: Path,
) -> dict[str, int]:
    with gzip.open(cis_instruments_path, "rt", encoding="utf-8", newline="") as handle:
        instrument_snps = {row["SNP"] for row in csv.DictReader(handle, delimiter="\t")}
    instrument_rows: list[dict[str, Any]] = []
    regional_rows: list[dict[str, Any]] = []
    with gzip.open(cad_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("CAD GWAS has no header")
        fields = list(reader.fieldnames)
        for row in reader:
            rsid = row.get("rsid") or row.get("variant_id")
            if rsid in instrument_snps:
                instrument_rows.append(row)
            chromosome = str(row.get("chromosome", "")).removeprefix("chr")
            try:
                position = int(row.get("base_pair_location", ""))
            except ValueError:
                continue
            if chromosome == "1" and SORT1_START <= position <= SORT1_END:
                regional_rows.append(row)
    instrument_count = _write_dict_rows(output / "cad_instrument_outcomes.tsv.gz", instrument_rows, fields)
    regional_count = _write_dict_rows(output / "cad_sort1_region.tsv.gz", regional_rows, fields)
    return {"cad_instrument_outcomes": instrument_count, "cad_sort1_region_rows": regional_count}


def fetch_liver_sort1(destination: Path, *, verify_tls: bool, existing: Path | None = None) -> int:
    if existing and existing.exists():
        shutil.copyfile(existing, destination)
        with gzip.open(destination, "rt", encoding="utf-8") as handle:
            return sum(1 for _ in handle) - 1
    try:
        import pysam
    except ImportError as error:  # pragma: no cover - exercised only by download workflow
        raise RuntimeError("pysam is required to fetch the remote liver QTL slice") from error
    previous = Path.cwd()
    rows: list[str] = []
    with tempfile.TemporaryDirectory(prefix="coregene-tabix-") as temporary:
        try:
            os.chdir(temporary)
            table = pysam.TabixFile(GTEX_LIVER_URL)
            for line in table.fetch("1", SORT1_START - 1, SORT1_END):
                fields = line.split("\t")
                if len(fields) == len(QTL_COLUMNS) and fields[16] in SORT1_GENE_IDS:
                    rows.append(line)
            table.close()
        finally:
            os.chdir(previous)
    with gzip.open(destination, "wt", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(QTL_COLUMNS) + "\n")
        handle.writelines(f"{row}\n" for row in rows)
    return len(rows)


def write_candidates(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("Cannot write an empty candidate table")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_real_snapshot(
    *,
    output: Path,
    cache: Path,
    endpoint: str,
    l2g_threshold: float,
    p_threshold: float,
    bayesv2g_directory: Path,
    existing_liver_snapshot: Path | None,
    verify_tls: bool,
) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    candidates, locus_rows, candidate_metadata = fetch_ldl_candidates(
        endpoint,
        l2g_threshold=l2g_threshold,
        p_threshold=p_threshold,
        bayesv2g_directory=bayesv2g_directory,
        verify_tls=verify_tls,
    )
    write_candidates(output / "candidates.csv", candidates)
    write_candidates(output / "candidate_locus_rows.csv", locus_rows)

    cis = download(EQTLGEN_CIS_URL, cache / "eqtlgen_cis_significant.tsv.gz", verify_tls=verify_tls)
    trans = download(EQTLGEN_TRANS_URL, cache / "eqtlgen_trans_full.tsv.gz", verify_tls=verify_tls)
    cad = download(CAD_URL, cache / "GCST90132314.h.tsv.gz", verify_tls=verify_tls)
    raw_sources = {
        "eqtlgen_cis": _source_fingerprint(cis, EQTLGEN_CIS_URL),
        "eqtlgen_trans": _source_fingerprint(trans, EQTLGEN_TRANS_URL),
        "cad_gwas": _source_fingerprint(cad, CAD_URL),
    }
    candidate_sha256 = sha256_file(output / "candidates.csv")
    eqtl_state_inputs = {
        "candidate_sha256": candidate_sha256,
        "eqtlgen_cis": raw_sources["eqtlgen_cis"],
        "eqtlgen_trans": raw_sources["eqtlgen_trans"],
    }
    eqtl_state_path = output / "eqtl_compaction_state.json"
    eqtl_counts = _reusable_compaction(
        eqtl_state_path,
        expected=eqtl_state_inputs,
        output=output,
        required_files=EQTL_COMPACT_FILES,
    )
    if eqtl_counts is None:
        eqtl_counts = compact_eqtlgen(
            cis_path=cis, trans_path=trans, candidates=candidates, output=output
        )
        write_json({"inputs": eqtl_state_inputs, "counts": eqtl_counts}, eqtl_state_path)

    cad_state_inputs = {
        "cis_instruments_sha256": sha256_file(
            output / "eqtlgen_cis_instrument_candidates.tsv.gz"
        ),
        "cad_gwas": raw_sources["cad_gwas"],
    }
    cad_state_path = output / "cad_compaction_state.json"
    cad_counts = _reusable_compaction(
        cad_state_path,
        expected=cad_state_inputs,
        output=output,
        required_files=CAD_COMPACT_FILES,
    )
    if cad_counts is None:
        cad_counts = compact_cad(
            cad_path=cad,
            cis_instruments_path=output / "eqtlgen_cis_instrument_candidates.tsv.gz",
            output=output,
        )
        write_json({"inputs": cad_state_inputs, "counts": cad_counts}, cad_state_path)
    liver_rows = fetch_liver_sort1(
        output / "gtex_v8_liver_sort1_genes.tsv.gz",
        verify_tls=verify_tls,
        existing=existing_liver_snapshot,
    )
    (output / "cad_metadata.yaml").write_bytes(
        _request_bytes(CAD_METADATA_URL, verify_tls=verify_tls, accept="text/yaml")
    )
    panelapp = json.loads(_request_bytes(PANELAPP_URL, verify_tls=verify_tls))
    write_json(panelapp, output / "panelapp_fh_v2_6.json")

    generated = sorted(path for path in output.iterdir() if path.is_file() and path.name != "source_manifest.json")
    manifest_payload = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "tls_verified": verify_tls,
        "candidate_selection": {
            "source_trait": "LDL cholesterol",
            "source_study": "GCST90019512",
            "gwas_p_threshold": p_threshold,
            "l2g_score_threshold": l2g_threshold,
            "bayesv2g_rank1_seeds_forced_in": True,
            **candidate_metadata,
        },
        "analysis_sources": {
            "disease": "coronary artery disease",
            "disease_gwas": "GCST90132314",
            "network_qtl": "eQTLGen phase I whole blood, n up to 31,684",
            "locus_qtl": "eQTL Catalogue r7 QTD000266 / GTEx v8 liver, n=208",
            "genome_build_note": (
                "GWAS and liver QTL coordinates are GRCh38; eQTLGen coordinates are GRCh37. "
                "Cross-source MR matching uses rsID plus signed alleles, never raw coordinate equality."
            ),
        },
        "counts": {**eqtl_counts, **cad_counts, "liver_sort1_rows": liver_rows},
        "raw_sources": raw_sources,
        "urls": {
            "open_targets": endpoint,
            "eqtlgen_cis": EQTLGEN_CIS_URL,
            "eqtlgen_trans": EQTLGEN_TRANS_URL,
            "cad_gwas": CAD_URL,
            "cad_metadata": CAD_METADATA_URL,
            "gtex_liver": GTEX_LIVER_URL,
            "panelapp": PANELAPP_URL,
        },
        "files": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in generated
        },
    }
    write_json(manifest_payload, output / "source_manifest.json")
    return manifest_payload
