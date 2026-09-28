"""Portable figures and HTML report."""

from __future__ import annotations

import html
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from coregene_bayes.network import top_edges_for_plot


COLOR = "#28536B"
ACCENT = "#D1495B"
MUTED = "#6B7280"


def plot_core_rankings(core: pd.DataFrame, path: Path, *, limit: int = 25) -> None:
    frame = core.dropna(subset=["p_core_m2"]).nlargest(limit, "p_core_m2").sort_values("p_core_m2")
    fig, axis = plt.subplots(figsize=(9, max(5, 0.28 * len(frame))))
    lower = frame["p_core_m2"] - frame["p_core_m2_sensitivity_min"]
    upper = frame["p_core_m2_sensitivity_max"] - frame["p_core_m2"]
    axis.errorbar(
        frame["p_core_m2"],
        frame["gene_symbol"],
        xerr=np.vstack([lower.clip(lower=0), upper.clip(lower=0)]),
        fmt="o",
        color=COLOR,
        ecolor="#93A8B3",
        capsize=2,
    )
    axis.set_xlabel("P(non-null CAD ratio evidence and ≥2 incoming loci | model/data)")
    axis.set_ylabel("")
    axis.set_xlim(left=0)
    axis.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_convergence(core: pd.DataFrame, path: Path) -> None:
    frame = core.dropna(subset=["disease_probability", "expected_incoming_loci"])
    fig, axis = plt.subplots(figsize=(8, 6))
    sizes = 30 + 300 * frame["p_core_m2"].fillna(0)
    axis.scatter(
        frame["expected_incoming_loci"],
        frame["disease_probability"],
        s=sizes,
        c=frame["p_core_m2"],
        cmap="viridis",
        alpha=0.78,
        edgecolor="white",
        linewidth=0.5,
    )
    offsets = ((6, 9), (6, -14), (6, 9), (6, -14), (6, 9), (6, -14))
    for row, offset in zip(frame.nlargest(6, "p_core_m2").itertuples(), offsets, strict=True):
        axis.annotate(
            row.gene_symbol,
            (row.expected_incoming_loci, row.disease_probability),
            xytext=offset,
            textcoords="offset points",
            fontsize=8,
        )
    axis.set_xlabel("Expected number of incoming physical source loci")
    axis.set_ylabel("Model probability of a non-null gene→CAD ratio signal")
    axis.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_sort1_coloc(summary: pd.DataFrame, path: Path) -> None:
    frame = summary.sort_values("pp_h4")
    fig, axis = plt.subplots(figsize=(7, 3.5))
    axis.barh(frame["gene_symbol"], frame["pp_h4"], color=COLOR)
    for index, value in enumerate(frame["pp_h4"]):
        axis.text(0.985, index, f"{value:.6f}", ha="right", va="center", color="white", fontsize=8)
    axis.set_xlim(0, 1)
    axis.set_xlabel("Colocalization PP(H4), CAD and liver expression")
    axis.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_directed_network(edges: pd.DataFrame, path: Path, *, limit: int = 30) -> None:
    """Plot the strongest exploratory directed edges as a source-target matrix.

    Rows encode instrumented source genes and columns encode target genes, so a
    populated cell always reads ``row → column``.  Colour represents the
    spike-and-slab non-null probability; the annotation also reports the sign
    of the Wald-ratio estimate.  Missing cells are unselected/untested in this
    top-edge view and must not be read as evidence for a null edge.
    """

    if limit < 1:
        raise ValueError("limit must be positive")
    required = {
        "source_gene_symbol",
        "target_gene_symbol",
        "edge_probability",
        "mr_beta",
        "mr_p_value",
    }
    if "harmonization_retained" not in edges and "harmonization_status" not in edges:
        required.add("harmonization_retained")
    missing = required.difference(edges.columns)
    if missing:
        raise ValueError(f"Edge table lacks plotting fields: {sorted(missing)}")

    # Ask for all eligible rows before de-duplicating source-target pairs.  This
    # avoids repeated instruments for one pair crowding distinct edges out of
    # the displayed top-N set.
    selected = top_edges_for_plot(edges, limit=max(len(edges), limit)).copy()
    selected["source_label"] = selected["source_gene_symbol"].astype("string")
    selected["target_label"] = selected["target_gene_symbol"].astype("string")
    selected = selected.loc[
        selected["source_label"].notna() & selected["target_label"].notna()
    ]
    selected = (
        selected.sort_values(
            ["edge_probability", "mr_p_value"], ascending=[False, True]
        )
        .drop_duplicates(["source_label", "target_label"])
        .head(limit)
    )

    if selected.empty:
        fig, axis = plt.subplots(figsize=(8, 3.5))
        axis.text(
            0.5,
            0.5,
            "No harmonized gene→gene edges were available for plotting",
            ha="center",
            va="center",
            color=MUTED,
        )
        axis.axis("off")
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        return

    source_order = (
        selected.groupby("source_label", observed=True)["edge_probability"]
        .max()
        .sort_values(ascending=False)
        .index.tolist()
    )
    target_order = (
        selected.groupby("target_label", observed=True)["edge_probability"]
        .max()
        .sort_values(ascending=False)
        .index.tolist()
    )
    probability = selected.pivot_table(
        index="source_label",
        columns="target_label",
        values="edge_probability",
        aggfunc="max",
    ).reindex(index=source_order, columns=target_order)
    beta = selected.pivot_table(
        index="source_label",
        columns="target_label",
        values="mr_beta",
        aggfunc="first",
    ).reindex(index=source_order, columns=target_order)

    width = max(7.0, 3.0 + 0.48 * len(target_order))
    height = max(4.5, 2.2 + 0.38 * len(source_order))
    fig, axis = plt.subplots(figsize=(width, height))
    colour_map = plt.colormaps["viridis"].with_extremes(bad="#F1F5F7")
    image = axis.imshow(
        np.ma.masked_invalid(probability.to_numpy(float)),
        cmap=colour_map,
        vmin=0,
        vmax=1,
        aspect="auto",
    )
    axis.set_xticks(np.arange(len(target_order)), labels=target_order)
    axis.set_yticks(np.arange(len(source_order)), labels=source_order)
    axis.tick_params(axis="x", labelrotation=55, labelsize=8)
    axis.tick_params(axis="y", labelsize=8)
    for label in axis.get_xticklabels():
        label.set_ha("right")
        label.set_rotation_mode("anchor")
    axis.set_xlabel("Target gene")
    axis.set_ylabel("Instrumented source gene  →")
    axis.set_title(
        f"Strongest {len(selected)} exploratory gene→gene edges\n"
        "cell = P(non-null); +/− = Wald-ratio sign",
        pad=12,
    )

    axis.set_xticks(np.arange(-0.5, len(target_order), 1), minor=True)
    axis.set_yticks(np.arange(-0.5, len(source_order), 1), minor=True)
    axis.grid(which="minor", color="white", linewidth=1.2)
    axis.tick_params(which="minor", bottom=False, left=False)
    for row_index in range(len(source_order)):
        for column_index in range(len(target_order)):
            value = probability.iat[row_index, column_index]
            if not np.isfinite(value):
                continue
            estimate = beta.iat[row_index, column_index]
            sign = "+" if np.isfinite(estimate) and estimate >= 0 else "−"
            text_colour = "white" if value >= 0.48 else "#17202A"
            axis.text(
                column_index,
                row_index,
                f"{value:.2f}\n{sign}",
                ha="center",
                va="center",
                fontsize=7,
                color=text_colour,
            )
    colour_bar = fig.colorbar(image, ax=axis, pad=0.02, shrink=0.85)
    colour_bar.set_label("Model P(non-null edge)")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _table(frame: pd.DataFrame, columns: list[str], *, rows: int = 20) -> str:
    present = [column for column in columns if column in frame.columns]
    return frame[present].head(rows).to_html(index=False, border=0, float_format=lambda x: f"{x:.4g}")


def write_html_report(
    path: Path,
    *,
    core: pd.DataFrame,
    edges: pd.DataFrame,
    disease_mr: pd.DataFrame,
    sort1: pd.DataFrame,
    attrition: pd.DataFrame,
    validation_metrics: pd.DataFrame,
    metadata: dict,
) -> None:
    eligible = core["p_core_m2"].notna().sum()
    disease_estimable = int(disease_mr["mr_beta"].notna().sum())
    network_sources = int(edges["source_gene_id"].nunique())
    model_supported_edges = edges["edge_qc"].isin(
        ["supported_pleiotropy_ambiguous", "model_supported_not_fdr"]
    ).sum()
    fdr_supported_edges = edges["edge_qc"].eq("supported_pleiotropy_ambiguous").sum()
    title = html.escape(metadata.get("project", "CoreGene-Bayes"))
    rankable = core.dropna(subset=["p_core_m2"])
    if rankable.empty:
        observed_result = "No gene met the joint disease/network coverage criteria."
    else:
        top_m2 = rankable.nlargest(1, "p_core_m2").iloc[0]
        top_m3 = rankable.nlargest(1, "p_core_m3").iloc[0]
        top_m5 = rankable.nlargest(1, "p_core_m5").iloc[0]
        observed_result = (
            f"{html.escape(str(top_m2['gene_symbol']))} ranks first at m=2 "
            f"(P={float(top_m2['p_core_m2']):.3f}; prior/slab range "
            f"{float(top_m2['p_core_m2_sensitivity_min']):.3f}–"
            f"{float(top_m2['p_core_m2_sensitivity_max']):.3f}). "
            f"{html.escape(str(top_m3['gene_symbol']))} ranks first at m=3 "
            f"(P={float(top_m3['p_core_m3']):.3f}) and "
            f"{html.escape(str(top_m5['gene_symbol']))} ranks first at m=5 "
            f"(P={float(top_m5['p_core_m5']):.3f})."
        )
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title>
<style>
body{{font-family:Inter,ui-sans-serif,system-ui,sans-serif;max-width:1120px;margin:0 auto;padding:32px;color:#17202a;line-height:1.55}}
h1,h2{{color:#17324d}} .hero{{background:#eef5f8;padding:24px;border-left:5px solid #28536B;border-radius:6px}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px;margin:18px 0}}
.card{{border:1px solid #d8e0e5;border-radius:7px;padding:15px}} .number{{font-size:1.7rem;font-weight:700;color:#28536B}}
.warning{{background:#fff7e6;border-left:5px solid #d99000;padding:16px}} table{{border-collapse:collapse;width:100%;font-size:.86rem;display:block;overflow-x:auto}}
th,td{{padding:7px 9px;border-bottom:1px solid #e5e7eb;text-align:left;white-space:nowrap}} th{{background:#f4f6f8}}
img{{max-width:100%;height:auto}} code{{background:#f2f4f5;padding:2px 4px;border-radius:3px}}
</style></head><body>
<section class="hero"><h1>{title}</h1>
<p>A real-data research prototype asking which LDL-locus candidate genes have both genetically anchored expression-to-CAD ratio evidence and convergent incoming whole-blood expression associations.</p></section>
<div class="cards"><div class="card"><div class="number">{len(core)}</div>candidate genes</div>
<div class="card"><div class="number">{disease_estimable}</div>genes with CAD ratio estimates</div>
<div class="card"><div class="number">{network_sources}</div>network source genes</div>
<div class="card"><div class="number">{eligible}</div>genes meeting coverage/QC</div>
<div class="card"><div class="number">{len(edges):,}</div>tested directed pairs</div>
<div class="card"><div class="number">{model_supported_edges:,}</div>model-supported edges</div>
<div class="card"><div class="number">{fdr_supported_edges:,}</div>also passing FDR</div></div>
<div class="warning"><strong>Interpretation boundary.</strong> These are exploratory model-based core-evidence probabilities, not probabilities of biological truth. Genome-wide regional colocalization was unavailable, and the gene→gene layer uses single-instrument, same-cohort whole-blood summary data. Horizontal pleiotropy, linkage, sample-overlap bias, and source-gene ambiguity remain unresolved. “Directed influence” does not mean a direct molecular regulatory edge.</div>
<h2>Analysis population and attrition</h2><p>Non-estimability in whole blood is a data-coverage result, not negative biological evidence.</p>
{_table(attrition, list(attrition.columns), rows=20)}
<h2>Primary result</h2><p>The prespecified definition is <code>C_j(m)=D_j × 1(K_j≥m)</code>, where D is non-null gene→CAD Wald-ratio evidence and K counts probabilistic incoming physical source loci. The primary threshold is m=2; m=3 and m=5 are reported as sensitivity definitions. With {int(metadata.get('source_locus_count', 0))} source loci and the configured priors alone, before observing effect estimates, P(K≥2) is {float(metadata.get('prior_only_convergence_probability', {}).get(2, float('nan'))):.3f} and P(C(m=2)) is {float(metadata.get('prior_only_core_probability', {}).get(2, float('nan'))):.3f}; rankings should be interpreted as evidence updates from that baseline.</p>
<p><strong>Observed:</strong> {observed_result} Wide sensitivity ranges are part of the result, not uncertainty bars for a calibrated probability of biological truth.</p>
<img src="figures/core_rankings.png" alt="Core-gene evidence ranking"><img src="figures/convergence.png" alt="Disease evidence versus convergence">
{_table(core.sort_values('p_core_m2',ascending=False), ['core_rank','gene_symbol','p_core_m2','p_core_m3','p_core_m5','disease_probability','expected_incoming_loci','network_coverage','p_core_m2_sensitivity_min','p_core_m2_sensitivity_max'], rows=30)}
<h2>Gene→CAD Wald-ratio screen</h2><p>Each estimate uses the strongest available non-palindromic eQTLGen cis instrument that could be allele-aligned to the separate CAD GWAS. Real CAD allele frequencies put the outcome beta and SE on the same genotype-SD scale as the approximate eQTL effect. This is an exploratory screening layer without genome-wide regional colocalization, not LD-aware causal confirmation.</p>
{_table(disease_mr.sort_values('disease_probability',ascending=False), ['gene_symbol','instrument_snp','mr_beta','mr_se','mr_p_value','mr_fdr','instrument_f','disease_probability','qc_status'], rows=25)}
<h2>Directed influence graph</h2><p>Source variants were selected from the fixed eQTLGen trans-tested SNP universe; complete target coefficients, including null results, were retained. Co-located source genes are collapsed into physical source loci before convergence is counted.</p>
<img src="figures/directed_network.png" alt="Heatmap of strongest exploratory directed gene-to-gene edges">
<p>The matrix displays the strongest harmonized pairs by model probability. Read each populated cell as row source → column target; +/− is the signed Wald-ratio estimate. Blank cells are outside this top-edge display or unavailable, not confirmed null relationships.</p>
{_table(edges.sort_values('edge_probability',ascending=False), ['source_gene_symbol','target_gene_symbol','instrument_snp','mr_beta','mr_se','mr_p_value','mr_fdr','edge_probability','edge_qc'], rows=30)}
<h2>Liver SORT1-region audit</h2><p>This independent locus analysis uses GTEx v8 liver cis-eQTLs (QTD000266 release 7) and CAD GCST90132314. Classic coloc assumes at most one causal signal per trait, so PP(H4) is shown with prior sensitivity and is not equated with mediation.</p>
<img src="figures/sort1_coloc.png" alt="SORT1-region colocalization">
{_table(sort1.sort_values('pp_h4',ascending=False), ['gene_symbol','shared_variant_count','pp_h3','pp_h4','pp_h4_given_h3_h4','top_shared_variant','mr_instrument','mr_beta','mr_se','mr_p_value'], rows=10)}
<h2>External proxy annotation audit</h2><p>Genomics England PanelApp familial-hypercholesterolaemia green genes were held out from fitting. Retrieval and enrichment metrics treat unannotated candidates as negatives for calculation only. They are diagnostic disease genes—not verified negatives/positives for CAD core-gene status—so this is face validity, not independent predictive validation.</p>
{_table(validation_metrics, list(validation_metrics.columns), rows=20)}
<h2>What the probability means</h2><p>Normal spike-and-slab Bayes factors convert each signed MR estimate and standard error into a non-null probability. Incoming probabilities are collapsed per physical source locus, then a Poisson-binomial calculation gives P(K≥m). Disease and network probabilities are combined only at the final latent-event definition; colocalization is reported separately to avoid multiplying dependent evidence based on the same variants.</p>
<h2>Major limitations</h2><ul><li>One selected cis instrument per exposure; no ancestry-matched regional LD matrix and no MR-link-2 in this compact run.</li><li>eQTLGen whole blood is a cross-tissue network layer for lipid-mediated CAD, while liver is audited only at SORT1.</li><li>eQTLGen source and target effects share samples, and trans tests cover a preselected trait-associated SNP universe.</li><li>Physical 1 Mb source-locus collapse is conservative bookkeeping, not demonstrated LD independence.</li><li>Core probabilities are conditional on candidate ascertainment, tested pairs, priors, and the spike-and-slab model.</li></ul>
<p>See <code>run_manifest.json</code>, <code>METHODS.md</code>, <code>DATA_DICTIONARY.md</code>, and <code>LIMITATIONS.md</code> for reproducibility and full assumptions.</p>
</body></html>"""
    path.write_text(document, encoding="utf-8")
