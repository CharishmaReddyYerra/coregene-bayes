const state = {
  data: null,
  threshold: "m2",
  selectedGene: new URL(window.location.href).searchParams.get("gene")?.toUpperCase() || "LDLR",
  visibleGenes: 10,
  edgeFilter: "supported"
};

const elements = {
  ranking: document.querySelector("#ranking-chart"),
  expand: document.querySelector("#ranking-expand"),
  dossier: document.querySelector("#gene-dossier"),
  search: document.querySelector("#gene-search"),
  options: document.querySelector("#gene-options"),
  network: document.querySelector("#network-svg"),
  edgeBody: document.querySelector("#edge-table-body"),
  edgeCount: document.querySelector("#edge-count"),
  locus: document.querySelector("#locus-chart"),
  sensitivity: document.querySelector("#sensitivity-cases"),
  attrition: document.querySelector("#attrition-list"),
  sources: document.querySelector("#source-table"),
  tooltip: document.querySelector("#tooltip")
};

const baseline = { m2: 0.025314764851508187, m3: 0.007394003418227538, m5: 0.0003010924506458796 };

const escapeHtml = value => String(value ?? "")
  .replaceAll("&", "&amp;")
  .replaceAll("<", "&lt;")
  .replaceAll(">", "&gt;")
  .replaceAll('"', "&quot;")
  .replaceAll("'", "&#039;");

const clamp = value => Math.max(0, Math.min(1, Number(value) || 0));
const probability = value => Number.isFinite(value) ? value.toFixed(3) : "—";
const estimate = value => Number.isFinite(value) ? value.toFixed(3) : "—";
const scientific = value => {
  if (!Number.isFinite(value)) return "—";
  if (value === 0) return "0";
  return value < 0.001 ? value.toExponential(1).replace("e-", "e−") : value.toFixed(3);
};
const integer = value => Number.isFinite(value) ? Math.round(value).toLocaleString() : "—";
const bytes = value => {
  if (!Number.isFinite(value)) return "—";
  const units = ["B", "KB", "MB", "GB"];
  let size = value;
  let unit = 0;
  while (size >= 1024 && unit < units.length - 1) { size /= 1024; unit += 1; }
  return `${size.toFixed(unit > 1 ? 2 : 0)} ${units[unit]}`;
};

function selectedRecord() {
  return state.data?.genes.find(item => item.gene === state.selectedGene) || state.data?.genes.find(item => item.gene === "LDLR");
}

function setSelectedGene(gene, options = {}) {
  const record = state.data.genes.find(item => item.gene === gene);
  if (!record) return;
  state.selectedGene = record.gene;
  elements.search.value = record.gene;
  const url = new URL(window.location.href);
  url.searchParams.set("gene", record.gene);
  history.replaceState({}, "", url);
  renderRanking();
  renderDossier();
  if (state.edgeFilter === "selected") renderNetwork();
  if (options.scroll) elements.dossier.scrollIntoView({ behavior: "smooth", block: "nearest" });
}

function renderRanking() {
  const eligible = state.data.genes
    .filter(item => Number.isFinite(item[state.threshold]))
    .sort((a, b) => b[state.threshold] - a[state.threshold]);
  const visible = eligible.slice(0, state.visibleGenes);
  const selected = selectedRecord();
  if (selected && Number.isFinite(selected[state.threshold]) && !visible.some(item => item.gene === selected.gene)) visible.push(selected);
  const baselineValue = baseline[state.threshold];
  const baselinePct = clamp(baselineValue) * 100;
  const showRange = state.threshold === "m2";
  const rangeKey = document.querySelector(".chart-key .line")?.parentElement;
  if (rangeKey) rangeKey.hidden = !showRange;

  elements.ranking.innerHTML = `
    <div class="axis"><span>GENE</span><span class="axis-track"><i>0</i><i>.25</i><i>.50</i><i>.75</i><i>1.0</i></span><span>P(CORE)</span></div>
    ${visible.map((item, index) => {
      const pinned = index >= state.visibleGenes;
      const selectedClass = item.gene === state.selectedGene ? " selected" : "";
      const range = showRange && Number.isFinite(item.sensitivityMin) && Number.isFinite(item.sensitivityMax)
        ? `<i class="range" style="left:${clamp(item.sensitivityMin) * 100}%;width:${Math.max(0, clamp(item.sensitivityMax) - clamp(item.sensitivityMin)) * 100}%"></i>`
        : "";
      return `<div class="gene-row${selectedClass}" data-gene="${escapeHtml(item.gene)}">
        <button class="gene-name" type="button" aria-label="Inspect ${escapeHtml(item.gene)}">${escapeHtml(item.gene)}${item.panelAppPositive ? "<small>◆</small>" : ""}${pinned ? "<small>selected</small>" : ""}</button>
        <span class="track" tabindex="0" aria-label="${escapeHtml(item.gene)}: core probability for ${state.threshold.replace("m", "m equals ")} is ${probability(item[state.threshold])}${showRange ? `; sensitivity range ${probability(item.sensitivityMin)} to ${probability(item.sensitivityMax)}` : ""}" data-tooltip="${escapeHtml(item.gene)} · P(core | ${state.threshold.replace("m", "m=")}) ${probability(item[state.threshold])}${showRange ? ` · sensitivity ${probability(item.sensitivityMin)}–${probability(item.sensitivityMax)}` : ""}">
          <i class="baseline" style="left:${baselinePct}%"></i>${range}
          <i class="point" style="left:${clamp(item[state.threshold]) * 100}%"></i>
        </span>
        <span class="value">${probability(item[state.threshold])}</span>
      </div>`;
    }).join("")}
  `;

  elements.ranking.querySelectorAll(".gene-name").forEach(button => {
    button.addEventListener("click", () => setSelectedGene(button.closest(".gene-row").dataset.gene));
  });
  elements.ranking.querySelectorAll("[data-tooltip]").forEach(bindTooltip);
  elements.expand.textContent = state.visibleGenes === 10 ? "Show top 25" : state.visibleGenes === 25 ? `Show all ${eligible.length}` : "Return to top 10";
}

function renderDossier() {
  const item = selectedRecord();
  if (!item) return;
  const ciLow = Number.isFinite(item.mrBeta) && Number.isFinite(item.mrSe) ? item.mrBeta - 1.96 * item.mrSe : null;
  const ciHigh = Number.isFinite(item.mrBeta) && Number.isFinite(item.mrSe) ? item.mrBeta + 1.96 * item.mrSe : null;
  const geneEdges = state.data.edges.filter(edge => edge.source === item.gene || edge.target === item.gene);
  const tags = [
    item.panelAppPositive ? '<span class="positive">PanelApp FH proxy ◆</span>' : "",
    item.bayesv2gSeed ? "<span>BayesV2G seed</span>" : "",
    Number.isFinite(item.rank) ? `<span>eligible rank ${integer(item.rank)}</span>` : "<span>not primary-rank eligible</span>"
  ].filter(Boolean).join("");

  elements.dossier.innerHTML = `
    <p class="eyebrow">Selected candidate</p>
    <h3>${escapeHtml(item.gene)}</h3>
    <p class="dossier-name">${escapeHtml(item.name)}</p>
    <div class="dossier-meta">
      <span>CHR <b>${escapeHtml(item.chromosome || "—")}</b></span>
      <span>L2G <b>${probability(item.l2g)}</b></span>
      <span>INSTRUMENT <b>${escapeHtml(item.instrument || "—")}</b></span>
    </div>
    <div class="evidence-list">
      ${evidenceRow("Primary core probability", probability(item.m2), item.m2, "m=2")}
      ${evidenceRow("Stricter convergence", `${probability(item.m3)} / ${probability(item.m5)}`, item.m3, "m=3 / m=5")}
      ${evidenceRow("Disease evidence", probability(item.diseaseProbability), item.diseaseProbability, "posterior")}
      ${evidenceRow("Expected incoming loci", Number.isFinite(item.expectedIncoming) ? item.expectedIncoming.toFixed(2) : "—", Math.min((item.expectedIncoming || 0) / 5, 1), `${integer(item.sourceLociTested)} loci tested`)}
      ${evidenceRow("Network coverage", Number.isFinite(item.networkCoverage) ? `${(item.networkCoverage * 100).toFixed(1)}%` : "—", item.networkCoverage, `${integer(item.sourceLociTested)} / ${integer(item.sourceLocusUniverse)}`)}
    </div>
    <div class="mr-record">
      <span>Gene → CAD Wald ratio</span>
      <b>β ${estimate(item.mrBeta)} · 95% CI ${estimate(ciLow)} to ${estimate(ciHigh)}</b>
      <span>SE ${estimate(item.mrSe)} · p ${scientific(item.mrP)} · F ${Number.isFinite(item.instrumentF) ? item.instrumentF.toFixed(1) : "—"}</span>
    </div>
    <div class="dossier-tags">${tags}</div>
    <p class="dossier-note">${geneEdges.length ? `${geneEdges.length} displayed high-PIP directed ${geneEdges.length === 1 ? "edge involves" : "edges involve"} ${escapeHtml(item.gene)}.` : `No PIP ≥ 0.50 directed edge involving ${escapeHtml(item.gene)} is displayed.`} Ratio evidence is compatible with an effect; it is not causal confirmation.</p>
  `;
}

function evidenceRow(label, value, width, note) {
  return `<div class="evidence-item"><header><span>${escapeHtml(label)} · ${escapeHtml(note)}</span><strong>${escapeHtml(value)}</strong></header><div class="meter"><i style="width:${clamp(width) * 100}%"></i></div></div>`;
}

function filteredEdges() {
  if (state.edgeFilter === "supported") return state.data.edges.filter(edge => Number.isFinite(edge.fdr) && edge.fdr <= 0.05);
  if (state.edgeFilter === "selected") return state.data.edges.filter(edge => edge.source === state.selectedGene || edge.target === state.selectedGene);
  return state.data.edges;
}

function renderNetwork() {
  const edges = filteredEdges();
  elements.edgeCount.textContent = `${edges.length.toLocaleString()} / ${state.data.summary.directed_pair_count.toLocaleString()}`;
  elements.edgeBody.innerHTML = edges.length ? edges.map(edge => `
    <tr>
      <td><button class="edge-link" type="button" data-gene="${escapeHtml(edge.target)}" aria-label="Inspect target gene ${escapeHtml(edge.target)}">${escapeHtml(edge.source)} <span class="edge-arrow">→</span> ${escapeHtml(edge.target)}</button></td>
      <td>${probability(edge.probability)}</td>
      <td>${scientific(edge.fdr)}</td>
      <td>${estimate(edge.beta)}</td>
    </tr>`).join("") : `<tr><td colspan="4">No PIP ≥ 0.50 edge involves ${escapeHtml(state.selectedGene)}.</td></tr>`;

  elements.edgeBody.querySelectorAll(".edge-link").forEach(button => {
    button.addEventListener("click", () => setSelectedGene(button.dataset.gene));
  });

  if (!edges.length) {
    elements.network.innerHTML = `<text x="360" y="205" text-anchor="middle" class="network-empty">No high-PIP edge for ${escapeHtml(state.selectedGene)}</text><text x="360" y="235" text-anchor="middle" class="network-empty">Try “All PIP ≥ .50” to restore the network.</text>`;
    return;
  }

  const sources = [...new Set(edges.map(edge => edge.source))].sort();
  const targets = [...new Set(edges.map(edge => edge.target))].sort();
  const sourceY = positions(sources, 42, 390);
  const targetY = positions(targets, 42, 390);
  const paths = edges.map((edge, index) => {
    const y1 = sourceY.get(edge.source);
    const y2 = targetY.get(edge.target);
    const selected = edge.source === state.selectedGene || edge.target === state.selectedGene;
    return `<path class="network-edge${edge.beta < 0 ? " negative" : ""}${selected ? " selected-edge" : ""}" tabindex="0" data-index="${index}" d="M 130 ${y1} C 285 ${y1}, 435 ${y2}, 590 ${y2}" stroke-width="${(1.2 + edge.probability * 3).toFixed(2)}" marker-end="url(#arrow)"><title>${escapeHtml(edge.source)} to ${escapeHtml(edge.target)}, PIP ${probability(edge.probability)}, beta ${estimate(edge.beta)}</title></path>`;
  }).join("");
  const sourceNodes = sources.map(gene => nodeMarkup(gene, 130, sourceY.get(gene), "source")).join("");
  const targetNodes = targets.map(gene => nodeMarkup(gene, 590, targetY.get(gene), "target")).join("");
  elements.network.innerHTML = `
    <defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#66706d"></path></marker></defs>
    ${paths}${sourceNodes}${targetNodes}
  `;
  elements.network.querySelectorAll("[data-gene]").forEach(node => {
    node.addEventListener("click", () => setSelectedGene(node.dataset.gene));
    node.addEventListener("keydown", event => { if (event.key === "Enter" || event.key === " ") { event.preventDefault(); setSelectedGene(node.dataset.gene); } });
  });
  elements.network.querySelectorAll("[data-index]").forEach(path => {
    const edge = edges[Number(path.dataset.index)];
    bindTooltip(path, `${edge.source} → ${edge.target} · PIP ${probability(edge.probability)} · FDR ${scientific(edge.fdr)} · β ${estimate(edge.beta)}`);
  });
}

function positions(items, start, end) {
  const map = new Map();
  const step = items.length > 1 ? (end - start) / (items.length - 1) : 0;
  items.forEach((item, index) => map.set(item, items.length === 1 ? (start + end) / 2 : start + index * step));
  return map;
}

function nodeMarkup(gene, x, y, side) {
  const selected = gene === state.selectedGene;
  const labelX = side === "source" ? x - 14 : x + 14;
  const anchor = side === "source" ? "end" : "start";
  return `<g tabindex="0" role="button" aria-label="Inspect ${escapeHtml(gene)}" data-gene="${escapeHtml(gene)}">
    <circle cx="${x}" cy="${y}" r="6" class="network-node${selected ? " selected" : ""}"></circle>
    <text x="${labelX}" y="${y + 4}" text-anchor="${anchor}" class="network-label${selected ? " selected" : ""}">${escapeHtml(gene)}</text>
  </g>`;
}

function renderLocus() {
  const order = ["SORT1", "PSRC1", "CELSR2"];
  const rows = [...state.data.coloc].sort((a, b) => order.indexOf(a.gene) - order.indexOf(b.gene));
  elements.locus.innerHTML = `
    <p class="scroll-cue">Scroll horizontally to inspect the full locus figure →</p>
    <div class="locus-scale"><span>GENE</span><div><i>.995</i><i>.9975</i><i>1.000</i></div><span>PP(H4)</span></div>
    ${rows.map(item => {
      const zoomed = clamp((item.ppH4 - 0.995) / 0.005) * 100;
      return `<div class="locus-row">
        <h3>${escapeHtml(item.gene)}</h3>
        <div class="coloc-track" tabindex="0" aria-label="${escapeHtml(item.gene)}: posterior probability of a shared signal ${item.ppH4.toFixed(6)}; prior sensitivity ${item.priorLow.toFixed(6)} to ${item.priorHigh.toFixed(6)}" data-tooltip="${escapeHtml(item.gene)} · PP(H4) ${item.ppH4.toFixed(6)} · prior sensitivity ${item.priorLow.toFixed(6)}–${item.priorHigh.toFixed(6)}">
          <div class="coloc-bar" style="width:${zoomed}%"></div><i class="coloc-pin" style="left:${zoomed}%"></i>
        </div>
        <div class="coloc-value">${item.ppH4.toFixed(6)}</div>
        <div class="locus-detail"><span>β ${estimate(item.beta)} ± ${estimate(item.se)}</span><span>p ${scientific(item.p)}</span><span>F ${item.instrumentF.toFixed(1)}</span></div>
      </div>`;
    }).join("")}
    <div class="locus-legend"><span><i></i>zoomed PP(H4) scale</span><span><i></i>observed value</span></div>
  `;
  elements.locus.querySelectorAll("[data-tooltip]").forEach(bindTooltip);
}

function renderSensitivity() {
  const cases = [
    ["HP", "Primary m=2 leader; broad prior dependence."],
    ["LDLR", "Stricter-threshold leader; narrower configured range."],
    ["SORT1", "Strong disease evidence; prior-sensitive convergence."]
  ];
  elements.sensitivity.innerHTML = cases.map(([gene, caption]) => {
    const item = state.data.genes.find(row => row.gene === gene);
    const left = clamp(item.sensitivityMin) * 100;
    const width = (clamp(item.sensitivityMax) - clamp(item.sensitivityMin)) * 100;
    return `<article class="sensitivity-card">
      <header><h3>${gene}</h3><span>P(core) ${probability(item.m2)}</span></header>
      <div class="sensitivity-track">
        <i class="sensitivity-range" style="left:${left}%;width:${width}%"></i>
        <i class="sensitivity-point" style="left:${clamp(item.m2) * 100}%"></i>
      </div>
      <div class="sensitivity-values"><span>${probability(item.sensitivityMin)}</span><span>configured range</span><span>${probability(item.sensitivityMax)}</span></div>
      <p class="sensitivity-caption">${caption}</p>
    </article>`;
  }).join("");
}

function renderAttrition() {
  const max = Math.max(...state.data.attrition.map(item => item.count));
  elements.attrition.innerHTML = state.data.attrition.map(item => `
    <div class="attrition-row" data-tooltip="${escapeHtml(item.note)}">
      <b>${escapeHtml(item.stage)}</b>
      <span class="attrition-bar"><i style="width:${(item.count / max) * 100}%"></i></span>
      <strong>${integer(item.count)}</strong>
      <span class="sr-only">${escapeHtml(item.note)}</span>
    </div>`).join("");
  elements.attrition.querySelectorAll("[data-tooltip]").forEach(bindTooltip);
}

function renderSources() {
  const { sources } = state.data;
  const raw = sources.raw;
  const rows = [
    { name: "Open Targets Genetics", detail: `LDL candidate ascertainment · release ${sources.selection.open_targets_release.year}.${sources.selection.open_targets_release.month} · ${sources.selection.source_study}`, meta: `${sources.selection.significant_locus_count} loci · L2G ≥ ${sources.selection.l2g_score_threshold}`, url: sources.urls.open_targets },
    { name: "eQTLGen Phase I", detail: "Whole-blood cis instruments and complete trans-tested coefficients", meta: `n up to 31,684 · cis ${bytes(raw.eqtlgen_cis.bytes)} · trans ${bytes(raw.eqtlgen_trans.bytes)}`, url: sources.urls.eqtlgen_cis },
    { name: "CAD GWAS Catalog", detail: `${sources.analysis.disease_gwas} · harmonized summary statistics`, meta: `${bytes(raw.cad_gwas.bytes)} · SHA ${raw.cad_gwas.sha256.slice(0, 12)}…`, url: sources.urls.cad_gwas },
    { name: "GTEx v8 liver", detail: "eQTL Catalogue r7 QTD000266 · locus-specific audit", meta: `n=208 · ${sources.counts.liver_sort1_rows.toLocaleString()} compact rows`, url: sources.urls.gtex_liver },
    { name: "PanelApp", detail: "Familial hypercholesterolaemia panel 772 · v2.6", meta: "External proxy annotation · not predictive validation", url: sources.urls.panelapp }
  ];
  elements.sources.innerHTML = rows.map(row => `<div class="source-row"><h3>${escapeHtml(row.name)}</h3><p>${escapeHtml(row.detail)}</p><code>${escapeHtml(row.meta)}</code><a href="${escapeHtml(row.url)}" target="_blank" rel="noreferrer">source ↗</a></div>`).join("") + `<p class="provenance-fine">Raw-file SHA-256 checksums and byte counts are preserved in the downloadable source manifest. Acquisition used checksum verification; transport TLS verification was disabled for the recorded fetch.</p>`;
}

function bindTooltip(node, explicitText) {
  const text = explicitText || node.dataset.tooltip;
  const show = event => {
    elements.tooltip.textContent = text;
    elements.tooltip.classList.add("visible");
    moveTooltip(event);
  };
  node.addEventListener("mouseenter", show);
  node.addEventListener("mousemove", moveTooltip);
  node.addEventListener("mouseleave", hideTooltip);
  node.addEventListener("focus", () => {
    const box = node.getBoundingClientRect();
    show({ clientX: box.left + box.width / 2, clientY: box.top });
  });
  node.addEventListener("blur", hideTooltip);
}

function moveTooltip(event) {
  const x = Math.min(window.innerWidth - 300, Math.max(10, event.clientX + 12));
  const y = Math.min(window.innerHeight - 100, Math.max(10, event.clientY + 12));
  elements.tooltip.style.left = `${x}px`;
  elements.tooltip.style.top = `${y}px`;
}
function hideTooltip() { elements.tooltip.classList.remove("visible"); }

function wireControls() {
  document.querySelectorAll("[data-threshold]").forEach(button => button.addEventListener("click", () => {
    state.threshold = button.dataset.threshold;
    document.querySelectorAll("[data-threshold]").forEach(item => {
      const active = item === button;
      item.classList.toggle("active", active);
      item.setAttribute("aria-pressed", String(active));
    });
    renderRanking();
  }));

  document.querySelectorAll("[data-edge-filter]").forEach(button => button.addEventListener("click", () => {
    state.edgeFilter = button.dataset.edgeFilter;
    document.querySelectorAll("[data-edge-filter]").forEach(item => {
      const active = item === button;
      item.classList.toggle("active", active);
      item.setAttribute("aria-pressed", String(active));
    });
    renderNetwork();
  }));

  elements.expand.addEventListener("click", () => {
    const eligibleCount = state.data.genes.filter(item => Number.isFinite(item[state.threshold])).length;
    state.visibleGenes = state.visibleGenes === 10 ? Math.min(25, eligibleCount) : state.visibleGenes < eligibleCount ? eligibleCount : 10;
    renderRanking();
  });

  const chooseSearch = () => {
    const query = elements.search.value.trim().toLowerCase();
    if (!query) return;
    const record = state.data.genes.find(item => item.gene.toLowerCase() === query)
      || state.data.genes.find(item => item.name.toLowerCase() === query)
      || state.data.genes.find(item => item.gene.toLowerCase().startsWith(query) || item.name.toLowerCase().includes(query));
    if (record) setSelectedGene(record.gene, { scroll: true });
  };
  elements.search.addEventListener("change", chooseSearch);
  elements.search.addEventListener("keydown", event => { if (event.key === "Enter") { event.preventDefault(); chooseSearch(); } });
}

async function init() {
  elements.ranking.innerHTML = '<p class="dossier-loading">Loading ranked candidates…</p>';
  try {
    const response = await fetch("data.json");
    if (!response.ok) throw new Error(`Data request returned ${response.status}`);
    state.data = await response.json();
    if (!state.data.genes.some(item => item.gene === state.selectedGene)) state.selectedGene = "LDLR";
    elements.options.innerHTML = [...state.data.genes]
      .sort((a, b) => a.gene.localeCompare(b.gene))
      .map(item => `<option value="${escapeHtml(item.gene)}">${escapeHtml(item.name)}</option>`).join("");
    elements.search.value = state.selectedGene;
    wireControls();
    renderRanking();
    renderDossier();
    renderNetwork();
    renderLocus();
    renderSensitivity();
    renderAttrition();
    renderSources();
    if (window.location.hash) {
      const restoreSection = () => document.querySelector(window.location.hash)?.scrollIntoView();
      requestAnimationFrame(() => requestAnimationFrame(restoreSection));
      window.setTimeout(restoreSection, 120);
    }
  } catch (error) {
    console.error(error);
    const message = '<div class="error-state">The analytical snapshot could not be loaded. Reload the page or use the downloadable source tables below.</div>';
    elements.ranking.innerHTML = message;
    elements.dossier.innerHTML = message;
    elements.edgeBody.innerHTML = `<tr><td colspan="4">Analytical snapshot unavailable.</td></tr>`;
    elements.locus.innerHTML = message;
  }
}

init();
