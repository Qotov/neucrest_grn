/* Neural Crest & Ectoderm GRN Explorer
   Fully client-side. Loads TSV edge lists from GRNs/, builds ego-networks
   around a focus gene, renders them with Cytoscape.js, and mirrors the
   filtered edges in a sortable/downloadable table.
   Kotov et al., PNAS 2024 (doi:10.1073/pnas.2311685121). */

'use strict';

// ---------------------------------------------------------------- config
const NETWORKS = {
  nc:       { label: 'Neural crest',   file: 'GRNs/Global_NC_GRN.tsv',       kind: 'global',    note: '16,979 predicted TF→target edges' },
  ectoderm: { label: 'Ectoderm',       file: 'GRNs/Global_Ectoderm_GRN.tsv', kind: 'global',    note: '10,086 predicted TF→target edges' },
  pax3:     { label: 'Pax3',           file: 'GRNs/pax3_GRN.tsv',            kind: 'validated', hub: 'PAX3',   note: 'GRNboost2 + morpholino + ChIP-seq' },
  tfap2e:   { label: 'TFAP2e',         file: 'GRNs/tfap2e_GRN.tsv',          kind: 'validated', hub: 'TFAP2E', note: 'GRNboost2 + morpholino + ChIP-seq' },
};

const state = {
  net: 'nc',
  gene: '',
  dir: 'both',
  minWeight: 0,
  maxEdges: 60,
  validatedOnly: false,
  sort: { key: null, asc: false },
  visibleEdges: [],
};

const cache = {};      // netKey -> parsed dataset
let cy = null;

const cssVar = (name) => getComputedStyle(document.body).getPropertyValue(name).trim();

// ---------------------------------------------------------------- parsing
function parseTSV(text, kind) {
  const lines = text.split(/\r?\n/).filter((l) => l.length);
  const header = lines[0].split('\t');
  const edges = [];
  for (let i = 1; i < lines.length; i++) {
    const c = lines[i].split('\t');
    if (kind === 'global') {
      edges.push({
        tf: c[0], target: c[1],
        weight: +c[2],
        tfExpr: c[3] === '' ? null : +c[3],
        targetExpr: c[4] === '' ? null : +c[4],
      });
    } else {
      const moChange = (c[6] || '').trim();
      const chip = c[7] === undefined || c[7].trim() === '' ? null : +c[7];
      const moL2fc = c[4] === undefined || c[4].trim() === '' ? null : +c[4];
      edges.push({
        tf: c[0], target: c[1],
        targetInMO: (c[2] || '').trim(),
        weight: +c[3],
        moLog2FC: moL2fc,
        moAbsDiff: c[5] === undefined || c[5].trim() === '' ? null : +c[5],
        moChange: moChange || 'not tested',
        chipScore: chip,
      });
    }
  }
  return { header, edges };
}

// group validated homeolog rows (foo.l / foo.s) into one edge per target
function aggregateValidated(edges) {
  const byTarget = new Map();
  for (const e of edges) {
    if (!byTarget.has(e.target)) {
      byTarget.set(e.target, {
        tf: e.tf, target: e.target, weight: e.weight,
        moChange: 'not tested', chipScore: e.chipScore,
        homeologs: [],
      });
    }
    const g = byTarget.get(e.target);
    g.weight = Math.max(g.weight, e.weight);
    if (e.chipScore != null) g.chipScore = Math.max(g.chipScore || 0, e.chipScore);
    if (e.moChange === 'decreased' || e.moChange === 'increased') {
      // any real change promotes the aggregate; keep the strongest signal
      if (g.moChange === 'not tested' || g.moChange === 'not changed') g.moChange = e.moChange;
    } else if (e.moChange === 'not changed' && g.moChange === 'not tested') {
      g.moChange = 'not changed';
    }
    if (e.targetInMO) g.homeologs.push({ id: e.targetInMO, log2fc: e.moLog2FC, change: e.moChange, absDiff: e.moAbsDiff });
  }
  return [...byTarget.values()];
}

function validationTier(e) {
  const chip = e.chipScore != null && e.chipScore > 0;
  const mo = e.moChange === 'decreased' || e.moChange === 'increased';
  if (chip && mo) return 'triple';
  if (mo) return 'mo';
  if (chip) return 'chip';
  return 'pred';
}

async function loadNetwork(key) {
  if (cache[key]) return cache[key];
  const cfg = NETWORKS[key];
  showStatus(`Loading ${cfg.label} network…`);
  let text;
  try {
    const res = await fetch(cfg.file);
    if (!res.ok) throw new Error(res.status);
    text = await res.text();
  } catch (err) {
    hideStatus();
    throw new Error(
      `Could not load ${cfg.file}. If you opened this file directly (file://), ` +
      `browsers block local data loads — use the live GitHub Pages site or run a ` +
      `local server (e.g. \`python3 -m http.server\`).`
    );
  }
  const parsed = parseTSV(text, cfg.kind);
  let edges = parsed.edges;
  if (cfg.kind === 'validated') edges = aggregateValidated(edges);

  const byTF = new Map(), byTarget = new Map(), tfSet = new Set(), geneSet = new Set();
  let maxWeight = 0;
  for (const e of edges) {
    e.tier = cfg.kind === 'validated' ? validationTier(e) : null;
    if (!byTF.has(e.tf)) byTF.set(e.tf, []);
    if (!byTarget.has(e.target)) byTarget.set(e.target, []);
    byTF.get(e.tf).push(e);
    byTarget.get(e.target).push(e);
    tfSet.add(e.tf);
    geneSet.add(e.tf); geneSet.add(e.target);
    if (e.weight > maxWeight) maxWeight = e.weight;
  }
  const ds = {
    kind: cfg.kind, hub: cfg.hub, rawEdges: parsed.edges, edges,
    byTF, byTarget, tfSet, genes: [...geneSet].sort(), maxWeight,
  };
  cache[key] = ds;
  hideStatus();
  return ds;
}

// ---------------------------------------------------------------- ego network
function buildEgo(ds) {
  const { gene, dir, minWeight, maxEdges, validatedOnly } = state;
  let pool = [];

  if (!gene) {
    // no focus: for validated nets center on the hub; for global show top edges
    if (ds.kind === 'validated') {
      pool = ds.edges.slice();
    } else {
      pool = ds.edges.slice().sort((a, b) => b.weight - a.weight);
    }
  } else {
    const asTargetOf = ds.byTarget.get(gene) || []; // regulators of gene
    const asTFof = ds.byTF.get(gene) || [];          // targets of gene
    if (dir === 'targets') pool = asTFof.slice();
    else if (dir === 'regulators') pool = asTargetOf.slice();
    else pool = asTFof.concat(asTargetOf);
  }

  pool = pool.filter((e) => e.weight >= minWeight);
  if (validatedOnly && ds.kind === 'validated') pool = pool.filter((e) => e.tier !== 'pred');
  pool.sort((a, b) => b.weight - a.weight);
  pool = pool.slice(0, maxEdges);
  return pool;
}

function renderGraph(ds, edges) {
  const emptyEl = document.getElementById('graph-empty');
  if (!edges.length) {
    if (cy) cy.elements().remove();
    emptyEl.hidden = false;
    emptyEl.textContent = state.gene
      ? `No edges for "${state.gene}" in the ${NETWORKS[state.net].label} network with the current filters.`
      : 'No edges match the current filters.';
    return;
  }
  emptyEl.hidden = true;

  const focus = state.gene ? state.gene.toUpperCase() : (ds.hub || null);
  const nodeMap = new Map();
  const addNode = (id) => {
    if (nodeMap.has(id)) return;
    const isFocus = focus && id === focus;
    const isTF = ds.tfSet.has(id);
    nodeMap.set(id, {
      data: { id, label: id, role: isFocus ? 'focus' : (isTF ? 'tf' : 'target') },
    });
  };

  const els = [];
  for (const e of edges) {
    addNode(e.tf); addNode(e.target);
    els.push({
      data: {
        id: `${e.tf}->${e.target}`,
        source: e.tf, target: e.target,
        weight: e.weight,
        tier: e.tier || 'global',
        moChange: e.moChange || '',
      },
    });
  }
  for (const n of nodeMap.values()) els.push(n);

  const wRange = [Math.min(...edges.map((e) => e.weight)), Math.max(...edges.map((e) => e.weight))];

  if (cy) cy.destroy();
  cy = cytoscape({
    container: document.getElementById('cy'),
    elements: els,
    minZoom: 0.15, maxZoom: 3,
    wheelSensitivity: 0.2,
    style: cyStyle(wRange),
    layout: layoutFor(nodeMap.size, focus),
  });

  cy.on('tap', 'node', (evt) => selectGene(evt.target.id()));
  cy.on('tap', 'edge', (evt) => showEdgeDetails(ds, evt.target.data()));
}

function layoutFor(nNodes, focus) {
  if (focus && nNodes > 2) {
    return {
      name: 'concentric',
      concentric: (n) => (n.data('role') === 'focus' ? 3 : (n.data('role') === 'tf' ? 2 : 1)),
      levelWidth: () => 1,
      minNodeSpacing: 22,
      padding: 30,
    };
  }
  return { name: 'cose', animate: false, padding: 30, nodeRepulsion: 9000, idealEdgeLength: 90 };
}

function cyStyle(wRange) {
  const [wMin, wMax] = wRange[0] === wRange[1] ? [0, wRange[1] || 1] : wRange;
  return [
    { selector: 'node', style: {
        'label': 'data(label)', 'font-size': 11, 'font-family': cssVar('--font') || 'sans-serif',
        'color': cssVar('--ink'), 'text-valign': 'center', 'text-halign': 'center',
        'text-outline-width': 2, 'text-outline-color': cssVar('--surface'),
        'width': 20, 'height': 20, 'background-color': cssVar('--target'),
        'border-width': 0,
    }},
    { selector: 'node[role="tf"]', style: { 'background-color': cssVar('--accent'), 'width': 26, 'height': 26 } },
    { selector: 'node[role="focus"]', style: {
        'background-color': cssVar('--accent-2'), 'width': 44, 'height': 44,
        'font-size': 13, 'font-weight': 'bold', 'border-width': 3, 'border-color': cssVar('--surface'),
    }},
    { selector: 'edge', style: {
        'curve-style': 'bezier', 'target-arrow-shape': 'triangle',
        'width': `mapData(weight, ${wMin}, ${wMax}, 1, 6)`,
        'line-color': cssVar('--val-pred'), 'target-arrow-color': cssVar('--val-pred'),
        'arrow-scale': 0.9, 'opacity': 0.85,
    }},
    // validated tiers
    { selector: 'edge[tier="triple"]', style: { 'line-color': cssVar('--val-triple'), 'target-arrow-color': cssVar('--val-triple'), 'opacity': 1, 'z-index': 10 } },
    { selector: 'edge[tier="mo"]',     style: { 'line-color': cssVar('--val-mo'),     'target-arrow-color': cssVar('--val-mo') } },
    { selector: 'edge[tier="chip"]',   style: { 'line-color': cssVar('--val-chip'),   'target-arrow-color': cssVar('--val-chip') } },
    { selector: 'edge[tier="global"]', style: { 'line-color': cssVar('--accent'), 'target-arrow-color': cssVar('--accent'), 'opacity': 0.55 } },
    { selector: 'node:selected', style: { 'border-width': 3, 'border-color': cssVar('--up') } },
  ];
}

// ---------------------------------------------------------------- side panel
function showNodeDetails(ds, gene) {
  const panel = document.getElementById('side-panel');
  const targets = (ds.byTF.get(gene) || []).slice().sort((a, b) => b.weight - a.weight);
  const regulators = (ds.byTarget.get(gene) || []).slice().sort((a, b) => b.weight - a.weight);
  const isTF = ds.tfSet.has(gene);
  const roleColor = isTF ? cssVar('--accent') : cssVar('--target');
  const roleLabel = isTF ? 'Transcription factor' : 'Target gene';

  const mkList = (arr, keyField) => arr.slice(0, 12).map((e) => {
    const g = e[keyField];
    return `<li data-gene="${g}"><span>${g}</span><span class="w">${e.weight.toFixed(1)}</span></li>`;
  }).join('') || '<li class="w" style="cursor:default">— none —</li>';

  panel.innerHTML = `
    <div class="node-card">
      <h3>${gene} <span class="role-tag" style="background:${roleColor}">${isTF ? 'TF' : 'target'}</span></h3>
      <p class="hint" style="border:none;margin:0;padding:0">${roleLabel} in the ${NETWORKS[state.net].label} network</p>
      <div class="stat-grid">
        <div class="stat"><div class="k">Targets</div><div class="v">${targets.length}</div></div>
        <div class="stat"><div class="k">Regulators</div><div class="v">${regulators.length}</div></div>
      </div>
      <div class="section-h">Top targets ▸</div>
      <ul class="mini-list">${mkList(targets, 'target')}</ul>
      <div class="section-h">◂ Top regulators</div>
      <ul class="mini-list">${mkList(regulators, 'tf')}</ul>
    </div>`;

  panel.querySelectorAll('.mini-list li[data-gene]').forEach((li) =>
    li.addEventListener('click', () => selectGene(li.dataset.gene)));
}

function showEdgeDetails(ds, d) {
  const panel = document.getElementById('side-panel');
  const src = d.source, tgt = d.target;
  let extra = '';
  if (ds.kind === 'validated') {
    const e = (ds.byTF.get(src) || []).find((x) => x.target === tgt);
    const tiers = {
      triple: ['Triple-validated', cssVar('--val-triple')],
      mo:     ['Morpholino-confirmed', cssVar('--val-mo')],
      chip:   ['ChIP-seq peak', cssVar('--val-chip')],
      pred:   ['Predicted only', cssVar('--val-pred')],
    };
    const tags = [];
    tags.push(`<span class="vtag" style="background:${cssVar('--accent')}">GRNboost2</span>`);
    if (e && e.chipScore != null && e.chipScore > 0) tags.push(`<span class="vtag" style="background:${cssVar('--val-chip')}">ChIP peak ${e.chipScore.toFixed(0)}</span>`);
    if (e && (e.moChange === 'decreased' || e.moChange === 'increased'))
      tags.push(`<span class="vtag" style="background:${e.moChange === 'decreased' ? cssVar('--down') : cssVar('--up')}">MO ${e.moChange}</span>`);
    const homeologs = (e && e.homeologs || []).map((h) =>
      `<li style="cursor:default"><span>${h.id}</span><span class="w">${h.log2fc != null ? 'log2FC ' + h.log2fc.toFixed(2) : '—'} · ${h.change}</span></li>`).join('');
    extra = `
      <div class="tag-row">${tags.join('')}</div>
      ${homeologs ? `<div class="section-h">Morpholino homeologs</div><ul class="mini-list">${homeologs}</ul>` : ''}`;
  } else {
    extra = `<div class="stat-grid">
      <div class="stat"><div class="k">TF expr.</div><div class="v">${fmtExpr(d, ds, src, 'tf')}</div></div>
      <div class="stat"><div class="k">Target expr.</div><div class="v">${fmtExpr(d, ds, tgt, 'target')}</div></div>
    </div>`;
  }
  panel.innerHTML = `
    <div class="edge-card">
      <div class="section-h">Regulatory edge</div>
      <div class="flow"><span style="color:${cssVar('--accent')}">${src}</span><span class="arrow">▸</span><span style="color:${cssVar('--target')}">${tgt}</span></div>
      <div class="stat" style="margin-top:10px"><div class="k">GRNboost2 weight</div><div class="v">${d.weight.toFixed(2)}</div></div>
      ${extra}
      <p class="hint">Click either gene name in the graph to re-center the view on it.</p>
    </div>`;
}

function fmtExpr(d, ds, gene, which) {
  const e = ds.rawEdges.find((x) => x.tf === d.source && x.target === d.target);
  if (!e) return '—';
  const v = which === 'tf' ? e.tfExpr : e.targetExpr;
  return v == null ? '—' : v.toFixed(0);
}

// ---------------------------------------------------------------- table
function renderTable(ds, edges) {
  const head = document.getElementById('table-head');
  const body = document.getElementById('table-body');
  const isVal = ds.kind === 'validated';
  const cols = isVal
    ? [['tf', 'TF'], ['target', 'Target'], ['weight', 'Weight'], ['moChange', 'MO change'], ['chipScore', 'ChIP score'], ['tier', 'Validation']]
    : [['tf', 'TF'], ['target', 'Target'], ['weight', 'Weight'], ['tfExpr', 'TF expr.'], ['targetExpr', 'Target expr.']];

  head.innerHTML = cols.map(([k, label]) => {
    const sorted = state.sort.key === k ? ` sorted${state.sort.asc ? ' asc' : ''}` : '';
    return `<th data-key="${k}" class="${sorted}">${label}</th>`;
  }).join('');
  head.querySelectorAll('th').forEach((th) => th.addEventListener('click', () => {
    const k = th.dataset.key;
    if (state.sort.key === k) state.sort.asc = !state.sort.asc;
    else { state.sort.key = k; state.sort.asc = false; }
    renderTable(ds, state.visibleEdges);
  }));

  let rows = edges.slice();
  const filter = document.getElementById('table-filter').value.trim().toUpperCase();
  if (filter) rows = rows.filter((e) => (e.tf + ' ' + e.target).toUpperCase().includes(filter));
  if (state.sort.key) {
    const k = state.sort.key, asc = state.sort.asc ? 1 : -1;
    rows.sort((a, b) => {
      const av = a[k], bv = b[k];
      if (typeof av === 'number' || typeof bv === 'number') return ((av || 0) - (bv || 0)) * asc;
      return String(av || '').localeCompare(String(bv || '')) * asc;
    });
  }

  const tierColor = { triple: '--val-triple', mo: '--val-mo', chip: '--val-chip', pred: '--val-pred' };
  const tierLabel = { triple: 'Triple', mo: 'MO', chip: 'ChIP', pred: 'Predicted' };

  body.innerHTML = rows.slice(0, 500).map((e) => {
    const cells = isVal ? [
      `<td class="gene-cell">${e.tf}</td>`,
      `<td class="gene-cell" data-gene="${e.target}">${e.target}</td>`,
      `<td class="num">${e.weight.toFixed(1)}</td>`,
      `<td class="chg ${e.moChange === 'decreased' ? 'decreased' : e.moChange === 'increased' ? 'increased' : 'none'}">${e.moChange}</td>`,
      `<td class="num">${e.chipScore != null ? e.chipScore.toFixed(0) : '—'}</td>`,
      `<td><span class="dot" style="background:var(${tierColor[e.tier]})"></span>${tierLabel[e.tier]}</td>`,
    ] : [
      `<td class="gene-cell" data-gene="${e.tf}">${e.tf}</td>`,
      `<td class="gene-cell" data-gene="${e.target}">${e.target}</td>`,
      `<td class="num">${e.weight.toFixed(1)}</td>`,
      `<td class="num">${e.tfExpr != null ? e.tfExpr.toFixed(0) : '—'}</td>`,
      `<td class="num">${e.targetExpr != null ? e.targetExpr.toFixed(0) : '—'}</td>`,
    ];
    return `<tr>${cells.join('')}</tr>`;
  }).join('');

  body.querySelectorAll('[data-gene]').forEach((td) =>
    td.addEventListener('click', (ev) => { ev.stopPropagation(); selectGene(td.dataset.gene); }));
  body.querySelectorAll('tr').forEach((tr, i) =>
    tr.addEventListener('click', () => selectGene(rows[i].tf)));

  document.getElementById('row-count').textContent =
    `${rows.length} edge${rows.length === 1 ? '' : 's'}${rows.length > 500 ? ' (showing 500)' : ''}`;
  document.getElementById('table-title').textContent =
    state.gene ? `Edges around ${state.gene}` : `Top edges — ${NETWORKS[state.net].label}`;
}

// ---------------------------------------------------------------- orchestration
async function refresh() {
  let ds;
  try {
    ds = await loadNetwork(state.net);
  } catch (err) {
    document.getElementById('graph-empty').hidden = false;
    document.getElementById('graph-empty').textContent = err.message;
    return;
  }
  // sync weight slider bounds to this network
  const wr = document.getElementById('weight-range');
  wr.max = Math.ceil(ds.maxWeight);
  if (state.minWeight > ds.maxWeight) { state.minWeight = 0; wr.value = 0; document.getElementById('weight-val').textContent = 0; }

  const edges = buildEgo(ds);
  state.visibleEdges = edges;
  renderGraph(ds, edges);
  renderTable(ds, edges);
  updateGraphStats(ds, edges);
  updateLegend(ds);
  if (state.gene) showNodeDetails(ds, state.gene);
  syncURL();
}

function updateGraphStats(ds, edges) {
  const nodes = new Set();
  edges.forEach((e) => { nodes.add(e.tf); nodes.add(e.target); });
  let txt = `${nodes.size} genes · ${edges.length} edges shown`;
  if (ds.kind === 'validated') {
    const v = edges.filter((e) => e.tier !== 'pred').length;
    txt += ` · ${v} validated`;
  }
  document.getElementById('graph-stats').textContent = txt;
}

function updateLegend(ds) {
  const el = document.getElementById('legend');
  if (ds.kind === 'validated') {
    el.innerHTML = `
      <div class="row"><span class="swatch line" style="background:var(--val-triple)"></span>Triple-validated (MO + ChIP)</div>
      <div class="row"><span class="swatch line" style="background:var(--val-mo)"></span>Morpholino-confirmed</div>
      <div class="row"><span class="swatch line" style="background:var(--val-chip)"></span>ChIP-seq peak</div>
      <div class="row"><span class="swatch line" style="background:var(--val-pred)"></span>Predicted only</div>
      <div class="row"><span class="swatch" style="background:var(--accent-2)"></span>Focus gene</div>`;
  } else {
    el.innerHTML = `
      <div class="row"><span class="swatch" style="background:var(--accent-2)"></span>Focus gene</div>
      <div class="row"><span class="swatch" style="background:var(--accent)"></span>Transcription factor</div>
      <div class="row"><span class="swatch" style="background:var(--target)"></span>Target-only gene</div>
      <div class="row"><span class="swatch line" style="background:var(--accent)"></span>Edge width ∝ GRNboost2 weight</div>`;
  }
}

function selectGene(gene) {
  state.gene = gene;
  document.getElementById('gene-input').value = gene;
  refresh();
}

// ---------------------------------------------------------------- URL state
function syncURL() {
  const p = new URLSearchParams();
  p.set('net', state.net);
  if (state.gene) p.set('gene', state.gene);
  if (state.dir !== 'both') p.set('dir', state.dir);
  history.replaceState(null, '', '?' + p.toString());
}

function readURL() {
  const p = new URLSearchParams(location.search);
  if (p.get('net') && NETWORKS[p.get('net')]) state.net = p.get('net');
  if (p.get('gene')) state.gene = p.get('gene').toUpperCase();
  if (p.get('dir')) state.dir = p.get('dir');
}

// ---------------------------------------------------------------- UI wiring
function buildNetworkChips() {
  const wrap = document.getElementById('network-chips');
  wrap.innerHTML = Object.entries(NETWORKS).map(([k, cfg]) => {
    const tag = cfg.kind === 'validated' ? '<span class="badge">✓ validated</span>' : '';
    return `<button class="chip${k === state.net ? ' active' : ''}" data-net="${k}" title="${cfg.note}">${cfg.label}${tag}</button>`;
  }).join('');
  wrap.querySelectorAll('.chip').forEach((c) => c.addEventListener('click', () => {
    state.net = c.dataset.net;
    state.gene = '';
    document.getElementById('gene-input').value = '';
    state.sort = { key: null, asc: false };
    wrap.querySelectorAll('.chip').forEach((x) => x.classList.toggle('active', x === c));
    document.querySelector('.validation-only').hidden = NETWORKS[state.net].kind !== 'validated';
    loadNetwork(state.net).then(fillGeneList).then(refresh);
  }));
}

async function fillGeneList() {
  const ds = cache[state.net];
  if (!ds) return;
  document.getElementById('gene-list').innerHTML =
    ds.genes.map((g) => `<option value="${g}"></option>`).join('');
}

function debounce(fn, ms) {
  let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); };
}

function wireControls() {
  const geneInput = document.getElementById('gene-input');
  const applyGene = () => {
    const v = geneInput.value.trim().toUpperCase();
    const ds = cache[state.net];
    if (v && ds && !ds.byTF.has(v) && !ds.byTarget.has(v)) return; // ignore unknown until valid
    state.gene = v;
    refresh();
  };
  geneInput.addEventListener('change', applyGene);
  geneInput.addEventListener('input', debounce(() => {
    const v = geneInput.value.trim().toUpperCase();
    const ds = cache[state.net];
    if (v && ds && (ds.byTF.has(v) || ds.byTarget.has(v))) { state.gene = v; refresh(); }
  }, 350));

  document.getElementById('gene-clear').addEventListener('click', () => {
    geneInput.value = ''; state.gene = ''; refresh();
  });

  document.querySelectorAll('#direction-seg button').forEach((b) =>
    b.addEventListener('click', () => {
      state.dir = b.dataset.dir;
      document.querySelectorAll('#direction-seg button').forEach((x) => x.classList.toggle('active', x === b));
      refresh();
    }));

  const wr = document.getElementById('weight-range');
  wr.addEventListener('input', debounce(() => {
    state.minWeight = +wr.value;
    document.getElementById('weight-val').textContent = wr.value;
    refresh();
  }, 200));

  const mr = document.getElementById('maxnodes-range');
  mr.addEventListener('input', debounce(() => {
    state.maxEdges = +mr.value;
    document.getElementById('maxnodes-val').textContent = mr.value;
    refresh();
  }, 200));

  document.getElementById('validated-only').addEventListener('change', (e) => {
    state.validatedOnly = e.target.checked; refresh();
  });

  document.getElementById('table-filter').addEventListener('input', debounce(() =>
    renderTable(cache[state.net], state.visibleEdges), 200));

  document.getElementById('download-btn').addEventListener('click', downloadCSV);
  document.getElementById('fit-btn').addEventListener('click', () => cy && cy.fit(null, 40));
  document.getElementById('png-btn').addEventListener('click', downloadPNG);

  document.getElementById('theme-toggle').addEventListener('click', toggleTheme);
}

function downloadCSV() {
  const ds = cache[state.net];
  const isVal = ds.kind === 'validated';
  const cols = isVal
    ? ['tf', 'target', 'weight', 'moChange', 'chipScore', 'tier']
    : ['tf', 'target', 'weight', 'tfExpr', 'targetExpr'];
  const head = isVal
    ? ['TF', 'target', 'GRNboost2_weight', 'MO_change', 'MACS2_chip_score', 'validation_tier']
    : ['TF', 'target', 'weight', 'tf_expression', 'target_expression'];
  const lines = [head.join(',')];
  for (const e of state.visibleEdges) {
    lines.push(cols.map((c) => {
      const v = e[c];
      return v == null ? '' : (typeof v === 'number' ? v : `"${String(v).replace(/"/g, '""')}"`);
    }).join(','));
  }
  const blob = new Blob([lines.join('\n')], { type: 'text/csv' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `${state.net}_GRN${state.gene ? '_' + state.gene : ''}_filtered.csv`;
  a.click();
  URL.revokeObjectURL(a.href);
}

function downloadPNG() {
  if (!cy) return;
  const png = cy.png({ full: true, scale: 2, bg: cssVar('--plane') });
  const a = document.createElement('a');
  a.href = png;
  a.download = `${state.net}_GRN${state.gene ? '_' + state.gene : ''}.png`;
  a.click();
}

// ---------------------------------------------------------------- theme
function toggleTheme() {
  const cur = document.documentElement.getAttribute('data-theme');
  const next = cur === 'dark' ? 'light' : (cur === 'light' ? 'dark' : (matchMedia('(prefers-color-scheme: dark)').matches ? 'light' : 'dark'));
  document.documentElement.setAttribute('data-theme', next);
  try { localStorage.setItem('grn-theme', next); } catch (e) {}
  if (cache[state.net]) refresh(); // restyle graph to new theme vars
}

function initTheme() {
  let t;
  try { t = localStorage.getItem('grn-theme'); } catch (e) {}
  if (!t) t = matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
  document.documentElement.setAttribute('data-theme', t);
}

// ---------------------------------------------------------------- status toast
let statusTimer;
function showStatus(msg) {
  const el = document.getElementById('load-status');
  el.textContent = msg; el.hidden = false;
  clearTimeout(statusTimer);
}
function hideStatus() {
  const el = document.getElementById('load-status');
  statusTimer = setTimeout(() => { el.hidden = true; }, 250);
}

// ---------------------------------------------------------------- boot
function boot() {
  if (typeof cytoscape === 'undefined') {
    document.getElementById('graph-empty').hidden = false;
    document.getElementById('graph-empty').textContent =
      'The graph library (cytoscape.js) failed to load — check your network connection and reload.';
    return;
  }
  initTheme();
  readURL();
  buildNetworkChips();
  // sync UI to state read from URL
  if (state.dir !== 'both') {
    document.querySelectorAll('#direction-seg button').forEach((x) =>
      x.classList.toggle('active', x.dataset.dir === state.dir));
  }
  document.querySelector('.validation-only').hidden = NETWORKS[state.net].kind !== 'validated';
  if (state.gene) document.getElementById('gene-input').value = state.gene;
  wireControls();
  loadNetwork(state.net).then(fillGeneList).then(refresh).catch((err) => {
    document.getElementById('graph-empty').hidden = false;
    document.getElementById('graph-empty').textContent = err.message;
  });
}

document.addEventListener('DOMContentLoaded', boot);
