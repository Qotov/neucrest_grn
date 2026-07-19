# Neural Crest &amp; Ectoderm Gene Regulatory Networks

Gene regulatory networks (GRNs) for *Xenopus* neural crest and ectoderm development,
predicted with **GRNboost2** from single-cell transcriptomes and experimentally
validated by **morpholino knockdown** (RNA-seq) and **transcription-factor ChIP-seq**.

Companion data & tools for:

> Kotov A, Seal S, Alkobtawi M, Kappès V, Medina Ruiz S, Arbès H, Harland RM,
> Peshkin L, Monsoro-Burq AH. **A time-resolved single-cell roadmap of the logic
> driving anterior neural crest diversification from neural border to migration
> stages.** *PNAS* 121, e2311685121 (2024).
> [doi:10.1073/pnas.2311685121](https://www.pnas.org/doi/10.1073/pnas.2311685121)

## 🔎 Interactive explorer

**→ [Launch the GRN Explorer](https://qotov.github.io/neucrest_grn/)** (GitHub Pages)

Search any gene to see its regulators and targets, colored by how each interaction
was validated. Runs entirely in your browser — no install, no server.

- **Gene-centered views** — type a gene (e.g. `PAX3`, `SOX9`, `ZIC1`) to build its
  ego-network of upstream regulators and downstream targets.
- **Four networks** — the global **Neural crest** and **Ectoderm** connectomes, plus
  the **Pax3** and **TFAP2e** sub-networks validated by morpholino + ChIP-seq.
- **Cranial / Vagal subviews** — on the NC network, filter to the neighborhood of
  curated cranial (hox-negative) or vagal (hox-positive + *mafb/egr2*) lineage
  regulators from Fig. 2–3. These are a transparent filter on the single global NC
  GRN, not separately inferred sub-networks (the gene panels are an editable constant
  in `assets/js/app.js`).
- **Validation at a glance** — in the validated networks, edges are colored
  predicted-only · morpholino-confirmed · ChIP-seq peak · **triple-validated**.
- **Filter, sort, export** — weight threshold, edge cap, sortable table, and one-click
  CSV download of the current subnetwork. Views are deep-linkable, e.g.
  `?net=nc&gene=SOX9`.

> Replaces the previous Heroku apps (`ectoapp` / `ncapp`), which went offline when
> Heroku retired its free tier.

## 📊 Data tables (`GRNs/`)

All tables are tab-separated (`.tsv`), directed **TF → target**.

### Global connectomes

`Global_NC_GRN.tsv` (16,978 edges) · `Global_Ectoderm_GRN.tsv` (10,085 edges)

| column | description |
|---|---|
| `TF` | regulator (transcription factor) |
| `target` | predicted target gene |
| `weight` | GRNboost2 importance score for the TF→target link |
| `tf_expression` | summed TF expression in the relevant cell population |
| `target_expression` | summed target expression in the relevant cell population |

### Validated sub-networks

`pax3_GRN.tsv` · `tfap2e_GRN.tsv` — one row per target *homeolog* (`.l` / `.s`
*X. laevis* subgenomes).

| column | description |
|---|---|
| `TF` | the validated regulator (`PAX3` or `TFAP2E`) |
| `target` | predicted target gene |
| `target_in_MO` | homeolog measured in the morpholino RNA-seq |
| `GRNboost2_weight` | GRNboost2 importance score |
| `MO_log2FC` | log2 fold-change of the target after TF knockdown |
| `MO_abs_diff` | absolute expression difference in the morphant |
| `MO_change` | `decreased` / `increased` / `not changed` after knockdown |
| `MACS2_score (chip)` | MACS2 peak score if the TF binds near the target (ChIP-seq); empty if no peak |

An interaction is **triple-validated** when it is GRNboost2-predicted, changed upon
knockdown, **and** carries a ChIP-seq peak.

## 🧬 How the networks were built

Single-cell transcriptomes from eight consecutive *X. tropicalis* developmental
stages (Briggs et al. dataset, re-sequenced and re-annotated) yielded 6,135 neural
crest and 17,138 early ectoderm cells. GRNboost2 (Arboreto) was run against a curated
list of 1,417 *Xenopus* transcription factors to link TFs to candidate targets from
gene co-expression dynamics. Pax3 and TFAP2e — major predicted hubs — were then
functionally validated in vivo with antisense-morpholino RNA-seq and TF ChIP-seq.
See the paper's *Materials and Methods* for full detail.

## 🛠 Run the explorer locally

Because browsers block local file reads, serve the folder over HTTP:

```bash
git clone https://github.com/Qotov/neucrest_grn.git
cd neucrest_grn
python3 -m http.server 8000
# open http://localhost:8000
```

Pure static site — [Cytoscape.js](https://js.cytoscape.org/) for the graph, no build
step. Deployed via **GitHub Pages** from the repository root.

## 📄 Citation & links

- Paper: [PNAS 121, e2311685121 (2024)](https://www.pnas.org/doi/10.1073/pnas.2311685121)
- Lab: [Monsoro-Burq lab](https://curie.fr/equipe/monsoro-burq), Institut Curie / Université Paris-Saclay
- GRNboost2 / Arboreto: [Moerman et al., *Bioinformatics* 2019](https://doi.org/10.1093/bioinformatics/bty916)

If you use these networks, please cite the paper above.
