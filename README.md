# Single-cell Analysis of Tumor Ecosystem Heterogeneity and Tumor–Stromal Communication Across Breast Cancer Subtypes

## Overview

This project presents a complete single-cell RNA-seq (scRNA-seq) analysis pipeline investigating how breast cancer molecular subtypes differ in tumor ecosystem composition and cell-cell communication. Using a landmark 26-tumor dataset, I characterized 93,235 cells across nine major cell populations and identified subtype-specific patterns of immune infiltration, stromal programming, and intercellular signaling.

This work was developed as part of my PhD application portfolio to demonstrate proficiency in large-scale single-cell analysis, tumor microenvironment biology, and computational cancer research.

---

## Central Hypothesis

> TNBC (Triple Negative Breast Cancer) tumors harbor a more immunosuppressive microenvironment driven by CAF–tumor crosstalk, while ER+ tumors show fibroblast-dominant communication networks that facilitate tumor progression through distinct stromal signaling axes.

---

## Dataset

| Field | Details |
|---|---|
| Source | Gene Expression Omnibus |
| Accession | GSE176078 |
| Paper | Wu et al., 2021, *Nature Genetics* |
| Title | "A single-cell and spatially resolved atlas of human breast cancers" |
| Tumors | 26 primary breast tumors |
| Technology | 10x Genomics scRNA-seq |
| Subtypes | ER+, HER2+, TNBC |
| Cells after QC | 93,235 |

---

## Repository Structure

```
├── breast_cancer_scrna_pipeline.R   # Complete analysis pipeline
├── README.md                        # This file
├── figures/
│   ├── Figure1_global_UMAP.png
│   ├── Figure1b_UMAP_by_subtype.png
│   ├── Figure2_marker_genes.png
│   ├── Figure3_composition.png
│   ├── Figure4_tumor_subclusters.png
│   ├── Figure4_tumor_states_annotated.png
│   ├── Figure4b_states_by_subtype.png
│   ├── Figure5_CAF_annotated.png
│   ├── Figure5b_CAF_by_subtype.png
│   ├── Figure5c_CAF_composition.png
│   ├── Figure6_immune_annotated.png
│   ├── Figure6b_immune_by_subtype.png
│   ├── Figure6c_immune_composition.png
│   ├── Figure8_communication_heatmap.png
│   ├── Figure8b_subtype_communication.png
│   ├── Figure8c_key_LR_pairs.png
│   └── Figure9_biological_model.png
└── results/
    ├── immune_markers_all.csv
    ├── communication_results_all.csv
    ├── communication_subtype_results.csv
    └── key_LR_pairs.csv
```

---

## Analysis Pipeline

### Phase 1 — Data Loading
- Loaded all 26 processed scRNA-seq samples from GEO (GSE176078)
- Custom loading function built for the non-standard `count_matrix_*` file naming convention
- Attached author-provided metadata including subtype labels and cell type annotations

### Phase 2 — Quality Control
- Calculated mitochondrial percentage, gene counts, and UMI counts per cell
- Filtering thresholds: nFeature_RNA 200–8,000; percent.mito < 20%; nCount_RNA < 80,000
- Doublet detection with **scDblFinder** — removed 6,650 doublets (~6.6% doublet rate)
- Final dataset: **93,235 high-quality singlet cells**

### Phase 3 — Normalization and Dimensionality Reduction
- Log normalization (LogNormalize, scale factor 10,000)
- 3,000 highly variable features selected by VST
- ScaleData with mitochondrial percentage regression
- PCA (50 PCs) — elbow plot identified 30 PCs as optimal

### Phase 4 — Batch Correction and Clustering
- **Harmony** batch correction across 26 tumors — samples mix cleanly confirming successful integration
- UMAP visualization using Harmony-corrected embeddings (dims 1:30)
- Leiden clustering at resolution 0.5 — 30 global clusters resolved

### Phase 5 — Global Cell Type Annotation (Figures 1 & 2)
- Canonical marker genes used for validation: EPCAM, CD3D, CD68, CD79A, COL1A1, PECAM1, PTPRC
- Author metadata used as ground truth reference annotation
- Nine major populations identified: Cancer Epithelial, T-cells, Myeloid, B-cells, CAFs, Endothelial, Normal Epithelial, Plasmablasts, PVL

### Phase 6 — Cell Type Composition (Figure 3)
- Quantified proportions of each cell type per subtype
- Key finding: TNBC is immune-enriched; ER+ is epithelial-dominant; HER2+ shows unique plasmablast enrichment

### Phase 7 — Tumor Epithelial Subclustering (Figure 4)
- Isolated 33,000+ cancer epithelial cells and reclustered independently
- **13 tumor cell states** identified via data-driven marker discovery (Wilcoxon test, one cluster at a time)
- States include: Proliferative, EMT-like, Stem-like, ER+ Luminal, Luminal-like, Stress-response, Inflammatory-secretory, Basal-differentiated, TNBC-genomic-instability, Neuroendocrine-like
- Key finding: TNBC tumors harbor Proliferative, EMT-like, and Stem-like states; ER+ tumors are dominated by differentiated Luminal states — confirming subtype-specific tumor cell plasticity

### Phase 8 — CAF Heterogeneity (Figure 5)
- Isolated and reclustered cancer-associated fibroblasts
- **7 CAF states** identified: Myofibroblastic, Resting, Transitional, Interferon-stimulated, Inflammatory, Perivascular, ECM-remodeling
- Key finding: Myofibroblastic CAFs dominate ER+ and HER2+; Inflammatory and Interferon-stimulated CAFs are enriched in TNBC

### Phase 9 — Immune Microenvironment (Figure 6)
- Isolated 48,457 immune cells and reclustered
- **15 immune states** identified across T cells, macrophages, B cells, NK cells, dendritic cells
- Key finding: TNBC has higher M2-Macrophage and Classical-Monocyte proportions — a classic immunosuppressive myeloid signature; ER+ has more resting Naive-Memory T cells

### Phase 10 — Cell-Cell Communication (Figure 8)
- Ligand-receptor analysis using the **NicheNet human LR database** (Browaeys et al., 2020)
- 4,389 expressed LR pairs quantified across all cell type pairs
- Communication score = mean ligand expression in sender × mean receptor expression in receiver
- Key findings:
  - Myeloid cells are the central communication hub across all subtypes
  - **MIF → CD74**: dominant tumor-to-myeloid immunosuppressive signal
  - **LGALS1 → CD69**: CAF-mediated T cell suppression via Galectin-1
  - **CXCL9/10/11 → CXCR4**: inflammatory CAF chemokine signaling
  - **HLA-DRA → CD9**: myeloid antigen presentation feedback to tumor cells
  - Cancer Epithelial → Myeloid communication is stronger in TNBC than ER+

### Phase 11 — Biological Model (Figure 9)
- Integrated all findings into a proposed mechanistic schematic
- Model shows how TNBC tumor cells (via MIF→CD74) and CAFs (via LGALS1→CD69) cooperate to suppress immune activity while maintaining active stromal-immune crosstalk

---

## Key Biological Findings

1. **Subtype shapes ecosystem structure** — TNBC, ER+, and HER2+ tumors differ fundamentally in cell type composition, not just tumor cell gene expression

2. **TNBC tumor cells are transcriptionally plastic** — enriched for Proliferative, EMT-like, and Stem-like states compared to differentiated Luminal states in ER+

3. **CAF states are subtype-specific** — Inflammatory CAFs dominate TNBC while Myofibroblastic CAFs dominate ER+ and HER2+

4. **TNBC has an immunosuppressive myeloid-rich microenvironment** — higher M2-Macrophage proportion and stronger tumor-to-myeloid MIF signaling

5. **LGALS1→CD69 is the dominant CAF-T cell suppressive axis** — Galectin-1 signaling from CAFs to T cells represents a potential therapeutic target

6. **Myeloid cells are the universal communication hub** — the strongest cell-cell communication signal across all three subtypes is Myeloid → B-cell

---

## Tools and Methods

| Tool | Version | Purpose |
|---|---|---|
| Seurat | v5 | Main scRNA-seq analysis framework |
| Harmony | — | Batch correction across 26 tumors |
| scDblFinder | Bioconductor | Doublet detection and removal |
| NicheNet LR database | 2021-12-21 | Ligand-receptor communication analysis |
| ComplexHeatmap | Bioconductor | Communication strength heatmap |
| ggplot2 | — | All composition and barplots |
| ggraph / igraph | — | Biological model network visualization |
| R | 4.6.0 | Analysis environment |

---

## Computational Notes

**Hardware:** All analyses were run on a local machine with 8GB RAM running R 4.6 on Windows with WSL (Ubuntu 24.04).

**Memory management:** SCTransform was evaluated but exceeded memory limits on 8GB RAM. LogNormalize with ScaleData was used as a biologically equivalent alternative at this analysis scale. FindAllMarkers was replaced with a one-cluster-at-a-time loop to prevent memory crashes while producing identical results.

**Trajectory analysis — future direction:** Monocle3 pseudotime trajectory analysis was planned to trace tumor cell state transitions from Proliferative → EMT-like → Stem-like programs. Installation via GitHub and conda was unsuccessful due to package compatibility constraints between Monocle3, igraph 1.3+, and R 4.6. This analysis is identified as the immediate next phase of this project and will be completed on a high-performance computing cluster. The biological question — how tumor cell states transition across subtypes — remains a key open direction.

---

## How to Run

1. Download processed data from GEO accession GSE176078
2. Extract all 26 sample folders
3. Update `data_dir` path in `breast_cancer_scrna_pipeline.R`
4. Install required packages (installation lines included at top of script)
5. Run script sections sequentially — each section saves intermediate `.rds` objects so analysis can be resumed without rerunning from scratch

**Estimated runtime:** 3–5 hours on 8GB RAM; 45–90 minutes on 16GB+ RAM

---

## References

Wu SZ, Al-Eryani G, Roden DL, et al. A single-cell and spatially resolved atlas of human breast cancers. *Nature Genetics*. 2021;53(9):1334-1347. doi:10.1038/s41588-021-00911-1

Browaeys R, Saelens W, Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. *Nature Methods*. 2020;17(2):159-162.

Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*. 2019;16(12):1289-1296.

---

## Contact

Joseph Imhanbor
PhD Applicant — Computational Cancer Biology
