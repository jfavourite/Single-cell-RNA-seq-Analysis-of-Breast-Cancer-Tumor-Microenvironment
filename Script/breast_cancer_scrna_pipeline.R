# =============================================================================
# Single-cell RNA-seq Analysis of Breast Cancer Tumor Microenvironment
# =============================================================================
# Project: Tumor ecosystem heterogeneity and tumor-stromal communication
#          across breast cancer molecular subtypes
#
# Dataset: GSE176078 — Wu et al., 2021, Nature Genetics
#          26 human breast tumors (ER+, HER2+, TNBC)
#          93,235 high-quality cells after QC

# =============================================================================
# SECTION 1 — PACKAGE INSTALLATION AND LOADING
# =============================================================================

# Run install lines only once — comment out after first run
install.packages("Seurat")
install.packages("Matrix")
install.packages("tidyverse")
install.packages("harmony")
install.packages("NMF")
install.packages("circlize")
install.packages("ggraph")
install.packages("igraph")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scDblFinder")
BiocManager::install("ComplexHeatmap")
BiocManager::install("BiocNeighbors")

# Load all libraries
library(Seurat)
library(Matrix)
library(tidyverse)
library(ggplot2)
library(harmony)
library(scDblFinder)
library(ComplexHeatmap)
library(circlize)
library(ggraph)
library(igraph)


# =============================================================================
# SECTION 2 — DATA LOADING
# =============================================================================
# Dataset: GSE176078 — 26 processed scRNA-seq samples
# Each sample folder contains:
#   count_matrix_sparse.mtx  — count matrix
#   count_matrix_barcodes.tsv — cell barcodes
#   count_matrix_genes.tsv   — gene names
#   metadata.csv             — cell metadata including subtype and cell type

# Set path to data directory — update this to your local path
data_dir <- "C:/Users/User/OneDrive/Documents/scRNA files"

# List all 26 sample folders
sample_folders <- list.dirs(data_dir, recursive = FALSE)
sample_names   <- basename(sample_folders)

# Function to load one sample into a Seurat object
load_sample <- function(folder_path) {
  mat      <- readMM(file.path(folder_path, "count_matrix_sparse.mtx"))
  barcodes <- read.table(file.path(folder_path, "count_matrix_barcodes.tsv"),
                         header = FALSE, stringsAsFactors = FALSE)$V1
  genes    <- read.table(file.path(folder_path, "count_matrix_genes.tsv"),
                         header = FALSE, stringsAsFactors = FALSE)$V1
  rownames(mat) <- genes
  colnames(mat) <- barcodes
  seurat_obj <- CreateSeuratObject(counts = mat,
                                   project = basename(folder_path),
                                   min.cells = 3,
                                   min.features = 200)
  return(seurat_obj)
}

# Load all 26 samples with progress messages
seurat_list <- list()
for (i in seq_along(sample_folders)) {
  cat("Loading sample", i, "of", length(sample_folders), ":", sample_names[i], "\n")
  seurat_list[[sample_names[i]]] <- load_sample(sample_folders[i])
}


# =============================================================================
# SECTION 3 — METADATA ATTACHMENT AND MERGING
# =============================================================================
# Metadata contains: subtype, celltype_major, celltype_minor,
#                    celltype_subset, percent.mito, nCount_RNA, nFeature_RNA

for (i in seq_along(seurat_list)) {
  meta         <- read.csv(file.path(sample_folders[i], "metadata.csv"), row.names = 1)
  common_cells <- intersect(rownames(meta), colnames(seurat_list[[i]]))
  seurat_list[[i]] <- seurat_list[[i]][, common_cells]
  seurat_list[[i]] <- AddMetaData(seurat_list[[i]], meta[common_cells, ])
}

# Merge all 26 samples into one Seurat object
merged_seurat <- merge(seurat_list[[1]],
                       y       = seurat_list[-1],
                       add.cell.ids = sample_names,
                       project = "BreastCancer_TME")
merged_seurat


# =============================================================================
# SECTION 4 — QUALITY CONTROL
# =============================================================================

# Calculate mitochondrial percentage
merged_seurat[["percent.mito.check"]] <- PercentageFeatureSet(merged_seurat,
                                                               pattern = "^MT-")

# Visualize QC metrics
VlnPlot(merged_seurat,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
        ncol = 3, pt.size = 0)

# Summary statistics to guide filtering thresholds
summary(merged_seurat$nFeature_RNA)
summary(merged_seurat$nCount_RNA)
summary(merged_seurat$percent.mito)

# Filter cells:
# — nFeature_RNA > 200: remove empty droplets
# — nFeature_RNA < 8000: remove likely doublets
# — percent.mito < 20: remove dying cells
# — nCount_RNA < 80000: remove multiplets
merged_seurat <- subset(merged_seurat,
                        subset = nFeature_RNA > 200 &
                                 nFeature_RNA < 8000 &
                                 percent.mito < 20   &
                                 nCount_RNA < 80000)

# Doublet detection with scDblFinder
# Join layers before conversion (required for Seurat v5)
merged_seurat <- JoinLayers(merged_seurat)

sce                      <- as.SingleCellExperiment(merged_seurat)
sce                      <- scDblFinder(sce, samples = "orig.ident")
merged_seurat$scDblFinder.class <- sce$scDblFinder.class
merged_seurat$scDblFinder.score <- sce$scDblFinder.score

# Check doublet counts — expect ~5-7% doublet rate
table(merged_seurat$scDblFinder.class)

# Remove doublets — retained 93,235 singlet cells
merged_seurat <- subset(merged_seurat,
                        subset = scDblFinder.class == "singlet")


# =============================================================================
# SECTION 5 — NORMALIZATION AND DIMENSIONALITY REDUCTION
# =============================================================================
# Note: SCTransform was attempted but exceeded memory limits on 8GB RAM.
# LogNormalize with ScaleData is used as a memory-efficient alternative.
# Results are biologically equivalent for this analysis scale.

# Normalize
merged_seurat <- NormalizeData(merged_seurat,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000)

# Identify highly variable features
merged_seurat <- FindVariableFeatures(merged_seurat,
                                      selection.method = "vst",
                                      nfeatures = 3000)

# Scale data — regress out mitochondrial percentage
merged_seurat <- ScaleData(merged_seurat,
                           vars.to.regress = "percent.mito",
                           verbose = TRUE)

# PCA
merged_seurat <- RunPCA(merged_seurat, npcs = 50, verbose = TRUE)

# Elbow plot to determine optimal number of PCs
ElbowPlot(merged_seurat, ndims = 50)
# Elbow levels off around PC30 — use 30 PCs for downstream analysis


# =============================================================================
# SECTION 6 — BATCH CORRECTION AND CLUSTERING
# =============================================================================
# Harmony removes inter-tumor batch effects while preserving biological signal

merged_seurat <- RunHarmony(merged_seurat,
                            group.by.vars  = "orig.ident",
                            reduction      = "pca",
                            reduction.save = "harmony",
                            verbose        = TRUE)

# UMAP using Harmony-corrected embeddings
merged_seurat <- RunUMAP(merged_seurat,
                         reduction = "harmony",
                         dims      = 1:30,
                         verbose   = TRUE)

# Neighbor graph and clustering
merged_seurat <- FindNeighbors(merged_seurat,
                               reduction = "harmony",
                               dims      = 1:30,
                               verbose   = TRUE)

merged_seurat <- FindClusters(merged_seurat,
                              resolution = 0.5,
                              verbose    = TRUE)

# Set cell type identity using author annotations
Idents(merged_seurat) <- "celltype_major"

# Save merged object — critical checkpoint
saveRDS(merged_seurat, "merged_seurat.rds")


# =============================================================================
# SECTION 7 — GLOBAL CELL TYPE ANNOTATION (FIGURES 1 & 2)
# =============================================================================

# Figure 1a — Global UMAP colored by cell type
p <- DimPlot(merged_seurat,
             reduction  = "umap",
             group.by   = "celltype_major",
             label      = TRUE,
             repel      = TRUE,
             pt.size    = 0.3) +
     ggtitle("Breast Cancer Tumor Ecosystem — 93,235 cells") +
     theme(plot.title = element_text(size = 16, face = "bold"))

ggsave("Figure1_global_UMAP.png", plot = p, width = 12, height = 9, dpi = 200)

# Figure 1b — UMAP colored by breast cancer subtype
p2 <- DimPlot(merged_seurat,
              reduction = "umap",
              group.by  = "subtype",
              pt.size   = 0.3) +
      ggtitle("Cells colored by breast cancer subtype")

ggsave("Figure1b_UMAP_by_subtype.png", plot = p2, width = 12, height = 9, dpi = 200)

# Figure 2 — Canonical marker gene feature plots
# Used to validate cluster identities against known cell type markers
marker_genes <- c(
  "EPCAM",  # Tumor epithelial cells
  "CD3D",   # T cells
  "CD68",   # Macrophages
  "CD79A",  # B cells
  "COL1A1", # Fibroblasts / CAFs
  "PECAM1", # Endothelial cells
  "PTPRC"   # All immune cells (CD45)
)

p <- FeaturePlot(merged_seurat,
                 features  = marker_genes,
                 ncol      = 3,
                 reduction = "umap",
                 pt.size   = 0.1)

ggsave("Figure2_marker_genes.png", plot = p, width = 18, height = 12, dpi = 150)


# =============================================================================
# SECTION 8 — CELL TYPE COMPOSITION ACROSS SUBTYPES (FIGURE 3)
# =============================================================================

composition_data <- merged_seurat@meta.data %>%
  group_by(subtype, celltype_major) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(subtype) %>%
  mutate(proportion = count / sum(count) * 100)

p <- ggplot(composition_data,
            aes(x = subtype, y = proportion, fill = celltype_major)) +
     geom_bar(stat = "identity") +
     scale_fill_brewer(palette = "Set2") +
     labs(title = "Cell type composition across breast cancer subtypes",
          x = "Subtype", y = "Proportion (%)", fill = "Cell type") +
     theme_classic() +
     theme(plot.title  = element_text(size = 14, face = "bold"),
           axis.text   = element_text(size = 12))

ggsave("Figure3_composition.png", plot = p, width = 10, height = 7, dpi = 200)


# =============================================================================
# SECTION 9 — TUMOR EPITHELIAL SUBCLUSTERING (FIGURE 4)
# =============================================================================

# Subset cancer epithelial cells
cancer_epi <- subset(merged_seurat, subset = celltype_major == "Cancer Epithelial")

# Recluster within tumor epithelial compartment
cancer_epi <- FindVariableFeatures(cancer_epi, nfeatures = 3000)
cancer_epi <- ScaleData(cancer_epi)
cancer_epi <- RunPCA(cancer_epi, npcs = 30)
cancer_epi <- RunHarmony(cancer_epi,
                         group.by.vars  = "orig.ident",
                         reduction      = "pca",
                         reduction.save = "harmony")
cancer_epi <- RunUMAP(cancer_epi,    reduction = "harmony", dims = 1:20)
cancer_epi <- FindNeighbors(cancer_epi, reduction = "harmony", dims = 1:20)
cancer_epi <- FindClusters(cancer_epi,  resolution = 0.4)

# Figure 4a — Tumor clusters by subtype and cluster number
p1 <- DimPlot(cancer_epi, reduction = "umap",
              group.by = "subtype", pt.size = 0.3) +
      ggtitle("Tumor cells by subtype")

p2 <- DimPlot(cancer_epi, reduction = "umap",
              label = TRUE, pt.size = 0.3) +
      ggtitle("Tumor cell clusters")

p1 + p2
ggsave("Figure4_tumor_subclusters.png", width = 14, height = 6, dpi = 200)

# Marker discovery — one cluster at a time to manage memory
clusters    <- levels(Idents(cancer_epi))
all_markers <- list()

for (clust in clusters) {
  cat("Finding markers for cluster", clust, "\n")
  markers           <- FindMarkers(cancer_epi,
                                   ident.1        = clust,
                                   only.pos       = TRUE,
                                   min.pct        = 0.25,
                                   logfc.threshold = 0.5,
                                   test.use       = "wilcox")
  markers$cluster   <- clust
  markers$gene      <- rownames(markers)
  all_markers[[clust]] <- markers
}

tumor_markers <- do.call(rbind, all_markers)
top_markers   <- tumor_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print(top_markers, n = 65)

# Assign biological state labels based on marker gene interpretation
# Labels derived from top marker genes per cluster (data-driven annotation)
tumor_labels <- c(
  "0"  = "Luminal-like",              # GJA1, SERPINA5
  "1"  = "Stress-response",           # BTG2, GDF15, EGR3
  "2"  = "ER+ Luminal",               # ESR1, SCGB2A2, ANKRD30A
  "3"  = "Inflammatory-secretory",    # SAA1, SAA2, MFGE8
  "4"  = "Proliferative",             # RRM2, PLK1, DLGAP5 (cell cycle genes)
  "5"  = "Stem-like",                 # LIN28B (cancer stem cell marker)
  "6"  = "Immune-contamination",      # IGKV genes — not epithelial
  "7"  = "EMT-like",                  # SERPINE1, NCAM1, SLIT2
  "8"  = "Stromal-contamination",     # COL3A1, PDGFRB, THY1 — fibroblast markers
  "9"  = "Basal-differentiated",      # KRT4, ALPL
  "10" = "TNBC-genomic-instability",  # HORMAD1 (cancer-testis antigen)
  "11" = "EMT-like-2",                # SERPINE1, FGF5
  "12" = "Neuroendocrine-like"        # NEUROD2, PAX6 — rare population
)

cancer_epi$tumor_state <- unname(tumor_labels[as.character(Idents(cancer_epi))])

# Remove contaminating non-epithelial cells
cancer_epi_clean <- subset(cancer_epi,
                            subset = tumor_state != "Immune-contamination" &
                                     tumor_state != "Stromal-contamination")

# Figure 4b — Annotated tumor states
p <- DimPlot(cancer_epi_clean,
             reduction = "umap",
             group.by  = "tumor_state",
             label     = TRUE, repel = TRUE,
             pt.size   = 0.3) +
     ggtitle("Tumor cell states across breast cancer") +
     theme(plot.title = element_text(size = 14, face = "bold"))

ggsave("Figure4_tumor_states_annotated.png", plot = p, width = 14, height = 9, dpi = 200)

# Figure 4c — Tumor states split by subtype
p2 <- DimPlot(cancer_epi_clean,
              reduction = "umap",
              group.by  = "tumor_state",
              split.by  = "subtype",
              pt.size   = 0.3, ncol = 3) +
      ggtitle("Tumor states by subtype")

ggsave("Figure4b_states_by_subtype.png", plot = p2, width = 18, height = 6, dpi = 200)

# Save clean tumor epithelial object
saveRDS(cancer_epi_clean, "cancer_epi_clean.rds")


# =============================================================================
# SECTION 10 — CAF SUBCLUSTERING (FIGURE 5)
# =============================================================================

cafs <- subset(merged_seurat, subset = celltype_major == "CAFs")

cafs <- FindVariableFeatures(cafs, nfeatures = 2000)
cafs <- ScaleData(cafs)
cafs <- RunPCA(cafs, npcs = 20)
cafs <- RunHarmony(cafs,
                   group.by.vars  = "orig.ident",
                   reduction      = "pca",
                   reduction.save = "harmony")
cafs <- RunUMAP(cafs,     reduction = "harmony", dims = 1:15)
cafs <- FindNeighbors(cafs, reduction = "harmony", dims = 1:15)
cafs <- FindClusters(cafs,  resolution = 0.3)

# Marker discovery for CAF clusters
caf_clusters    <- levels(Idents(cafs))
caf_all_markers <- list()

for (clust in caf_clusters) {
  cat("Finding markers for CAF cluster", clust, "\n")
  markers           <- FindMarkers(cafs,
                                   ident.1        = clust,
                                   only.pos       = TRUE,
                                   min.pct        = 0.25,
                                   logfc.threshold = 0.5,
                                   test.use       = "wilcox")
  markers$cluster   <- clust
  markers$gene      <- rownames(markers)
  caf_all_markers[[clust]] <- markers
}

caf_markers     <- do.call(rbind, caf_all_markers)
top_caf_markers <- caf_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print(top_caf_markers, n = 35)

# CAF state labels derived from marker gene interpretation
caf_labels <- c(
  "0" = "Myofibroblastic-CAF",       # LRRC15, COL11A1
  "1" = "Resting-CAF",               # PI16, LEPR — universal fibroblast markers
  "2" = "Transitional-CAF",          # Weak markers — intermediate state
  "3" = "Interferon-stimulated-CAF", # IRF1, KLF4
  "4" = "Inflammatory-CAF",          # CXCL9, CXCL10, CXCL11, IDO1
  "5" = "Perivascular-CAF",          # RGS5, MCAM
  "6" = "ECM-remodeling-CAF"         # CXCL14, THBS4, SBSPON
)

cafs$caf_state <- unname(caf_labels[as.character(Idents(cafs))])

# Figure 5a — Annotated CAF states
p1 <- DimPlot(cafs, reduction = "umap",
              group.by = "caf_state",
              label = TRUE, repel = TRUE, pt.size = 0.5) +
      ggtitle("CAF states") +
      theme(plot.title = element_text(size = 14, face = "bold"))

ggsave("Figure5_CAF_annotated.png", plot = p1, width = 12, height = 8, dpi = 200)

# Figure 5b — CAF states split by subtype
p2 <- DimPlot(cafs, reduction = "umap",
              group.by = "caf_state",
              split.by = "subtype",
              pt.size  = 0.3, ncol = 3) +
      ggtitle("CAF states by subtype")

ggsave("Figure5b_CAF_by_subtype.png", plot = p2, width = 18, height = 6, dpi = 200)

# Figure 5c — CAF state composition barplot
caf_composition <- cafs@meta.data %>%
  group_by(subtype, caf_state) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(subtype) %>%
  mutate(proportion = count / sum(count) * 100)

p3 <- ggplot(caf_composition,
             aes(x = subtype, y = proportion, fill = caf_state)) +
      geom_bar(stat = "identity") +
      labs(title = "CAF state composition across subtypes",
           x = "Subtype", y = "Proportion (%)", fill = "CAF state") +
      theme_classic() +
      theme(plot.title = element_text(size = 14, face = "bold"))

ggsave("Figure5c_CAF_composition.png", plot = p3, width = 10, height = 7, dpi = 200)


# =============================================================================
# SECTION 11 — IMMUNE MICROENVIRONMENT SUBCLUSTERING (FIGURE 6)
# =============================================================================

immune <- subset(merged_seurat,
                 subset = celltype_major %in% c("T-cells", "Myeloid",
                                                 "B-cells", "Plasmablasts"))

immune <- FindVariableFeatures(immune, nfeatures = 2000)
immune <- ScaleData(immune)
immune <- RunPCA(immune, npcs = 30)
immune <- RunHarmony(immune,
                     group.by.vars  = "orig.ident",
                     reduction      = "pca",
                     reduction.save = "harmony")
immune <- RunUMAP(immune,     reduction = "harmony", dims = 1:20)
immune <- FindNeighbors(immune, reduction = "harmony", dims = 1:20)
immune <- FindClusters(immune,  resolution = 0.4)

saveRDS(immune, "immune.rds")

# Marker discovery for immune clusters
immune_clusters    <- levels(Idents(immune))
immune_all_markers <- list()

for (clust in immune_clusters) {
  cat("Finding markers for immune cluster", clust, "\n")
  markers           <- FindMarkers(immune,
                                   ident.1        = clust,
                                   only.pos       = TRUE,
                                   min.pct        = 0.25,
                                   logfc.threshold = 0.5,
                                   test.use       = "wilcox")
  markers$cluster   <- clust
  markers$gene      <- rownames(markers)
  immune_all_markers[[clust]] <- markers
}

immune_markers     <- do.call(rbind, immune_all_markers)
top_immune_markers <- immune_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print(top_immune_markers, n = 80)
write.csv(immune_markers, "immune_markers_all.csv", row.names = TRUE)

# Immune state labels derived from marker gene interpretation
immune_labels <- c(
  "0"  = "Naive-Memory-Tcell",    # IL7R, CCR7
  "1"  = "Cytotoxic-CD8-Tcell",  # CD8A, CD8B, IFNG
  "2"  = "Regulatory-Tcell",     # FOXP3, CTLA4, IL2RA
  "3"  = "M2-Macrophage",        # APOE, CD209, SEPP1
  "4"  = "NK-cell",              # GNLY, KLRD1, KLRF1
  "5"  = "B-cell",               # MS4A1, BANK1, FCRLA
  "6"  = "Classical-Monocyte",   # FCN1, VCAN, CD300E
  "7"  = "Plasma-cell",          # IGKV immunoglobulin genes
  "8"  = "Proliferating-immune", # RRM2, PLK1, DLGAP5
  "9"  = "Plasma-cell-lambda",   # IGLV immunoglobulin genes
  "10" = "Interferon-stimulated",# IFIT1, IFIT2, IFIT3, RSAD2
  "11" = "Naive-Bcell",          # LTB, TCF7, CD24
  "12" = "Mature-DC",            # CLEC9A, LAMP3, FLT3
  "13" = "Plasmacytoid-DC",      # LILRA4, CLEC4C
  "14" = "Proliferating-immune-2",# UHRF1, ANLN, DEPDC1
  "15" = "Epithelial-contamination" # KRT19, AGR2 — not immune
)

immune$immune_state <- unname(immune_labels[as.character(Idents(immune))])

# Remove epithelial contamination
immune_clean <- subset(immune,
                       subset = immune_state != "Epithelial-contamination")

# Figure 6a — Annotated immune states
p1 <- DimPlot(immune_clean, reduction = "umap",
              group.by = "immune_state",
              label = TRUE, repel = TRUE, pt.size = 0.3) +
      ggtitle("Immune cell states") +
      theme(plot.title = element_text(size = 14, face = "bold"))

ggsave("Figure6_immune_annotated.png", plot = p1, width = 14, height = 9, dpi = 200)

# Figure 6b — Immune states split by subtype
p2 <- DimPlot(immune_clean, reduction = "umap",
              group.by = "immune_state",
              split.by = "subtype",
              pt.size  = 0.3, ncol = 3) +
      ggtitle("Immune states by subtype")

ggsave("Figure6b_immune_by_subtype.png", plot = p2, width = 18, height = 6, dpi = 200)

# Figure 6c — Immune composition barplot
immune_composition <- immune_clean@meta.data %>%
  group_by(subtype, immune_state) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(subtype) %>%
  mutate(proportion = count / sum(count) * 100)

p3 <- ggplot(immune_composition,
             aes(x = subtype, y = proportion, fill = immune_state)) +
      geom_bar(stat = "identity") +
      labs(title = "Immune composition across subtypes",
           x = "Subtype", y = "Proportion (%)", fill = "Immune state") +
      theme_classic() +
      theme(plot.title = element_text(size = 14, face = "bold"))

ggsave("Figure6c_immune_composition.png", plot = p3, width = 12, height = 7, dpi = 200)

saveRDS(immune_clean, "immune_clean.rds")


# =============================================================================
# SECTION 12 — CELL-CELL COMMUNICATION ANALYSIS 
# =============================================================================
# Method: NicheNet ligand-receptor database (Browaeys et al., 2020)
# Communication score = mean ligand expression in sender x
#                       mean receptor expression in receiver
# Database: lr_network_human_21122021 (4,986 LR pairs, 4,389 expressed in data)

# Load NicheNet ligand-receptor database
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))

# Get average expression per cell type
Idents(merged_seurat) <- "celltype_major"
avg_expr <- AverageExpression(merged_seurat,
                              group.by = "celltype_major",
                              layer    = "data")$RNA

# Filter to expressed LR pairs only
expressed_genes <- rownames(merged_seurat)
lr_expressed    <- lr_network %>%
  filter(from %in% expressed_genes & to %in% expressed_genes)

cat("Total LR pairs in database:", nrow(lr_network), "\n")
cat("Expressed LR pairs in data:", nrow(lr_expressed), "\n")

cell_types <- colnames(avg_expr)

# Calculate pairwise communication scores
results <- list()
for (sender in cell_types) {
  for (receiver in cell_types) {
    if (sender != receiver) {
      ligand_expr   <- avg_expr[lr_expressed$from, sender]
      receptor_expr <- avg_expr[lr_expressed$to,   receiver]
      scores        <- ligand_expr * receptor_expr
      top_idx       <- order(scores, decreasing = TRUE)[1:20]
      results[[paste(sender, receiver, sep = "->")]] <- data.frame(
        sender   = sender,
        receiver = receiver,
        ligand   = lr_expressed$from[top_idx],
        receptor = lr_expressed$to[top_idx],
        score    = scores[top_idx]
      )
    }
  }
}

comm_results <- do.call(rbind, results)
rownames(comm_results) <- NULL

# Figure 8a — Global communication heatmap
comm_summary <- comm_results %>%
  group_by(sender, receiver) %>%
  summarise(total_score = sum(score), .groups = "drop")

comm_matrix <- comm_summary %>%
  pivot_wider(names_from  = receiver,
              values_from = total_score,
              values_fill = 0) %>%
  column_to_rownames("sender")

col_fun <- colorRamp2(c(0,
                         max(as.matrix(comm_matrix)) / 2,
                         max(as.matrix(comm_matrix))),
                       c("white", "orange", "red"))

png("Figure8_communication_heatmap.png", width = 2400, height = 2000, res = 200)
Heatmap(as.matrix(comm_matrix),
        name             = "Communication\nstrength",
        col              = col_fun,
        row_title        = "Sender",
        column_title     = "Receiver",
        column_title_gp  = gpar(fontsize = 14, fontface = "bold"),
        row_title_gp     = gpar(fontsize = 14, fontface = "bold"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", as.matrix(comm_matrix)[i, j]),
                    x, y, gp = gpar(fontsize = 8))
        })
dev.off()

# Figure 8b — Top 20 interactions per subtype
subtypes     <- c("ER+", "HER2+", "TNBC")
subtype_comm <- list()

for (st in subtypes) {
  sub_obj    <- subset(merged_seurat, subset = subtype == st)
  avg_sub    <- AverageExpression(sub_obj,
                                  group.by = "celltype_major",
                                  layer    = "data")$RNA
  ct_present <- colnames(avg_sub)
  lr_sub     <- lr_expressed %>%
    filter(from %in% rownames(avg_sub) & to %in% rownames(avg_sub))

  sub_results <- list()
  for (sender in ct_present) {
    for (receiver in ct_present) {
      if (sender != receiver) {
        ligand_expr   <- avg_sub[lr_sub$from, sender]
        receptor_expr <- avg_sub[lr_sub$to,   receiver]
        scores        <- ligand_expr * receptor_expr
        total         <- sum(scores, na.rm = TRUE)
        sub_results[[paste(sender, receiver, sep = "->")]] <- data.frame(
          subtype     = st,
          sender      = sender,
          receiver    = receiver,
          total_score = total
        )
      }
    }
  }
  subtype_comm[[st]] <- do.call(rbind, sub_results)
}

all_subtype_comm <- do.call(rbind, subtype_comm)
rownames(all_subtype_comm) <- NULL

top_interactions <- all_subtype_comm %>%
  group_by(subtype) %>%
  top_n(n = 20, wt = total_score) %>%
  mutate(interaction = paste(sender, "->", receiver))

p <- ggplot(top_interactions,
            aes(x = reorder(interaction, total_score),
                y = total_score, fill = subtype)) +
     geom_bar(stat = "identity") +
     facet_wrap(~subtype, scales = "free_y", ncol = 1) +
     coord_flip() +
     labs(title = "Top 20 cell-cell interactions per subtype",
          x = "Interaction", y = "Communication strength") +
     theme_classic() +
     theme(plot.title  = element_text(size = 14, face = "bold"),
           strip.text  = element_text(size = 12, face = "bold"))

ggsave("Figure8b_subtype_communication.png", plot = p, width = 12, height = 14, dpi = 200)

# Figure 8c — Key ligand-receptor pairs in tumor-immune crosstalk
key_interactions <- comm_results %>%
  filter(
    (sender == "Cancer Epithelial" & receiver == "Myeloid") |
    (sender == "CAFs"              & receiver == "T-cells") |
    (sender == "Myeloid"           & receiver == "Cancer Epithelial")
  ) %>%
  group_by(sender, receiver) %>%
  top_n(n = 10, wt = score) %>%
  mutate(interaction = paste(ligand, "->", receptor))

p <- ggplot(key_interactions,
            aes(x    = reorder(interaction, score),
                y    = score,
                fill = paste(sender, "->", receiver))) +
     geom_bar(stat = "identity") +
     facet_wrap(~paste(sender, "->", receiver), scales = "free", ncol = 1) +
     coord_flip() +
     labs(title = "Key ligand-receptor pairs in tumor-immune crosstalk",
          x = "Ligand -> Receptor", y = "Communication score") +
     theme_classic() +
     theme(plot.title      = element_text(size = 14, face = "bold"),
           legend.position = "none",
           strip.text      = element_text(size = 11, face = "bold"))

ggsave("Figure8c_key_LR_pairs.png", plot = p, width = 12, height = 12, dpi = 200)

# Save communication results
write.csv(comm_results,      "communication_results_all.csv",      row.names = FALSE)
write.csv(all_subtype_comm,  "communication_subtype_results.csv",  row.names = FALSE)
write.csv(key_interactions,  "key_LR_pairs.csv",                   row.names = FALSE)


# =============================================================================
# SECTION 13 — BIOLOGICAL MODEL SCHEMATIC 
# =============================================================================
# Summarizes key signaling axes identified from communication analysis:
# — MIF → CD74: tumor cells polarize myeloid cells toward M2/immunosuppressive state
# — LGALS1 → CD69: CAFs suppress T cell activation via Galectin-1
# — HLA-DRA → CD9: myeloid antigen presentation feedback to tumor cells
# — COL1A1 → CD44: CAF ECM signals to T cells
# — S100A8 → CD68: tumor-derived alarmins recruit macrophages

nodes <- data.frame(
  name = c("Cancer\nEpithelial", "CAFs", "Myeloid\n(M2)",
           "T-cells\n(suppressed)", "NK cells"),
  type = c("Tumor", "Stromal", "Immune", "Immune", "Immune"),
  x    = c(0, -2, 2, -1.5, 1.5),
  y    = c(0,  1, 1, -1.5, -1.5)
)

edges <- data.frame(
  from     = c("Cancer\nEpithelial", "CAFs", "CAFs",
               "Cancer\nEpithelial", "Myeloid\n(M2)", "Myeloid\n(M2)"),
  to       = c("Myeloid\n(M2)", "Myeloid\n(M2)", "T-cells\n(suppressed)",
               "T-cells\n(suppressed)", "Cancer\nEpithelial", "T-cells\n(suppressed)"),
  label    = c("MIF→CD74", "LGALS1→CD69",
               "COL1A1→CD44", "S100A8→CD68",
               "HLA-DRA→CD9", "Immune suppression"),
  strength = c(3500, 480, 320, 200, 450, 300)
)

g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

set.seed(42)
p <- ggraph(g, layout = data.frame(x = nodes$x, y = nodes$y)) +
     geom_edge_arc(aes(label = label, width = strength),
                   arrow       = arrow(length = unit(4, "mm"), type = "closed"),
                   start_cap   = circle(8, "mm"),
                   end_cap     = circle(8, "mm"),
                   label_size  = 3,
                   angle_calc  = "along",
                   label_dodge = unit(3, "mm"),
                   color       = "gray40",
                   strength    = 0.2) +
     geom_node_point(aes(color = type), size = 18) +
     geom_node_text(aes(label = name), size = 3.5,
                    fontface = "bold", color = "white") +
     scale_edge_width(range = c(0.5, 2.5)) +
     scale_color_manual(values = c("Tumor"   = "#E64B35",
                                   "Stromal" = "#F5A623",
                                   "Immune"  = "#4DBBD5")) +
     labs(title    = "Proposed biological model — breast cancer tumor ecosystem",
          subtitle = "Key signaling axes driving subtype-specific microenvironment",
          color    = "Cell type") +
     theme_void() +
     theme(plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
           plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "gray40"),
           legend.position  = "bottom",
           plot.background  = element_rect(fill = "white", color = NA),
           panel.background = element_rect(fill = "white", color = NA))

ggsave("Figure9_biological_model.png", plot = p, width = 12, height = 10, dpi = 200)


# =============================================================================
# SECTION 14 — SESSION INFO
# =============================================================================

sessionInfo()

# =============================================================================
# THE END
# =============================================================================
