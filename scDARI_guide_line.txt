########################################################################################################
########################################################################################################
########################################################################################################
#####                                                                                              #####
#####     A Guide to Constructing eGRNs Using scDERI with Paired scRNA-seq and scATAC-seq Data     #####
#####                                                                                              #####
########################################################################################################
########################################################################################################
########################################################################################################

### Use R v4.4

### Environment Setup
Sys.setenv(LANGUAGE = "en") # Set system output language to English.
options(stringsAsFactors = FALSE) # Disable automatic string-to-factor conversion.
rm(list = ls()) # Clear all objects from the global environment.
gc() # Free unused memory.

### Path Setup
# Working directory for all subsequent code execution.
Working_directory = "/mnt/i/scDERI" # WSL2 path (recommended)
# Working_directory = "I:/scDERI" # Windows path (not recommended)
# Working_directory = "/home/scDERI" # Linux path (available)

setwd(Working_directory)





##############################################################################
##############################################################################
###                                                                        ###
###   Preparations: Installing R Packages and Setting Up the Environment   ###
###                                                                        ###
##############################################################################
##############################################################################
### Environment Configuration
# Most packages are compatible with 2024–2026 releases. Listed versions are for reference only.


### Install core dependency packages.
# Install igraph 2.0.3 and dplyr 1.1.4 from the CRAN archive.
if (!requireNamespace("igraph", quietly = TRUE) || packageVersion("igraph") != "2.0.3") {
  install.packages("https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz", repos = NULL, type = "source")
}
packageVersion("igraph") # Version: 2.0.3
if (!requireNamespace("dplyr", quietly = TRUE) || packageVersion("dplyr") != "1.1.4") {
  install.packages("https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_1.1.4.tar.gz", repos = NULL, type = "source")
}
packageVersion("dplyr") # Version: 1.1.4

# Install CRAN packages.
install.packages(c("ggplot2", "Signac", "SeuratObject", "Seurat", "VGAM", "ggrepel", "minpack.lm",
                   "data.table", "tibble", "tidyr", "pROC", "DDRTree", "GA", "psych"))
packageVersion("ggplot2") # Version: 4.0.0
packageVersion("Signac") # Version: 1.15.0
packageVersion("SeuratObject") # Version: 5.2.0
packageVersion("Seurat") # Version: 5.3.0
packageVersion("VGAM") # Version: 1.1.13
packageVersion("ggrepel") # Version: 0.9.6
packageVersion("minpack.lm") # Version: 1.2.4
packageVersion("data.table") # Version: 1.17.8
packageVersion("tibble") # Version: 3.3.0
packageVersion("tidyr") # Version: 1.3.1
packageVersion("pROC") # Version: 1.19.0.1
packageVersion("DDRTree") # Version: 0.1.5
packageVersion("GA") # Version: 3.2.4
packageVersion("psych") # Version: 2.5.6

# Install Bioconductor packages.
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomeInfoDb", "TFBSTools", "JASPAR2024", "IRanges", "ChIPseeker",
                       "Biobase", "BiocGenerics", "monocle", "Biostrings"))
packageVersion("GenomeInfoDb") # Version: 1.42.3
packageVersion("TFBSTools") # Version: 1.44.0
packageVersion("JASPAR2024") # Version: 0.99.6
packageVersion("IRanges") # Version: 2.40.1
packageVersion("ChIPseeker") # Version: 1.42.1
packageVersion("Biobase") # Version: 2.66.0
packageVersion("BiocGenerics") # Version: 0.52.0
packageVersion("monocle") # Version: 2.34.0
packageVersion("Biostrings") # Version: 2.74.1


### Install genome annotation packages corresponding to the analysis data.
# Install transcript database (TxDb).
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene") # hg38 (GRCh38, 2013)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene") # hg19 (GRCh37, 2010)
packageVersion("TxDb.Hsapiens.UCSC.hg38.knownGene") # Version: 3.20.0
packageVersion("TxDb.Hsapiens.UCSC.hg19.knownGene") # Version: 3.2.2

# Install Ensembl genome annotation package.
BiocManager::install("EnsDb.Hsapiens.v86") # hg38 (GRCh38, 2013)
# BiocManager::install("EnsDb.Hsapiens.v75") # hg19 (GRCh37, 2010)
packageVersion("EnsDb.Hsapiens.v86") # Version: 2.99.0
packageVersion("EnsDb.Hsapiens.v75") # Version: 2.99.0

# Install BSgenome reference genome package.
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38") # hg38 (GRCh38, 2013)
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19") # hg19 (GRCh37, 2010)
packageVersion("BSgenome.Hsapiens.UCSC.hg38") # Version: 1.4.5
packageVersion("BSgenome.Hsapiens.UCSC.hg19") # Version: 1.4.3


### Install scDERI package from GitHub repository.
devtools::install_github("hth-buaa/scDERI")
packageVersion("scDERI") # Version: 2.0


### Load required packages.
library(scDERI)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v86)
# library(EnsDb.Hsapiens.v75)
library(dplyr)
library(DDRTree)
library(BSgenome.Hsapiens.UCSC.hg38)
# library(BSgenome.Hsapiens.UCSC.hg19)


### Usage Notes
# Memory management notes:
# - scDERI uses parallel computing; memory may not be fully released after parallel steps.
# - Restart the R session to free memory if needed.
# - All scDERI functions are modular and independent; save intermediate results for modular execution.
# This tutorial demonstrates a modular workflow, loading required data and packages before each step.





################################
################################
###                          ###
###   Part 0: Loading Data   ###
###                          ###
################################
################################




########################################################
#                                                      #
# 0.1 Loading the Seurat Object of the Example Dataset #
#                                                      #
########################################################
setwd(Working_directory)




##############################################################
# 0.1.1 Loading the Raw Seurat Object of the Example Dataset #
##############################################################
### Download data in R.
# Install SeuratData (development version 0.2.2.9002 used in this tutorial).
devtools::install_github('satijalab/seurat-data')
packageVersion("SeuratData") # Version: 0.2.2.9002

options(timeout = 999999999) # Extend download timeout for large data files.
SeuratData::InstallData("pbmcMultiome")
scRNAseq = SeuratData::LoadData("pbmcMultiome", "pbmc.rna")
scATACseq = SeuratData::LoadData("pbmcMultiome", "pbmc.atac")


### Alternative: manual download and local installation (used in this tutorial).
# Use manual download if automated download fails repeatedly.
# Download pbmcMultiome.SeuratData_0.1.4.tar.gz to the working directory first.
# Data source: http://seurat.nygenome.org/src/contrib/pbmcMultiome.SeuratData_0.1.4.tar.gz (267.5 MB)
install.packages("./pbmcMultiome.SeuratData_0.1.4.tar.gz", repos = NULL, type = "source")
packageVersion("pbmcMultiome.SeuratData") # Version: 0.1.4
scRNAseq = pbmcMultiome.SeuratData::pbmc.rna
scATACseq = pbmcMultiome.SeuratData::pbmc.atac

# Inspect scRNA-seq Seurat object.
base::print(scRNAseq)
# Expected output: 36601 features × 11909 cells, RNA assay with counts/data layers.

# Inspect scATAC-seq Seurat object.
base::print(scATACseq)
# Expected output: 108377 features × 11909 cells, ATAC assay with counts/data layers.

# Verify matched cell barcodes between scRNA-seq and scATAC-seq data.
base::print(base::table(base::sort(base::colnames(scRNAseq)) == base::sort(base::colnames(scATACseq))))
# Expected output: 11909 TRUE (all cells match).


### For unmatched multi-omics data:
# - Match cells across modalities and discard unmatched cells.
# - Unify cell identifiers before proceeding.
# Cell name unification code is not provided in this tutorial.



##############################################################
# 0.1.2 Upgrading the Seurat v4 Object to a Seurat v5 Object #
##############################################################
### Convert Seurat v4 Assay to v5 Assay5 format for scDERI compatibility.
if (class(scRNAseq[["RNA"]]) == "Assay") {
  scRNAseq[["RNA"]] = SeuratObject::CreateAssay5Object(counts = scRNAseq[["RNA"]]@counts)
}
base::print(class(scRNAseq[["RNA"]])) # Expected output: Assay5 (SeuratObject)
base::print(class(scATACseq[["ATAC"]])) # Expected output: ChromatinAssay (Signac)



#################################################
# 0.1.3 Removing Data from Abnormal Chromosomes #
#################################################
### Inspect chromosome names in ATAC peak data.
base::print(utils::head(base::rownames(scATACseq)))
# ↓↓↓ Expected output ↓↓↓ #
# [1] "chr1-10109-10357"   "chr1-180730-181630" "chr1-191491-191736" "chr1-267816-268196" "chr1-586028-586373" "chr1-629721-630172"
# ↑↑↑ Expected output ↑↑↑ #

base::print(base::as.data.frame(base::table(base::sapply(base::strsplit(base::rownames(scATACseq), "-"), `[`, 1))))
# ↓↓↓ Expected output ↓↓↓ #
#          Var1  Freq
# 1  GL000194.1     4
# 2  GL000195.1     2
# 3  GL000205.2     6
# 4  GL000219.1    10
# 5  KI270713.1     6
# 6  KI270721.1     1
# 7  KI270726.1     2
# 8  KI270727.1     1
# 9  KI270734.1     1
# 10       chr1 10538
# 11      chr10  4967
# 12      chr11  5249
# 13      chr12  5670
# 14      chr13  2410
# 15      chr14  3734
# 16      chr15  3513
# 17      chr16  4033
# 18      chr17  5766
# 19      chr18  2013
# 20      chr19  4881
# 21       chr2  8556
# 22      chr20  3153
# 23      chr21  1370
# 24      chr22  2445
# 25       chr3  6808
# 26       chr4  4268
# 27       chr5  5332
# 28       chr6  6652
# 29       chr7  5281
# 30       chr8  4416
# 31       chr9  4414
# 32       chrX  2872
# 33       chrY     3
# ↑↑↑ Expected output ↑↑↑ #


### Retain only standard chromosomes; remove non-canonical contigs.
peaks_keep = as.character(Signac::seqnames(Signac::granges(scATACseq))) %in%
  GenomeInfoDb::standardChromosomes(Signac::granges(scATACseq))
base::print(base::table(peaks_keep))
# Expected output: 33 non-standard peaks removed, 108344 peaks retained.

atac_assay = scATACseq[["ATAC"]]
atac_assay@counts = atac_assay@counts[peaks_keep, ]
atac_assay@data = atac_assay@data[peaks_keep, ]
if (base::nrow(atac_assay@scale.data) > 0) {atac_assay@scale.data = atac_assay@scale.data[peaks_keep, ]}
atac_assay@ranges = atac_assay@ranges[peaks_keep]
atac_assay@meta.features = atac_assay@meta.features[peaks_keep, , drop = FALSE]
if (length(atac_assay@bias) > 0) {atac_assay@bias = atac_assay@bias[peaks_keep]}
if (length(atac_assay@motifs) > 0) {atac_assay@motifs@data = atac_assay@motifs@data[peaks_keep, ]}
scATACseq[["ATAC"]] = atac_assay; rm(peaks_keep, atac_assay)



############################
# 0.1.4 Adding Cell Labels #
############################
### Requirement: identical cell names and matching `seurat_annotations` metadata column in both Seurat objects.
# The `seurat_annotations` column stores cell type labels for each cell.

base::print(base::as.data.frame(base::table(scRNAseq@meta.data$seurat_annotations)))
base::print(base::as.data.frame(base::table(scATACseq@meta.data$seurat_annotations)))
# ↓↓↓ Expected output ↓↓↓ #
#              Var1 Freq
# 1       CD14 Mono 2812
# 2       CD16 Mono  514
# 3       CD4 Naive 1419
# 4         CD4 TCM 1149
# 5         CD4 TEM  298
# 6       CD8 Naive 1410
# 7       CD8 TEM_1  325
# 8       CD8 TEM_2  358
# 9            HSPC   26
# 10 Intermediate B  353
# 11           MAIT  137
# 12       Memory B  371
# 13             NK  468
# 14        Naive B  142
# 15         Plasma   18
# 16           Treg  162
# 17            cDC  198
# 18       filtered 1497
# 19            gdT  146
# 20            pDC  106
# ↑↑↑ Expected output ↑↑↑ #


### If `seurat_annotations` is missing or does not contain cell labels:
# - Manually add cell type labels to the `seurat_annotations` metadata column (character vector).



############################
# 0.1.5 Sampling the Cells #
############################
### Subset cells to retain only cell types of interest.

# Example differentiation pathways inferred from cell labels:
# 1. Myeloid lineage: monocyte differentiation (CD14 Mono → CD16 Mono), dendritic cell subsets (cDC, pDC)
# 2. Lymphoid lineage: B cell (Naive B → Intermediate B → Plasma/Memory B) and T cell subsets (CD4+/CD8+ memory, Treg)


# Example: retain CD14 Mono and CD16 Mono for downstream analysis.
BiocManager::install("S4Vectors")
packageVersion("S4Vectors") # Version: 0.44.0
scRNAseq = S4Vectors::subset(scRNAseq, subset = seurat_annotations %in% c("CD14 Mono", "CD16 Mono"))

cells_keep = scATACseq$seurat_annotations %in% c("CD14 Mono", "CD16 Mono")
atac_assay = scATACseq[["ATAC"]]
atac_assay@counts = atac_assay@counts[, cells_keep]
atac_assay@data   = atac_assay@data[, cells_keep]
if (base::nrow(atac_assay@scale.data) > 0) {atac_assay@scale.data = atac_assay@scale.data[, cells_keep]}
scATACseq@meta.data = scATACseq@meta.data[cells_keep, , drop = FALSE]
scATACseq@active.ident = scATACseq@active.ident[cells_keep]
if (length(scATACseq@reductions) > 0) {
  for (red_name in base::names(scATACseq@reductions)) {
    scATACseq@reductions[[red_name]]@cell.embeddings = scATACseq@reductions[[red_name]]@cell.embeddings[cells_keep, , drop = FALSE]
  }
}
if (length(scATACseq@graphs) > 0) {
  for (graph_name in base::names(scATACseq@graphs)) {
    scATACseq@graphs[[graph_name]] = scATACseq@graphs[[graph_name]][cells_keep, cells_keep]
  }
}
scATACseq[["ATAC"]] = atac_assay; rm(cells_keep, atac_assay)

# Verify matched cells after subsetting.
base::print(base::table(base::sort(base::colnames(scRNAseq)) == base::sort(base::colnames(scATACseq))))
# Expected output: 3326 matched cells.

# Store original labels, then set `seurat_annotations` to the unified group name.
scRNAseq@meta.data$original_seurat_annotations = scRNAseq@meta.data$seurat_annotations
scATACseq@meta.data$original_seurat_annotations = scATACseq@meta.data$seurat_annotations
scRNAseq@meta.data$seurat_annotations = "Mono"
scATACseq@meta.data$seurat_annotations = "Mono"



######################################
# 0.1.6 Downsampling (Optional Step) #
######################################
### Usage Notes:
# Optional downsampling to reduce runtime.
# Downsampling speeds up computation but may lose weak biological signals.
# This tutorial applies downsampling for workflow demonstration.


### Cell downsampling strategy: retain top 3000 cells ranked by combined non-zero frequency (RNA + ATAC).
# Note: excessive cell downsampling reduces metacell groups and gene coverage, impairing downstream analysis.
rna_counts = SeuratObject::GetAssayData(scRNAseq, assay = "RNA", layer = "counts")
rna_nonzero_freq = Matrix::colSums(rna_counts > 0) / base::nrow(rna_counts)
base::print(utils::head(rna_nonzero_freq))
# ↓↓↓ Expected output ↓↓↓ #
# AAACAGCCATCCAGGT-1 AAACCAACACAATGCC-1 AAACCAACAGGAACTG-1 AAACCAACATAATCCG-1 AAACCAACATTGTGCA-1 AAACCGAAGCTGGACC-1
#         0.08272998         0.06691074         0.06650092         0.09196470         0.06576323         0.07546242
# ↑↑↑ Expected output ↑↑↑ #

atac_counts = SeuratObject::GetAssayData(scATACseq, assay = "ATAC", layer = "counts")
atac_nonzero_freq = Matrix::colSums(atac_counts > 0) / base::nrow(atac_counts)
base::print(utils::head(atac_nonzero_freq))
# ↓↓↓ Expected output ↓↓↓ #
# AAACAGCCATCCAGGT-1 AAACCAACACAATGCC-1 AAACCAACAGGAACTG-1 AAACCAACATAATCCG-1 AAACCAACATTGTGCA-1 AAACCGAAGCTGGACC-1
#         0.07938603         0.06072325         0.08798826         0.12506461         0.07789079         0.07502031
# ↑↑↑ Expected output ↑↑↑ #

# Verify cell barcode order match.
base::print(base::table(base::names(rna_nonzero_freq) == base::names(atac_nonzero_freq)))
# Expected output: all 3326 cells match.

combined_freq = rna_nonzero_freq + atac_nonzero_freq
top_cells = base::names(combined_freq)[base::order(combined_freq, decreasing = TRUE)][1:3000]

scRNAseq = S4Vectors::subset(scRNAseq, cells = top_cells)

atac_assay = scATACseq[["ATAC"]]
atac_assay@counts = atac_assay@counts[, top_cells]
atac_assay@data   = atac_assay@data[, top_cells]
if (base::nrow(atac_assay@scale.data) > 0) {atac_assay@scale.data = atac_assay@scale.data[, top_cells]}
scATACseq@meta.data = scATACseq@meta.data[top_cells, , drop = FALSE]
scATACseq@active.ident = scATACseq@active.ident[top_cells]
if (length(scATACseq@reductions) > 0) {
  for (red_name in base::names(scATACseq@reductions)) {
    scATACseq@reductions[[red_name]]@cell.embeddings = scATACseq@reductions[[red_name]]@cell.embeddings[top_cells, , drop = FALSE]
  }
}
if (length(scATACseq@graphs) > 0) {
  for (graph_name in base::names(scATACseq@graphs)) {
    scATACseq@graphs[[graph_name]] = scATACseq@graphs[[graph_name]][top_cells, top_cells]
  }
}
scATACseq[["ATAC"]] = atac_assay; rm(rna_counts, rna_nonzero_freq, atac_counts, atac_nonzero_freq, combined_freq, top_cells, atac_assay)

# Check cell type distribution after downsampling.
base::print(base::table(scRNAseq@meta.data[["original_seurat_annotations"]]))
base::print(base::table(scATACseq@meta.data[["original_seurat_annotations"]]))
# Expected output: ~2501 CD14 Mono, ~499 CD16 Mono.


### Gene downsampling strategy: filter non-coding RNAs and mitochondrial genes, then keep top 10,000 highly variable genes.
genes = base::setdiff(base::rownames(scRNAseq), base::grep("(^AC|^AL|^LINC)[0-9]+", base::rownames(scRNAseq), value = TRUE))
genes = base::grep("MT-", genes, invert = TRUE, value = TRUE)
scRNAseq = S4Vectors::subset(scRNAseq, features = genes)

scRNAseq = Seurat::NormalizeData(scRNAseq)
scRNAseq = Seurat::FindVariableFeatures(scRNAseq, selection.method = "vst", nfeatures = base::nrow(scRNAseq))

RNA_var_meta = scRNAseq@assays[["RNA"]]@meta.data
base::rownames(RNA_var_meta) = base::rownames(scRNAseq)
RNA_var_meta = RNA_var_meta[base::order(RNA_var_meta$vf_vst_counts_variance.standardized, decreasing = TRUE),]
genes = base::rownames(RNA_var_meta)[1:10000]
scRNAseq = S4Vectors::subset(scRNAseq, features = genes)
rm(RNA_var_meta, genes)


### Peak downsampling strategy: keep top 10,000 most variable chromatin peaks.
scATACseq = Signac::RunTFIDF(scATACseq)
scATACseq = Signac::FindTopFeatures(scATACseq, min.cutoff = "q0")

ATAC_var_meta = scATACseq@assays[["ATAC"]]@meta.features
ATAC_var_meta = ATAC_var_meta[base::order(ATAC_var_meta$percentile, decreasing = TRUE), ]
peaks = base::rownames(ATAC_var_meta)[1:10000]

peaks_keep = base::rownames(scATACseq) %in% peaks
base::print(base::table(peaks_keep))
# Expected output: 98344 peaks removed, 10000 peaks retained.

atac_assay = scATACseq[["ATAC"]]
atac_assay@counts = atac_assay@counts[peaks_keep, ]
atac_assay@data = atac_assay@data[peaks_keep, ]
if (base::nrow(atac_assay@scale.data) > 0) {atac_assay@scale.data = atac_assay@scale.data[peaks_keep, ]}
atac_assay@ranges = atac_assay@ranges[peaks_keep]
atac_assay@meta.features = atac_assay@meta.features[peaks_keep, , drop = FALSE]
if (length(atac_assay@bias) > 0) {atac_assay@bias = atac_assay@bias[peaks_keep]}
if (length(atac_assay@motifs) > 0) {atac_assay@motifs@data = atac_assay@motifs@data[peaks_keep, ]}
atac_assay@var.features = character(0)
scATACseq[["ATAC"]] = atac_assay; rm(peaks, atac_assay, ATAC_var_meta)



####################################
# 0.1.7 Saving the Example Dataset #
####################################
# Save processed Seurat objects.
# Default I/O: base R saveRDS / readRDS.
# saveRDS(scRNAseq, file = "./scRNAseq.rds")
# saveRDS(scATACseq, file = "./scATACseq.rds")

# Alternative: use qs package for faster I/O with large datasets.
install.packages("qs")
packageVersion("qs") # Version: 0.27.3
qs::qsave(scRNAseq, file = "./scRNAseq.qsave")
qs::qsave(scATACseq, file = "./scATACseq.qsave")




################################################
#                                              #
# 0.2 Downloading the scATAC-seq Fragment File #
#                                              #
################################################
### Fragment file (.tsv.gz) and tabix index (.tbi) are required for Signac-based scATAC-seq analysis.
# Choose workflow based on available files:
# - Both .tsv.gz and .tbi exist: skip this section.
# - Only .tsv.gz exists: generate the .tbi index.
# - Neither exists: generate both files from the count matrix.
# Generated fragment files are comparable to official source files.



##################################
# 0.2.1 Obtain _fragments.tsv.gz #
##################################
### Data import
setwd(Working_directory)
scATACseq = qs::qread("./scATACseq.qsave")


### Option 1: download official fragment file
# File URL is stored in the Seurat object's fragment slot.
base::print(scATACseq@assays[["ATAC"]]@fragments[[1]]@path)
# File URL: https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
# Download this file to the working directory.
# This tutorial uses the the fragment file above.


### Option 2: generate fragment file from count matrix
# Use this method when no official fragment file is provided.

# Extract count matrix from ATAC assay.
count_matrix = scATACseq@assays[["ATAC"]]@counts
base::print(class(count_matrix)) # Expected output: dgCMatrix

fragments = base::rownames(count_matrix)
base::print(utils::head(fragments))
# ↓↓↓ Expected output ↓↓↓ #
# [1] "chr1-777634-779926"   "chr1-958518-963388"   "chr1-994735-1006799"  "chr1-1012606-1015854" "chr1-1016432-1028059" "chr1-1031406-1044350"
# ↑↑↑ Expected output ↑↑↑ #

cells = base::colnames(count_matrix)
base::print(utils::head(cells))
# ↓↓↓ Expected output ↓↓↓ #
# [1] "GCAATGTTCTGTTCAT-1" "CATAGACTCCCTGGAA-1" "GTGCATTAGTTAGACC-1" "TTACCTGTCTCCTCTT-1" "GACTATTCAAACTGTT-1" "ACTAACTCAGGAACCA-1"
# ↑↑↑ Expected output ↑↑↑ #

# Peak name format: chr-start-end (hyphen-separated).
# Extract chromosome, start, and end coordinates.
# Adjust parsing if your data uses a different peak naming format.
fragments_split = base::strsplit(fragments, "-")
chromosome = base::sapply(fragments_split, function(x) x[1])
start = as.integer(base::sapply(fragments_split, function(x) x[2]))
end = as.integer(base::sapply(fragments_split, function(x) x[3]))
base::print(base::as.data.frame(base::table(chromosome)))
# ↓↓↓ Expected output ↓↓↓ #
#    chromosome Freq
# 1        chr1  992
# 2       chr10  425
# 3       chr11  546
# 4       chr12  544
# 5       chr13  188
# 6       chr14  332
# 7       chr15  318
# 8       chr16  477
# 9       chr17  668
# 10      chr18  142
# 11      chr19  725
# 12       chr2  665
# 13      chr20  318
# 14      chr21  117
# 15      chr22  270
# 16       chr3  556
# 17       chr4  322
# 18       chr5  423
# 19       chr6  573
# 20       chr7  431
# 21       chr8  351
# 22       chr9  423
# 23       chrX  194
# ↑↑↑ Expected output ↑↑↑ #

# Convert sparse matrix to triplet format and build fragment dataframe.
triplet = DOSE::summary(count_matrix)
result_df = data.frame(
  chromosome = chromosome[triplet$i],
  start = start[triplet$i],
  end = end[triplet$i],
  cell = cells[triplet$j],
  counts = triplet$x,
  stringsAsFactors = FALSE
)

# Sort by genomic coordinates and cell barcodes; retain only standard chromosomes.
result_df = result_df[base::order(result_df$chromosome, result_df$start, result_df$end, result_df$cell), ]
result_df = result_df[base::grepl("chr", result_df$chromosome), ] # Remove non‑standard chromosomes.
utils::head(result_df, n = 3)
# ↓↓↓ Expected output ↓↓↓ #
#          chromosome  start    end               cell counts
# 3890295        chr1 777634 779926 AAACAGCCATCCAGGT-1      2
# 10964721       chr1 777634 779926 AAACCAACACAATGCC-1      2
# 1509806        chr1 777634 779926 AAACCAACATAATCCG-1      2
# ↑↑↑ Expected output ↑↑↑ #

# Export uncompressed fragment TSV file.
write.table(result_df, file = "./pbmc_atac_fragments.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
rm(count_matrix, fragments, cells, fragments_split, chromosome, start, end, triplet, result_df)



########################################
# 0.2.2 Generate _fragments.tsv.gz.tbi #
########################################
### Linux environment required for the following commands.

# Install tabix (tool for indexing and fast retrieval of genomic interval files).
sudo apt-get install tabix

# Navigate to working directory.
cd /mnt/i/scDERI

# Create backup of uncompressed fragment file.
cp pbmc_atac_fragments.tsv pbmc_atac_fragments_bak.tsv

# Compress with bgzip (creates gzip file compatible with tabix indexing).
# -i flag generates an index during compression.
bgzip -i pbmc_atac_fragments.tsv

# Create tabix index for BED-format file.
# Output file: pbmc_atac_fragments.tsv.gz.tbi
tabix -p bed pbmc_atac_fragments.tsv.gz





######################################
######################################
###                                ###
###   Part 1: Data Preprocessing   ###
###                                ###
######################################
######################################




#######################
#                     #
# 1.1 Identifying TFs #
#                     #
#######################
### Function
# Integrates transcription factor (TF) candidates from JASPAR2024 and ChIP-Atlas.
# Applies multi-layer expression filtering based on scRNA-seq data and validates with ChIP-seq prior evidence.
# Returns a curated TF list for downstream gene regulatory network construction.

# See function documentation for input, logic and output details.
?get_TFs_from_priors

# Extend timeout for remote file downloads.
options(timeout = 999999999)


### Data import
setwd(Working_directory)

# scRNA-seq Seurat object with counts layer for expression quantification and TF screening
scRNAseq = qs::qread("./scRNAseq.qsave")

# Species for JASPAR2024 PWM retrieval.
# Accepts Latin name (e.g. "Homo sapiens") or NCBI tax ID (e.g. "9606").
species = "Homo sapiens"

# Install rbioapi to list valid JASPAR species options.
install.packages("rbioapi")
base::print(utils::head(rbioapi::rba_jaspar_species(release = 2024)$species))
# ↓↓↓ Expected output ↓↓↓ #
# Retrieving a list of species available in JASPAR release 2024.
# [1] "Antirrhinum majus"                "Arabidopsis lyrata subsp. lyrata" "Arabidopsis thaliana"
# [4] "Aspergillus nidulans"             "Bos taurus"                       "Caenorhabditis briggsae"
# ↑↑↑ Expected output ↑↑↑ #

# JASPAR matrix collections to include
# Available options: "CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED"
collection = c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED")

# Genome assembly version for ChIP-Atlas data retrieval
# Supported assemblies: "hg38", "hg19", "mm10", "mm9", "rn6", "dm6", "dm3", "ce11", "ce10", "sacCer3"
genome_assembly = "hg38"
# genome_assembly = "hg19" # Use for hg19-aligned datasets.

# Cell type keywords to exclude from ChIP-Atlas samples
excluded_cell_types_in_priors = "NotExclude" # Set to "NotExclude" to disable filtering.
# excluded_cell_types_in_priors = c("pbmc", "Peripheral", "blood", "mononuclear") # Exclude PBMC-related ChIP-seq priors.


### Hyperparameter configuration
# Minimum non-zero cell expression frequency threshold (0–1). Defaults to 0.05.
min_nonzero_freq = 0.05

# Minimum mean expression in non-zero cells. Defaults to 0 (no filter).
min_nonzero_mean = 0

# Minimum overall mean expression across all cells. Defaults to 0 (no filter).
min_overall_mean = 0

# Minimum genome-wide quantile of non-zero cell frequency (0–1). Defaults to 0 (no filter).
min_nonzero_freq_quantile = 0

# Minimum genome-wide quantile of non-zero cell mean expression (0–1). Defaults to 0 (no filter).
min_nonzero_mean_quantile = 0

# Minimum genome-wide quantile of overall mean expression (0–1). Defaults to 0 (no filter).
min_overall_mean_quantile = 0

# Distance cutoff from TSS for ChIP-Atlas target gene annotation (kilobases).
# Allowed values: 1, 5, 10. Defaults to 10.
distance_from_TSS = 10

# Minimum number of valid target genes to retain a TF prior.
# TFs with fewer targets are discarded after filtering. Defaults to 2000.
prior_min_TG_number = 2000

# Maximum number of top target genes to retain per TF. Defaults to 5000.
prior_top_n_keep = 5000

# Proportional cap of top target genes to retain per TF (0–1).
# Final count is the smaller value of prior_top_n_keep and nrow(prior) * prior_top_q_keep. Defaults to 0.7.
prior_top_q_keep = 0.7


### Run get_TFs_from_priors with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  TFs = get_TFs_from_priors(
    scRNAseq = scRNAseq,
    species = species, collection = collection, genome_assembly = genome_assembly,
    min_nonzero_freq = min_nonzero_freq, min_nonzero_mean = min_nonzero_mean, min_overall_mean = min_overall_mean,
    min_nonzero_freq_quantile = min_nonzero_freq_quantile, min_nonzero_mean_quantile = min_nonzero_mean_quantile, min_overall_mean_quantile = min_overall_mean_quantile,
    distance_from_TSS = distance_from_TSS,
    excluded_cell_types_in_priors = excluded_cell_types_in_priors,
    prior_min_TG_number = prior_min_TG_number, prior_top_n_keep = prior_top_n_keep, prior_top_q_keep = prior_top_q_keep
  )
})
message("Get ", length(TFs), " TFs, ", print_run_stats(monitor))
# Example output: Get 80 TFs, running: 2 hours 14 mins 13 secs, peak RAM: 0 GB 408 MB


### Check per-TF download status.
download_information = readRDS("./1.1 Data Preprocessing - Get TFs From Database/ChIP_Atlas_Priors/download_information.rds")
base::print(download_information[download_information$status == "URL not found or download error", ])
# ↓↓↓ Expected output ↓↓↓ #
#                                                           url
# ESR1 https://chip-atlas.dbcls.jp/data/hg38/target/ESR1.10.tsv
#                                                                                  destination
# ESR1 /mnt/i/scDERI/1.1 Data Preprocessing - Get TFs From Database/ChIP_Atlas_Priors/ESR1.tsv
#                               status          time
# ESR1 URL not found or download error 30.34024 mins
# ↑↑↑ Expected output ↑↑↑ #

# Manually download failed files. URLs are stored in download_information$url; save paths are specified in download_information$destination.


### Re-run get_TFs_from_priors until all TF prior files are downloaded successfully.
monitor = peakRAM::peakRAM({
  TFs = get_TFs_from_priors(
    scRNAseq = scRNAseq,
    species = species, collection = collection, genome_assembly = genome_assembly,
    min_nonzero_freq = min_nonzero_freq, min_nonzero_mean = min_nonzero_mean, min_overall_mean = min_overall_mean,
    min_nonzero_freq_quantile = min_nonzero_freq_quantile, min_nonzero_mean_quantile = min_nonzero_mean_quantile, min_overall_mean_quantile = min_overall_mean_quantile,
    distance_from_TSS = distance_from_TSS,
    excluded_cell_types_in_priors = excluded_cell_types_in_priors,
    prior_min_TG_number = prior_min_TG_number, prior_top_n_keep = prior_top_n_keep, prior_top_q_keep = prior_top_q_keep
  )
})
message("Get ", length(TFs), " TFs, ", print_run_stats(monitor))
# Example output: Get 81 TFs, running: 0 hours 2 mins 5 secs, peak RAM: 0 GB 604 MB


### Inspect results.
# Print sorted TF gene symbols passing all filtering steps.
base::print(TFs)
# ↓↓↓ Expected output ↓↓↓ #
#  [1] "ATF3"    "ATF7"    "BACH2"   "BATF"    "BATF3"   "BCL11A"  "BCL6"    "BHLHE40" "CEBPA"   "CEBPB"   "CEBPD"   "CLOCK"   "CREM"    "EGR1"
# [15] "ESR1"    "ETS1"    "FOSL1"   "FOXK1"   "HES1"    "HIF1A"   "HSF2"    "IKZF1"   "IRF8"    "JDP2"    "JUN"     "JUNB"    "KLF10"   "KLF11"
# [29] "KLF12"   "KLF3"    "KLF9"    "MAFF"    "MEF2C"   "MGA"     "MITF"    "MXI1"    "MYC"     "NFATC3"  "NFIA"    "NFIL3"   "NFYA"    "NR2C1"
# [43] "NR2C2"   "NR4A1"   "NRL"     "PBX3"    "PKNOX1"  "POU2F1"  "POU2F2"  "PRDM1"   "RELB"    "REST"    "RFX3"    "RORA"    "RUNX2"   "RUNX3"
# [57] "SATB1"   "SMAD3"   "SMAD5"   "SP4"     "SREBF2"  "TCF4"    "TCF7"    "TCF7L2"  "TFCP2"   "TFDP1"   "TGIF2"   "TP53"    "TRPS1"   "XBP1"
# [71] "ZBED4"   "ZBTB11"  "ZEB1"    "ZNF143"  "ZNF148"  "ZNF331"  "ZNF341"  "ZNF524"  "ZNF652"  "ZNF76"   "ZNF766"
# ↑↑↑ Expected output ↑↑↑ #

### Save output.
qs::qsave(TFs, file = "./1.1 Data Preprocessing - Get TFs From Database/TFs.qsave")




##########################################################
#                                                        #
# 1.2 Annotating Chromatin Fragments and Identifying TGs #
#                                                        #
##########################################################



############################################################
# 1.2.1 Annotating Chromatin Fragments and Identifying TGs #
############################################################
### Function
# Performs genomic annotation on scATAC-seq chromatin accessibility peaks.
# Identifies target genes associated with open promoter regions.
# Retains only genes detected in matched scRNA-seq data and generates annotation tables.

# See function documentation for input, logic and output details.
?identify_TGs


### Data import
setwd(Working_directory)

# scRNA-seq Seurat object with gene names as rownames, used for target gene filtering.
scRNAseq = qs::qread("./scRNAseq.qsave")

# scATAC-seq Seurat object with peak ranges for genomic annotation.
# Requires matched _fragments.tsv.gz and .tbi index files in the Signac framework.
scATACseq = qs::qread("./scATACseq.qsave")

# TxDb object for genome annotation.
genome_for_anno = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# genome_for_anno = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene # Use for hg19-aligned datasets.


# OrgDb object for organism-level gene ID conversion.
annoDb = "org.Hs.eg.db"


### Hyperparameter configuration
# Annotation range around TSS in base pairs, used to define promoter-proximal regions.
# Smaller values reduce runtime but may exclude distal regulatory elements. Defaults to 500000.
fragment_annotation_range = 500000

# Distance (bp) extended around each peak to search for flanking genes.
# See ?ChIPseeker::annotatePeak for implementation details. Defaults to 10000.
flankDistance = 10000


### Run identify_TGs with runtime and peak memory monitoring.
library(dplyr)
monitor = peakRAM::peakRAM({
  results_identify_TGs_all = identify_TGs(
    scRNAseq = scRNAseq, scATACseq = scATACseq,
    genome_for_anno = genome_for_anno, annoDb = annoDb, flankDistance = flankDistance
  )
})
message(print_run_stats(monitor))
# Example output: running: 0 hours 0 mins 58 secs, peak RAM: 0 GB 725 MB


### Inspect results.
# Preview unique target gene symbol vector.
base::print(utils::head(results_identify_TGs_all[["TGs"]])) # Example output: "ISG15"    "AGRN"     "C1orf159" "FNDC10"   "SKI"      "MMEL1"
base::print(length(results_identify_TGs_all[["TGs"]])) # Total identified target genes: 2731

# Annotated promoter peak table with flanking gene symbols; all genes are present in scRNA-seq.
base::print(base::colnames(results_identify_TGs_all[["peak_anno_promoter"]]))
# View(results_identify_TGs_all[["peak_anno_promoter"]])
# ↓↓↓ Expected output ↓↓↓ #
#  [1] "seqnames"             "start"                "end"                  "width"                "strand"
#  [6] "annotation"           "geneChr"              "geneStart"            "geneEnd"              "geneLength"
# [11] "geneStrand"           "geneId"               "transcriptId"         "distanceToTSS"        "ENSEMBL"
# [16] "SYMBOL"               "GENENAME"             "flank_txIds"          "flank_geneIds"        "flank_gene_distances"
# [21] "flank_geneSymbols"
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(results_identify_TGs_all, file = "./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs_all.qsave")



######################################
# 1.2.2 Filtering Annotation Results #
######################################
### Function
# Applies distance-adaptive chromatin accessibility quantile thresholds to refine peak-gene annotations.
# Integrates gene expression quality adjustment and supports an optional per-gene fallback rule.
# Built-in disk caching accelerates repeated runs.

# See function documentation for input, logic and output details.
?filter_peak_annotations


### Data import
setwd(Working_directory)

# Output list from the identify_TGs function.
results_identify_TGs_all = qs::qread("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs_all.qsave")

# scRNA-seq Seurat object with counts layer for gene quality score calculation.
scRNAseq = qs::qread("./scRNAseq.qsave")

# scATAC-seq Seurat object with counts layer for peak accessibility calculation.
scATACseq = qs::qread("./scATACseq.qsave")


### Hyperparameter configuration
# Genomic distance anchor points (bp) to TSS for threshold curve fitting.
# Must be strictly increasing and match the length of min_quantile_points. Defaults to c(200, 5000, 500000).
distance_points = c(1, 500, 500000)

# Minimum quantile thresholds corresponding to each distance anchor.
# Must be strictly increasing, with values in [0, 1]. Defaults to c(0, 0.6, 1).
min_quantile_points = c(0.35, 0.55, 1.55)

# Starting values for constrained nonlinear least squares curve fitting (named list: a, b, c).
# Default values are optimized for the default anchor points.
distance_quantile_fit_start = list(a = 0.006667, b = 1.659167, c = 0.6003)

# If TRUE, guarantee at least one retained peak per target gene.
# Single-peak genes are kept unconditionally; multi-peak genes with zero passing peaks retain the most accessible one. Defaults to TRUE.
one_peak_at_least_per_gene = TRUE

# Minimum scaling factor (0–1) for high-quality genes (smaller value = more relaxed thresholds). Defaults to 0.8.
gene_quality_factor_min = 0.9

# Weight of non-zero cell frequency in the composite gene quality score (0–1).
weight_gene_nonzero_freq = 0.7

# Weight of non-zero mean expression in the composite gene quality score (0–1).
weight_gene_nonzero_mean = 0.2

# Weight of overall mean expression in the composite gene quality score (0–1).
# The three weights must sum to 1. Defaults: 0.2, 0.2, 0.6.
weight_gene_overall_mean = 0.1

# Scaling factor for peak non-zero frequency quantile threshold (smaller = more relaxed) (>0). Default: 0.7.
min_peak_nonzero_freq_quantile_scale = 0.4

# Scaling factor for peak non-zero mean accessibility quantile threshold (smaller = more relaxed) (>0). Default: 0.7.
min_peak_nonzero_mean_scale = 1.2

# Scaling factor for peak overall mean accessibility quantile threshold (smaller = more relaxed) (>0). Default: 1.
min_peak_overall_mean_quantile_scale = 0.1


### Run filter_peak_annotations with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  results_identify_TGs = filter_peak_annotations(
    results_identify_TGs = results_identify_TGs_all, scRNAseq = scRNAseq, scATACseq = scATACseq,
    distance_points = distance_points, min_quantile_points = min_quantile_points, distance_quantile_fit_start = distance_quantile_fit_start,
    one_peak_at_least_per_gene = one_peak_at_least_per_gene,
    gene_quality_factor_min = gene_quality_factor_min,
    weight_gene_nonzero_freq = weight_gene_nonzero_freq, weight_gene_nonzero_mean = weight_gene_nonzero_mean, weight_gene_overall_mean = weight_gene_overall_mean,
    min_peak_nonzero_freq_quantile_scale = min_peak_nonzero_freq_quantile_scale, min_peak_nonzero_mean_scale = min_peak_nonzero_mean_scale, min_peak_overall_mean_quantile_scale = min_peak_overall_mean_quantile_scale
  )
})
message(base::nrow(results_identify_TGs[["peak_anno_promoter_all"]]), "/", length(base::unique(results_identify_TGs[["peak_anno_promoter"]][["peak"]])), " peaks, ",
        length(base::unique(results_identify_TGs[["peak_anno_promoter_all"]][["SYMBOL"]])), "/", length(results_identify_TGs[["TGs"]]), " TGs, ",
        base::nrow(results_identify_TGs[["peak_anno_promoter_all"]]), "/", base::nrow(base::unique(results_identify_TGs[[1]][,c("peak","SYMBOL")])), " peak-TG, peaks/TGs = ",
        round(length(base::unique(results_identify_TGs[["peak_anno_promoter"]][["peak"]]))/length(results_identify_TGs[["TGs"]]), digits = 2), ", ",
        print_run_stats(monitor))
# Example output: 3473/2724 peaks, 2731/2868 TGs, 3473/2911 peak-TG, peaks/TGs = 0.95, running: 0 hours 0 mins 6 secs, peak RAM: 0 GB 597 MB


### Inspect results.
# Preview filtered peak-gene annotation table with quality metrics.
View(results_identify_TGs_all[["peak_anno_promoter"]])


### Save output.
qs::qsave(results_identify_TGs, file = "./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.qsave")




##########################################################################
#                                                                        #
# 1.3 Generating Counts Matrices for TGs/TFs Expression and TGs Activity #
#                                                                        #
##########################################################################
### Function
# Extracts expression matrices of target genes (TGs) and transcription factors (TFs) from multi-omics data.
# Computes a gene activity matrix from scATAC-seq data for downstream GRN construction.

# See function documentation for input, logic and output details.
?get_expression_and_activity_matrix


### Data import
setwd(Working_directory)

# Character vector of target gene symbols.
TGs = qs::qread("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.qsave")$TGs # Character vector of target gene symbols.

# Character vector of transcription factor symbols from get_TFs_from_priors output.
TFs = qs::qread("./1.1 Data Preprocessing - Get TFs From Database/TFs.qsave")

# scRNA-seq Seurat object.
scRNAseq = qs::qread("./scRNAseq.qsave")

# scATAC-seq Seurat object.
scATACseq = qs::qread("./scATACseq.qsave")

# EnsDb genomic annotation object passed to Signac::GetGRangesFromEnsDb.
# Must match the species and genome build of input data.
Ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
# Ensdb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75 # Use for hg19-aligned datasets.

# Absolute path to scATAC-seq fragment file (.tsv.gz). A matching .tbi index file must exist in the same directory.
path = paste0(Working_directory, "/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")
# path = paste0(Working_directory, "/pbmc_atac_fragments.tsv.gz") # Use self-generated file as alternative.


### Hyperparameter configuration
# Base pairs to extend for gene activity calculation (>0).
# Applied both upstream of TSS and downstream of TTS, passed to Signac::GeneActivity. Defaults to 10000.
flanking_region = 10000


### Run get_expression_and_activity_matrix with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  Basic_Info = get_expression_and_activity_matrix(
    TGs = TGs, TFs = TFs, scRNAseq = scRNAseq, scATACseq = scATACseq, Ensdb = Ensdb,
    flanking_region = flanking_region, path = path
  )
})
message(length(Basic_Info[["TFs"]]), " TFs, ", length(Basic_Info[["TGs"]]), " TGs, ",
        length(Basic_Info[["scRNAseq_for_GRN"]]@meta.data$seurat_annotations), " cells, ", print_run_stats(monitor))
# Example output: 79 TFs, 2657 TGs, 3000 cells, running: 0 hours 19 mins 24 secs, peak RAM: 3 GB 254 MB


### Inspect results.
# scRNA-seq Seurat object subsetted to genes shared with the scATAC-seq dataset.
base::print(Basic_Info[["scRNAseq_for_GRN"]])
# ↓↓↓ Expected output ↓↓↓ #
# An object of class Seurat
# 2691 features across 3000 samples within 1 assay
# Active assay: RNA (2691 features, 2691 variable features)
#  2 layers present: counts, data
# ↑↑↑ Expected output ↑↑↑ #

# scATAC-seq Seurat object with ACTIVITY assay, subsetted to shared genes; ACTIVITY is the default assay.
base::print(Basic_Info[["scATACseq_for_GRN"]])
# ↓↓↓ Expected output ↓↓↓ #
# An object of class Seurat
# 2691 features across 3000 samples within 1 assay
# Active assay: ACTIVITY (2691 features, 0 variable features)
#  2 layers present: counts, data
# ↑↑↑ Expected output ↑↑↑ #

# Genes present in both scRNA-seq and scATAC-seq datasets.
base::print(utils::head(Basic_Info[["gene_for_GRN"]])) # Example output: "ISG15"    "AGRN"     "C1orf159" "FNDC10"   "SKI"      "MMEL1"

# Filtered target genes retained in both omics datasets.
base::print(utils::head(Basic_Info[["TGs"]])) # Example output: "ISG15"    "AGRN"     "C1orf159" "FNDC10"   "SKI"      "MMEL1"

# Filtered transcription factors retained in both omics datasets.
base::print(utils::head(Basic_Info[["TFs"]])) # Example output: "ATF3"   "ATF7"   "BACH2"  "BATF"   "BATF3"  "BCL11A"


### Save output.
qs::qsave(Basic_Info, file = "./1.3 Data Preprocessing - Get Expression And Activity Matrix/Basic_Info.qsave")





#######################################################################
#######################################################################
###                                                                 ###
###   Part 2: Inferring Pseudotime, Branching, and Grouping Cells   ###
###                                                                 ###
#######################################################################
#######################################################################




########################################################
#                                                      #
# 2.1 Pseudotime Analysis and Cell Branch Partitioning #
#                                                      #
########################################################



###########################################################
# 2.1.1 Obtaining Relevant Data of Cell Types of Interest #
###########################################################
### Function
# Extracts scRNA-seq and scATAC-seq data for specified cell types.
# Performs quality control by filtering genes and cells based on zero-expression rates.
# Aligns the two omics datasets for downstream GRN construction and pseudotime analysis.

# See function documentation for input, logic and output details.
?get_interest_cell_type_data


### Data import
setwd(Working_directory)

# Output list from get_expression_and_activity_matrix.
# Contains scRNAseq_for_GRN, scATACseq_for_GRN, gene_for_GRN, TGs and TFs elements.
Basic_Info = qs::qread("./1.3 Data Preprocessing - Get Expression And Activity Matrix/Basic_Info.qsave")

# Target cell type names to extract.
# Must match values in the seurat_annotations metadata column.
interest_cell_type = base::unique(Basic_Info[["scRNAseq_for_GRN"]]@meta.data[["seurat_annotations"]]) # Target cell type names to extract
base::print(interest_cell_type) # Example output: Mono

# Number of cores for parallel computation via mclapply.
# Parallel processing only works on Unix-like systems; set to 1 on Windows.
# Adjust based on available RAM; reduce core count if memory is exhausted.
install.packages("parallel")
packageVersion("parallel") # Version: 4.4.1
ncores = parallel::detectCores() - 1 # Linux setting; 79 cores detected during development
# ncores = 1 # Windows setting
base::print(ncores) # Example output: 78


### Run get_interest_cell_type_data with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_data = get_interest_cell_type_data(
    interest_cell_type = interest_cell_type, Basic_Info = Basic_Info, ncores = 1
  )
})
message(print_run_stats(monitor))
# Example output: running: 0 hours 0 mins 3 secs, peak RAM: 0 GB 567 MB


### Inspect results.
# Filtered scRNA-seq Seurat object.
base::print(interest_cell_type_data[[1]][["interest_scRNAseq_for_GRN"]])
# ↓↓↓ Expected output ↓↓↓ #
# An object of class Seurat
# 2691 features across 3000 samples within 1 assay
# Active assay: RNA (2691 features, 2691 variable features)
#  2 layers present: counts, data
# ↑↑↑ Expected output ↑↑↑ #

# Filtered scATAC-seq Seurat object, fully aligned with scRNA-seq data.
base::print(Basic_Info[["scATACseq_for_GRN"]])
# ↓↓↓ Expected output ↓↓↓ #
# An object of class Seurat
# 2691 features across 3000 samples within 1 assay
# Active assay: ACTIVITY (2691 features, 0 variable features)
#  2 layers present: counts, data
# ↑↑↑ Expected output ↑↑↑ #

# Gene set for GRN construction (union of filtered TFs and TGs).
base::print(utils::head(interest_cell_type_data[[1]][["interest_gene_for_GRN"]])) # Example output: "A2M"   "AAK1"  "AAMDC" "AANAT" "AATK"  "ABCA1"

# Filtered target genes retained in the current cell type.
base::print(utils::head(interest_cell_type_data[[1]][["interest_TGs"]])) # Example output: "A2M"   "AAK1"  "AAMDC" "AANAT" "AATK"  "ABCA1"

# Filtered transcription factors retained in the current cell type.
base::print(utils::head(interest_cell_type_data[[1]][["interest_TFs"]])) # Example output: "ATF3"   "ATF7"   "BACH2"  "BATF"   "BATF3"  "BCL11A"


### Save output.
qs::qsave(interest_cell_type_data, file = "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.qsave")




##########################################################################################################################
# 2.1.2 Performing Pseudotemporal Ordering and Branching Analysis of Cell Types of Interest Based on Transcriptomic Data #
##########################################################################################################################
### Function
# Builds developmental trajectories and performs pseudotime ordering using Monocle2 (DDRTree).
# Supports both automatic and manual parameter tuning for trajectory construction.
# Partitions cells into branches via interactive user input based on visualization results.

# See function documentation for input, logic and output details.
?order_pseudotime_and_divide_branches


### Data import
setwd(Working_directory)

# Output list from get_interest_cell_type_data.
# Contains filtered scRNA-seq objects and gene sets for each cell type.
interest_cell_type_data = qs::qread("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.qsave")


### Hyperparameter configuration
# Whether to enable automatic parameter selection for DDRTree.
# Passed to monocle::reduceDimension; uses default DDRTree parameters when TRUE.
auto_param_selection = TRUE # Used in this tutorial.

# Manual DDRTree tuning (used when auto_param_selection = FALSE).
# maxIter: Maximum iterations for DDRTree. Defaults to 20.
# sigma: Bandwidth parameter for DDRTree. Defaults to 0.001.
# lambda: Regularization parameter for DDRTree. Defaults to NULL.
# ncenter: Number of centroids for DDRTree. Defaults to NULL.
# param.gamma: Gamma parameter for DDRTree. Defaults to 10.
# tol: Convergence tolerance for DDRTree. Defaults to 0.001.

# Configuration notes:
# - Use auto_param_selection = TRUE for default behavior.
# - Switch to FALSE if you are unsatisfied with automatic branching results.
# - Parameters are passed to DDRTree::DDRTree; see ?DDRTree for details.
# - Adjust ncenter alone to control branch count; ncenter = 4 often produces 1 branch (empirical suggestion).


### Run order_pseudotime_and_divide_branches with runtime and peak memory monitoring.
library(DDRTree)
monitor = peakRAM::peakRAM({
  interest_cell_type_Branches = order_pseudotime_and_divide_branches(
    interest_cell_type_data = interest_cell_type_data,
    auto_param_selection = TRUE
    # auto_param_selection = FALSE, ncenter = 4 # Force 1 branch; not demonstrated in this tutorial
  )
})
# ↓↓↓ Running Logs ↓↓↓ #
# ...
# Interactively input the number of branches (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).
# Please enter the number of branches based on the pseudotime 2D visualization (e.g., 3) (e.g., 5):
# ↑↑↑ Running Logs ↑↑↑ #

# The function pauses here to wait for user input.
# 1. Open the output folder to view trajectory plots.
# 2. Identify the trajectory start (darkest color) from interest_Pseudotime_in_trajectory.png.
# 3. Count the total number of branches.
# For this dataset, enter 2 in the console and press Enter.

# ↓↓↓ Running Logs ↓↓↓ #
# ...
# Interactively input the number of branches (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).
# Please enter the number of branches based on the pseudotime 2D visualization (e.g., 3) (e.g., 5): 2
# Interactively input the integer vector of cell state order for each branch (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).
# Please enter the integer vector of cell state order for branch 1 (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5):
# ↑↑↑ Running Logs ↑↑↑ #

# The function pauses a second time for cell state order input.
# Use interest_AllCellState_in_trajectory.png as the primary reference.
# For this dataset, the order for branch 1 is 1 → 2; enter 1 2 and press Enter.

# ↓↓↓ Running Logs ↓↓↓ #
# ...
# Interactively input the number of branches (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).
# Please enter the number of branches based on the pseudotime 2D visualization (e.g., 3) (e.g., 5): 2
# Interactively input the integer vector of cell state order for each branch (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).
# Please enter the integer vector of cell state order for branch 1 (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5):
#   1 2
# Please enter the integer vector of cell state order for branch 2 (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5):
# ↑↑↑ Running Logs ↑↑↑ #

# Enter 1 3 for branch 2 and press Enter.

# ↓↓↓ Running Logs ↓↓↓ #
# ...
# Interactively input the number of branches (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).
# Please enter the number of branches based on the pseudotime 2D visualization (e.g., 3) (e.g., 5): 2
# Interactively input the integer vector of cell state order for each branch (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).
# Please enter the integer vector of cell state order for branch 1 (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5):
#   1 2
# Please enter the integer vector of cell state order for branch 2 (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5):
#   1 3
# Extracting pseudotime and cell states for each branch.
# ...
# ↑↑↑ Running Logs ↑↑↑ #

# Interactive input is complete; wait for the function to finish.

message(print_run_stats(monitor))
# Example output: running: 0 hours 3 mins 55 secs, peak RAM: 1 GB 164 MB


### Inspect results.
# Pseudotime and state annotations for cells in each branch.
base::print(utils::head(interest_cell_type_Branches[[1]][[1]], n = 3)) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#                    Pseudotime State
# GATAATCGTTGCACGG-1 0.00000000     1
# GAAGTATAGAGAAGGG-1 0.02918536     1
# CTTGAACAGGCTGTCA-1 0.03536132     1
# ↑↑↑ Expected output ↑↑↑ #

base::print(utils::head(interest_cell_type_Branches[[1]][[2]], n = 3)) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#                    Pseudotime State
# GATAATCGTTGCACGG-1 0.00000000     1
# GAAGTATAGAGAAGGG-1 0.02918536     1
# CTTGAACAGGCTGTCA-1 0.03536132     1
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(interest_cell_type_Branches, file = "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_Branches.qsave")




##############################################################
#                                                            #
# 2.2 Cell Grouping and Group-Wise Information Specification #
#                                                            #
##############################################################



#############################################################################
# 2.2.1 Retrieving Pseudotemporal Expression/Activity Vectors for Each Gene #
#############################################################################
### Function
# Calculates gene importance scores and filters genes and cells for each cell type and branch.
# Extracts pseudotime-series matrices of TF expression, TG expression and chromatin activity.
# Prepares data for downstream GRN inference.

# See function documentation for input, logic and output details.
?get_genes_pseudotime_info


### Data import
setwd(Working_directory)

# Output list from order_pseudotime_and_divide_branches.
interest_cell_type_Branches = qs::qread("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_Branches.qsave")

# Output list from get_interest_cell_type_data.
interest_cell_type_data = qs::qread("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.qsave")

# Number of cores for parallel computation. Defaults to 1.
# Parallel processing only works on Unix-like systems; set to 1 on Windows.
ncores = parallel::detectCores() - 1 # Linux setting
# ncores = 1 # Windows setting


### Hyperparameter configuration
# Quantile threshold for TF filtering (0–1). Defaults to 0.63.
# TFs with importance scores above this quantile are retained.
alpha_TF = 0

# When alpha_TF = 0, no TF filtering is performed by importance score.
# Weight parameters for TF importance score (TFIS) calculation (range: -0.2 to 1):
# w_TFNER: Weight for TF non-zero expression ratio
# w_TFNERR: Weight for normalized rank of TF non-zero expression ratio
# w_TFNEA: Weight for normalized TF non-zero expression average
# w_TFNEAR: Weight for normalized rank of TF non-zero expression average
# w_TFNAR: Weight for TF non-zero activity ratio
# w_TFNARR: Weight for normalized rank of TF non-zero activity ratio
# w_TFNAA: Weight for normalized TF non-zero activity average
# w_TFNAAR: Weight for normalized rank of TF non-zero activity average

# Quantile threshold for TG filtering (0–1). Defaults to 0.69.
# TGs with importance scores above this quantile are retained.
alpha_TG = 0

# When alpha_TG = 0, no TG filtering is performed by importance score.
# Weight parameters for TG importance score (TGIS) calculation (range: -0.2 to 1):
# w_TGNER: Weight for TG non-zero expression ratio
# w_TGNERR: Weight for normalized rank of TG non-zero expression ratio
# w_TGNEA: Weight for normalized TG non-zero expression average
# w_TGNEAR: Weight for normalized rank of TG non-zero expression average
# w_TGNAR: Weight for TG non-zero activity ratio
# w_TGNARR: Weight for normalized rank of TG non-zero activity ratio
# w_TGNAA: Weight for normalized TG non-zero activity average
# w_TGNAAR: Weight for normalized rank of TG non-zero activity average

# Zero-rate threshold for cell filtering (0–1). Defaults to 1 (no filtering).
# Cells with combined zero-expression rate ≤ this value are retained.
alpha_cell = 1
# No cell filtering is applied here, as downsampling already removed sparse cells.

# Optional post-filtering by non-zero rate. NULL by default.
# Top-level names must match cell types; each element is a data frame with rows TFE, TGA, TGE and columns for branches.
# Values range from 0 to 1. Invalid formats are reset to NULL with a message.
sparsity_cut = list()
sparsity_cut[["Mono"]] = base::as.data.frame(base::matrix(
  0, nrow = 3, ncol = length(interest_cell_type_Branches[["Mono"]]),
  dimnames = list(c("TFE", "TGA", "TGE"), paste0("branch", base::seq_along(interest_cell_type_Branches[["Mono"]])))
))
base::print(sparsity_cut)
# ↓↓↓ Expected output ↓↓↓ #
# $Mono
#     branch1 branch2
# TFE       0       0
# TGA       0       0
# TGE       0       0
# ↑↑↑ Expected output ↑↑↑ #

sparsity_cut[["Mono"]]["TFE", "branch1"] = 0.06
sparsity_cut[["Mono"]]["TFE", "branch2"] = 0.06
sparsity_cut[["Mono"]]["TGA", "branch1"] = 0.10
sparsity_cut[["Mono"]]["TGA", "branch2"] = 0.10
sparsity_cut[["Mono"]]["TGE", "branch1"] = 0.11
sparsity_cut[["Mono"]]["TGE", "branch2"] = 0.11
base::print(sparsity_cut)
# ↓↓↓ Expected output ↓↓↓ #
# $Mono
#     branch1 branch2
# TFE    0.06    0.06
# TGA    0.10    0.10
# TGE    0.11    0.11
# ↑↑↑ Expected output ↑↑↑ #


### Run get_genes_pseudotime_info with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_genes_pseudotime_info = get_genes_pseudotime_info(
    interest_cell_type_Branches = interest_cell_type_Branches, interest_cell_type_data = interest_cell_type_data,
    alpha_cell = alpha_cell, alpha_TF = alpha_TF, alpha_TG = alpha_TG,
    ncores = ncores, sparsity_cut = sparsity_cut
  )
})
message(print_run_stats(monitor))
# Example output: running: 0 hours 0 mins 12 secs, peak RAM: 0 GB 1011 MB


### Inspect results.
# TF expression matrices (genes × cells ordered by pseudotime, normalized values).
base::print(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TFE"]][[1]][1:3, 1:10]) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#             1      2       3        4        5        6        7 8        9       10
# ATF3 2.584367 3.5928 2.94028 3.316381 2.767914 3.962961 3.665382 0 0.000000 0.000000
# ATF7 0.000000 0.0000 0.00000 0.000000 0.000000 0.000000 0.000000 0 2.745812 0.000000
# BATF 0.000000 0.0000 0.00000 0.000000 0.000000 0.000000 0.000000 0 0.000000 2.651496
# ↑↑↑ Expected output ↑↑↑ #

base::print(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TFE"]][[2]][1:3, 1:10]) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#              1      2       3        4        5        6        7 8        9 10
# ATF3  2.584367 3.5928 2.94028 3.316381 2.767914 3.962961 3.664146 0 0.000000  0
# ATF7  0.000000 0.0000 0.00000 0.000000 0.000000 0.000000 0.000000 0 2.745812  0
# BACH2 0.000000 0.0000 0.00000 0.000000 0.000000 0.000000 0.000000 0 0.000000  0
# ↑↑↑ Expected output ↑↑↑ #

# TG chromatin activity matrices (genes × cells ordered by pseudotime, normalized values).
base::print(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TGA"]][[1]][1:3, 1:10]) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#              1        2        3        4         5        6        7        8        9       10
# AAK1   2.197047 2.538513 1.248494 2.360564 2.5018032 1.874685 1.804158 2.963056 2.208864 1.904237
# ABCA1  0.000000 1.465953 2.134790 1.576661 0.0000000 0.000000 1.263400 2.880730 2.409798 2.519998
# ABCB10 2.197047 2.233146 1.786777 1.910750 0.8757021 1.874685 1.804158 1.720847 1.619920 1.570598
# ↑↑↑ Expected output ↑↑↑ #

base::print(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TGA"]][[2]][1:3, 1:10]) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#               1        2        3        4         5        6        7        8        9       10
# AAK1   2.197225 2.538360 1.248494 2.360564 2.5011606 1.874685 1.803946 2.963056 2.208864 1.904237
# ABCA1  0.000000 1.465824 2.134790 1.576661 0.0000000 0.000000 1.263218 2.880730 2.409798 2.519998
# ABCB10 2.197225 2.232997 1.786777 1.910750 0.8752938 1.874685 1.803946 1.720847 1.619920 1.570598
# ↑↑↑ Expected output ↑↑↑ #

# TG expression matrices (genes × cells ordered by pseudotime, normalized values).
base::print(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TGE"]][[1]][1:3, 1:10]) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#               1        2 3 4        5        6        7        8 9       10
# AAK1   2.584367 0.000000 0 0 3.429161 3.288643 0.000000 2.815324 0 3.308733
# ABCA1  0.000000 0.000000 0 0 0.000000 0.000000 2.616692 0.000000 0 0.000000
# ABCB10 0.000000 2.926802 0 0 0.000000 0.000000 0.000000 0.000000 0 0.000000
# ↑↑↑ Expected output ↑↑↑ #

base::print(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TGE"]][[2]][1:3, 1:10]) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#               1        2 3 4        5        6        7        8 9       10
# AAK1   2.584367 0.000000 0 0 3.429161 3.288643 0.000000 2.815324 0 3.308733
# ABCA1  0.000000 0.000000 0 0 0.000000 0.000000 2.615516 0.000000 0 0.000000
# ABCB10 0.000000 2.926802 0 0 0.000000 0.000000 0.000000 0.000000 0 0.000000
# ↑↑↑ Expected output ↑↑↑ #

# Number of retained cells per branch.
base::print(base::unlist(interest_cell_type_genes_pseudotime_info[[1]][["Branches_n_cell"]])) # Example output: 2620 2194

# Full branch annotation with cell-pseudotime index mapping.
base::print(utils::head(interest_cell_type_genes_pseudotime_info[[1]][["Branches"]][[1]], n = 3)) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#                    Pseudotime State               cell k~
# GATAATCGTTGCACGG-1 0.00000000     1 GATAATCGTTGCACGG-1  1
# GAAGTATAGAGAAGGG-1 0.02918536     1 GAAGTATAGAGAAGGG-1  2
# CTTGAACAGGCTGTCA-1 0.03536132     1 CTTGAACAGGCTGTCA-1  3
# ↑↑↑ Expected output ↑↑↑ #

base::print(utils::head(interest_cell_type_genes_pseudotime_info[[1]][["Branches"]][[2]], n = 3)) # Branch 3
# ↓↓↓ Expected output ↓↓↓ #
#                    Pseudotime State               cell k~
# GATAATCGTTGCACGG-1 0.00000000     1 GATAATCGTTGCACGG-1  1
# GAAGTATAGAGAAGGG-1 0.02918536     1 GAAGTATAGAGAAGGG-1  2
# CTTGAACAGGCTGTCA-1 0.03536132     1 CTTGAACAGGCTGTCA-1  3
# ↑↑↑ Expected output ↑↑↑ #

# Retained TF names per branch.
base::print(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TFs"]][[1]]) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#  [1] "ATF3"    "ATF7"    "BATF"    "BATF3"   "BCL11A"  "BCL6"    "BHLHE40" "CEBPA"   "CEBPB"   "CEBPD"   "CLOCK"   "CREM"    "EGR1"    "FOSL1"   "FOXK1"
# [16] "HIF1A"   "HSF2"    "IKZF1"   "IRF8"    "JDP2"    "JUN"     "JUNB"    "KLF10"   "KLF11"   "KLF12"   "KLF3"    "KLF9"    "MAFF"    "MEF2C"   "MGA"
# [31] "MITF"    "MXI1"    "NFATC3"  "NFIL3"   "NFYA"    "NR2C1"   "NR2C2"   "NR4A1"   "PBX3"    "PKNOX1"  "POU2F1"  "POU2F2"  "RELB"    "REST"    "RFX3"
# [46] "RUNX2"   "RUNX3"   "SATB1"   "SMAD3"   "SMAD5"   "SP4"     "SREBF2"  "TCF4"    "TCF7L2"  "TFCP2"   "TFDP1"   "TGIF2"   "TP53"    "TRPS1"   "XBP1"
# [61] "ZBED4"   "ZBTB11"  "ZNF143"  "ZNF148"  "ZNF331"  "ZNF341"  "ZNF524"  "ZNF652"  "ZNF76"   "ZNF766"
# ↑↑↑ Expected output ↑↑↑ #

base::print(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TFs"]][[2]]) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#  [1] "ATF3"    "ATF7"    "BACH2"   "BATF"    "BCL11A"  "BCL6"    "BHLHE40" "CEBPA"   "CEBPB"   "CEBPD"   "CLOCK"   "CREM"    "EGR1"    "ESR1"    "ETS1"
# [16] "FOSL1"   "FOXK1"   "HIF1A"   "HSF2"    "IKZF1"   "IRF8"    "JDP2"    "JUN"     "JUNB"    "KLF10"   "KLF11"   "KLF12"   "KLF3"    "KLF9"    "MAFF"
# [31] "MEF2C"   "MGA"     "MITF"    "MXI1"    "NFATC3"  "NFIL3"   "NFYA"    "NR2C1"   "NR2C2"   "NR4A1"   "PBX3"    "PKNOX1"  "POU2F1"  "POU2F2"  "RELB"
# [46] "REST"    "RFX3"    "RUNX2"   "RUNX3"   "SATB1"   "SMAD3"   "SMAD5"   "SP4"     "SREBF2"  "TCF4"    "TCF7"    "TCF7L2"  "TFCP2"   "TFDP1"   "TGIF2"
# [61] "TP53"    "TRPS1"   "XBP1"    "ZBED4"   "ZBTB11"  "ZEB1"    "ZNF143"  "ZNF148"  "ZNF331"  "ZNF652"  "ZNF76"   "ZNF766"
# ↑↑↑ Expected output ↑↑↑ #

# Retained TG names per branch.
base::print(utils::head(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TGs"]][[1]])) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
# [1] "AAK1"   "ABCA1"  "ABCB10" "ABCD3"  "ABCF1"  "ABHD5"
# ↑↑↑ Expected output ↑↑↑ #

base::print(utils::head(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TGs"]][[2]])) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
# [1] "AAK1"   "ABCA1"  "ABCB10" "ABCD3"  "ABCF1"  "ABHD5"
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(interest_cell_type_genes_pseudotime_info, file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")



###############################################################################################################################
# 2.2.2 Grouping Cells Based on Pseudotime Analysis and Calculating Valid (Non-Zero) Expression and Activity Values per Group #
###############################################################################################################################
### Function
# Performs adaptive cell grouping along pseudotime trajectories for each cell type.
# Fits a logarithmic model to set minimum cells per group based on total cell count.
# Iteratively adjusts group sizes to ensure sufficient non-zero values per gene.
# Returns aggregated gene profiles and pseudotime endpoints for dynamic modeling.

# See function documentation for input, logic and output details.
?cell_grouping


### Data import
setwd(Working_directory)

# Output list from get_genes_pseudotime_info.
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")

# Number of cores for parallel computation. Defaults to 1.
# Parallel processing only works on Unix-like systems; set to 1 on Windows.
ncores = parallel::detectCores() - 1 # Linux setting
# ncores = 1 # Windows setting


### Hyperparameter configuration
# Three manually defined points to fit a logarithmic model y = c * ln(ax + b).
# x: total cell count; y: minimum cells per group.

# Total cell count anchor points (x-axis), strictly increasing positive integers.
n_cell = round(base::mean(base::unlist(interest_cell_type_genes_pseudotime_info[[1]][["Branches_n_cell"]])))
fit_scale_factor = 1.05 # generally >0.8, <1.2
points_x_for_fitting_nc_and_nmin = c(
  round(n_cell * 0.5 * fit_scale_factor),
  round(n_cell * 1.05),
  round(n_cell * 1.5 * fit_scale_factor)
)
base::print(points_x_for_fitting_nc_and_nmin) # Example output: 1264 2527 3791

# Corresponding minimum cells per group (y-axis), strictly increasing positive integers.
# Defaults to c(5, 10, 15).
n_min_expected = 4 # Target minimum cells per group, around 4.
points_y_for_fitting_nc_and_nmin = c(n_min_expected - 1, n_min_expected, n_min_expected + 1)
base::print(points_y_for_fitting_nc_and_nmin) # Example output: 3 4 5

# Fitting notes:
# - For arithmetic sequence y-values, ensure 2 * median(x) < sum of first and last x.
# - To force a specific n_min for a cell count, set the corresponding (x, y) pair directly.
# - Example: for 1234 cells and n_min = 12, set the i-th x to 1234*factor (factor is in [1, 1.1]) and i-th y to 12.


### Run cell_grouping with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_group = cell_grouping(
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    points_x_for_fitting_nc_and_nmin = points_x_for_fitting_nc_and_nmin,
    points_y_for_fitting_nc_and_nmin = points_y_for_fitting_nc_and_nmin,
    ncores = ncores
  )
})
message(print_run_stats(monitor))
# Example output: running: 0 hours 1 mins 3 secs, peak RAM: 0 GB 910 MB


### Inspect results.
# Number of valid groups per branch; 0 if grouping fails.
base::print(base::unlist(interest_cell_type_group[[1]][["Branches_n_group"]])) # Example output: 8 8

# Exact cell count of each group per branch.
base::print(interest_cell_type_group[[1]][["Branches_n_k"]])
# ↓↓↓ Expected output ↓↓↓ #
# [[1]]
# [1] 583 452 309 340 135 159  94 548
# [[2]]
# [1] 738 479 412 242  80  65  88  90
# ↑↑↑ Expected output ↑↑↑ #

# Fitted minimum cell threshold per branch.
base::print(interest_cell_type_group[[1]][["n_min"]]) # Example output: 5 4

# Pseudotime value at the endpoint of each group.
base::print(interest_cell_type_group[[1]][["Branches_T_k"]])
# ↓↓↓ Expected output ↓↓↓ #
# [[1]]
# [1]  1.834408  3.535982  4.685209  7.771460  9.469274 11.265137 11.886710 15.484866
# [[2]]
# [1]  2.504493  4.344136  7.384443 11.549756 13.594401 13.760324 13.944713 14.597503
# ↑↑↑ Expected output ↑↑↑ #

# Length of pseudotime interval for each group.
base::print(interest_cell_type_group[[1]][["Branches_Tslot_k"]])
# ↓↓↓ Expected output ↓↓↓ #
# [[1]]
# [1] 1.8344078 1.7015745 1.1492267 3.0862507 1.6978140 1.7958636 0.6215733 3.5981557
# [[2]]
# [1] 2.5044931 1.8396424 3.0403072 4.1653131 2.0446453 0.1659231 0.1843885 0.6527900
# ↑↑↑ Expected output ↑↑↑ #

# Aggregated mean non-zero expression for TFs (genes × groups).
base::print(utils::head(interest_cell_type_group[[1]][["Branches_TFE_T"]][[1]], n = 3)) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#             1        2        3        4        5        6        7        8
# ATF3 3.072106 3.067088 3.019299 2.958824 3.003423 2.910599 2.963062 2.482744
# ATF7 2.825813 2.735095 2.738140 2.666481 2.535333 2.601107 2.483336 2.385111
# BATF 2.580502 2.555763 2.515973 2.456288 2.688231 2.354480 2.179854 2.138373
# ↑↑↑ Expected output ↑↑↑ #

base::print(utils::head(interest_cell_type_group[[1]][["Branches_TFE_T"]][[2]], n = 3)) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#              1        2        3        4        5        6        7        8
# ATF3  3.071233 3.081425 2.916346 2.956322 2.424782 2.506254 2.397610 2.378085
# ATF7  2.810332 2.744660 2.673462 2.555244 2.146415 2.264963 2.265312 2.170351
# BACH2 2.391223 2.426087 2.454012 2.775427 3.025674 2.789505 3.152494 3.361768
# ↑↑↑ Expected output ↑↑↑ #

# Aggregated mean non-zero activity for TGs (genes × groups).
base::print(utils::head(interest_cell_type_group[[1]][["Branches_TGA_T"]][[1]], n = 3)) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#               1        2        3        4        5        6        7        8
# AAK1   2.120084 2.081187 2.160225 2.150473 2.179560 2.153994 2.114437 2.206777
# ABCA1  1.993366 2.034464 2.005056 2.010421 1.966451 2.044633 2.087603 1.998765
# ABCB10 1.564771 1.584545 1.594070 1.692407 1.657556 1.721049 1.698226 1.710255
# ↑↑↑ Expected output ↑↑↑ #

base::print(utils::head(interest_cell_type_group[[1]][["Branches_TGA_T"]][[2]], n = 3)) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#               1        2        3        4        5        6        7        8
# AAK1   2.104856 2.114029 2.164139 2.196094 2.153543 2.068728 2.149999 2.143485
# ABCA1  1.994128 2.040264 2.010113 2.014280 1.601724 1.484141 1.535953 1.554293
# ABCB10 1.576125 1.560878 1.676047 1.666620 1.498931 1.503316 1.438212 1.505404
# ↑↑↑ Expected output ↑↑↑ #

# Aggregated mean non-zero expression for TGs (genes × groups).
base::print(utils::head(interest_cell_type_group[[1]][["Branches_TGE_T"]][[1]], n = 3)) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#               1        2        3        4        5        6        7        8
# AAK1   3.093741 2.998561 2.979363 2.906093 2.864815 2.801649 2.837255 2.489660
# ABCA1  2.820616 2.851578 2.698928 2.653027 2.628227 2.644240 2.522406 2.456225
# ABCB10 2.581650 2.606795 2.587192 2.521982 2.456509 2.432524 2.418175 2.212099
# ↑↑↑ Expected output ↑↑↑ #

base::print(utils::head(interest_cell_type_group[[1]][["Branches_TGE_T"]][[2]], n = 3)) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#               1        2        3        4        5        6        7        8
# AAK1   3.084084 2.988308 2.922666 2.859432 2.801591 2.956251 3.040081 3.063276
# ABCA1  2.827137 2.817448 2.611040 2.608745 2.424004 1.972459 2.170264 2.057759
# ABCB10 2.587728 2.616320 2.530707 2.459481 1.957653 2.056608 1.912462 1.963751
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(interest_cell_type_group, file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")


### Notes
# Verify group count per branch and minimum non-zero counts per group.
# Adjust fitting parameters if grouping results are unsatisfactory.

# Troubleshooting for too few groups or failed grouping:
# 1. Relax the minimum cell requirement by adjusting points_x_for_fitting_nc_and_nmin and points_y_for_fitting_nc_and_nmin.
# 2. Tune parameters in get_genes_pseudotime_info:
#    - Increase non-zero ratio weights to reduce gene sparsity.
#    - Raise alpha_TF and alpha_TG to reduce the number of selected genes.

# Troubleshooting for too many groups:
# - Apply the reverse adjustments as above.

# Recommendations:
# - A minimum of 8 groups is recommended to avoid overfitting in downstream model training (23 parameters total).
# - 5–7 groups are acceptable if 8 cannot be achieved; the model includes regularization and structural constraints.
# - Avoid more than ~20 groups, as small subgroups reduce statistical power and increase computational cost.





############################################################
############################################################
###                                                      ###
###   Part 3: Constructing Initial GRNs from TFBS PWMs   ###
###                                                      ###
############################################################
############################################################




#################################################################
#                                                               #
# 3.1 Constructing an Initial GRN for Comprehensive TG-TF Pairs #
#                                                               #
#################################################################



##########################################################################################################
# 3.1.1  Constructing an Initial GRN for Comprehensive TG-TF Pairs Based on Base PWM of TF Binding Sites #
##########################################################################################################
### Function
# Scans promoter-associated accessible peaks with JASPAR2024 TF binding site PWMs for each cell type.
# Constructs an initial gene regulatory network (iGRN) matrix.
# Rows are target genes; columns are transcription factors; values are normalized PWM matching scores.

# See function documentation for input, logic and output details.
?get_iGRN_by_TFBS_pwm_by_JASPAR2024


### Data import
setwd(Working_directory)

# Output list from get_genes_pseudotime_info. Provides target TF and TG sets.
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")

# Output list from identify_TGs or filter_peak_annotations.
# Contains peak_anno_promoter element with promoter-associated ATAC peak annotations.
results_identify_TGs = qs::qread("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.qsave")

# BSgenome reference genome object for sequence extraction.
genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
# genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 # Use for hg19-aligned datasets.

# Species name for JASPAR2024 database query. See ?get_TFs_from_priors for details.
species = "Homo sapiens"

# JASPAR matrix collections to use. See ?get_TFs_from_priors for details.
collection = c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED")

# Whether to include predicted TFBS sequences in the output CSV.
# Defaults to FALSE to reduce runtime and file size, exporting only gene, TF, max score and peak name.
output_predicted_TFBS = FALSE

# Number of cores for parallel computation via mclapply.
# Parallel processing only works on Unix-like systems; set to 1 on Windows.
# Adjust based on available RAM; reduce core count if memory is exhausted.
ncores = parallel::detectCores() - 1 # Linux setting
# ncores = 1 # Windows setting


### Hyperparameter configuration
# Minimum score threshold for PWM matching, passed to Biostrings::matchPWM.
# Higher thresholds produce sparser networks with higher specificity.
# Use "0%" or low values to retain potential TF–TG regulatory interactions.
min_score_for_matchPWM = "0%"


### Run get_iGRN_by_TFBS_pwm_by_JASPAR2024 with runtime and peak memory monitoring.
library(dplyr)
monitor = peakRAM::peakRAM({
  interest_cell_type_iGRN_all = get_iGRN_by_TFBS_pwm_by_JASPAR2024(
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    results_identify_TGs = results_identify_TGs, genome = genome,
    min_score_for_matchPWM = min_score_for_matchPWM, species = species, collection =collection, ncores = ncores
  )
})
message(base::ncol(interest_cell_type_iGRN_all[[1]]), " TFs, ", base::nrow(interest_cell_type_iGRN_all[[1]]), " TGs, ",
        sum(interest_cell_type_iGRN_all[[1]] > 0), " pairs, ", print_run_stats(monitor))
# Example output: 75 TFs, 1379 TGs, 103425 pairs, running: 0 hours 12 mins 20 secs, peak RAM: 0 GB 128 MB


### Inspect results.
# Named list with one element per cell type. Each element is a TG × TF association matrix.
# Values are normalized PWM matching scores ranging from 0 to ~0.95.
# A value of 0 indicates no predicted TF–TG regulatory interaction.
base::print(interest_cell_type_iGRN_all[[1]][1:3, 1:3])
# ↓↓↓ Expected output ↓↓↓ #
#             ATF3      ATF7     BACH2
# AAK1   0.4869373 0.5010867 0.6282499
# ABCA1  0.5550473 0.5685793 0.6187872
# ABCB10 0.4778124 0.4657027 0.5650279
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(interest_cell_type_iGRN_all, file = "./3 get iGRN/interest_cell_type_iGRN_all.qsave")



#########################################################
# 3.1.2 Filtering  the Initial GRN by ChIP-Atlas Priors #
#########################################################
### Function
# Refines cell-type-specific iGRNs using preprocessed ChIP-Atlas ChIP-seq prior data.
# Performs prior-weighted edge correction, dual-threshold QC with target gene rescue, and zero-preserving transformation.
# Balances network sparsity and numerical stability for downstream dynamic modeling.

# See function documentation for input, logic and output details.
?filter_iGRN_by_ChIP_Atlas_priors


### Data import
setwd(Working_directory)

# Output list from get_iGRN_by_TFBS_pwm_by_JASPAR2024. Each element is a TG × TF association matrix.
interest_cell_type_iGRN_all = qs::qread("./3 get iGRN/interest_cell_type_iGRN_all.qsave")

# Calculate peak-to-TG ratio for network size calibration.
results_identify_TGs = qs::qread("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.qsave")
peak_anno_promoter = results_identify_TGs[["peak_anno_promoter"]]; rm(results_identify_TGs)
peak_anno_promoter = peak_anno_promoter[peak_anno_promoter$SYMBOL %in% rownames(interest_cell_type_iGRN_all[[1]]), ]
peak_TG_ratio = length(base::unique(peak_anno_promoter[["peak"]])) / length(base::unique(peak_anno_promoter[["SYMBOL"]]))
message("peak_TG_ratio = ", round(peak_TG_ratio, digits = 4))
# Example output: peak_TG_ratio = 1.0015


### Hyperparameter configuration
# Expected number of peaks per TF for controlling iGRN network scale.
expected_peak_number_per_TF = 50

# Global quantile threshold for edge retention.
# Edges with corrected scores in the top proportion are retained. Defaults to 0.15 (top 15%).
theta_i_top_quantile = expected_peak_number_per_TF / peak_TG_ratio / base::nrow(interest_cell_type_iGRN_all[[1]])
base::print(theta_i_top_quantile) # Example output: 0.03620565

# Per-TF top-N threshold. Top N targets per TF are always retained regardless of global scores. Defaults to 50.
theta_i_top_n = round(expected_peak_number_per_TF/peak_TG_ratio)
base::print(theta_i_top_n) # Example output: 50

# Absolute score floor. Corrected values below this threshold are set to zero before re-normalization. Defaults to 0.05.
theta_i_absolute_min = 0.05

# Minimum number of regulators guaranteed per target gene after QC.
# TGs with fewer non-zero edges will have their top pre-QC edges restored. Set to 0 to disable rescue. Defaults to 1.
TG_rescue_number = 5


### Run filter_iGRN_by_ChIP_Atlas_priors with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_iGRN = filter_iGRN_by_ChIP_Atlas_priors(
    interest_cell_type_iGRN = interest_cell_type_iGRN_all,
    theta_i_top_quantile = theta_i_top_quantile, theta_i_top_n = theta_i_top_n, TG_rescue_number = TG_rescue_number
  )
})
message("top_q=", round(theta_i_top_quantile, digits = 4),
        ", top_n=", theta_i_top_n, ", ",
        base::ncol(interest_cell_type_iGRN[[1]]), " TFs, ", base::nrow(interest_cell_type_iGRN[[1]]), " TGs, ",
        sum(interest_cell_type_iGRN[[1]] > 0), " pairs, ",
        round(sum(interest_cell_type_iGRN[[1]] > 0) / base::ncol(interest_cell_type_iGRN[[1]]), digits = 2), " TGs per TF, ", print_run_stats(monitor))
# Example output: top_q=0.0362, top_n=50, 75 TFs, 1358 TGs, 7628 pairs, 101.71 TGs per TF, running: 0 hours 0 mins 6 secs, peak RAM: 0 GB 269 MB


### Inspect results.
# Named list of refined iGRN matrices with the same structure as input.
# Matrices are prior-weighted, quality-controlled and value-transformed.
# Non-zero edges range from 0.1 to 1.0; zero edges are preserved to maintain sparsity.
# All-zero rows (TGs) and columns (TFs) are automatically removed.
base::print(interest_cell_type_iGRN[["Mono"]][1:3, 1:3])
# ↓↓↓ Expected output ↓↓↓ #
#        ATF3 ATF7 BACH2
# AAK1      0    0     0
# ABCA1     0    0     0
# ABCB10    0    0     0
# ↑↑↑ Expected output ↑↑↑ #

# Distribution of TG counts per TF.
res = base::t(base::apply(interest_cell_type_iGRN[[1]], 2, function(x) {
  nz = x[x != 0]
  c(nonzero_freq = length(nz) / length(x),
    nonzero_count = length(nz),
    nonzero_min = if(length(nz)) min(nz) else NA)
}))
summary(res[, "nonzero_count"])
# ↓↓↓ Expected output ↓↓↓ #
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 50.0    61.0    87.0   101.7   127.0   280.0
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(interest_cell_type_iGRN, file = "./3 get iGRN/interest_cell_type_iGRN.qsave")




############################################
#                                          #
# 3.2 Setting the Regulatory Threshold tau #
#                                          #
############################################
### Function
# Computes a cell-type-specific binarization threshold (τ) for each initial GRN.
# Derived from the minimum non-zero interaction strength with a fixed lower bound to avoid over-sensitivity.
# Pairs with scores above τ are considered regulatory; pairs below are not.

# See function documentation for input, logic and output details.
?get_cell_type_tau


### Data import
setwd(Working_directory)

# Output list from get_iGRN_by_TFBS_pwm_by_JASPAR2024 or filter_iGRN_by_ChIP_Atlas_priors.
# Each element is a TG × TF association matrix.
interest_cell_type_iGRN = qs::qread("./3 get iGRN/interest_cell_type_iGRN.qsave")


### Run get_cell_type_tau with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_tau = get_cell_type_tau(
    interest_cell_type_iGRN = interest_cell_type_iGRN
  )
})
message("tau is ", round(interest_cell_type_tau[1], digits = 4), ", ", print_run_stats(monitor))
# Example output: tau is 0.1433, running: 0 hours 0 mins 1 secs, peak RAM: 0 GB 3 MB


### Save output.
qs::qsave(interest_cell_type_tau, file = "./3 get iGRN/interest_cell_type_tau.qsave")




#############################################
#                                           #
# 3.3 Constructing the iGRN for Each Branch #
#                                           #
#############################################
### Function
# Extracts branch-specific initial GRNs from cell-type-level iGRNs.
# Creates subnetworks for each branch by filtering the global matrix to branch-specific TFs and TGs.

# See function documentation for input, logic and output details.
?get_branch_iGRN


### Data import
setwd(Working_directory)

# Output from get_iGRN_by_TFBS_pwm_by_JASPAR2024 or filter_iGRN_by_ChIP_Atlas_priors.
interest_cell_type_iGRN = qs::qread("./3 get iGRN/interest_cell_type_iGRN.qsave")

# Output from the cell_grouping function.
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")

# Number of cores for parallel computation via mclapply.
# Parallel processing only works on Unix-like systems; set to 1 on Windows.
# Adjust based on available RAM; reduce core count if memory is exhausted.
ncores = parallel::detectCores() - 1 # Linux setting
# ncores = 1 # Windows setting


### Run get_branch_iGRN with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_iGRN = get_branch_iGRN(
    interest_cell_type_iGRN = interest_cell_type_iGRN,
    interest_cell_type_group = interest_cell_type_group,
    ncores = ncores
  )
})
message(base::paste(base::unlist(base::lapply(interest_cell_type_branch_iGRN[[1]], base::ncol)), collapse = "/"), " TFs, ",
        base::paste(base::unlist(base::lapply(interest_cell_type_branch_iGRN[[1]], base::nrow)), collapse = "/"), " TGs, ",
        base::paste(base::sapply(interest_cell_type_branch_iGRN[[1]], function(df) sum(df > 0, na.rm = TRUE)), collapse = "/"), " pairs, ", print_run_stats(monitor))
# Example output: 70/72 TFs, 1307/1325 TGs, 6764/7153 pairs, running: 0 hours 0 mins 19 secs, peak RAM: 0 GB 868 MB


### Inspect results.
# Nested list structure:
# - Level 1: named by cell type
# - Level 2: list of branch-specific iGRNs per cell type
# Each branch iGRN is a TG × TF data frame with normalized regulatory strengths.
# All-zero rows and columns are automatically removed.
base::print(interest_cell_type_branch_iGRN[[1]][[1]][1:3, 1:3]) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#        ATF3 ATF7 BATF
# AAK1      0    0    0
# ABCA1     0    0    0
# ABCB10    0    0    0
# ↑↑↑ Expected output ↑↑↑ #

base::print(interest_cell_type_branch_iGRN[[1]][[2]][1:4, 1:4]) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#        ATF3 ATF7 BACH2
# AAK1      0    0     0
# ABCA1     0    0     0
# ABCB10    0    0     0
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(interest_cell_type_branch_iGRN, file = "./3 get iGRN/interest_cell_type_branch_iGRN.qsave")





##########################################################################################
##########################################################################################
###                                                                                    ###
###   Part 4: Developing DERI-Based GRN Regression Predictive Models for Each Branch   ###
###                                                                                    ###
##########################################################################################
##########################################################################################




##################################
#                                #
# 4.1 Partitioning Training Sets #
#                                #
##################################
### Function
# Extracts high-confidence TF-TG interactions from branch-level initial GRNs to construct training sets.
# Implements multi-layer threshold filtering with optional TF coverage guarantee.
# Appends pseudotime-resolved expression and activity features for each branch.

# See function documentation for input, logic and output details.
?select_training_set


### Data import
setwd(Working_directory)

# Branch-level iGRN matrices for each cell type, output from get_branch_iGRN.
interest_cell_type_branch_iGRN = qs::qread("./3 get iGRN/interest_cell_type_branch_iGRN.qsave")

# Cell-type-specific binarization thresholds, output from get_cell_type_tau.
interest_cell_type_tau = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")

# Cell grouping results with aggregated expression/activity matrices, output from cell_grouping.
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")

# Number of cores for parallel computation via mclapply.
# Parallel processing only works on Unix-like systems; set to 1 on Windows.
# Adjust based on available RAM; reduce core count if memory is exhausted.
ncores = parallel::detectCores() - 1 # Linux setting
# ncores = 1 # Windows setting


### Hyperparameter configuration
# Number of top interactions retained per TF and per TG; their union forms the initial training set. Defaults to 1.
# Set to 1 in most cases to avoid excessive computational load.
top_n_for_max_theta_in_each_regulon = 1

# Quantile retention ratio of regulatory strength within the initial training set (0–1). Defaults to 1 (no filtering).
max_theta_retention_rate_in_tops = 1

# Quantile retention ratio of regulatory strength across all candidate interactions (0–1). Defaults to 1 (no filtering).
max_theta_retention_rate_in_all = 1

# Absolute number of top-ranked interactions by regulatory strength. Defaults to 1000.
retention_number_of_tops = 9999999 # large number: no filtering

# Whether to perform mutual-top refinement on the training set.
# Defaults to TRUE: further refines the training set by requiring mutual top membership (TF in top regulons of TG, and TG in top regulons of TF) and retains a top ratio.
refine = TRUE

# Whether to ensure every TF has at least one interaction retained after refinement. Only works when refine = TRUE. Defaults to TRUE.
keep_all_TF = TRUE

# Quantile retention ratio of regulatory strength after the refinement step (0–1). Only works when refine = TRUE. Defaults to 1 (no filtering).
top_ratio_in_refine = 1


### Run select_training_set with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_training_set = select_training_set(
    interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
    interest_cell_type_tau  = interest_cell_type_tau, interest_cell_type_group = interest_cell_type_group,
    top_n_for_max_theta_in_each_regulon = top_n_for_max_theta_in_each_regulon,
    max_theta_retention_rate_in_tops = max_theta_retention_rate_in_tops,
    max_theta_retention_rate_in_all = max_theta_retention_rate_in_all,
    retention_number_of_tops = retention_number_of_tops,
    refine = refine, keep_all_TF = keep_all_TF, top_ratio_in_refine = top_ratio_in_refine, ncores = ncores
  )
})
message(base::paste(base::lapply(interest_cell_type_branch_training_set[[1]], base::nrow), collapse = ","), " training TF-TG pairs, ", print_run_stats(monitor))
# Example output: 63,72,75 training TF-TG pairs, running: 0 hours 0 mins 10 secs, peak RAM: 0 GB 660 MB


### Inspect results.
# Training set data frame with TF-TG pairs as rows.
# Columns include TF, TG, theta_i, and time-point-specific TFE_K, TGA_K, TGE_K features.
base::print(interest_cell_type_branch_training_set[[1]][["branch1"]][1:3, 1:10]) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#                 TF    TG   theta_i    TFE_1    TGA_1    TGE_1    TFE_2    TGA_2    TGE_2    TFE_3
# ATF3_to_RPS29 ATF3 RPS29 0.9048870 3.072106 1.493751 3.489330 3.067088 1.527613 3.466311 3.019299
# ATF7_to_PUM3  ATF7  PUM3 0.9062229 2.825813 1.411290 2.722523 2.735095 1.453676 2.730181 2.738140
# BATF_to_ATF3  BATF  ATF3 0.9075643 2.580502 2.601344 3.072106 2.555763 2.681091 3.067088 2.515973
# ↑↑↑ Expected output ↑↑↑ #

base::print(interest_cell_type_branch_training_set[[1]][["branch2"]][1:3, 1:10]) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#                 TF    TG   theta_i    TFE_1    TGA_1    TGE_1    TFE_2    TGA_2    TGE_2    TFE_3
# ATF3_to_RPS29 ATF3 RPS29 0.9048870 3.071233 1.498308 3.480479 3.081425 1.534615 3.480738 2.916346
# ATF7_to_PUM3  ATF7  PUM3 0.9062229 2.810332 1.416562 2.724441 2.744660 1.465871 2.732649 2.673462
# BATF_to_ATF3  BATF  ATF3 0.9075643 2.583434 2.618875 3.071233 2.593358 2.691502 3.081425 2.434905
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(interest_cell_type_branch_training_set, file = "./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.qsave")


### Tips
# Of course, users can manually add some training sets and remove any training sets themselves.




##############################################
#                                            #
# 4.2 Build Predictive Model for Each Branch #
#                                            #
##############################################



###############################################################################
# 4.2.1 Setting Initial Values, Lower Bounds, and Upper Bounds for Parameters #
###############################################################################
### Function
# Generates initialized parameter structures for the scDERI GRN model for each TF-TG pair.
# Assigns initial values, upper bounds and lower bounds for subsequent optimization.
# Adaptively initializes baseline expression parameters from training set means.


# See function documentation for input, logic and output details.
?set_init_params


### Data import
setwd(Working_directory)

# Branch-level training sets with expression and activity features, output from select_training_set.
interest_cell_type_branch_training_set = qs::qread("./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.qsave")


### Run set_init_params with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_init_params = set_init_params(
    interest_cell_type_branch_training_set = interest_cell_type_branch_training_set
  )
})
message(print_run_stats(monitor))
# Example output: running: 0 hours 0 mins 2 secs, peak RAM: 0 GB 14 MB


### Inspect results.
# List model parameter names.
base::print(base::unlist(base::lapply(interest_cell_type_branch_init_params[[1]][["branch1"]][[1]], base::names), use.names = FALSE))
# ↓↓↓ Expected output ↓↓↓ #
#  [1] "r2"     "alpha1" "alpha2" "beta1"  "beta2"  "beta3"  "s1"     "s2"     "u11"    "u21"    "u311"   "u321"   "u331"
# [14] "u12"    "u22"    "u312"   "u322"   "u332"   "E~_TG"  "E~_TF"  "v1"     "v2"     "v3"
# ↑↑↑ Expected output ↑↑↑ #

# Total number of model parameters.
base::print(length(base::unlist(base::lapply(interest_cell_type_branch_init_params[[1]][["branch1"]][[1]], base::names), use.names = FALSE))) # 23

# Preview parameter initialization values and bounds.
base::print(base::unlist(interest_cell_type_branch_init_params[[1]][[1]][[1]]))
# ↓↓↓ Expected output ↓↓↓ #
#  r.r2.initialization                  r.r2.upper                  r.r2.lower alpha.alpha1.initialization          alpha.alpha1.upper          alpha.alpha1.lower alpha.alpha2.initialization
#             1.000000                    1.113383                    0.370506                    1.000000                   50.000000                    0.100000                    1.000000
#   alpha.alpha2.upper          alpha.alpha2.lower   beta.beta1.initialization            beta.beta1.upper            beta.beta1.lower   beta.beta2.initialization            beta.beta2.upper
#            50.000000                    0.100000                    0.000000                   10.000000                  -10.000000                    0.000000                   10.000000
#     beta.beta2.lower   beta.beta3.initialization            beta.beta3.upper            beta.beta3.lower         s.s1.initialization                  s.s1.upper                  s.s1.lower
#           -10.000000                    0.000000                   10.000000                  -10.000000                    1.000000                    4.900000                    0.100000
#  s.s2.initialization                  s.s2.upper                  s.s2.lower        u.u11.initialization                 u.u11.upper                 u.u11.lower        u.u21.initialization
#             1.000000                    4.900000                    0.100000                    1.000000               777600.000000                    0.100000                    1.000000
#          u.u21.upper                 u.u21.lower       u.u311.initialization                u.u311.upper                u.u311.lower       u.u321.initialization                u.u321.upper
#        777600.000000                    0.100000                    1.000000               777600.000000                    0.100000                    1.000000               777600.000000
#         u.u321.lower       u.u331.initialization                u.u331.upper                u.u331.lower        u.u12.initialization                 u.u12.upper                 u.u12.lower
#             0.100000                    1.000000               777600.000000                    0.100000                    1.000000                    4.900000                    0.100000
# u.u22.initialization                 u.u22.upper                 u.u22.lower       u.u312.initialization                u.u312.upper                u.u312.lower       u.u322.initialization
#             1.000000                    4.900000                    0.100000                    1.000000                    4.900000                    0.100000                    1.000000
#         u.u322.upper                u.u322.lower       u.u332.initialization                u.u332.upper                u.u332.lower     E~.E~_TG.initialization              E~.E~_TG.upper
#             4.900000                    0.100000                    1.000000                    4.900000                    0.100000                    3.507137                    9.000000
#       E~.E~_TG.lower     E~.E~_TF.initialization              E~.E~_TF.upper              E~.E~_TF.lower         v.v1.initialization                  v.v1.upper                  v.v1.lower
#             1.000000                    2.934643                    9.000000                    1.000000                    1.000000                  200.000000                    0.100000
#  v.v2.initialization                  v.v2.upper                  v.v2.lower         v.v3.initialization                  v.v3.upper                  v.v3.lower
#             1.000000                  200.000000                    0.100000                    1.000000                  200.000000                    0.100000
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(interest_cell_type_branch_init_params, file = "./4.2 BUild Prediction Model/interest_cell_type_branch_init_params.qsave")



#########################
# 4.2.2 Training models #
#########################
### Function
# Trains GRN models for each TF-TG pair using a hybrid optimization strategy.
# Combines genetic algorithm global search and gradient ascent local refinement.
# Supports multiple early stopping criteria and regulon-level parallel computing.

# See function documentation for input, logic and output details.
?model_train


### Data import
setwd(Working_directory)

# Branch-level training data with expression and activity features, output from select_training_set.
interest_cell_type_branch_training_set = qs::qread("./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.qsave")

# Initialized parameters and bounds for each TF-TG pair, output from set_init_params.
interest_cell_type_branch_init_params = qs::qread("./4.2 BUild Prediction Model/interest_cell_type_branch_init_params.qsave")

# Pseudotime interval lengths per branch, output from cell_grouping.
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")

# Cell-type-specific binarization thresholds, output from get_cell_type_tau.
interest_cell_type_tau = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")

# Cell counts per branch for weighted result integration, output from get_genes_pseudotime_info.
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")

# Number of cores for parallel computation via mclapply.
# Parallel processing only works on Unix-like systems; set to 1 on Windows.
# Adjust based on available RAM; reduce core count if memory is exhausted.
ncores = parallel::detectCores() - 1 # Linux setting
# ncores = 1 # Windows setting


### Hyperparameter configuration
# Maximum number of training epochs. Defaults to 1000.
# Use a moderate value first to estimate required epochs.
# Training is primarily controlled by early stopping rather than max epochs.
max_epochs = 300

# Number of consecutive epochs without loss improvement to trigger early stopping. Defaults to 10.
# Higher values increase runtime and may improve performance but raise overfitting risk.
early_stop = 10

# Convergence threshold for predicted vs. observed regulatory strength. Defaults to 1e-3.
eps_theta = 1e-3

# Threshold for negligible loss improvement. Defaults to 1e-5.
eps_loss = 1e-5

# Genetic algorithm parameters
# Population size for genetic algorithm. Defaults to 50. See popSize in GA::ga.
# Higher values increase runtime but improve local optima escape.
popSize_ga = 50

# Maximum generations per genetic algorithm run. Defaults to 150. See maxiter in GA::ga.
# Higher values increase runtime but improve local optima escape.
maxiter_ga = 150

# Crossover probability for genetic algorithm (0–1). Defaults to 0.9. See pcrossover in GA::ga.
# Modification is not recommended; GA serves as a corrective mechanism for gradient descent.
pcrossover_ga = 0.5

# Mutation probability for genetic algorithm (0–1). Defaults to 0.5. See pmutation in GA::ga.
# Modification is not recommended; reasoning matches the crossover parameter above.
pmutation_ga = 0.5

# Gradient ascent parameters
# Number of gradient ascent steps per epoch. Defaults to 30.
# Higher values increase runtime with diminishing marginal returns.
iterations_grad = 30


### Run model_train with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_model_train_all = model_train(
    interest_cell_type_branch_training_set = interest_cell_type_branch_training_set,
    interest_cell_type_branch_init_params = interest_cell_type_branch_init_params,
    interest_cell_type_group = interest_cell_type_group, interest_cell_type_tau = interest_cell_type_tau,
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    max_epochs = max_epochs, early_stop = early_stop, eps_theta = eps_theta, eps_loss = eps_loss,
    popSize_ga = popSize_ga, maxiter_ga = maxiter_ga, pcrossover_ga = pcrossover_ga, pmutation_ga = pmutation_ga,
    iterations_grad = iterations_grad,
    ncores = ncores
  )
})
message(print_run_stats(monitor))
# Example output: running: 0 hours 46 mins 19 secs, peak RAM: 0 GB 32 MB


### Save output.
qs::qsave(interest_cell_type_branch_model_train_all, file = "./4.2 BUild Prediction Model/interest_cell_type_branch_model_train_all.qsave")


### Remove low-quality training sets
interest_cell_type_branch_model_train = interest_cell_type_branch_model_train_all
for (cell_type in base::names(interest_cell_type_branch_model_train)) {
  df = interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]]
  # bad_flag = with(df, theta_i < 0.011 | theta_i > 0.989 | theta_difference >= 0.1 | fit_loss >= 0.2) # Set threshold by yourself
  bad_flag = with(df, theta_i < 0.011) # Set threshold by yourself
  rownames_filtered = base::rownames(df)[bad_flag]
  interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]] = df[!bad_flag, ]
  for (n in 1:(length(interest_cell_type_branch_model_train[[cell_type]]) - 1)) {
    lst = interest_cell_type_branch_model_train[[cell_type]][[paste0("branch", n)]]
    interest_cell_type_branch_model_train[[cell_type]][[paste0("branch", n)]] = lst[!base::names(lst) %in% rownames_filtered]
  }
}


### Save output.
qs::qsave(interest_cell_type_branch_model_train, file = "./4.2 BUild Prediction Model/interest_cell_type_branch_model_train.qsave")


### Inspect results.
# Summarize parameter distribution across all trained models.
interest_cell_type_branch_training_set = qs::qread("./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.qsave")
pa_init = round(base::unlist(interest_cell_type_branch_init_params[[1]][["branch1"]][[1]]), digits = 3)
res_list = list()
for (pa in base::names(interest_cell_type_branch_model_train[[1]][["branch1"]][[1]][["best_params"]])) {
  for (n in 1:(length(interest_cell_type_branch_model_train[[1]]) - 1)) {
    pa_values = round(as.numeric(base::sapply(interest_cell_type_branch_model_train[[1]][[paste0("branch", n)]], \(x) x[["best_params"]][[pa]])), digits = 3)
    pa_quants = round(stats::quantile(pa_values), digits = 3)
    row_name = paste0(pa, "-Branch", n)
    res_list[[row_name]] = c(
      `Lower Bound`    = as.numeric(pa_init[paste0(pa, ".lower")]),
      Min              = as.numeric(pa_quants["0%"]),
      Q1               = as.numeric(pa_quants["25%"]),
      Median           = as.numeric(pa_quants["50%"]),
      Q3               = as.numeric(pa_quants["75%"]),
      Max              = as.numeric(pa_quants["100%"]),
      `Upper Bound`    = as.numeric(pa_init[paste0(pa, ".upper")]),
      `Range Coverage` = as.numeric(round((pa_quants["100%"] - pa_quants["0%"]) / (pa_init[paste0(pa, ".upper")] - pa_init[paste0(pa, ".lower")]), digits = 3)),
      Var              = round(stats::var(pa_values), digits = 3)
    )
  }
}
result_df = base::as.data.frame(base::do.call(base::rbind, res_list))
options(scipen = 999)
base::print(result_df)
# ↓↓↓ Expected output ↓↓↓ #
#                      Lower Bound    Min     Q1 Median     Q3    Max Upper Bound Range Coverage    Var
# r.r2-Branch1               0.371  0.506  0.732  0.801  0.843  0.996       1.113          0.660  0.011
# r.r2-Branch2               0.371  0.546  0.751  0.798  0.865  0.980       1.113          0.585  0.008
# alpha.alpha1-Branch1       0.100  0.366  0.460  0.622  0.990  1.268      50.000          0.018  0.068
# alpha.alpha1-Branch2       0.100  0.222  0.701  0.960  1.016  1.599      50.000          0.028  0.081
# alpha.alpha2-Branch1       0.100  0.572  0.782  0.832  0.849  3.582      50.000          0.060  0.135
# alpha.alpha2-Branch2       0.100  0.350  0.621  0.775  0.839  1.380      50.000          0.021  0.037
# beta.beta1-Branch1       -10.000 -0.119 -0.051 -0.045 -0.036  0.056      10.000          0.009  0.001
# beta.beta1-Branch2       -10.000 -0.085 -0.059 -0.051 -0.037  0.092      10.000          0.009  0.001
# beta.beta2-Branch1       -10.000 -0.011  0.006  0.012  0.016  0.038      10.000          0.002  0.000
# beta.beta2-Branch2       -10.000 -0.055 -0.019 -0.005  0.006  0.064      10.000          0.006  0.000
# beta.beta3-Branch1       -10.000 -0.032 -0.014 -0.001  0.010  0.037      10.000          0.003  0.000
# beta.beta3-Branch2       -10.000 -0.041 -0.021 -0.012 -0.005  0.040      10.000          0.004  0.000
# s.s1-Branch1               0.100  0.393  0.500  0.680  0.996  2.353       4.900          0.408  0.225
# s.s1-Branch2               0.100  0.402  0.814  1.360  1.946  3.006       4.900          0.542  0.491
# s.s2-Branch1               0.100  0.449  0.611  0.778  0.853  1.982       4.900          0.319  0.070
# s.s2-Branch2               0.100  0.335  0.453  0.568  0.818  1.699       4.900          0.284  0.058
# u.u11-Branch1              0.100  0.876  1.000  1.009  1.032  1.079  777600.000          0.000  0.001
# u.u11-Branch2              0.100  0.945  0.999  1.013  1.048 24.559  777600.000          0.000  8.160
# u.u21-Branch1              0.100  0.978  1.003  1.011  1.032  1.084  777600.000          0.000  0.001
# u.u21-Branch2              0.100  0.997  1.011  1.035  1.074  1.187  777600.000          0.000  0.002
# u.u311-Branch1             0.100  0.771  0.952  0.985  0.996  1.447  777600.000          0.000  0.006
# u.u311-Branch2             0.100  0.299  0.915  0.981  0.997  1.079  777600.000          0.000  0.016
# u.u321-Branch1             0.100  0.958  1.003  1.020  1.060 12.210  777600.000          0.000  2.360
# u.u321-Branch2             0.100  0.911  1.006  1.023  1.067 18.206  777600.000          0.000  4.525
# u.u331-Branch1             0.100  0.989  1.003  1.017  1.046  3.461  777600.000          0.000  0.095
# u.u331-Branch2             0.100  0.990  1.006  1.030  1.049  6.529  777600.000          0.000  0.675
# u.u12-Branch1              0.100  0.384  0.567  0.624  1.046  2.067       4.900          0.351  0.130
# u.u12-Branch2              0.100  0.378  0.838  1.026  1.114  4.013       4.900          0.757  0.338
# u.u22-Branch1              0.100  0.377  0.492  0.543  0.720  1.223       4.900          0.176  0.047
# u.u22-Branch2              0.100  0.328  0.430  0.474  0.619  0.994       4.900          0.139  0.027
# u.u312-Branch1             0.100  0.848  1.333  1.781  2.282  3.049       4.900          0.459  0.319
# u.u312-Branch2             0.100  0.524  1.066  1.563  1.970  2.808       4.900          0.476  0.330
# u.u322-Branch1             0.100  0.465  0.694  0.843  0.977  1.515       4.900          0.219  0.042
# u.u322-Branch2             0.100  0.417  0.817  0.936  1.050  2.311       4.900          0.395  0.088
# u.u332-Branch1             0.100  0.360  0.677  0.781  0.902  1.645       4.900          0.268  0.053
# u.u332-Branch2             0.100  0.398  0.762  0.882  0.963  1.408       4.900          0.210  0.044
# E~.E~_TG-Branch1           1.000  2.398  2.486  2.598  2.745  3.506       9.000          0.138  0.069
# E~.E~_TG-Branch2           1.000  2.206  2.293  2.369  2.576  3.785       9.000          0.197  0.102
# E~.E~_TF-Branch1           1.000  2.302  2.490  2.617  2.766  3.747       9.000          0.181  0.098
# E~.E~_TF-Branch2           1.000  2.213  2.339  2.458  2.614  3.700       9.000          0.186  0.092
# v.v1-Branch1               0.100  0.979  1.175  1.296  1.439 39.002     200.000          0.190 22.598
# v.v1-Branch2               0.100  0.992  1.139  1.327  1.600 10.325     200.000          0.047  1.332
# v.v2-Branch1               0.100  0.766  0.911  0.964  0.989  1.985     200.000          0.006  0.021
# v.v2-Branch2               0.100  0.735  0.885  0.953  0.983 11.414     200.000          0.053  1.724
# v.v3-Branch1               0.100  0.774  0.868  0.938  0.990 31.354     200.000          0.153 14.697
# v.v3-Branch2               0.100  0.652  0.859  0.948  0.981  1.260     200.000          0.003  0.010
# ↑↑↑ Expected output ↑↑↑ #

# Helper function to print distribution statistics.
print_stats = function(name, x, digits = 6) {
  q = round(stats::quantile(x), digits)
  message(sprintf("[%s]  Min: %s,  Q1: %s,  Median: %s,  Q3: %s,  Max: %s,  Mean: %s,  Var: %s",
                  name, q[1], q[2], q[3], q[4], q[5],
                  round(base::mean(x), digits), round(stats::var(x), digits)))
}

# Overall theta difference distribution.
print_stats("theta_difference", interest_cell_type_branch_model_train[[1]][["best_pred_result"]][["theta_difference"]])
# ↓↓↓ Expected output ↓↓↓ #
# [theta_difference]  Min: 0.00027,  Q1: 0.010527,  Median: 0.015839,  Q3: 0.039182,  Max: 0.492467,  Mean: 0.048572,  Var: 0.009118
# ↑↑↑ Expected output ↑↑↑ #

# Overall fitting loss distribution.
print_stats("fit_loss", interest_cell_type_branch_model_train[[1]][["best_pred_result"]][["fit_loss"]])
# ↓↓↓ Expected output ↓↓↓ #
# [fit_loss]  Min: 0.021013,  Q1: 0.047189,  Median: 0.062355,  Q3: 0.083531,  Max: 0.174243,  Mean: 0.068949,  Var: 0.000876
# ↑↑↑ Expected output ↑↑↑ #

# Training epoch distribution per branch.
all_epoch_valid = c()
epoch_res_list = list()
for (n in 1:(length(interest_cell_type_branch_model_train[[1]]) - 1)) {
  best_epoch_branch = interest_cell_type_branch_model_train[[1]][["best_pred_result"]][[paste0("best_epoch_branch", n)]]
  epoch_valid = best_epoch_branch[best_epoch_branch > 0]
  all_epoch_valid = c(all_epoch_valid, epoch_valid)
  epoch_quants = round(stats::quantile(epoch_valid), 3)
  epoch_res_list[[paste0("Branch ", n)]] = c(
    Min    = as.numeric(epoch_quants["0%"]),
    Q1     = as.numeric(epoch_quants["25%"]),
    Median = as.numeric(epoch_quants["50%"]),
    Q3     = as.numeric(epoch_quants["75%"]),
    Max    = as.numeric(epoch_quants["100%"]),
    Mean   = round(base::mean(epoch_valid), 3),
    Var    = round(stats::var(epoch_valid), 3)
  )
}
epoch_quants_all = round(stats::quantile(all_epoch_valid), 3)
epoch_res_list[["Overall"]] = c(
  Min    = as.numeric(epoch_quants_all["0%"]),
  Q1     = as.numeric(epoch_quants_all["25%"]),
  Median = as.numeric(epoch_quants_all["50%"]),
  Q3     = as.numeric(epoch_quants_all["75%"]),
  Max    = as.numeric(epoch_quants_all["100%"]),
  Mean   = round(base::mean(all_epoch_valid), 3),
  Var    = round(stats::var(all_epoch_valid), 3)
)
epoch_stats_df = base::as.data.frame(base::do.call(base::rbind, epoch_res_list), check.names = FALSE)
base::print(epoch_stats_df)
# ↓↓↓ Expected output ↓↓↓ #
#          Min   Q1 Median   Q3 Max   Mean      Var
# Branch 1   2  5.5     16 33.5  96 22.222  396.176
# Branch 2   1 13.5     35 76.5 300 55.819 4637.051
# Overall    1 10.0     26 47.5 300 40.141 2923.286
# ↑↑↑ Expected output ↑↑↑ #

# Training time distribution per branch.
parse_time_to_mins = function(time_str) {
  time_str = trimws(time_str)
  num_val = as.numeric(sub("^([0-9.]+).*", "\\1", time_str))
  unit = sub("^[0-9.]+\\s*(.*)$", "\\1", time_str)
  num_val[unit == "secs"]  = num_val[unit == "secs"] / 60
  num_val[unit == "hours"] = num_val[unit == "hours"] * 60
  num_val[unit == "days"]  = num_val[unit == "days"] * 1440
  return(num_val)
}
all_time_valid = c()
res_list = list()
for (n in 1:(length(interest_cell_type_branch_model_train[[1]]) - 1)) {
  time_raw = interest_cell_type_branch_model_train[[1]][["best_pred_result"]][[paste0("train_time_branch", n)]]
  time_mins = parse_time_to_mins(time_raw)
  time_valid = time_mins[time_mins > 0]
  all_time_valid = c(all_time_valid, time_valid)
  time_quants = round(stats::quantile(time_valid), 3)
  res_list[[paste0("Branch ", n)]] = c(
    Min    = as.numeric(time_quants["0%"]),
    Q1     = as.numeric(time_quants["25%"]),
    Median = as.numeric(time_quants["50%"]),
    Q3     = as.numeric(time_quants["75%"]),
    Max    = as.numeric(time_quants["100%"]),
    Mean   = round(base::mean(time_valid), 3),
    Var    = round(stats::var(time_valid), 3)
  )
}
time_quants_all = round(stats::quantile(all_time_valid), 3)
res_list[["Overall"]] = c(
  Min    = as.numeric(time_quants_all["0%"]),
  Q1     = as.numeric(time_quants_all["25%"]),
  Median = as.numeric(time_quants_all["50%"]),
  Q3     = as.numeric(time_quants_all["75%"]),
  Max    = as.numeric(time_quants_all["100%"]),
  Mean   = round(base::mean(all_time_valid), 3),
  Var    = round(stats::var(all_time_valid), 3)
)
time_stats_df = base::as.data.frame(base::do.call(base::rbind, res_list), check.names = FALSE)
base::print(time_stats_df)
# ↓↓↓ Expected output ↓↓↓ #
#           Min    Q1 Median     Q3   Max   Mean    Var
# Branch 1 3.20 4.555  5.990  8.065 12.78  6.332  4.481
# Branch 2 4.25 7.060 11.595 16.462 33.39 12.582 48.628
# Overall  3.20 5.230  7.800 12.380 33.39  9.665 37.633
# ↑↑↑ Expected output ↑↑↑ #

# Training stop cause distribution per branch.
all_cause_vec = c()
cause_res_list = list()
cause_order = c("loss improvement is below the threshold", "no loss improvement for consecutive early_stop epochs", "predicted theta_p converges to observed theta_i")
for (n in 1:(length(interest_cell_type_branch_model_train[[1]]) - 1)) {
  cause_vec = interest_cell_type_branch_model_train[[1]][["best_pred_result"]][[paste0("stop_cause_branch", n)]]
  all_cause_vec = c(all_cause_vec, cause_vec)
  count_tb = base::table(cause_vec)
  freq_tb = round(prop.table(count_tb), 3)
  row_vec = c()
  for (cause in cause_order) {
    cnt = base::ifelse(cause %in% base::names(count_tb), count_tb[cause], 0)
    frq = base::ifelse(cause %in% base::names(freq_tb), freq_tb[cause], 0)
    row_vec[paste0(cause, " - count")] = as.numeric(cnt)
    row_vec[paste0(cause, " - frequency")] = as.numeric(frq)
  }
  cause_res_list[[paste0("Branch ", n)]] = row_vec
}
count_tb_all = base::table(all_cause_vec)
freq_tb_all = round(prop.table(count_tb_all), 3)
row_vec_all = c()
for (cause in cause_order) {
  cnt = base::ifelse(cause %in% base::names(count_tb_all), count_tb_all[cause], 0)
  frq = base::ifelse(cause %in% base::names(freq_tb_all), freq_tb_all[cause], 0)
  row_vec_all[paste0(cause, " - count")] = as.numeric(cnt)
  row_vec_all[paste0(cause, " - frequency")] = as.numeric(frq)
}
cause_res_list[["Overall"]] = row_vec_all
cause_stats_df = base::as.data.frame(base::do.call(base::rbind, cause_res_list), check.names = FALSE)
base::print(base::t(cause_stats_df))
# ↓↓↓ Expected output ↓↓↓ #
#                                                                   Branch 1 Branch 2 Overall
# loss improvement is below the threshold - count                     37.000   46.000  83.000
# loss improvement is below the threshold - frequency                  0.587    0.639   0.615
# no loss improvement for consecutive early_stop epochs - count       23.000   15.000  38.000
# no loss improvement for consecutive early_stop epochs - frequency    0.365    0.208   0.281
# predicted theta_p converges to observed theta_i - count              3.000    8.000  11.000
# predicted theta_p converges to observed theta_i - frequency          0.048    0.111   0.081
# ↑↑↑ Expected output ↑↑↑ #





##################################
##################################
###                            ###
###   Part 5: Predicting GRN   ###
###                            ###
##################################
##################################




######################################
#                                    #
# 5.1 Getting Data for Inferring GRN #
#                                    #
######################################
### Function
# Generates inference-ready datasets by combining branch-level iGRN pairs with time-series cell state features.
# Supports partition-based parallel processing for large datasets.
# Outputs both branch-level full feature sets and a global deduplicated regulatory pair set.

# See function documentation for input, logic and output details.
?get_data_for_inferring_GRN


### Data import
setwd(Working_directory)

# Branch-specific iGRN matrices (TGs as rows, TFs as columns) for each cell type, output from get_branch_iGRN.
interest_cell_type_branch_iGRN = qs::qread("./3 get iGRN/interest_cell_type_branch_iGRN.qsave")

# Aggregated TFE/TGA/TGE matrices per branch, output from cell_grouping.
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")

# Number of cores for parallel computation via mclapply.
# Parallel processing only works on Unix-like systems; set to 1 on Windows.
# Adjust based on available RAM; reduce core count if memory is exhausted.
ncores = parallel::detectCores() - 1 # Linux setting
# ncores = 1 # Windows setting

# Number of data partitions for parallel feature binding. Defaults to 1.
# Recommended to match the number of cores used in execution.
n_part = ncores

### Run model_train with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_data_for_GRN_infer = get_data_for_inferring_GRN(
    interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
    interest_cell_type_group = interest_cell_type_group,
    n_part = n_part, ncores = ncores
  )
})
message(print_run_stats(monitor))
# Example output: running: 0 hours 1 mins 57 secs, peak RAM: 0 GB 210 MB


### Inspect results.
# Branch-specific TF-TG pair data frames with time-point features.
# Columns include TG, TF, theta_i, and TFE_K, TGA_K, TGE_K for each time point. Pairs with theta_i = 0 are filtered out.
base::print(interest_cell_type_branch_data_for_GRN_infer[[1]][["branch1"]][1:3, 1:10]) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#                  TG   TF   theta_i    TFE_1    TGA_1    TGE_1    TFE_2    TGA_2    TGE_2    TFE_3
# ATF3_to_ADNP2 ADNP2 ATF3 0.7573849 3.072106 1.508963 2.719392 3.067088 1.538563 2.596248 3.019299
# ATF3_to_APTX   APTX ATF3 0.7646319 3.072106 1.909924 2.715656 3.067088 1.876471 2.654134 3.019299
# ATF3_to_BCL6   BCL6 ATF3 0.6101782 3.072106 1.845628 3.061429 3.067088 1.915393 2.892772 3.019299
# ↑↑↑ Expected output ↑↑↑ #

base::print(interest_cell_type_branch_data_for_GRN_infer[[1]][["branch2"]][1:3, 1:10]) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#                  TG   TF   theta_i    TFE_1    TGA_1    TGE_1    TFE_2    TGA_2    TGE_2    TFE_3
# ATF3_to_ADNP2 ADNP2 ATF3 0.7573849 3.071233 1.526293 2.676563 3.081425 1.530290 2.599428 2.916346
# ATF3_to_APTX   APTX ATF3 0.7646319 3.071233 1.911295 2.688297 3.081425 1.871386 2.689494 2.916346
# ATF3_to_BCL6   BCL6 ATF3 0.6101782 3.071233 1.854588 3.034467 3.081425 1.921161 2.874683 2.916346
# ↑↑↑ Expected output ↑↑↑ #

# Global deduplicated TF-TG pair data frame across all branches.
# Contains only TG, TF and theta_i columns; pairs with theta_i = 0 are filtered out.
base::print(utils::head(interest_cell_type_branch_data_for_GRN_infer[[1]][["GRN"]], n = 3)) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#                  TG   TF   theta_i
# ATF3_to_ADNP2 ADNP2 ATF3 0.7573849
# ATF3_to_APTX   APTX ATF3 0.7646319
# ATF3_to_BCL6   BCL6 ATF3 0.6101782
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(interest_cell_type_branch_data_for_GRN_infer, file = "./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.qsave")




#####################################
#                                   #
# 5.2 Inferring GRN for Each Branch #
#                                   #
#####################################
### Function
# Predicts optimal regulatory strengths and prediction losses for all candidate TF-TG pairs using pre-trained scDERI models.
# Supports subset processing by cell type and branch, chunk-based parallel computation.
# Automatically saves branch-level results.

# See function documentation for input, logic and output details.
?infer_GRN


### Data import
setwd(Working_directory)

# Branch-level TF-TG pairs with time-series expression features, output from get_data_for_inferring_GRN.
interest_cell_type_branch_data_for_GRN_infer = qs::qread("./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.qsave")

# Cell-type-specific binarization thresholds, output from get_cell_type_tau.
interest_cell_type_tau = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")

# Trained model parameters for each branch, output from model_train.
interest_cell_type_branch_model_train = qs::qread("./4.2 BUild Prediction Model/interest_cell_type_branch_model_train.qsave")

# Pseudotime interval lengths per branch, output from cell_grouping.
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")

# Number of cores for parallel computation via mclapply.
# Parallel processing only works on Unix-like systems; set to 1 on Windows.
# Adjust based on available RAM; reduce core count if memory is exhausted.
ncores = parallel::detectCores() - 1 # Linux setting
# ncores = 1 # Windows setting


### Way 1: One-step inference (for sufficient memory)
# Maximum number of regulatory pairs per parallel chunk. Adjust to control peak memory usage. Defaults to 1000000.
# Recommended to distribute data load evenly across cores.
each = ceiling((max(base::sapply(
  interest_cell_type_branch_data_for_GRN_infer[[1]][1:(length(interest_cell_type_branch_data_for_GRN_infer[[1]]) - 1)],
  base::nrow
))) / ncores)
base::print(each) # Example output: 92

# Run model_train with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_pGRN = infer_GRN(
    interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
    interest_cell_type_tau = interest_cell_type_tau,
    interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
    interest_cell_type_group = interest_cell_type_group,
    select_cell_type = NULL, select_branch = NULL, each = each, ncores = ncores
  )
})
message(print_run_stats(monitor))
# Example output: running: 0 hours 10 mins 12 secs, peak RAM: 0 GB 195 MB

# This method may cause out-of-memory errors for large datasets with many branches.
# Use Way 2 below for more stable execution.


### Way 2: Branch-by-branch inference (stable but more steps)
# Process each branch separately to avoid gradual memory accumulation.
# Run branches one by one manually if the loop still risks memory overflow.
interest_cell_type_pGRN = list()
cell_type = base::names(interest_cell_type_branch_data_for_GRN_infer)[1] # Mono
interest_cell_type_pGRN[[cell_type]] = list()

for (n in 1:(length(interest_cell_type_branch_data_for_GRN_infer[[cell_type]]) - 1)) {
  message("Infer cell type ", cell_type, " branch ", n, ".")
  each = ceiling(base::nrow(interest_cell_type_branch_data_for_GRN_infer[[1]][[paste0("branch", n)]]) / ncores)
  base::print(each) # Example output: 87 92
  gc()

  # Run model_train with runtime and peak memory monitoring.
  monitor = peakRAM::peakRAM({
    interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]] = infer_GRN(
      interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
      interest_cell_type_tau = interest_cell_type_tau,
      interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
      interest_cell_type_group = interest_cell_type_group,
      select_cell_type = cell_type, select_branch = n,
      each = each, ncores = ncores
    )[[cell_type]][[paste0("branch", n)]]
  })
  message("Cell type ", cell_type, " - Branch ", n, ": ",print_run_stats(monitor))
  # ↓↓↓ Expected output ↓↓↓ #
  # Cell type Mono - Branch 1: running: 0 hours 4 mins 17 secs, peak RAM: 0 GB 138 MB
  # Cell type Mono - Branch 2: running: 0 hours 4 mins 5 secs, peak RAM: 0 GB 166 MB
  # ↑↑↑ Expected output ↑↑↑ #

  save.image(paste0("./5 Infer GRN/5_Infer_GRN_branch", n, "middle.RData")) # Check point
  # load(paste0("./5 Infer GRN/5_Infer_GRN_branch", n, "middle.RData"))

  gc()
}


### Inspect results.
# Predicted regulatory strength matrix. Rows: candidate TF-TG pairs; columns: trained models.
base::print(interest_cell_type_pGRN[[1]][["branch1"]][["theta_p"]][1:3, 1:3]) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#               ATF3_to_RPS29 ATF7_to_PUM3 BATF3_to_APTX
# ATF3_to_ADNP2    0.01000001   0.01000001     0.1313364
# ATF3_to_APTX     0.01000001   0.98999997     0.7185931
# ATF3_to_BCL6     0.01000001   0.45137411     0.4257530
# ↑↑↑ Expected output ↑↑↑ #

base::print(interest_cell_type_pGRN[[1]][["branch2"]][["theta_p"]][1:3, 1:3]) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#               ATF3_to_RPS29 ATF7_to_PUM3 BACH2_to_GPI
# ATF3_to_ADNP2    0.01000001    0.9838706     0.794008
# ATF3_to_APTX     0.01000001    0.9134817     0.990000
# ATF3_to_BCL6     0.01000001    0.5783428     0.990000
# ↑↑↑ Expected output ↑↑↑ #

# Minimum prediction loss matrix corresponding to optimal theta values.
base::print(interest_cell_type_pGRN[[1]][["branch1"]][["loss"]][1:3, 1:3]) # Branch 1
# ↓↓↓ Expected output ↓↓↓ #
#               ATF3_to_RPS29 ATF7_to_PUM3 BATF3_to_APTX
# ATF3_to_ADNP2    0.07270821   0.07570539    0.08204194
# ATF3_to_APTX     0.07919798   0.07902704    0.07706172
# ATF3_to_BCL6     0.07255131   0.07439087    0.07609968
# ↑↑↑ Expected output ↑↑↑ #

base::print(interest_cell_type_pGRN[[1]][["branch2"]][["loss"]][1:3, 1:3]) # Branch 2
# ↓↓↓ Expected output ↓↓↓ #
#               ATF3_to_RPS29 ATF7_to_PUM3 BACH2_to_GPI
# ATF3_to_ADNP2     0.1601121    0.1386336    0.2182068
# ATF3_to_APTX      0.1431057    0.1250938    0.1944710
# ATF3_to_BCL6      0.1197168    0.1047312    0.1697980
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.qsave")




################################################
#                                              #
# 5.3 Selecting Best Prediction for each TF-TG #
#                                              #
################################################
### Function
# Post-processes inferred GRN results to select high-confidence regulatory predictions.
# Applies four-level hierarchical filtering to select the optimal prediction from multiple trained models per pair.
# Integrates branch-level results into cell-type-level summaries using cell-count weighted averaging.

# See function documentation for input, logic and output details.
?select_best_prediction


### Data import
setwd(Working_directory)

# Inferred GRN results, output from infer_GRN.
interest_cell_type_pGRN = qs::qread("./5 Infer GRN/interest_cell_type_pGRN.qsave")

# Cell-type-specific binarization thresholds, output from get_cell_type_tau.
interest_cell_type_tau = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")

# Cell counts per branch for weighted integration, output from get_genes_pseudotime_info.
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")

# Number of cores for parallel computation via mclapply.
# Parallel processing only works on Unix-like systems; set to 1 on Windows.
# Adjust based on available RAM; reduce core count if memory is exhausted.
ncores = parallel::detectCores() - 1 # Linux setting
# ncores = 1 # Windows setting

# Maximum number of TF-TG pairs per chunk for parallel processing. Defaults to 1,000,000.
# Recommended to distribute data load evenly across cores.
interest_cell_type_branch_data_for_GRN_infer = qs::qread("./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.qsave")
each = ceiling((max(base::sapply(
  interest_cell_type_branch_data_for_GRN_infer[[1]][1:(length(interest_cell_type_branch_data_for_GRN_infer[[1]]) - 1)],
  base::nrow
))) / ncores)
base::print(each) # Example output: 92


### Hyperparameter configuration
# Quantile threshold between 0 and 1 for filtering predictions based on loss. Defaults to 0.2.
fit_loss_quantile_threshold = 0.2

# Quantile threshold between 0 and 1 for regulatory strength filtering. Defaults to 0.8.
theta_p_quantile_threshold = 0.8


### Run select_best_prediction with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_pGRN = select_best_prediction(
    interest_cell_type_pGRN = interest_cell_type_pGRN, interest_cell_type_tau = interest_cell_type_tau,
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    fit_loss_quantile_threshold = fit_loss_quantile_threshold, theta_p_quantile_threshold = theta_p_quantile_threshold,
    each = each, ncores = ncores
  )
})
message(base::nrow(interest_cell_type_pGRN[[1]][["pred_summary"]]), " pairs remaining, ", print_run_stats(monitor))
# Example output: 7560 pairs remaining, running: 0 hours 0 mins 34 secs, peak RAM: 0 GB 58 MB


### Inspect results.
# Per-branch best prediction results.
# Each best_result data frame contains optimal regulatory strength (theta_p), fitting loss, and selected model name for each TF-TG pair.
base::print(utils::head(interest_cell_type_pGRN[[1]][["branch1"]][["best_result"]], n = 3))
# ↓↓↓ Expected output ↓↓↓ #
#               theta_p_branch1 loss_branch1   fit_branch1
# ATF3_to_ADNP2            0.99   0.07562459 CEBPD_to_EID1
# ATF3_to_APTX             0.99   0.06948342 RELB_to_BIRC3
# ATF3_to_BCL6             0.99   0.06583464 RELB_to_BIRC3
# ↑↑↑ Expected output ↑↑↑ #

base::print(utils::head(interest_cell_type_pGRN[[1]][["branch2"]][["best_result"]], n = 3))
# ↓↓↓ Expected output ↓↓↓ #
#               theta_p_branch2 loss_branch2     fit_branch2
# ATF3_to_ADNP2       0.9900000    0.1312688 NFYA_to_TMEM160
# ATF3_to_APTX        0.9660375    0.1179775  NR2C1_to_SNRPF
# ATF3_to_BCL6        0.8884499    0.1012106  CLOCK_to_NR1D2
# ↑↑↑ Expected output ↑↑↑ #

# Cell-type-level integrated summary with branch-specific values and weighted averages.
# Theta_p and fit_loss are weighted by cell count across branches.
base::print(utils::head(interest_cell_type_pGRN[[1]][["pred_summary"]], n = 3))
# ↓↓↓ Expected output ↓↓↓ #
#               theta_p_branch1 loss_branch1   fit_branch1 theta_p_branch2 loss_branch2     fit_branch2   theta_p   fit_loss
# ATF3_to_ADNP2            0.99   0.07562459 CEBPD_to_EID1       0.9900000    0.1312688 NFYA_to_TMEM160 0.9900000 0.10098468
# ATF3_to_APTX             0.99   0.06948342 RELB_to_BIRC3       0.9660375    0.1179775  NR2C1_to_SNRPF 0.9790790 0.09158478
# ATF3_to_BCL6             0.99   0.06583464 RELB_to_BIRC3       0.8884499    0.1012106  CLOCK_to_NR1D2 0.9437181 0.08195739
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.qsave")




#####################
#                   #
# 5.4 Filtering GRN #
#                   #
#####################
### Function
# Filters inferred GRN results at two stringency levels using dual thresholds.
# Combines quantile statistics and absolute bounds for fitting loss and regulatory strength.
# Outputs both loose and strict filtered high-confidence interaction sets.

# See function documentation for input, logic and output details.
?filter_GRN


### Data import
setwd(Working_directory)

# Inferred GRN results, output from infer_GRN.
interest_cell_type_pGRN = qs::qread("./5 Infer GRN/interest_cell_type_pGRN.qsave")


### Hyperparameter configuration
# Quantile threshold for fitting loss. Final cutoff = min(quantile value, fit_loss_upper). Lower values = stricter filtering. Defaults to 0.2.
fit_loss_quantile_threshold = 0.80 # Retain top 80%; filter bottom 20%.

# Quantile threshold for predicted regulatory strength. Final cutoff = max(quantile value, theta_p_lower). Higher values = stricter filtering. Defaults to 0.8.
theta_p_quantile_threshold = 0.15 # Retain top 85%; filter bottom 15%.

# Absolute upper bound of acceptable fitting loss as a safety cutoff. Defaults to 0.2.
fit_loss_upper = 0.2 # Filter out results with loss > 0.2.

# Absolute lower bound of acceptable regulatory strength as a safety cutoff. Defaults to 0.2.
theta_p_lower = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")["Mono"] # Filter out results with theta_p < tau.

# Final fit loss cutoff is the minimum of the quantile value and absolute upper bound.
# Final theta_p cutoff is the maximum of the quantile value and absolute lower bound.
# Quantile-based thresholds are generally preferred over absolute bounds for most datasets.


### Run filter_GRN with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_pGRN = filter_GRN(
    interest_cell_type_pGRN = interest_cell_type_pGRN,
    fit_loss_quantile_threshold = fit_loss_quantile_threshold, theta_p_quantile_threshold = theta_p_quantile_threshold,
    fit_loss_upper = fit_loss_upper, theta_p_lower = theta_p_lower
  )
})
message(base::nrow(interest_cell_type_pGRN[[1]][["pred_summary_loose"]]), " loose pairs, ", print_run_stats(monitor))
# Example output: 5133 loose pairs, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 3 MB


### Inspect results.
# Two GRN sets with different stringency levels.
# Loose set includes interactions present in any branch.
# Strict set retains only interactions present in all branches.
base::print(utils::head(interest_cell_type_pGRN[[1]][["pred_summary_loose"]], n = 3))
# ↓↓↓ Expected output ↓↓↓ #
#               theta_p_branch1 loss_branch1   fit_branch1 theta_p_branch2 loss_branch2    fit_branch2   theta_p   fit_loss
# ATF3_to_APTX             0.99   0.06948342 RELB_to_BIRC3       0.9660375    0.1179775 NR2C1_to_SNRPF 0.9790790 0.09158478
# ATF3_to_BCL6             0.99   0.06583464 RELB_to_BIRC3       0.8884499    0.1012106 CLOCK_to_NR1D2 0.9437181 0.08195739
# ATF3_to_BUD31            0.99   0.07886349 CEBPD_to_EID1              NA           NA           <NA> 0.9900000 0.07886349
# ↑↑↑ Expected output ↑↑↑ #

base::print(utils::head(interest_cell_type_pGRN[[1]][["pred_summary_strict"]], n = 3))
# ↓↓↓ Expected output ↓↓↓ #
#                theta_p_branch1 loss_branch1   fit_branch1 theta_p_branch2 loss_branch2    fit_branch2   theta_p   fit_loss
# ATF3_to_APTX              0.99   0.06948342 RELB_to_BIRC3       0.9660375   0.11797748 NR2C1_to_SNRPF 0.9790790 0.09158478
# ATF3_to_BCL6              0.99   0.06583464 RELB_to_BIRC3       0.8884499   0.10121062 CLOCK_to_NR1D2 0.9437181 0.08195739
# ATF3_to_CDC14A            0.99   0.06601661 RELB_to_BIRC3       0.9900000   0.09426061   TCF4_to_TCF7 0.9900000 0.07888893
# ↑↑↑ Expected output ↑↑↑ #

# The loose GRN preserves more comprehensive regulatory interactions.
# The strict GRN captures core regulatory relationships conserved across all branches.


### Save output.
qs::qsave(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.qsave")




#######################
#                     #
# 5.5 Inferring peGRN #
#                     #
#######################
### Function
# Integrates iGRNs, pGRNs and single-cell multi-omics features to construct peak-level enhancer-driven GRNs (peGRNs).
# Returns results at two confidence levels (loose and strict) for each cell type.

# See function documentation for input, logic and output details.
?infer_eGRN


### Data import
setwd(Working_directory)

# scRNA-seq Seurat object for gene-level expression feature calculation.
scRNAseq = qs::qread(file = "./scRNAseq.qsave")

# scATAC-seq Seurat object for peak-level accessibility feature calculation.
scATACseq = qs::qread("./scATACseq.qsave")

# Peak-to-gene promoter mapping results, output from identify_TGs or filter_peak_annotations.
# Must contain peak_anno_promoter element with peak and SYMBOL columns.
results_identify_TGs = qs::qread("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.qsave")

# Cell-type-level iGRN matrices (TGs × TFs), output from filter_iGRN_by_ChIP_Atlas_priors.
interest_cell_type_iGRN = qs::qread("./3 get iGRN/interest_cell_type_iGRN.qsave")

# Filtered pGRN results with loose and strict summaries, output from filter_GRN.
interest_cell_type_pGRN = qs::qread("./5 Infer GRN/interest_cell_type_pGRN.qsave")


### Run infer_eGRN with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_peGRN = infer_eGRN(
    scRNAseq = scRNAseq, scATACseq = scATACseq, results_identify_TGs = results_identify_TGs,
    interest_cell_type_iGRN = interest_cell_type_iGRN, interest_cell_type_pGRN = interest_cell_type_pGRN,
  )
})
message(print_run_stats(monitor))
# Example output: running: 0 hours 0 mins 16 secs, peak RAM: 0 GB 586 MB


### Inspect results.
# Nested list organized by cell type. Each cell type contains two TF-peak-TG edge tables: loose and strict.
# Core columns: TF, peak, TG; additional features include PWM scores, regulatory strengths and multi-omics metrics.
base::print(utils::head(interest_cell_type_peGRN[[1]][["loose"]], n = 3))
# ↓↓↓ Expected output ↓↓↓ #
#     TF                     peak       TG      FP-S    FP-S-R FG-theta_p       FG-L     FG-NL FG-theta_i FG-theta_i-R    F-NZER F-NZER-R   F-NZEA
# 1 ATF3 chr1-100350863-100355404   CDC14A 0.3915097 0.8024976   0.990000 0.07888893 0.8192691  0.5525195    0.2911641 0.5723333   0.9616 2.241002
# 2 ATF3 chr1-209782958-209785634 TRAF3IP3 0.4904210 0.9024786   0.965576 0.08138587 0.8491474  0.7107178    0.7924751 0.5723333   0.9616 2.241002
# 3 ATF3   chr1-28884498-28888553    EPB41 0.4899576 0.9020988   0.990000 0.08658159 0.9113191  0.6982575    0.7659937 0.5723333   0.9616 2.241002
#   F-NZEA-R   F-EA F-EA-R    P-NZAR P-NZAR-R    P-NZAA P-NZAA-R      P-AA P-AA-R    G-NZER G-NZER-R   G-NZEA G-NZEA-R      G-EA G-EA-R
# 1   0.9593 1.2826  0.964 0.3356667   0.2696 1.0155226   0.6789 0.3408771 0.1795 0.3953333   0.9192 2.060560   0.9026 0.8146082 0.9225
# 2   0.9593 1.2826  0.964 0.3716667   0.3891 0.9835651   0.6161 0.3655584 0.3993 0.4463333   0.9364 1.971479   0.8420 0.8799368 0.9312
# 3   0.9593 1.2826  0.964 0.4236667   0.5402 0.8857056   0.4326 0.3752439 0.4850 0.4426667   0.9348 2.115190   0.9255 0.9363243 0.9385
# ↑↑↑ Expected output ↑↑↑ #

base::print(utils::head(interest_cell_type_peGRN[[1]][["strict"]], n = 3))
# ↓↓↓ Expected output ↓↓↓ #
#     TF                     peak       TG      FP-S    FP-S-R FG-theta_p       FG-L     FG-NL FG-theta_i FG-theta_i-R    F-NZER F-NZER-R   F-NZEA
# 1 ATF3 chr1-100350863-100355404   CDC14A 0.3915097 0.8024976   0.990000 0.07888893 0.7967862  0.5525195    0.2911641 0.5723333   0.9616 2.241002
# 2 ATF3 chr1-209782958-209785634 TRAF3IP3 0.4904210 0.9024786   0.965576 0.08138587 0.8303813  0.7107178    0.7924751 0.5723333   0.9616 2.241002
# 3 ATF3   chr1-28884498-28888553    EPB41 0.4899576 0.9020988   0.990000 0.08658159 0.9002872  0.6982575    0.7659937 0.5723333   0.9616 2.241002
#   F-NZEA-R   F-EA F-EA-R    P-NZAR P-NZAR-R    P-NZAA P-NZAA-R      P-AA P-AA-R    G-NZER G-NZER-R   G-NZEA G-NZEA-R      G-EA G-EA-R
# 1   0.9593 1.2826  0.964 0.3356667   0.2696 1.0155226   0.6789 0.3408771 0.1795 0.3953333   0.9192 2.060560   0.9026 0.8146082 0.9225
# 2   0.9593 1.2826  0.964 0.3716667   0.3891 0.9835651   0.6161 0.3655584 0.3993 0.4463333   0.9364 1.971479   0.8420 0.8799368 0.9312
# 3   0.9593 1.2826  0.964 0.4236667   0.5402 0.8857056   0.4326 0.3752439 0.4850 0.4426667   0.9348 2.115190   0.9255 0.9363243 0.9385
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(interest_cell_type_peGRN, file = "./5 Infer GRN/interest_cell_type_peGRN.qsave")




#######################
#                     #
# 5.6 Filtering peGRN #
#                     #
#######################
### Function
# Filters cell-type-specific peGRNs using an adjusted quantile threshold strategy.
# Fits a base threshold curve inversely related to predicted TF-TG regulatory strength.
# Includes two sequential rescue mechanisms to ensure minimum network connectivity.

# See function documentation for input, logic and output details.
?filter_eGRN


### Data import
setwd(Working_directory)

# peGRN results, output from infer_eGRN.
interest_cell_type_peGRN = qs::qread("./5 Infer GRN/interest_cell_type_peGRN.qsave")

# Three anchor points of FG-theta_p (0–1) for threshold curve fitting.
# Must be sorted ascending and match length of initial_quantile_points. Default: c(0, 0.5, 1).
theta_p_points = c(0, 0.5, 1)

# Corresponding initial quantile thresholds at each theta_p anchor. Must match length of theta_p_points. Default: c(1.5, 0.9, 0.6).
initial_quantile_points = c(5.79, 4.66, 4.23)

# Starting values for inverse curve fitting via minpack.lm::nlsLM (named list: a, b, c). Parameter a is positive-constrained.
theta_p_quantile_fit_start = list(a = 0.646435, b = 0.397096, c = 3.271714)

# Whether to perform peak-TG pair rescue to guarantee minimum retained TFs per pair. Default: TRUE.
rescue_peakTG = TRUE

# Quality floor quantile of FP-S-R for peak-TG rescue candidates. Only candidates above this threshold are eligible. Defaults to 0.1.
peakTG_rescue_quantile_threshold = 0.25

# Minimum number of TFs to retain per peak-TG pair after rescue. Only effective when rescue_peakTG = TRUE. Default: 1.
min_TF_per_peakTG = 1

# Whether to perform TF rescue to guarantee minimum retained peak-TG pairs per TF. Default: TRUE.
rescue_TF = TRUE

# Quality floor quantile of FP-S-R for TF rescue candidates. Only candidates above this threshold are eligible. Defaults to 0.
TF_rescue_quantile_threshold = 0.0

# Minimum number of peak-TG pairs to retain per TF after rescue. Only effective when rescue_TF = TRUE. Default: 2.
min_peakTG_per_TF = 7

# Balance center values for each feature-specific adjustment factor. Controls the neutral position of adjustment.
factor_balance = list(`FG-NL` = 0.2, `FG-theta_i-R` = 0.2,
                      `F-NZER-R` = 0.2, `F-NZEA-R` = 0.2, `F-EA-R` = 0.2,
                      `P-NZAR-R` = 0.2, `P-NZAA-R` = 0.2, `P-AA-R` = 0.2,
                      `G-NZER-R` = 0.2, `G-NZEA-R` = 0.2, `G-EA-R` = 0.2)

# Adjustment amplitude for each feature. Larger values produce stronger penalty/relaxation effects.
factor_range = list(`FG-NL` = 0.43, `FG-theta_i-R` = 0.64,
                    `F-NZER-R` = 0.28, `F-NZEA-R` = 0.33, `F-EA-R` = 0.22,
                    `P-NZAR-R` = 0.19, `P-NZAA-R` = 0.67, `P-AA-R` = 0.13,
                    `G-NZER-R` = 0.21, `G-NZEA-R` = 0.12, `G-EA-R` = 0.26)


### Run filter_eGRN with runtime and peak memory monitoring.
monitor = peakRAM::peakRAM({
  interest_cell_type_peGRN = filter_eGRN(
    interest_cell_type_peGRN = interest_cell_type_peGRN,
    theta_p_points = theta_p_points, initial_quantile_points = initial_quantile_points, theta_p_quantile_fit_start = theta_p_quantile_fit_start,
    rescue_peakTG = rescue_peakTG, peakTG_rescue_quantile_threshold = peakTG_rescue_quantile_threshold, min_TF_per_peakTG = min_TF_per_peakTG,
    rescue_TF = rescue_TF, TF_rescue_quantile_threshold = TF_rescue_quantile_threshold, min_peakTG_per_TF = min_peakTG_per_TF,
    factor_balance = factor_balance, factor_range = factor_range
  )
})
message(print_run_stats(monitor), "\n",
        "[loose] ",
        length(base::unique(interest_cell_type_peGRN[[1]][["loose"]][["TF"]])), "/", length(base::unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][["TF"]])), " TFs, ",
        length(base::unique(interest_cell_type_peGRN[[1]][["loose"]][["peak"]])), "/", length(base::unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][["peak"]])), " peaks, ",
        length(base::unique(interest_cell_type_peGRN[[1]][["loose"]][["TG"]])), "/", length(base::unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][["TG"]])), " TGs, ",
        base::nrow(base::unique(interest_cell_type_peGRN[[1]][["loose"]][, c("TF", "peak")])), "/", base::nrow(base::unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][, c("TF", "peak")])), " TF-peak, ",
        base::nrow(base::unique(interest_cell_type_peGRN[[1]][["loose"]][, c("peak", "TG")])), "/", base::nrow(base::unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][, c("peak", "TG")])), " peak-TG, ",
        base::nrow(base::unique(interest_cell_type_peGRN[[1]][["loose"]][, c("TF", "TG")])), "/", base::nrow(base::unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][, c("TF", "TG")])), " TF-TG, ",
        base::nrow(base::unique(interest_cell_type_peGRN[[1]][["loose"]][, c("TF", "peak", "TG")])), "/", base::nrow(base::unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][, c("TF", "peak", "TG")])), " TF-peak-TG, rescue ratio ",
        round(1 - sum(interest_cell_type_peGRN[[1]][["loose_filtered"]][["filter_keep"]]) / base::nrow(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]]), 3), "\n",
        "[strict] ",
        length(base::unique(interest_cell_type_peGRN[[1]][["strict"]][["TF"]])), "/", length(base::unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][["TF"]])), " TFs, ",
        length(base::unique(interest_cell_type_peGRN[[1]][["strict"]][["peak"]])), "/", length(base::unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][["peak"]])), " peaks, ",
        length(base::unique(interest_cell_type_peGRN[[1]][["strict"]][["TG"]])), "/", length(base::unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][["TG"]])), " TGs, ",
        base::nrow(base::unique(interest_cell_type_peGRN[[1]][["strict"]][, c("TF", "peak")])), "/", base::nrow(base::unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][, c("TF", "peak")])), " TF-peak, ",
        base::nrow(base::unique(interest_cell_type_peGRN[[1]][["strict"]][, c("peak", "TG")])), "/", base::nrow(base::unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][, c("peak", "TG")])), " peak-TG, ",
        base::nrow(base::unique(interest_cell_type_peGRN[[1]][["strict"]][, c("TF", "TG")])), "/", base::nrow(base::unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][, c("TF", "TG")])), " TF-TG, ",
        base::nrow(base::unique(interest_cell_type_peGRN[[1]][["strict"]][, c("TF", "peak", "TG")])), "/", base::nrow(base::unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][, c("TF", "peak", "TG")])), " TF-peak-TG, rescue ratio ",
        round(1 - sum(interest_cell_type_peGRN[[1]][["strict_filtered"]][["filter_keep"]]) / base::nrow(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]]), 3))
# ↓↓↓ Expected output ↓↓↓ #
# running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 25 MB
# [loose] 71/71 TFs, 1287/997 peaks, 1280/992 TGs, 5214/1321 TF-peak, 1304/1009 peak-TG, 5133/1309 TF-TG, 5231/1325 TF-peak-TG, rescue ratio 0.566
# [strict] 64/64 TFs, 1187/923 peaks, 1176/916 TGs, 4547/1236 TF-peak, 1198/932 peak-TG, 4469/1222 TF-TG, 4555/1238 TF-peak-TG, rescue ratio 0.523
# ↑↑↑ Expected output ↑↑↑ #


### Inspect results.
# Full edge table with original columns plus computed adjustment factors, quality scores and retention flags.
base::print(utils::head(interest_cell_type_peGRN[[1]][["loose_filtered"]], n = 3))
# ↓↓↓ Expected output ↓↓↓ #
#     TF                     peak       TG      FP-S    FP-S-R FG-theta_p       FG-L     FG-NL FG-theta_i FG-theta_i-R    F-NZER F-NZER-R   F-NZEA F-NZEA-R   F-EA F-EA-R
# 1 ATF3 chr1-100350863-100355404   CDC14A 0.3915097 0.8024976   0.990000 0.07888893 0.8192691  0.5525195    0.2911641 0.5723333   0.9616 2.241002   0.9593 1.2826  0.964
# 2 ATF3 chr1-209782958-209785634 TRAF3IP3 0.4904210 0.9024786   0.965576 0.08138587 0.8491474  0.7107178    0.7924751 0.5723333   0.9616 2.241002   0.9593 1.2826  0.964
# 3 ATF3   chr1-28884498-28888553    EPB41 0.4899576 0.9020988   0.990000 0.08658159 0.9113191  0.6982575    0.7659937 0.5723333   0.9616 2.241002   0.9593 1.2826  0.964
#      P-NZAR P-NZAR-R    P-NZAA P-NZAA-R      P-AA P-AA-R    G-NZER G-NZER-R   G-NZEA G-NZEA-R      G-EA G-EA-R                              peak-TG       IQ  FG-NL-F
# 1 0.3356667   0.2696 1.0155226   0.6789 0.3408771 0.1795 0.3953333   0.9192 2.060560   0.9026 0.8146082 0.9225   chr1-100350863-100355404_to_CDC14A 4.235973 1.266286
# 2 0.3716667   0.3891 0.9835651   0.6161 0.3655584 0.3993 0.4463333   0.9364 1.971479   0.8420 0.8799368 0.9312 chr1-209782958-209785634_to_TRAF3IP3 4.250880 1.279133
# 3 0.4236667   0.5402 0.8857056   0.4326 0.3752439 0.4850 0.4426667   0.9348 2.115190   0.9255 0.9363243 0.9385      chr1-28884498-28888553_to_EPB41 4.235973 1.305867
#   FG-theta_i-R-F F-NZER-R-F F-NZEA-R-F F-EA-R-F P-NZAR-R-F P-NZAA-R-F P-AA-R-F G-NZER-R-F G-NZEA-R-F G-EA-R-F        AQ filter_keep rescue_peakTG_keep rescue_TF_keep keep
# 1      0.9416550   0.786752   0.749431  0.83192   0.986776   0.679137 1.002665   0.848968   0.915688 0.812150 1.0000000       FALSE               TRUE          FALSE TRUE
# 2      0.6208159   0.786752   0.749431  0.83192   0.964071   0.721213 0.974091   0.845356   0.922960 0.809888 0.7086442        TRUE              FALSE          FALSE TRUE
# 3      0.6377640   0.786752   0.749431  0.83192   0.935362   0.844158 0.962950   0.845692   0.912940 0.807990 0.8207880        TRUE              FALSE          FALSE TRUE
# ↑↑↑ Expected output ↑↑↑ #

base::print(utils::head(interest_cell_type_peGRN[[1]][["strict_filtered"]], n = 3))
# ↓↓↓ Expected output ↓↓↓ #
#     TF                     peak       TG      FP-S    FP-S-R FG-theta_p       FG-L     FG-NL FG-theta_i FG-theta_i-R    F-NZER F-NZER-R   F-NZEA F-NZEA-R   F-EA F-EA-R
# 1 ATF3 chr1-100350863-100355404   CDC14A 0.3915097 0.8024976   0.990000 0.07888893 0.7967862  0.5525195    0.2911641 0.5723333   0.9616 2.241002   0.9593 1.2826  0.964
# 2 ATF3 chr1-209782958-209785634 TRAF3IP3 0.4904210 0.9024786   0.965576 0.08138587 0.8303813  0.7107178    0.7924751 0.5723333   0.9616 2.241002   0.9593 1.2826  0.964
# 3 ATF3   chr1-28884498-28888553    EPB41 0.4899576 0.9020988   0.990000 0.08658159 0.9002872  0.6982575    0.7659937 0.5723333   0.9616 2.241002   0.9593 1.2826  0.964
#      P-NZAR P-NZAR-R    P-NZAA P-NZAA-R      P-AA P-AA-R    G-NZER G-NZER-R   G-NZEA G-NZEA-R      G-EA G-EA-R                              peak-TG       IQ  FG-NL-F
# 1 0.3356667   0.2696 1.0155226   0.6789 0.3408771 0.1795 0.3953333   0.9192 2.060560   0.9026 0.8146082 0.9225   chr1-100350863-100355404_to_CDC14A 4.235973 1.256618
# 2 0.3716667   0.3891 0.9835651   0.6161 0.3655584 0.3993 0.4463333   0.9364 1.971479   0.8420 0.8799368 0.9312 chr1-209782958-209785634_to_TRAF3IP3 4.250880 1.271064
# 3 0.4236667   0.5402 0.8857056   0.4326 0.3752439 0.4850 0.4426667   0.9348 2.115190   0.9255 0.9363243 0.9385      chr1-28884498-28888553_to_EPB41 4.235973 1.301124
#   FG-theta_i-R-F F-NZER-R-F F-NZEA-R-F F-EA-R-F P-NZAR-R-F P-NZAA-R-F P-AA-R-F G-NZER-R-F G-NZEA-R-F G-EA-R-F        AQ filter_keep rescue_peakTG_keep rescue_TF_keep keep
# 1      0.9416550   0.786752   0.749431  0.83192   0.986776   0.679137 1.002665   0.848968   0.915688 0.812150 1.0000000       FALSE               TRUE          FALSE TRUE
# 2      0.6208159   0.786752   0.749431  0.83192   0.964071   0.721213 0.974091   0.845356   0.922960 0.809888 0.7041737        TRUE              FALSE          FALSE TRUE
# 3      0.6377640   0.786752   0.749431  0.83192   0.935362   0.844158 0.962950   0.845692   0.912940 0.807990 0.8178063        TRUE              FALSE          FALSE TRUE
# ↑↑↑ Expected output ↑↑↑ #

# Simplified three-column edge list (TF, peak, TG) for direct downstream network analysis.
base::print(utils::head(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]], n = 3))
# ↓↓↓ Expected output ↓↓↓ #
#     TF                     peak       TG
# 1 ATF3 chr1-100350863-100355404   CDC14A
# 2 ATF3 chr1-209782958-209785634 TRAF3IP3
# 3 ATF3   chr1-28884498-28888553    EPB41
# ↑↑↑ Expected output ↑↑↑ #

base::print(utils::head(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]], n = 3))
# ↓↓↓ Expected output ↓↓↓ #
#     TF                     peak       TG
# 1 ATF3 chr1-100350863-100355404   CDC14A
# 2 ATF3 chr1-209782958-209785634 TRAF3IP3
# 3 ATF3   chr1-28884498-28888553    EPB41
# ↑↑↑ Expected output ↑↑↑ #


### Save output.
qs::qsave(interest_cell_type_peGRN, file = "./5 Infer GRN/interest_cell_type_peGRN.qsave")
