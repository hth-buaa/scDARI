################################################################################
################################################################################
################################################################################
#####                                                                      #####
#####     A Guide to Constructing GRNs with scDARI from scRNA-seq Data     #####
#####                                                                      #####
################################################################################
################################################################################
################################################################################


### Use R v4.4


### The pipeline is organized into six parts:
# Part 0 - Loading data: prepare a Seurat object with a counts layer and cell type labels.
# Part 1 - Data preprocessing: define the TF set, the target-gene (TG) set and the counts matrices used downstream.
# Part 2 - Pseudotime, branching and grouping: order cells along pseudotime, split them into branches, and aggregate cells into metacell groups.
# Part 3 - Initial GRNs: build the information-based initial GRN (iGRN), set the regulatory threshold tau, and derive branch iGRNs.
# Part 4 - Model building: train one DARI regression model per branch.
# Part 5 - GRN prediction: infer, select and filter the final GRN.

# Every scDARI function is self-contained and saves its result to disk,
# so each part can be re-run independently after loading the required intermediate files.


### Environment and path setup
Sys.setenv(LANGUAGE = "en") # Use English for system messages.
options(stringsAsFactors = FALSE) # Do not auto-convert strings to factors.
rm(list = ls()) # Clear the global environment.
gc() # Release unused memory.

# Working directory for all subsequent steps.
Working_directory = "/mnt/i/scDARI" # WSL2 path (recommended).
# Working_directory = "I:/scDARI" # Windows path (not recommended).
# Working_directory = "/home/scDARI" # Native Linux path (available).
setwd(Working_directory)





##############################################################################
##############################################################################
###                                                                        ###
###   Preparations: Installing R Packages and Setting Up the Environment   ###
###                                                                        ###
##############################################################################
##############################################################################

# Most dependencies work with 2024-2026 releases; pinned versions below are for reference only.

# Install igraph 2.0.3 and dplyr 1.1.4 from the CRAN archive.
if (!requireNamespace("igraph", quietly = TRUE) || packageVersion("igraph") != "2.0.3") {
  install.packages("https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz", repos = NULL, type = "source")
}
if (!requireNamespace("dplyr", quietly = TRUE) || packageVersion("dplyr") != "1.1.4") {
  install.packages("https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_1.1.4.tar.gz", repos = NULL, type = "source")
}

# Install CRAN packages.
CRAN_packages = c("Seurat", "Matrix", "utils", "SeuratObject", "parallel", 
                  "minpack.lm", "stats", "ggplot2", "ggrepel", "methods", 
                  "VGAM", "parmigene", "data.table", "GA", "DDRTree", 
                  "BiocManager", "devtools", "peakRAM", "qs")
for (pac in CRAN_packages) {
  if (!requireNamespace(pac, quietly = TRUE)) {
    install.packages(pac, dependencies = TRUE)
  }
  message("The version of CRAN package [", pac, "] is ", packageVersion(pac))
}
# Reference versions: Seurat 5.5.1, Matrix 1.7.6, SeuratObject 5.4.0, minpack.lm 1.2.4,
# ggplot2 4.0.3, VGAM 1.1.14, parmigene 1.1.1, data.table 1.18.6.1, GA 3.2.5,
# DDRTree 0.1.6, BiocManager 1.30.27, devtools 2.5.2, peakRAM 1.0.2 and qs 0.27.3 from CRAN.

# Install Bioconductor packages.
Bioconductor_packages = c("Biobase", "BiocGenerics", "monocle", "S4Vectors")
for (pac in Bioconductor_packages) {
  if (!requireNamespace(pac, quietly = TRUE)) {
    BiocManager::install(pac)
  }
  message("The version of Bioconductor package [", pac, "] is ", packageVersion(pac))
}
# Reference versions: Biobase 2.66.0, BiocGenerics 0.52.0, monocle 2.34.0, S4Vectors 0.44.0.

# Install scDARI package from GitHub repository.
if (!requireNamespace("scDARI", quietly = TRUE)) {
  devtools::install_github("hth-buaa/scDARI")
}
message("The version of CRAN package [scDARI] is ", packageVersion("scDARI")) # 1.0

# Load required packages.
library(scDARI)
library(DDRTree)

# Memory note: parallel steps may not release memory immediately.
# Restart the R session if needed; the modular design lets you resume from saved files.





################################
################################
###                          ###
###   Part 0: Loading Data   ###
###                          ###
################################
################################

# The example is the PBMC scRNA-seq dataset from SeuratData, held in a Seurat object.
# A raw count matrix (genes x cells) is recommended;
# if only an expression matrix is available, store it as the counts layer and the rest of the workflow remains valid.




############################################################
#                                                          #
# 0.1 Loading the Raw Seurat Object of the Example Dataset #
#                                                          #
############################################################
setwd(Working_directory)

# Download pbmcMultiome.SeuratData_0.1.4.tar.gz (about 267.5 MB) first:
# http://seurat.nygenome.org/src/contrib/pbmcMultiome.SeuratData_0.1.4.tar.gz
install.packages("./pbmcMultiome.SeuratData_0.1.4.tar.gz", repos = NULL, type = "source")

scRNAseq = pbmcMultiome.SeuratData::pbmc.rna

# Inspect scRNA-seq Seurat object.
base::print(scRNAseq) # Expected output: 36601 features × 11909 cells, RNA assay with counts/data layers.




############################################################
#                                                          #
# 0.2 Upgrading the Seurat v4 Object to a Seurat v5 Object #
#                                                          #
############################################################
setwd(Working_directory)

# scDARI expects the Seurat v5 Assay5 format; convert a legacy v4 Assay if present.
if (base::class(scRNAseq[["RNA"]]) == "Assay") {
  scRNAseq[["RNA"]] = SeuratObject::CreateAssay5Object(counts = scRNAseq[["RNA"]]@counts)
}
base::print(base::class(scRNAseq[["RNA"]])) # Expected output: Assay5 (SeuratObject)




##########################
#                        #
# 0.3 Adding Cell Labels #
#                        #
##########################
setwd(Working_directory)

# Cell type labels are read from the seurat_annotations metadata column.
base::print(base::as.data.frame(base::table(scRNAseq@meta.data$seurat_annotations)))
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

# If the column is missing, add cell labels manually as a character vector:
# scRNAseq@meta.data$seurat_annotations = your_cell_type_labels.




##########################
#                        #
# 0.4 Sampling the Cells #
#                        #
##########################
setwd(Working_directory)

# Keep only the cell types relevant to the trajectory under study.
# In this monocyte example we retain CD14 Mono and CD16 Mono.
scRNAseq = S4Vectors::subset(scRNAseq, subset = seurat_annotations %in% c("CD14 Mono", "CD16 Mono"))
message(base::ncol(scRNAseq), " cells.") # Expected output: 3326 cells.

# Preserve the original labels, then assign a unified group name to the cells.
scRNAseq@meta.data$original_seurat_annotations = scRNAseq@meta.data$seurat_annotations
scRNAseq@meta.data$seurat_annotations = "Mono"




####################################
#                                  #
# 0.5 Downsampling (Optional Step) #
#                                  #
####################################
setwd(Working_directory)

# Downsampling speeds up the demo but can weaken biological signal.


### Cell downsampling.
# scDARI does not recommend aggressive cell downsampling,
# because cells are later aggregated into metacells (Section 2.2.2) and too few cells harm that step.

# Here we keep 3000 cells for demonstration.
n_cell = 3000

# Rank cells by the fraction of detected genes and keep the top n_cell.
rna_counts = SeuratObject::GetAssayData(scRNAseq, assay = "RNA", layer = "counts")
rna_nonzero_freq = Matrix::colSums(rna_counts > 0) / base::nrow(rna_counts)
top_cells = base::names(rna_nonzero_freq)[base::order(rna_nonzero_freq, decreasing = TRUE)][base::seq_len(n_cell)]
scRNAseq = S4Vectors::subset(scRNAseq, cells = top_cells)
base::print(base::table(scRNAseq@meta.data[["original_seurat_annotations"]]))
# 2493 CD14 Mono and 507 CD16 Mono.


### Gene downsampling.
# Strategy: drop non-coding RNAs and mitochondrial genes, then keep the most variable genes.

# Fewer genes means faster runs; use about 10000 genes for a full analysis, while this demo keeps 1000 for speed.
n_gene = 1000

genes = base::setdiff(
  base::rownames(scRNAseq),
  base::grep("(^AC|^AL|^LINC)[0-9]+", base::rownames(scRNAseq), value = TRUE)
)
genes = base::grep("MT-", genes, invert = TRUE, value = TRUE)
scRNAseq = S4Vectors::subset(scRNAseq, features = genes)

scRNAseq = Seurat::NormalizeData(scRNAseq, verbose = FALSE)
scRNAseq = Seurat::FindVariableFeatures(scRNAseq, selection.method = "vst", nfeatures = base::nrow(scRNAseq), verbose = FALSE)

RNA_var_meta = scRNAseq@assays[["RNA"]]@meta.data
base::rownames(RNA_var_meta) = base::rownames(scRNAseq)
RNA_var_meta = RNA_var_meta[base::order(RNA_var_meta$vf_vst_counts_variance.standardized, decreasing = TRUE), ]
genes = base::rownames(RNA_var_meta)[base::seq_len(n_gene)]
scRNAseq = S4Vectors::subset(scRNAseq, features = genes)
base::rm(RNA_var_meta, genes)

base::print(scRNAseq) # Expected output: 1000 features x 3000 cells.




##################################
#                                #
# 0.6 Saving the Example Dataset #
#                                #
##################################
setwd(Working_directory)

# Save the prepared object. qs is faster than saveRDS for large objects.
# base::saveRDS(scRNAseq, file = "./scRNAseq.rds")
qs::qsave(scRNAseq, file = "./scRNAseq.qsave")





######################################
######################################
###                                ###
###   Part 1: Data Preprocessing   ###
###                                ###
######################################
######################################

# This part takes the prepared Seurat object from Part 0 (counts layer and cell labels) and defines the gene universe for GRN inference:
# - TFs (regulators): must be expressed in the data and supported by ChIP-seq evidence downloaded from ChIP-Atlas.
# - TGs (targets): only need to pass scRNA-seq expression-quality thresholds.
# The two sets are merged into the counts matrices stored in Basic_Info, which is the input to Part 2.




#######################
#                     #
# 1.1 Identifying TFs #
#                     #
#######################

# [get_TFs_from_priors] runs five stages in order:
#  1) start from the built-in ChIP-Atlas TF list for the chosen assembly and keep TFs present in the scRNA-seq data;
#  2) filter those candidates by expression quality;
#  3) download one ChIP-seq target-gene table per qualified TF (cached);
#  4) clean the tables and reject TFs with too few valid targets;
#  5) simplify the retained tables and return the TF names.
?get_TFs_from_priors


### Load inputs.
setwd(Working_directory)
scRNAseq = qs::qread("./scRNAseq.qsave") # Must contain a counts layer.

# Genome assembly: selects the built-in TF list and builds the download URLs.
# Supported: hg38, hg19, mm10, mm9, rn6, dm6, dm3, ce11, ce10, sacCer3.
genome_assembly = "hg38"

# Keywords used to drop irrelevant ChIP-Atlas sample columns by cell type.
excluded_cell_types_in_priors = "NotExclude" # "NotExclude" disables keyword filtering.
# excluded_cell_types_in_priors = base::c("pbmc", "Peripheral", "blood", "mononuclear")


### Expression-quality thresholds.
# Three metrics are computed per gene from raw counts:
#   nonzero_freq  : fraction of cells with count > 0;
#   nonzero_mean  : mean expression within cells with count > 0 only;
#   overall_mean  : mean expression over all cells, zeros included.
# Each metric also has a genome-wide percentile rank (0 = lowest, 1 = highest).
# A candidate TF is kept only when ALL six conditions below hold at once.
# Set any threshold to 0 to disable that condition.
min_TF_nonzero_freq = 0.05 # Absolute minimum detected-cell fraction.
min_TF_nonzero_mean = 0 # Absolute minimum conditional mean.
min_TF_overall_mean = 0 # Absolute minimum overall mean.
min_TF_nonzero_freq_quantile = 0 # Minimum percentile rank of nonzero_freq.
min_TF_nonzero_mean_quantile = 0 # Minimum percentile rank of nonzero_mean.
min_TF_overall_mean_quantile = 0 # Minimum percentile rank of overall_mean.


### Download and prior-size settings.
timeout = 60 # Seconds allowed for one file before the download is marked failed.
distance_from_TSS = 10 # Target genes within this distance (kb) of the TSS; 1/5/10.
prior_min_TG_number = 2000 # A TF needs at least this many non-zero targets to be kept.
prior_top_n_keep = 5000 # Absolute cap on targets retained per TF.
prior_top_q_keep = 0.7 # Proportional cap; the smaller of the two caps applies.


### Run (first pass): some downloads may fail on a slow network.
monitor = peakRAM::peakRAM({
  TFs = get_TFs_from_priors(
    scRNAseq = scRNAseq,
    genome_assembly = genome_assembly,
    min_TF_nonzero_freq = min_TF_nonzero_freq,
    min_TF_nonzero_mean = min_TF_nonzero_mean,
    min_TF_overall_mean = min_TF_overall_mean,
    min_TF_nonzero_freq_quantile = min_TF_nonzero_freq_quantile,
    min_TF_nonzero_mean_quantile = min_TF_nonzero_mean_quantile,
    min_TF_overall_mean_quantile = min_TF_overall_mean_quantile,
    timeout = timeout,
    distance_from_TSS = distance_from_TSS,
    delay_time = 0.1,
    excluded_cell_types_in_priors = excluded_cell_types_in_priors,
    prior_min_TG_number = prior_min_TG_number,
    prior_top_n_keep = prior_top_n_keep,
    prior_top_q_keep = prior_top_q_keep
  )
})
message("Get ", base::length(TFs), " TFs, ", print_run_stats(monitor))
# Example: Get 29 TFs, running: 0 hours 10 mins 29 secs, peak RAM: 0 GB 79 MB.


### Check failed downloads.
download_information = readRDS("./1.1 Data Preprocessing - Get TFs From Database/ChIP_Atlas_Priors/download_information.rds")
base::print(download_information[download_information$status == "URL not found or download error", ])
# ↓↓↓ Expected output ↓↓↓ #
#                                                           url                                                                             destination                          status          time
# EGR1 https://chip-atlas.dbcls.jp/data/hg38/target/EGR1.10.tsv /mnt/i/scDARI/1.1 Data Preprocessing - Get TFs From Database/ChIP_Atlas_Priors/EGR1.tsv URL not found or download error 1.000402 mins
# ESR1 https://chip-atlas.dbcls.jp/data/hg38/target/ESR1.10.tsv /mnt/i/scDARI/1.1 Data Preprocessing - Get TFs From Database/ChIP_Atlas_Priors/ESR1.tsv URL not found or download error  1.00037 mins
# JUN   https://chip-atlas.dbcls.jp/data/hg38/target/JUN.10.tsv  /mnt/i/scDARI/1.1 Data Preprocessing - Get TFs From Database/ChIP_Atlas_Priors/JUN.tsv URL not found or download error 1.000377 mins
# ↑↑↑ Expected output ↑↑↑ #

# For each failed TF, manually download the file listed in $url to the path in $destination (e.g. EGR1, ESR1, JUN may fail on the first pass).


### Re-run: cached files are reused, so failed TFs are recovered.
monitor = peakRAM::peakRAM({
  TFs = get_TFs_from_priors(
    scRNAseq = scRNAseq,
    genome_assembly = genome_assembly,
    min_TF_nonzero_freq = min_TF_nonzero_freq,
    min_TF_nonzero_mean = min_TF_nonzero_mean,
    min_TF_overall_mean = min_TF_overall_mean,
    min_TF_nonzero_freq_quantile = min_TF_nonzero_freq_quantile,
    min_TF_nonzero_mean_quantile = min_TF_nonzero_mean_quantile,
    min_TF_overall_mean_quantile = min_TF_overall_mean_quantile,
    timeout = timeout,
    distance_from_TSS = distance_from_TSS,
    delay_time = 0.1,
    excluded_cell_types_in_priors = excluded_cell_types_in_priors,
    prior_min_TG_number = prior_min_TG_number,
    prior_top_n_keep = prior_top_n_keep,
    prior_top_q_keep = prior_top_q_keep
  )
})
message("Get ", base::length(TFs), " TFs, ", print_run_stats(monitor))
# Example: Get 32 TFs, running: 0 hours 0 mins 37 secs, peak RAM: 0 GB 257 MB.

### Inspect and save.
base::print(utils::head(TFs))
# AHR, ARID5B, BACH2, BCL11B, CDK6, CREB5

qs::qsave(TFs, file = "./1.1 Data Preprocessing - Get TFs From Database/TFs.qsave")




#####################
#                   #
# 1.2 Filtering TGs #
#                   #
#####################

# [filter_TGs] computes the same three metrics and percentile ranks as 1.1 over all genes,
# keeps genes passing all thresholds, and caches the feature table.
?filter_TGs


### Load inputs and set thresholds.
setwd(Working_directory)
scRNAseq = qs::qread("./scRNAseq.qsave") # Must contain a counts layer.

# Same metric definitions as in 1.1; all six conditions must hold, and 0
# disables a condition.
min_TG_nonzero_freq = 0.05
min_TG_nonzero_mean = 0
min_TG_overall_mean = 0
min_TG_nonzero_freq_quantile = 0
min_TG_nonzero_mean_quantile = 0
min_TG_overall_mean_quantile = 0

# [forced_feature_calculation]:
#   TRUE  -> recompute the metrics and overwrite the cached table;
#   FALSE -> load the cached TG_features_all.rds when it exists.
forced_feature_calculation = TRUE


### Run, inspect and save.
monitor = peakRAM::peakRAM({
  TGs = filter_TGs(
    scRNAseq = scRNAseq,
    forced_feature_calculation = forced_feature_calculation,
    min_TG_nonzero_freq = min_TG_nonzero_freq,
    min_TG_nonzero_mean = min_TG_nonzero_mean,
    min_TG_overall_mean = min_TG_overall_mean,
    min_TG_nonzero_freq_quantile = min_TG_nonzero_freq_quantile,
    min_TG_nonzero_mean_quantile = min_TG_nonzero_mean_quantile,
    min_TG_overall_mean_quantile = min_TG_overall_mean_quantile
  )
})
message("Get ", base::length(TGs), " TGs, ", print_run_stats(monitor))
# Get 519 TGs, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 10 MB

base::print(utils::head(TGs)) # ABCA1, ABCC3, ABHD5, ABLIM1, ABTB2, ACSL1.
qs::qsave(TGs, file = "./1.2 Data Preprocessing - Filter TGs/TGs.qsave")





#########################################################
#                                                       #
# 1.3 Generating Counts Matrices for TFs/TGs Expression #
#                                                       #
#########################################################

# [get_expression_matrix] merges the TF and TG sets,
# subsets the Seurat object to that union,
# and returns the harmonized inputs for GRN construction.
?get_expression_matrix


### Load inputs.
setwd(Working_directory)
TFs = qs::qread("./1.1 Data Preprocessing - Get TFs From Database/TFs.qsave")
TGs = qs::qread("./1.2 Data Preprocessing - Filter TGs/TGs.qsave")
scRNAseq = qs::qread("./scRNAseq.qsave")


### Run.
monitor = peakRAM::peakRAM({
  Basic_Info = get_expression_matrix(
    TFs = TFs,
    TGs = TGs,
    scRNAseq = scRNAseq
  )
})
message(base::length(Basic_Info[["TFs"]]), " TFs, ",
        base::length(Basic_Info[["TGs"]]), " TGs, ",
        base::ncol(Basic_Info[["scRNAseq_for_GRN"]]), " cells, ",
        print_run_stats(monitor))
# 32 TFs, 519 TGs, 3000 cells, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 28 MB


### Inspect the four returned elements.
# scRNAseq_for_GRN: Seurat object restricted to the union of TFs and TGs.
base::print(Basic_Info[["scRNAseq_for_GRN"]])
# 519 features across 3000 samples, RNA assay with counts and data layers.

# gene_for_GRN: sorted de-duplicated union of TFs and TGs.
base::print(utils::head(Basic_Info[["gene_for_GRN"]]))

# TFs / TGs: the harmonized input gene lists (sorted, de-duplicated).
base::print(utils::head(Basic_Info[["TFs"]]))
base::print(utils::head(Basic_Info[["TGs"]]))

# Note: gene_for_GRN, TFs and TGs are derived from the inputs only and are not intersected with the object's features,
# whereas scRNAseq_for_GRN silently keeps only features available in the object.


### Save.
qs::qsave(Basic_Info, file = "./1.3 Data Preprocessing - Filter Expression Matrix/Basic_Info.qsave")





#######################################################################
#######################################################################
###                                                                 ###
###   Part 2: Inferring Pseudotime, Branching, and Grouping Cells   ###
###                                                                 ###
#######################################################################
#######################################################################

# Part 2 adds the time axis and turns it into compact pseudotime-series inputs for modeling:
#   2.1  For each cell type, infer a trajectory with Monocle2 (DDRTree), assign pseudotime,
#        and divide the cells into branches from the trajectory plot.
#   2.2  On every branch, order cells by pseudotime, extract TF expression (TFE) and TG expression (TGE),
#        then group consecutive cells adaptively and average their non-zero expression.




########################################################
#                                                      #
# 2.1 Pseudotime Analysis and Cell Branch Partitioning #
#                                                      #
########################################################



###########################################################
# 2.1.1 Obtaining Relevant Data of Cell Types of Interest #
###########################################################

# [get_interest_cell_type_data] processes each cell type independently:
# it subsets the GRN-ready object, removes all-zero genes/cells,
# optionally filters by zero-expression rate, and rebuilds the TF/TG sets for that cell type.
?get_interest_cell_type_data


### Load inputs.
setwd(Working_directory)
Basic_Info = qs::qread("./1.3 Data Preprocessing - Filter Expression Matrix/Basic_Info.qsave")

# Cell types to analyze; values must exist in seurat_annotations.
interest_cell_type = base::unique(
  Basic_Info[["scRNAseq_for_GRN"]]@meta.data[["seurat_annotations"]]
)
base::print(interest_cell_type) # Mono.

# Parallel cores for mclapply (Unix/macOS only); use 1 on Windows.
ncores = parallel::detectCores() - 1 # 78 in development
# ncores = 1 # Windows.


### Run, inspect and save.
monitor = peakRAM::peakRAM({
  interest_cell_type_data = get_interest_cell_type_data(
    Basic_Info = Basic_Info,
    interest_cell_type = interest_cell_type,
    ncores = ncores
  )
})
message(print_run_stats(monitor))
# running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 118 MB

# Each cell type returns 4 elements:
# interest_scRNAseq_for_GRN, interest_gene_for_GRN, interest_TFs, interest_TGs.
base::print(interest_cell_type_data[[1]][["interest_scRNAseq_for_GRN"]])
# 519 features across 3000 samples, RNA assay with counts and data layers.

base::print(utils::head(interest_cell_type_data[[1]][["interest_gene_for_GRN"]]))
# ABCA1, ABCC3, ABHD5, ABLIM1, ABTB2, ACSL1.

base::print(utils::head(interest_cell_type_data[[1]][["interest_TFs"]]))
# AHR, ARID5B, BACH2, BCL11B, CDK6, CREB5.

base::print(utils::head(interest_cell_type_data[[1]][["interest_TGs"]]))
# ABCA1, ABCC3, ABHD5, ABLIM1, ABTB2, ACSL1.

qs::qsave(interest_cell_type_data, file = "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.qsave")




##################################################################################
# 2.1.2 Pseudotemporal Ordering and Branching Analysis of Cell Types of Interest #
##################################################################################

# [order_pseudotime_and_divide_branches] builds a Monocle2 DDRTree trajectory per cell type,
# saves the CellDataSet (CDS) and trajectory plots, then asks the user to read the branch structure from those plots.
# This is an interactive function: run it in the console, not in a non-interactive Rscript job.
?order_pseudotime_and_divide_branches


### Load inputs.
setwd(Working_directory)
interest_cell_type_data = qs::qread("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.qsave")


### Choose DDRTree parameters.
# TRUE  -> use DDRTree defaults (demonstrated here);
# FALSE -> forward the manual arguments ([maxIter], [sigma], [lambda], [ncenter], [param.gamma], [tol]) to monocle::reduceDimension.
auto_param_selection = TRUE
# Tip: tune ncenter in manual mode to change the number of branches
# (e.g. ncenter = 4 often yields a single branch on small datasets).


### Run and answer the prompts using the saved plots.
library(DDRTree)
monitor = peakRAM::peakRAM({
  interest_cell_type_Branches = order_pseudotime_and_divide_branches(
    interest_cell_type_data = interest_cell_type_data,
    auto_param_selection = auto_param_selection
  )
  # Manual alternative:
  # interest_cell_type_Branches = order_pseudotime_and_divide_branches(
  #   interest_cell_type_data = interest_cell_type_data,
  #   auto_param_selection = FALSE, ncenter = 4
  # )
})

# Plots are saved in ".../2.1 .../<CellType>/":
#   interest_Pseudotime_in_trajectory.png   -> darkest color marks the trajectory start (lowest pseudotime);
#   interest_AllCellState_in_trajectory.png -> numbers identify every cell state.

# One console interaction per cell type:
#   1) enter the number of branches;
#   2) for each branch, enter the space-separated state IDs along that branch, from root to tip.

# A root state can be shared by several branches. For Mono enter:
#   Please enter the number of branches ... : 2
#   Please enter the state IDs for branch 1 : 1 2
#   Please enter the state IDs for branch 2 : 1 3

# State 1 is the shared root; states 2 and 3 are the two terminal branches.

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

message(print_run_stats(monitor))
# running: 0 hours 3 mins 0 secs, peak RAM: 1 GB 462 MB


### Inspect: a named list per cell type; each branch is a data frame of cells sorted by ascending pseudotime, with columns Pseudotime and State.
base::print(utils::head(interest_cell_type_Branches[[1]][[1]], n = 3)) # Branch 1
#                    Pseudotime State
# GACTATTCATGTTGGC-1 0.00000000     1
# GGACGGATCTCCTCTT-1 0.04467413     1
# TGCATCCTCCTTGCAC-1 0.05221350     1


### Save.
qs::qsave(interest_cell_type_Branches, file = "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_Branches.qsave")




##############################################################
#                                                            #
# 2.2 Cell Grouping and Group-Wise Information Specification #
#                                                            #
##############################################################



####################################################################
# 2.2.1 Retrieving Pseudotemporal Expression Vectors for Each Gene #
####################################################################

# [get_genes_pseudotime_info] scores genes on each branch from four expression features, filters genes/cells, orders cells by pseudotime,
# and returns the log-normalized TFE and TGE matrices (genes x cells, columns ordered by time).
?get_genes_pseudotime_info


### Load inputs.
setwd(Working_directory)
interest_cell_type_Branches = qs::qread("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_Branches.qsave")
interest_cell_type_data = qs::qread("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.qsave")
ncores = parallel::detectCores() - 1
# ncores = 1 # Windows.


### Gene-score quantile cutoffs.
# Each gene receives an importance score from four features:
# non-zero expression ratio (NER), its rank (NERR), mean in non-zero cells (NEA), and its rank (NEAR).
# A gene is kept with a STRICT comparison score > quantile(score, alpha).
# Consequences worth noting:
#   alpha = 0 keeps all but the single lowest-scoring gene (effectively no cut);
#   alpha = 1 keeps no gene.
# Here both are set to 0; sparsity_cut below is used as the main gene filter.
alpha_TF = 0 # Defaults to 0.63; set to 0 in this demo.
alpha_TG = 0 # Defaults to 0.69; set to 0 in this demo.


### Cell cutoff: drop cells whose zero fraction across retained genes exceeds alpha_cell. 1 disables cell filtering.
alpha_cell = 1


### Optional detection-frequency filter: [sparsity_cut].
# A nested list named by cell type.
# Each element is a data frame whose rows are TFE and TGE and whose columns are branch1, branch2, ....
# A value is the MINIMUM required non-zero-cell fraction for a gene on that branch;
# genes below it are removed after the TFE/TGE matrices are built.
sparsity_cut = base::list()
sparsity_cut[["Mono"]] = base::as.data.frame(base::matrix(
  0.15, nrow = 2,
  ncol = base::length(interest_cell_type_Branches[["Mono"]]),
  dimnames = base::list(
    base::c("TFE", "TGE"),
    base::paste0("branch", base::seq_along(interest_cell_type_Branches[["Mono"]]))
  )
))
# 0.15 is used only to speed up the demo and to show this filtering step;
# in a real analysis choose the cutoffs from the research question (0 disables it).

base::print(sparsity_cut[[1]])
#     branch1 branch2
# TFE    0.15    0.15
# TGE    0.15    0.15


### Run, inspect and save.
monitor = peakRAM::peakRAM({
  interest_cell_type_genes_pseudotime_info = get_genes_pseudotime_info(
    interest_cell_type_Branches = interest_cell_type_Branches,
    interest_cell_type_data = interest_cell_type_data,
    alpha_TF = alpha_TF, alpha_TG = alpha_TG, alpha_cell = alpha_cell,
    ncores = ncores,
    sparsity_cut = sparsity_cut
  )
})
message(print_run_stats(monitor))
# running: 0 hours 0 mins 3 secs, peak RAM: 0 GB 302 MB

# Cells retained per branch after filtering.
base::print(base::unlist(interest_cell_type_genes_pseudotime_info[[1]][["Branches_n_cell"]])) # 1625 2659.

# TFE / TGE previews: rows are genes, integer columns are pseudotime positions.
base::print(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TFE"]][[1]][1:3, 1:6]) # Branch 1
#        1        2 3        4        5        6
# AHR    0 0.000000 0 2.971364 0.000000 0.000000
# ARID5B 0 3.394199 0 3.359603 2.097141 2.161078
# BACH2  0 0.000000 0 0.000000 0.000000 0.000000

base::print(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TGE"]][[1]][1:3, 1:6]) # Branch 1
#       1 2        3 4 5        6
# ABCA1 0 0 2.873253 0 0 0.000000
# ABCC3 0 0 0.000000 0 0 2.161078
# ABHD5 0 0 0.000000 0 0 2.161078

# Full branch annotation with cell-pseudotime index mapping.
base::print(utils::head(interest_cell_type_genes_pseudotime_info[[1]][["Branches"]][[1]], n = 3)) # Branch 1
#                    Pseudotime State               cell k~
# GACTATTCATGTTGGC-1 0.00000000     1 GACTATTCATGTTGGC-1  1
# GGACGGATCTCCTCTT-1 0.04467413     1 GGACGGATCTCCTCTT-1  2
# TGCATCCTCCTTGCAC-1 0.05221350     1 TGCATCCTCCTTGCAC-1  3

# Retained TF names per branch.
base::print(utils::head(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TFs"]][[1]])) # Branch 1
# [1] "AHR"    "ARID5B" "BACH2"  "BCL11B" "CDK6"   "CREB5"

# Retained TG names per branch.
base::print(utils::head(utils::head(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TGs"]][[1]]))) # Branch 1
# [1] "ABCA1" "ABCC3" "ABHD5" "ACSL1" "ACTB"  "ACTG1"

qs::qsave(interest_cell_type_genes_pseudotime_info, file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")



#####################################################################################
# 2.2.2 Grouping Cells and Calculating Valid (Non-Zero) Expression Values per Group #
#####################################################################################

# [cell_grouping] divides pseudotime-ordered cells of each branch into groups and averages non-zero expression per group,
# producing compact time series for the dynamical model.
?cell_grouping


### Load inputs.
setwd(Working_directory)
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")

# Number of cores for parallel computation. Defaults to 1.
# Parallel processing only works on Unix-like systems; set to 1 on Windows.
ncores = parallel::detectCores() - 1
# ncores = 1 # Windows.


### Calibrate the minimum group size.
# Three (x, y) points fit y = c*ln(a*x+b), where x is a branch's total cell count and y is the minimum cells per group;
# the rounded-up prediction is n_min. Both point vectors must be strictly increasing.
n_cell = base::round(base::mean(base::unlist(
  interest_cell_type_genes_pseudotime_info[[1]][["Branches_n_cell"]]
)))
fit_scale_factor = 1.05 # generally >0.8, <1.2
points_x_for_fitting_nc_and_nmin = base::c(
  base::round(n_cell * 0.5 * fit_scale_factor),
  base::round(n_cell * 1.05),
  base::round(n_cell * 1.5 * fit_scale_factor)
)
base::print(points_x_for_fitting_nc_and_nmin) # 1125 2249 3374.

n_min_expected = 5 # Target minimum cells per group.
points_y_for_fitting_nc_and_nmin = base::c(
  n_min_expected - 1, n_min_expected, n_min_expected + 1
)
base::print(points_y_for_fitting_nc_and_nmin) # 4 5 6.
# To force a desired n_min at the current cell count, place an (x, y) anchor at
# x = n_cell (times 1 to 1.1) and y = that n_min.


### Run, inspect and save.
monitor = peakRAM::peakRAM({
  interest_cell_type_group = cell_grouping(
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    points_x_for_fitting_nc_and_nmin = points_x_for_fitting_nc_and_nmin,
    points_y_for_fitting_nc_and_nmin = points_y_for_fitting_nc_and_nmin,
    ncores = ncores
  )
})
message(print_run_stats(monitor))
# running: 0 hours 0 mins 20 secs, peak RAM: 0 GB 286 MB

# Number of groups per branch (0 means the branch failed the non-zero constraint).
base::print(base::unlist(interest_cell_type_group[[1]][["Branches_n_group"]])) # 7 9.

# Cells in each group; groups sum to the branch cell count.
base::print(interest_cell_type_group[[1]][["Branches_n_k"]])
# [[1]] 374 318 588 79 64 69 133
# [[2]] 450 163 179 307 247 188 207 178 740

# Fitted n_min per branch.
base::print(interest_cell_type_group[[1]][["n_min"]]) # 5 6.

# Group endpoints in pseudotime (T_k) and interval lengths (Tslot_k).
base::print(interest_cell_type_group[[1]][["Branches_Tslot_k"]])

# Aggregated non-zero mean expression per group (TFE_T / TGE_T); NA means the gene is zero throughout that group.
base::print(utils::head(interest_cell_type_group[[1]][["Branches_TGE_T"]][[1]], n = 3))
#              1        2        3        4        5        6        7
# ABCA1 2.807489 3.163883 3.384037 2.749738 2.303753 2.770932 2.486290
# ABCC3 2.852150 2.915726 3.350048 2.319853 2.330691 2.585778 2.217952
# ABHD5 2.712514 2.962555 3.449351 2.635150 2.825011 2.807445 2.775495

qs::qsave(interest_cell_type_group, file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")


### How to read and adjust the grouping.
# Aim for roughly 8-15 groups per branch: fewer groups lose temporal resolution,
# while groups with very few cells weaken the per-group averages.
# Too few groups, or a branch returning 0 groups:
#   - lower the y calibration points to reduce n_min;
#   - upstream, raise alpha_TF/alpha_TG to keep fewer genes or increase the non-zero-ratio weights to reduce gene sparsity.
# Too many groups: apply the opposite adjustments.





#############################################
#############################################
###                                       ###
###   Part 3: Constructing Initial GRNs   ###
###                                       ###
#############################################
#############################################

# Part 3 builds the initial GRN (iGRN):
# a statistical edge-strength network used BOTH to decide which TG-TF pairs enter dynamical modeling and to supply their initial strengths.
# The steps are:
#   3.1.1 estimate TG-TF association by mutual information, correct the background with CLR, and rank-normalize scores per TF;
#   3.1.2 weight the edges with ChIP-Atlas ChIP-seq priors and apply dual-threshold quality control with target-gene rescue;
#   3.2   derive the binarization threshold tau for each cell type;
#   3.3   extract one iGRN per branch (subsetting only, no recomputation).




#################################################################
#                                                               #
# 3.1 Constructing an Initial GRN for Comprehensive TG-TF Pairs #
#                                                               #
#################################################################



##############################################################
# 3.1.1  Constructing an Initial GRN for Comprehensive Pairs #
##############################################################

# calculate_iGRN for each cell type,
# estimates the pairwise dependence between target genes (rows) and TFs (columns) from log-normalized expression.
?calculate_iGRN


### Load inputs.
setwd(Working_directory)
interest_cell_type_data = qs::qread("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.qsave")
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")


### Estimator settings.
# Mutual information (MI) measures nonlinear dependence;
# the Kraskov k-nearest-neighbor estimator computes it directly from cells.
k = 3 # Nearest neighbors; must be smaller than the number of cells.

noise = 1e-12 # Tiny random perturbation added to break ties in distances.
# Use 0 if reproducibility is required. This typically has no impact on single-cell data.

# The estimator is run on several noise-perturbed copies and averaged,
# because random tie-breaking makes a single run noisy.
# Reduce [replicate_time] on very large datasets.
replicate_time = if (noise == 0) 1 else 10

seed = 123 # Integer seed for reproducibility, or NULL for an unseeded run.


### Run, inspect and save.
monitor = peakRAM::peakRAM({
  interest_cell_type_iGRN_all = calculate_iGRN(
    interest_cell_type_data = interest_cell_type_data,
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    replicate_time = replicate_time,
    k = k,
    noise = noise,
    seed = seed
  )
})
message(base::ncol(interest_cell_type_iGRN_all[[1]]), " TFs, ",
        base::nrow(interest_cell_type_iGRN_all[[1]]), " TGs, ",
        base::sum(interest_cell_type_iGRN_all[[1]] > 0), " pairs, ",
        print_run_stats(monitor))
# 23 TFs, 357 TGs, 4841 pairs, running: 0 hours 0 mins 30 secs, peak RAM: 0 GB 300 MB

# CLR first z-scores every edge against its row and column background and combines the two non-negative sides in quadrature,
# which removes broad MI background.
# Scores are then rank-normalized WITHIN each TF column to (0, 1],
# so every TF's strongest target is near 1 and 0 marks an unsupported pair.
# Absolute values are therefore not comparable across different TF columns.
base::print(interest_cell_type_iGRN_all[[1]][1:3, 1:3])
#              AHR    ARID5B     BACH2
# ABCA1 0.41107872 0.3525424 0.0000000
# ABCC3 0.09912536 0.0000000 0.0000000
# ABHD5 0.32361516 0.6881356 0.6343284

qs::qsave(interest_cell_type_iGRN_all, file = "./3 get iGRN/interest_cell_type_iGRN_all.qsave")



########################################################
# 3.1.2 Filtering the Initial GRN by ChIP-Atlas Priors #
########################################################

# filter_iGRN_by_ChIP_Atlas_priors multiplies the association scores by ChIP-seq prior intensity,
# sparsifies the network with a global and a per-TF threshold,
# and rescues target genes that would otherwise lose all regulators.
?filter_iGRN_by_ChIP_Atlas_priors


### Load inputs.
setwd(Working_directory)
interest_cell_type_iGRN_all = qs::qread("./3 get iGRN/interest_cell_type_iGRN_all.qsave")
ncores = parallel::detectCores() - 1
# ncores = 1 # Windows.


### Threshold settings.
# Global cutoff:
# keep edges whose corrected score is in the top [theta_i_top_quantile] of all edges (0.15 = global top 15%).
theta_i_top_quantile = 0.15

# Per-TF cutoff:
# the top [theta_i_top_n] targets of every TF are kept even if they miss the global cutoff,
# guaranteeing each TF a minimal target set.
theta_i_top_n = 50

# After a square-root transform,
# corrected values below [theta_i_absolute_min] are set to zero before re-normalization.
theta_i_absolute_min = 0.05

# Every target gene is guaranteed this many non-zero TF regulators after QC;
# TGs below it have their strongest pre-QC edges restored. Set 0 to disable.
TG_rescue_number = 5

# Choose the quantile/top-N from the expected final network size by yourself;
# smaller values give a sparser network.


### Run, inspect and save.
monitor = peakRAM::peakRAM({
  interest_cell_type_iGRN = filter_iGRN_by_ChIP_Atlas_priors(
    interest_cell_type_iGRN = interest_cell_type_iGRN_all,
    theta_i_top_quantile = theta_i_top_quantile,
    theta_i_top_n = theta_i_top_n,
    theta_i_absolute_min = theta_i_absolute_min,
    TG_rescue_number = TG_rescue_number,
    ncores
  )
})
message("top_q=", base::round(theta_i_top_quantile, 4),
        ", top_n=", theta_i_top_n, ", ",
        base::ncol(interest_cell_type_iGRN[[1]]), " TFs, ",
        base::nrow(interest_cell_type_iGRN[[1]]), " TGs, ",
        base::sum(interest_cell_type_iGRN[[1]] > 0), " pairs, ",
        print_run_stats(monitor))
# top_q=0.15, top_n=50, 23 TFs, 289 TGs, 1284 pairs, running: 0 hours 0 mins 1 secs, peak RAM: 0 GB 45 MB

# Non-zero edges are shifted above 0.1 and capped below 1; zero edges stay exactly zero.
# All-zero rows and columns are removed.
base::print(interest_cell_type_iGRN[[1]][1:3, 1:3])
#             AHR    ARID5B BACH2
# ABCA1 0.0000000 0.4903921     0
# ABCC3 0.2049325 0.0000000     0
# ABHD5 0.0000000 0.0000000     0

# Target counts per TF (how many non-zero targets each TF keeps).
res = base::t(base::apply(interest_cell_type_iGRN[[1]], 2, function(x) {
  nz = x[x != 0]
  base::c(nonzero_count = base::length(nz),
          nonzero_min = if (base::length(nz)) base::min(nz) else NA)
}))
base::summary(res[, "nonzero_count"])
# Min. 1st Qu.  Median  Mean  3rd Qu.   Max.
# 4.00  42.00    52.00 55.83   69.50  103.00

qs::qsave(interest_cell_type_iGRN, file = "./3 get iGRN/interest_cell_type_iGRN.qsave")




############################################
#                                          #
# 3.2 Setting the Regulatory Threshold tau #
#                                          #
############################################

# get_cell_type_tau derives one binarization threshold per cell type from the weakest non-zero edge in the refined iGRN.
?get_cell_type_tau


### Load inputs, run and save.
setwd(Working_directory)
interest_cell_type_iGRN = qs::qread("./3 get iGRN/interest_cell_type_iGRN.qsave")

monitor = peakRAM::peakRAM({
  interest_cell_type_tau = get_cell_type_tau(
    interest_cell_type_iGRN = interest_cell_type_iGRN
  )
})
message("tau is ", base::round(interest_cell_type_tau[1], 4), ", ",
        print_run_stats(monitor))
# tau is 0.0993, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 1 MB

qs::qsave(interest_cell_type_tau, file = "./3 get iGRN/interest_cell_type_tau.qsave")




####################################
#                                  #
# 3.3 Constructing the Branch iGRN #
#                                  #
####################################

# get_branch_iGRN maps the cell-type iGRN onto each branch.
# It only subsets the existing matrix to the genes present in that branch;
# it never recomputes any regulatory strength.
?get_branch_iGRN


### Load inputs.
setwd(Working_directory)
interest_cell_type_iGRN = qs::qread("./3 get iGRN/interest_cell_type_iGRN.qsave")
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")
ncores = parallel::detectCores() - 1
# ncores = 1 # Windows.


### Run, inspect and save.
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_iGRN = get_branch_iGRN(
    interest_cell_type_iGRN = interest_cell_type_iGRN,
    interest_cell_type_group = interest_cell_type_group,
    ncores = ncores
  )
})
message(
  base::paste(base::unlist(base::lapply(interest_cell_type_branch_iGRN[[1]], base::ncol)), collapse = "/"),
  " TFs, ",
  base::paste(base::unlist(base::lapply(interest_cell_type_branch_iGRN[[1]], base::nrow)), collapse = "/"),
  " TGs, ",
  base::paste(base::sapply(
    interest_cell_type_branch_iGRN[[1]],
    function(df) base::sum(df > 0, na.rm = TRUE)
  ), collapse = "/"),
  " pairs, ", print_run_stats(monitor)
)
# 22/15 TFs, 280/231 TGs, 1208/923 pairs, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 2 MB

# Structure: cell type -> list of branch iGRNs in branch order.
# Each branch matrix keeps the strengths inherited from the full network for matched pairs;
# pairs absent from the branch stay zero, and all-zero rows/columns are removed.
base::print(interest_cell_type_branch_iGRN[[1]][[1]][1:3, 1:3]) # Branch 1.
#             AHR    ARID5B BACH2
# ABCA1 0.0000000 0.4925925     0
# ABCC3 0.2047438 0.0000000     0
# ABHD5 0.0000000 0.0000000     0

qs::qsave(interest_cell_type_branch_iGRN, file = "./3 get iGRN/interest_cell_type_branch_iGRN.qsave")





##########################################################################################
##########################################################################################
###                                                                                    ###
###   Part 4: Developing DARI-Based GRN Regression Predictive Models for Each Branch   ###
###                                                                                    ###
##########################################################################################
##########################################################################################

# Part 4 turns the discrete dynamical model into fitted prediction models,
# one per high-confidence TF-TG pair:
#   4.1 select the training pairs and attach per-group TFE/TGE features;
#   4.2 build the shared parameter template (initial values and bounds);
#   4.3 train every pair with a genetic algorithm (global search) followed by gradient ascent (local refinement),
#       remove low-quality fits, and summarize training behavior.




##################################
#                                #
# 4.1 Partitioning Training Sets #
#                                #
##################################

# select_training_set picks high-confidence pairs from the branch iGRNs,
# and, for each branch, appends the group-level TFE/TG expression features needed for fitting.
?select_training_set


### Load inputs.
setwd(Working_directory)
interest_cell_type_branch_iGRN = qs::qread("./3 get iGRN/interest_cell_type_branch_iGRN.qsave")
interest_cell_type_tau = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")
ncores = parallel::detectCores() - 1
# ncores = 1 # Windows.


### Build the initial candidate set.
# Separately for each TF (its strongest targets) and each TG (its strongest regulators),
# keep the top N interactions; the union of the two selections is the initial training set.
# Use 1 in most cases to limit computation.
top_n_for_max_theta_in_each_regulon = 1


### Optional strength filters (defaults keep everything).
# Within the initial set, keep this fraction by strength (1 - rate quantile).
max_theta_retention_rate_in_tops = 1

# Over the global candidate set, keep this fraction by strength.
max_theta_retention_rate_in_all = 1

# Absolute cap: keep at most this many strongest interactions.
retention_number_of_tops = 9999999 # Effectively no cap.


### Optional mutual-top refinement.
# [refine] keeps only pairs that are mutually top-ranked
# (the TF is among the TG's strongest regulators AND the TG is among the TF's strongest targets).
refine = TRUE

# [keep_all_TF] adds back the strongest edge of any TF lost in refinement so each TF keeps at least one candidate.
keep_all_TF = TRUE

# Fraction retained by strength after refinement.
top_ratio_in_refine = 1


### Run, inspect and save.
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_training_set = select_training_set(
    interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
    interest_cell_type_tau = interest_cell_type_tau,
    interest_cell_type_group = interest_cell_type_group,
    top_n_for_max_theta_in_each_regulon = top_n_for_max_theta_in_each_regulon,
    max_theta_retention_rate_in_tops = max_theta_retention_rate_in_tops,
    max_theta_retention_rate_in_all = max_theta_retention_rate_in_all,
    retention_number_of_tops = retention_number_of_tops,
    refine = refine, keep_all_TF = keep_all_TF,
    top_ratio_in_refine = top_ratio_in_refine,
    ncores = ncores
  )
})
message(
  base::paste(base::lapply(
    interest_cell_type_branch_training_set[[1]], base::nrow
  ), collapse = ","),
  " training pairs, ", print_run_stats(monitor)
)
# 23,14,24 training pairs, running: 0 hours 0 mins 1 secs, peak RAM: 0 GB 23 MB

# Recommended training-set size: at least ~50 pairs to ensure a meaningful scale,
# and at most ~200 pairs to limit runtime.
# The demo set is small only because earlier steps filtered aggressively for speed.
# You can also add or remove training pairs manually before training.

# Each branch table has columns TF, TG, theta_i,
# followed by interleaved group features TFE_1, TGE_1, ..., TFE_K, TGE_K (K = number of groups in the branch;
# NA when a gene has no non-zero value in that group).
base::print(interest_cell_type_branch_training_set[[1]][["branch1"]][1:2, 1:8])
#                     TF    TG   theta_i    TFE_1    TGE_1    TFE_2    TGE_2    TFE_3
# AHR_to_PLD1        AHR  PLD1 0.9104743 2.856700 2.683767 3.555180 3.088929 3.837067
# ARID5B_to_ACSL1 ARID5B ACSL1 0.9126465 2.862271 2.858980 3.042618 3.558745 3.343183

qs::qsave(interest_cell_type_branch_training_set, file = "./4 Build Prediction Model/interest_cell_type_branch_training_set.qsave")




################################
#                              #
# 4.2 Initializeing Parameters #
#                              #
################################

# [initialize_params] builds the constant parameter template shared by every pair.
# It takes no input; each trainable scalar gets an initial value and bounds.
?initialize_params


### Run, inspect and save.
setwd(Working_directory)
monitor = peakRAM::peakRAM({initial_params = initialize_params()})
message(print_run_stats(monitor))
# running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 0 MB

options(scipen = 999)
base::print(initial_params)
#         initialization         upper     lower
# r2                1.00      1.113383  0.370506
# alpha1            1.00     50.000000  0.100000
# alpha2            1.00     50.000000  0.100000
# beta1             0.00      5.000000 -5.000000
# beta2             0.00      5.000000 -5.000000
# s                 1.01      4.900000  1.000010
# u1i               1.00 777600.000000  0.100000
# u1j               1.00 777600.000000  0.100000
# u2i               1.00      4.900000  0.100000
# u2j               1.00      4.900000  0.100000
# v_tilde           0.50      0.900000  0.100000

# Meaning of the 11 rows:
#   r2            : steepness of the effect-efficiency rise above tau;
#   alpha1/alpha2 : conversion factors of the net-growth and feedback terms;
#   beta1/beta2   : linear correction terms of the two update equations;
#   s             : shape parameter of the feedback effect magnitude;
#   u1i/u1j       : apparent dissociation constants of the TF (i) and target (j) Hill functions (half-saturation expression levels);
#   u2i/u2j       : Hill coefficients of the TF (i) and target (j);
#   v_tilde       : environmental demand saturation.
# Bounds are preset from biological interpretability and numerical stability
# and normally do not need adjustment.

qs::qsave(initial_params, file = "./4 Build Prediction Model/initial_params.qsave")




#######################################
#                                     #
# 4.3 Training Models for Each Branch #
#                                     #
#######################################

# [model_train] trains every pair independently.
# Within each epoch it runs two stages: a genetic algorithm for global exploration,
# then gradient ascent for local refinement from the GA's best point;
# the refined point seeds the next epoch's GA population.
?model_train


### Load inputs.
setwd(Working_directory)
interest_cell_type_branch_training_set = qs::qread("./4 Build Prediction Model/interest_cell_type_branch_training_set.qsave")
initial_params = qs::qread("./4 Build Prediction Model/initial_params.qsave")
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")
interest_cell_type_tau = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")
ncores = parallel::detectCores() - 1
# ncores = 1 # Windows.


### Epoch and early-stopping settings.
max_epochs = 300 # Hard epoch cap; early stopping usually ends training earlier.
early_stop = 10 # Stop after this many epochs without improvement (use 5 to speed up).
eps_theta = 1e-3 # Predicted strength is considered converged within this of theta_i.
eps_loss = 1e-5 # Smaller loss improvement than this is negligible.


### Genetic algorithm settings.
popSize_ga = 50 # Population size; larger helps escape local optima at higher cost.
maxiter_ga = 150 # Generations per GA run.
pcrossover_ga = 0.5 # Crossover probability.
pmutation_ga = 0.5 # Mutation probability.
parallel_ga = FALSE # Whether the GA itself runs in parallel.
seed_ga = 123 # GA seed for reproducibility.


### Gradient ascent settings.
iterations_grad = 30 # Gradient steps per epoch.
ncores_grad = 1 # Cores for parallel gradient evaluation.


### Run training.
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_model_train_all = model_train(
    interest_cell_type_branch_training_set = interest_cell_type_branch_training_set,
    initial_params = initial_params,
    interest_cell_type_group = interest_cell_type_group,
    interest_cell_type_tau = interest_cell_type_tau,
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    max_epochs = max_epochs, early_stop = early_stop,
    eps_theta = eps_theta, eps_loss = eps_loss,
    popSize_ga = popSize_ga, maxiter_ga = maxiter_ga,
    pcrossover_ga = pcrossover_ga, pmutation_ga = pmutation_ga,
    parallel_ga = parallel_ga, seed_ga = seed_ga,
    iterations_grad = iterations_grad, ncores_grad = ncores_grad,
    ncores = ncores
  )
})
message(print_run_stats(monitor))
# running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 0 MB

# Stopping conditions are checked in priority order:
# max epochs reached,
# no improvement for early_stop epochs,
# predicted strength converged to theta_i,
# and negligible loss improvement.

# Branch integration: binary edge calls use union logic,
# while continuous strengths and losses are averaged across branches weighted by branch cell counts.

qs::qsave(interest_cell_type_branch_model_train_all, file = "./4 Build Prediction Model/interest_cell_type_branch_model_train_all.qsave")


### Remove low-quality fits.
interest_cell_type_branch_model_train = interest_cell_type_branch_model_train_all
for (cell_type in base::names(interest_cell_type_branch_model_train)) {
  df = interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]]

  # Demo rule: drop pairs whose integrated strength is below tau.
  bad_flag = base::with(df, theta_p < interest_cell_type_tau[cell_type])

  # Alternative, stricter rule (uncomment and choose thresholds as needed):
  # bad_flag = base::with(df, theta_p < interest_cell_type_tau[cell_type] | theta_difference >= 0.1 | fit_loss >= 0.2)

  rownames_filtered = base::rownames(df)[bad_flag]
  interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]] = df[!bad_flag, ]
  for (n in base::seq_len(base::length(interest_cell_type_branch_model_train[[cell_type]]) - 1)) {
    lst = interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]]
    interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]] =
      lst[!base::names(lst) %in% rownames_filtered]
  }
}

qs::qsave(interest_cell_type_branch_model_train, file = "./4 Build Prediction Model/interest_cell_type_branch_model_train.qsave")


### Training diagnostics (optional but recommended before Part 5).
# Helper: print min/quartiles/median/max, mean and variance of a numeric vector.
print_stats = function(name, x, digits = 6) {
  q = base::round(stats::quantile(x, na.rm = TRUE), digits)
  message(base::sprintf(
    "[%s] Min: %s, Q1: %s, Median: %s, Q3: %s, Max: %s, Mean: %s, Var: %s",
    name, q[1], q[2], q[3], q[4], q[5],
    base::round(base::mean(x, na.rm = TRUE), digits),
    base::round(stats::var(x, na.rm = TRUE), digits)))
}

# 1) Distribution of each fitted parameter across regulons and branches.
pa_init = base::round(initial_params, digits = 3)
res_list = base::list()
for (pa in base::names(interest_cell_type_branch_model_train[[1]][["branch1"]][[1]][["best_params"]])) {
  for (n in base::seq_len(base::length(interest_cell_type_branch_model_train[[1]]) - 1)) {
    pa_values = base::round(base::as.numeric(base::sapply(
      interest_cell_type_branch_model_train[[1]][[base::paste0("branch", n)]],
      function(x) x[["best_params"]][[pa]]
    )), digits = 3)
    pa_quants = base::round(stats::quantile(pa_values), digits = 3)
    res_list[[base::paste0(pa, "-Branch", n)]] = base::c(
      `Lower` = pa_init[pa, "lower"],
      Min = pa_quants["0%"], Median = pa_quants["50%"], Max = pa_quants["100%"],
      `Upper` = pa_init[pa, "upper"],
      Var = base::round(stats::var(pa_values), digits = 3)
    )
  }
}
options(scipen = 999)
base::print(base::as.data.frame(base::do.call(base::rbind, res_list)))
#                  Lower     Min.0% Median.50%   Max.100%      Upper             Var
# r2-Branch1       0.371      0.372      0.618      1.064      1.113           0.048
# r2-Branch2       0.371      0.380      0.628      0.939      1.113           0.038
# alpha1-Branch1   0.100      0.278     40.496     49.890     50.000         215.324
# ...

# 2) Integrated prediction error and fit loss.
best_result = interest_cell_type_branch_model_train[[1]][["best_pred_result"]]
print_stats("theta_difference", best_result[["theta_difference"]])
# [theta_difference] Min: 0.000795, Q1: 0.015089, Median: 0.059218, Q3: 0.078402, Max: 0.110823, Mean: 0.050488, Var: 0.001117

print_stats("fit_loss", best_result[["fit_loss"]])
# [fit_loss] Min: 0.024442, Q1: 0.064304, Median: 0.100517, Q3: 0.151597, Max: 0.200868, Mean: 0.105471, Var: 0.002558

# 3) Epochs actually used per branch (and overall).
# 4) Training time per branch (minutes); parse "secs/minutes/hours" strings.
parse_time_to_mins = function(time_str) {
  num_val = base::as.numeric(base::sub("^([0-9.]+).*", "\\1", base::trimws(time_str)))
  unit = base::sub("^[0-9.]+\\s*(.*)$", "\\1", base::trimws(time_str))
  num_val[unit == "secs"] = num_val[unit == "secs"] / 60
  num_val[unit == "hours"] = num_val[unit == "hours"] * 60
  return(num_val)
}

# Build per-branch summaries of best_epoch and training time.
epoch_rows = base::list(); time_rows = base::list()
for (n in base::seq_len(base::length(interest_cell_type_branch_model_train[[1]]) - 1)) {
  epoch_v = best_result[[base::paste0("best_epoch_branch", n)]]
  epoch_v = epoch_v[epoch_v > 0]
  epoch_rows[[base::paste0("Branch", n)]] = base::c(
    Median = stats::median(epoch_v), Max = base::max(epoch_v),
    Mean = base::mean(epoch_v)
  )
  time_v = parse_time_to_mins(best_result[[base::paste0("train_time_branch", n)]])
  time_v = time_v[time_v > 0]
  time_rows[[base::paste0("Branch", n)]] = base::c(
    Median = stats::median(time_v), Max = base::max(time_v),
    Mean = base::mean(time_v)
  )
}
base::print(base::do.call(base::rbind, epoch_rows))
#         Median Max      Mean
# Branch1    5.0  36  8.217391
# Branch2   13.5  72 19.500000

base::print(base::do.call(base::rbind, time_rows))
#         Median  Max     Mean
# Branch1   1.10 2.89 1.240681
# Branch2   1.26 5.17 1.731786

# 5) Why training stopped (counts and frequencies per branch).
cause_order = base::c(
  "loss improvement is below the threshold",
  "no loss improvement for consecutive early_stop epochs",
  "predicted theta_p converges to observed theta_i"
)
for (n in base::seq_len(base::length(interest_cell_type_branch_model_train[[1]]) - 1)) {
  cause_v = best_result[[base::paste0("stop_cause_branch", n)]]
  base::print(base::round(base::prop.table(base::table(
    factor(cause_v, levels = cause_order))), 3))
}
# loss improvement is below the threshold no loss improvement for consecutive early_stop epochs       predicted theta_p converges to observed theta_i
#                                   0.304                                                 0.652                                                 0.043
# loss improvement is below the threshold no loss improvement for consecutive early_stop epochs       predicted theta_p converges to observed theta_i
#                                   0.643                                                 0.357                                                 0.000





##################################
##################################
###                            ###
###   Part 5: Predicting GRN   ###
###                            ###
##################################
##################################

# Part 5 applies those trained models to EVERY candidate TG-TF pair and turns the predictions into the final network:
#   5.1 assemble per-branch candidate pairs with TFE/TGE features;
#   5.2 predict each pair under every trained model (strength and loss);
#   5.3 combine the lowest-loss models per pair and normalize within each TF;
#   5.4 apply an adaptive, feature-adjusted threshold with TF rescue to obtain the final loose and strict GRNs.




######################################
#                                    #
# 5.1 Getting Data for Inferring GRN #
#                                    #
######################################

# [get_data_for_inferring_GRN] expands every TG-TF pair of the branch iGRNs,
# and attaches the same group-level TFE/TGE features used in training.
?get_data_for_inferring_GRN


### Load inputs.
setwd(Working_directory)
interest_cell_type_branch_iGRN = qs::qread("./3 get iGRN/interest_cell_type_branch_iGRN.qsave")
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")
ncores = parallel::detectCores() - 1
# ncores = 1 # Windows.

# Partitions used for parallel feature binding; set to 1 to run serially.
n_part = ncores


### Run, inspect and save.
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_data_for_GRN_infer = get_data_for_inferring_GRN(
    interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
    interest_cell_type_group = interest_cell_type_group,
    n_part = n_part, ncores = ncores
  )
})
message(print_run_stats(monitor))
# running: 0 hours 0 mins 11 secs, peak RAM: 0 GB 56 MB

# Branch tables: one row per pair, columns TG, TF, theta_i and TFE_K, TGE_K for each group.
# Pairs with theta_i = 0 (no iGRN support) are removed.
base::print(interest_cell_type_branch_data_for_GRN_infer[[1]][["branch1"]][1:2, 1:8])
#                 TG  TF   theta_i  TFE_1   TGE_1   TFE_2    TGE_2    TFE_3
# AHR_to_ABCC3 ABCC3 AHR 0.2049325 2.8567 2.85215 3.55518 2.915726 3.837067
# AHR_to_ACTB   ACTB AHR 0.4978267 2.8567 6.07083 3.55518 5.817856 3.837067

# The GRN element is the global, de-duplicated pair set across branches, with only TG, TF and theta_i (no group features).
base::print(utils::head(interest_cell_type_branch_data_for_GRN_infer[[1]][["GRN"]], 2))
#                 TG  TF   theta_i
# AHR_to_ABCC3 ABCC3 AHR 0.2047438
# AHR_to_ACTB   ACTB AHR 0.4971116

qs::qsave(interest_cell_type_branch_data_for_GRN_infer, file = "./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.qsave")




################################
#                              #
# 5.2 Inferring the Branch GRN #
#                              #
################################

# [infer_GRN] applies every trained model of a branch to each candidate pair and re-optimizes the regulatory strength for that pair,
# returning the strength and loss under each model.
?infer_GRN


### Load inputs.
setwd(Working_directory)
interest_cell_type_branch_data_for_GRN_infer = qs::qread("./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.qsave")
interest_cell_type_tau = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")
interest_cell_type_branch_model_train = qs::qread("./4 Build Prediction Model/interest_cell_type_branch_model_train.qsave")
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")
ncores = parallel::detectCores() - 1
# ncores = 1 # Windows.


### Way 1: one-step inference (simplest; use when memory is sufficient).
# Chunk size: candidate pairs are split into chunks of at most `each`,
# which controls peak memory. Size chunks to spread work evenly across cores.
each = base::ceiling(base::max(base::sapply(
  interest_cell_type_branch_data_for_GRN_infer[[1]][
    -base::length(interest_cell_type_branch_data_for_GRN_infer[[1]])],
  base::nrow
)) / ncores)
base::print(each) # Example: 16.

monitor = peakRAM::peakRAM({
  interest_cell_type_pGRN = infer_GRN(
    interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
    interest_cell_type_tau = interest_cell_type_tau,
    interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
    interest_cell_type_group = interest_cell_type_group,
    select_cell_type = NULL, select_branch = NULL,
    each = each, ncores = ncores
  )
})
message(print_run_stats(monitor))
# running: 0 hours 0 mins 26 secs, peak RAM: 0 GB 12 MB

# For large datasets with many branches this can exhaust memory; Way 2 is safer.


### Way 2: branch-by-branch inference (stable; more steps).
# Process one branch at a time to avoid memory accumulating across branches.
interest_cell_type_pGRN = base::list()
cell_type = base::names(interest_cell_type_branch_data_for_GRN_infer)[1] # Mono
interest_cell_type_pGRN[[cell_type]] = base::list()
for (n in base::seq_len(base::length(interest_cell_type_branch_data_for_GRN_infer[[cell_type]]) - 1)) {
  message("Infer cell type ", cell_type, " branch ", n, ".")
  each = base::ceiling(base::nrow(
    interest_cell_type_branch_data_for_GRN_infer[[cell_type]][[base::paste0("branch", n)]]
  ) / ncores)
  base::gc()
  monitor = peakRAM::peakRAM({
    interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]] = infer_GRN(
      interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
      interest_cell_type_tau = interest_cell_type_tau,
      interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
      interest_cell_type_group = interest_cell_type_group,
      select_cell_type = cell_type, select_branch = n,
      each = each, ncores = ncores
    )[[cell_type]][[base::paste0("branch", n)]]
  })
  message("Cell type ", cell_type, " - Branch ", n, ": ", print_run_stats(monitor))
  # Cell type Mono - Branch 1: running: 0 hours 0 mins 12 secs, peak RAM: 0 GB 9 MB
  # Cell type Mono - Branch 2: running: 0 hours 0 mins 8 secs, peak RAM: 0 GB 5 MB

  # base::save.image(base::paste0("./5 Infer GRN/5_Infer_GRN_branch", n, "middle.RData")) # Checkpoint.
  base::gc()
}

# Branch results are also saved automatically as RDS by infer_GRN.

# Each branch holds two frames, theta_p and loss:
# rows are candidate pairs and columns are trained models.
# theta_p is the strength re-optimized for that pair-model combination;
# loss is the corresponding minimum prediction loss.
base::print(interest_cell_type_pGRN[[1]][["branch1"]][["theta_p"]][1:2, 1:3])
#              AHR_to_PLD1 ARID5B_to_ACSL1 BACH2_to_BACH2
# AHR_to_ABCC3  0.98999997      0.98999997     0.01000001
# AHR_to_ACTB   0.09773147      0.09897411     0.01000001

base::print(interest_cell_type_pGRN[[1]][["branch1"]][["loss"]][1:2, 1:3])
#              AHR_to_PLD1 ARID5B_to_ACSL1 BACH2_to_BACH2
# AHR_to_ABCC3   0.4303630       0.4551814      0.4654984
# AHR_to_ACTB    0.2506862       0.2735410      0.2713641

qs::qsave(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.qsave")




################################################
#                                              #
# 5.3 Selecting Best Prediction for Each TF-TG #
#                                              #
################################################

# [select_best_prediction] collapses the many trained models into one strength and loss per pair,
# integrates branches, normalizes within each TF, and attaches the observed strength and gene expression features.
?select_best_prediction


### Load inputs.
setwd(Working_directory)
interest_cell_type_pGRN = qs::qread("./5 Infer GRN/interest_cell_type_pGRN.qsave")
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")
interest_cell_type_iGRN = qs::qread("./3 get iGRN/interest_cell_type_iGRN.qsave")
Basic_Info = qs::qread("./1.3 Data Preprocessing - Filter Expression Matrix/Basic_Info.qsave")
interest_cell_type_tau = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")
interest_cell_type_branch_data_for_GRN_infer = qs::qread("./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.qsave")
ncores = parallel::detectCores() - 1
# ncores = 1 # Windows.

each = base::ceiling(base::max(base::sapply(
  interest_cell_type_branch_data_for_GRN_infer[[1]][
    -base::length(interest_cell_type_branch_data_for_GRN_infer[[1]])],
  base::nrow
)) / ncores)


### Number of models averaged per pair.
# The count is round(n_train * topn_ratio), clamped to [topn_min, topn_max].
topn_ratio = 0.1
topn_min = 5
topn_max = 25


### Weight settings.
# [fix_factor] multiplies the weight of boundary strengths (below 0.0101 or above 0.9899, which reflect the optimization bounds);
# 1 applies no down-weighting.
fix_factor = 1.0

# [min_scale] multiplies the smallest positive value when setting the min-max lower bound;
# 0.99 keeps the weakest positive edge slightly above tau.
min_scale = 0.99


### Run, inspect and save.
monitor = peakRAM::peakRAM({
  interest_cell_type_pGRN = select_best_prediction(
    interest_cell_type_pGRN = interest_cell_type_pGRN,
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    interest_cell_type_iGRN = interest_cell_type_iGRN,
    Basic_Info = Basic_Info,
    interest_cell_type_tau = interest_cell_type_tau,
    topn_max = topn_max, topn_min = topn_min, topn_ratio = topn_ratio,
    fix_factor = fix_factor, min_scale = min_scale,
    each = each, ncores = ncores
  )
})
message(base::nrow(interest_cell_type_pGRN[[1]][["pred_summary"]]),
        " pairs remaining, ", print_run_stats(monitor))
# 1283 pairs remaining, running: 0 hours 0 mins 11 secs, peak RAM: 0 GB 29 MB

# Per pair, the selected models receive triangular weights from loss
# (1 for the smallest loss, 0 for the largest; equal weights if all losses tie),
# and their strengths and losses are weighted means.
# Branch strength/loss are then integrated by branch cell counts, and normalized within each TF:
# positive strengths are shifted to [tau, 1) while absent pairs reset to 0;
# losses use the same min-max scaling without the shift.
base::print(utils::head(interest_cell_type_pGRN[[1]][["branch1"]][["best_result"]], 2))
#              theta_p_branch1 loss_branch1
# AHR_to_ABCC3       0.9900000    0.4214795
# AHR_to_ACTB        0.2464091    0.2415191

# [pred_summary] is the integrated cell-type table. Key column groups:
#   theta_p_branchN / loss_branchN : per-branch values;
#   FG.theta_p / FG.Loss : cell-count weighted integration;
#   FG.theta_p_norm(.R) / FG.Loss_norm(.R) : within-TF normalized values/ranks;
#   FG.theta_i(.R) : observed iGRN strength and its rank (rank 0 when theta_i = 0);
#   F.* / G.* : TF/TG expression features (non-zero ratio NZER, non-zero average NZEA, overall average EA) and normalized ranks.
base::print(utils::head(interest_cell_type_pGRN[[1]][["pred_summary"]][, base::c(
  "TF", "TG", "FG.theta_p", "FG.theta_p_norm.R", "FG.Loss_norm.R", "FG.theta_i", "FG.theta_i.R"
)], 3))
#               TF    TG FG.theta_p FG.theta_p_norm.R FG.Loss_norm.R FG.theta_i FG.theta_i.R
# AHR_to_ABCC3 AHR ABCC3  0.9900000        0.96336711      1.0000000  0.2047438   0.05529595
# AHR_to_ACTB  AHR  ACTB  0.1551572        0.01169135      0.4731099  0.4971116   0.44158879
# AHR_to_AHR   AHR   AHR  0.7029097        0.59859704      0.7560405  0.6902820   0.80062305

qs::qsave(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.qsave")




#####################
#                   #
# 5.4 Filtering GRN #
#                   #
#####################

# [filter_GRN] converts the continuous predictions into the final edge set with an adaptive threshold,
# rather than one fixed global cutoff.
?filter_GRN


### Load inputs.
setwd(Working_directory)
interest_cell_type_pGRN = qs::qread("./5 Infer GRN/interest_cell_type_pGRN.qsave")


### Anchor points for the initial threshold curve.
# IT = 1/(a theta_p + b) + c is fitted from paired anchors (predicted strength, initial threshold);
# weaker predicted strengths face a higher threshold.
theta_p_points = base::c(0, 0.5, 1)
initial_quantile_points = base::c(1.5, 0.9, 0.6)
theta_p_quantile_fit_start = base::list(a = 0.333333, b = 0.333333, c = -0.5)


### Per-feature threshold adjustment.
# Each feature has a center fb (factor_balance) and strength fr (factor_range).
# "Relaxation" lowers the threshold for favorable feature values; "penalty" raises it.
# By default only [FG.Loss_norm.R] is a penalty feature.
factor_balance = base::list(
  FG.Loss_norm.R = 0.2, FG.theta_i.R = 0.2,
  F.NZER.R = 0.2, F.NZEA.R = 0.2, F.EA.R = 0.2,
  G.NZER.R = 0.2, G.NZEA.R = 0.2, G.EA.R = 0.2
)
factor_range = base::list(
  FG.Loss_norm.R = 0.25, FG.theta_i.R = 0.25,
  F.NZER.R = 0.05, F.NZEA.R = 0.05, F.EA.R = 0.05,
  G.NZER.R = 0.15, G.NZEA.R = 0.15, G.EA.R = 0.15
)


### TF rescue.
# Ensure every TF keeps at least [min_TG_per_TF] target(s);
# rescue candidates must still pass the pooled low-quantile quality cut of [FG.theta_p_norm.R].
rescue_TF = TRUE
min_TG_per_TF = 1
TF_rescue_quantile_threshold = 0.1


### Run, inspect and save.
monitor = peakRAM::peakRAM({
  interest_cell_type_pGRN = filter_GRN(
    interest_cell_type_pGRN = interest_cell_type_pGRN,
    theta_p_points = theta_p_points,
    initial_quantile_points = initial_quantile_points,
    theta_p_quantile_fit_start = theta_p_quantile_fit_start,
    factor_balance = factor_balance,
    factor_range = factor_range,
    rescue_TF = rescue_TF,
    min_TG_per_TF = min_TG_per_TF,
    TF_rescue_quantile_threshold = TF_rescue_quantile_threshold
  )
})
message(base::nrow(interest_cell_type_pGRN[[1]][["pred_summary_loose"]]),
        " loose pairs, ",
        base::nrow(interest_cell_type_pGRN[[1]][["pred_summary_strict"]]),
        " strict pairs, ", print_run_stats(monitor))
# 562 loose pairs, 349 strict pairs, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 2 MB

# A pair is kept when its normalized strength rank reaches the adjusted threshold AT (IT multiplied by all feature factors);
# TF rescue fills TFs that would otherwise have too few targets.
#   [pred_summary_loose] : edges kept by the threshold or rescue (includes edges supported by any branch; suitable for exploration);
#   [pred_summary_strict]: complete cases of the loose set,
#                          i.e. edges with no missing branch value (core, cross-branch-conserved relationships;
#                          suitable for high-confidence analysis).
base::print(utils::head(interest_cell_type_pGRN[[1]][["pred_summary_loose"]][, base::c(
  "TF", "TG", "FG.theta_p", "keep"
)], 3))
#               TF    TG FG.theta_p keep
# AHR_to_ABCC3 AHR ABCC3  0.9900000 TRUE
# AHR_to_AHR   AHR   AHR  0.7029097 TRUE
# AHR_to_AIG1  AHR  AIG1  0.9874184 TRUE

qs::qsave(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.qsave")

# The final deliverables are [pred_summary_loose] and [pred_summary_strict].
# Each row is one directed TF -> TG edge with its integrated strength, loss and supporting features;
# export the TF/TG/strength columns to build downstream network tables or visualizations.
