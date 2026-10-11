#' @title Retrieve and Filter Transcription Factors from ChIP-Atlas Priors
#'
#' @description Retrieve transcription factor (TF) candidates from the built-in ChIP-Atlas TF list,
#' apply multi-layer expression quality filtering based on scRNA-seq data,
#' and validate them with ChIP-seq prior evidence downloaded from ChIP-Atlas.
#' Returns a curated TF list for downstream gene regulatory network construction.
#'
#' @param scRNAseq Seurat object containing single-cell RNA-seq data with a \code{counts} layer,
#' used for expression quantification and TF screening.
#' @param genome_assembly Character string of genome assembly version,
#' determining the built-in ChIP-Atlas TF list and download URLs.
#' Supported: \code{"hg38"}, \code{"hg19"}, \code{"mm10"}, \code{"mm9"}, \code{"rn6"},
#' \code{"dm6"}, \code{"dm3"}, \code{"ce11"}, \code{"ce10"}, \code{"sacCer3"}. Defaults to \code{"hg38"}.
#' @param min_TF_nonzero_freq Numeric between 0 and 1. Minimum absolute non-zero cell frequency. Defaults to \code{0.05}.
#' @param min_TF_nonzero_mean Numeric. Minimum absolute mean expression in non-zero cells. Defaults to \code{0} (no filter).
#' @param min_TF_overall_mean Numeric. Minimum absolute mean expression across all cells. Defaults to \code{0} (no filter).
#' @param min_TF_nonzero_freq_quantile Numeric between 0 and 1. Minimum genome-wide quantile of non-zero cell frequency. Defaults to \code{0} (no filter).
#' @param min_TF_nonzero_mean_quantile Numeric between 0 and 1. Minimum genome-wide quantile of mean expression in non-zero cells. Defaults to \code{0} (no filter).
#' @param min_TF_overall_mean_quantile Numeric between 0 and 1. Minimum genome-wide quantile of overall mean expression. Defaults to \code{0} (no filter).
#' @param timeout Positive numeric, time in seconds allowed for a single ChIP-Atlas file download before it is treated as failed.
#' Temporarily passed to \code{\link[base]{options}} via \code{options(timeout = ...)} during downloading and restored to the original value afterwards.
#' Defaults to \code{60}.
#' @param distance_from_TSS Integer, distance cutoff (kilobases) from TSS for ChIP-Atlas target gene annotation.
#' Allowed values: \code{1}, \code{5}, \code{10}. Defaults to \code{10}.
#' @param delay_time Numeric, delay in seconds between consecutive downloads to avoid rate limiting.
#' Defaults to \code{0.5}.
#' @param excluded_cell_types_in_priors Character vector of cell type keywords to exclude from ChIP-Atlas sample columns.
#' Set to \code{"NotExclude"} to disable filtering (default).
#' @param ignore_case_when_filtering_priors Logical, whether to ignore case when matching exclusion keywords.
#' Defaults to \code{TRUE}.
#' @param intensity_scale_factor Positive numeric, exponent for power transformation of min-max normalized ChIP-seq intensity.
#' Smaller values compress the intensity distribution more strongly.
#' Defaults to \code{0.5} (square-root transform).
#' @param prior_min_TG_number Positive integer, minimum number of valid target genes required to retain a TF prior.
#' TFs with fewer targets after filtering are discarded. Defaults to \code{2000}.
#' @param prior_top_n_keep Positive integer, maximum number of top target genes to retain per TF in the simplified prior table.
#' Defaults to \code{5000}.
#' @param prior_top_q_keep Numeric between 0 and 1, proportional cap of top target genes to retain per TF.
#' The actual retained count is \code{max(1, round(min(prior_top_n_keep, nrow(prior) * prior_top_q_keep)))}.
#' Defaults to \code{0.7}.
#'
#' @return Sorted character vector of TF gene symbols that pass all screening layers.
#'
#' @details The function implements a five-stage TF screening pipeline:
#' \enumerate{
#'   \item \strong{ChIP-Atlas candidate generation}
#'   \itemize{
#'     \item Extract TFs from the built-in ChIP-Atlas TF list for the specified \code{genome_assembly},
#'       and keep only those present in the scRNA-seq count matrix.
#'     \item Use the retained set as the initial TF candidate pool.
#'   }
#'   \item \strong{Expression quality filtering}
#'   \itemize{
#'     \item Compute three per-gene metrics from raw counts: non-zero cell frequency,
#'       mean expression in non-zero cells, and overall mean expression.
#'     \item Calculate genome-wide percentile ranks for each metric (0 = lowest, 1 = highest).
#'     \item Subset the feature table to the candidate TFs and prefix all columns with \code{TF_}.
#'     \item Retain TFs satisfying all 3 absolute thresholds and 3 quantile thresholds simultaneously.
#'   }
#'   \item \strong{ChIP-Atlas prior download with caching}
#'   \itemize{
#'     \item Download target gene tables from ChIP-Atlas for all expression-qualified TFs.
#'     \item Reuse local TSV files if present; automatically discard TFs with failed or empty downloads.
#'   }
#'   \item \strong{Prior processing and quality validation}
#'   \itemize{
#'     \item Remove the STRING database column, retain only ChIP-seq evidence columns.
#'     \item Filter out sample columns matching exclusion keywords; discard TFs with fewer than 2 valid sample columns left.
#'     \item Compute row-wise mean ChIP-seq intensity and remove entries with zero intensity.
#'     \item Discard TFs whose valid target gene count falls below \code{prior_min_TG_number}.
#'     \item Apply min-max normalization followed by power transformation to the mean intensity.
#'   }
#'   \item \strong{Prior simplification and final output}
#'   \itemize{
#'     \item Trim each valid prior to top target genes under the dual cap rule.
#'     \item Save intermediate files (full processed priors, simplified priors, download records) under the \code{ChIP_Atlas_Priors} directory.
#'     \item Extract TF names with successfully processed priors as the final output.
#'   }
#' }
#'
#' @note
#' \itemize{
#'   \item The function temporarily switches the working directory to
#'     \code{1.1 Data Preprocessing - Get TFs From Database} and restores the original path after completion.
#'   \item The ChIP-Atlas TF list is built into the package, curated from the ChIP-Atlas official resource.
#'   \item Expression quantile ranks are computed against all genes in the count matrix,
#'     providing a genome-wide baseline for quality assessment.
#'   \item Internet access is required for the first run; subsequent runs reuse cached local files automatically.
#'   \item Both full processed priors and simplified TF-target intensity tables are saved to disk for direct reuse in downstream analyses.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load scRNA-seq data
#' scRNAseq = readRDS("./scRNAseq.rds")
#'
#' # Retrieve and filter TFs with default parameters
#' TFs = get_TFs_from_priors(
#'   scRNAseq = scRNAseq,
#'   genome_assembly = "hg38",
#'   min_TF_nonzero_freq = 0.045,
#'   distance_from_TSS = 10,
#'   excluded_cell_types_in_priors = "NotExclude",
#'   prior_min_TG_number = 5000,
#'   prior_top_n_keep = 3000,
#'   prior_top_q_keep = 0.5
#' )
#'
#' saveRDS(TFs, file = "./1.1 Data Preprocessing - Get TFs From Database/TFs.rds")
#' }
get_TFs_from_priors = function(
    scRNAseq,
    genome_assembly = "hg38",
    min_TF_nonzero_freq = 0.05,
    min_TF_nonzero_mean = 0,
    min_TF_overall_mean = 0,
    min_TF_nonzero_freq_quantile = 0,
    min_TF_nonzero_mean_quantile = 0,
    min_TF_overall_mean_quantile = 0,
    timeout = 60,
    distance_from_TSS = 10,
    delay_time = 0.5,
    excluded_cell_types_in_priors = "NotExclude",
    ignore_case_when_filtering_priors = TRUE,
    intensity_scale_factor = 0.5,
    prior_min_TG_number = 2000,
    prior_top_n_keep = 5000,
    prior_top_q_keep = 0.7
) {
  ### Start.
  t_start = base::Sys.time()
  original_dir = base::getwd()
  new_folder = "1.1 Data Preprocessing - Get TFs From Database"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  message("Run: Retrieve TFs from prior database. Time: ", t_start, ".")

  ### Get candidate TFs from the built-in ChIP-Atlas list.
  message("Get TFs from ChIP-Atlas.")
  if (genome_assembly == "hg38") {
    TFs_in_ChIP_Atlas = c(
      "AATF", "ABI1", "ACSS2", "ACTR5", "ADAR", "ADNP", "ADNP2", "AEBP2", "AFF1", "AFF4", "AGO1", "AGO2",
      "AHCTF1", "AHDC1", "AHR", "AHRR", "AICDA", "AIM2", "AKAP8", "AKNA", "ALKBH1", "ALKBH3", "ALOX5",
      "ALYREF", "ANAPC1", "APC", "APOBEC3B", "AR", "ARHGAP35", "ARID1A", "ARID1B", "ARID2", "ARID3A",
      "ARID3B", "ARID4A", "ARID4B", "ARID5B", "ARNT", "ARNT2", "ARNTL", "ARNTL2", "ARRB1", "ASCL1", "ASCL2",
      "ASF1A", "ASH1L", "ASH2L", "ASPSCR1", "ASXL1", "ASXL2", "ASXL3", "ATAD3A", "ATF1", "ATF2", "ATF3", "ATF4",
      "ATF5", "ATF6", "ATF7", "ATF7IP", "ATM", "ATOH1", "ATOH8", "ATR", "ATRX", "ATXN7L3", "AUTS2", "BACH1",
      "BACH2", "BAHCC1", "BAHD1", "BANF1", "BANP", "BAP1", "BARX1", "BATF", "BATF3", "BBX", "BCL10", "BCL11A",
      "BCL11B", "BCL3", "BCL6", "BCL6B", "BCL9", "BCLAF1", "BCOR", "BCORL1", "BDP1", "BHLHE40", "BHLHE41",
      "BMI1", "BNC2", "BOD1L1", "BPTF", "BRCA1", "BRCA2", "BRD1", "BRD2", "BRD3", "BRD4", "BRD7", "BRD8",
      "BRD9", "BRDT", "BRF1", "BRF2", "BRPF1", "BRPF3", "BRWD1", "BUB3", "CAMTA1", "CAMTA2", "CARM1",
      "CASP8AP2", "CASZ1", "CBFA2T2", "CBFA2T3", "CBFB", "CBX1", "CBX2", "CBX3", "CBX4", "CBX5", "CBX6", "CBX7",
      "CBX8", "CC2D1A", "CCAR2", "CCND2", "CCNE1", "CCNT1", "CCNT2", "CD74", "CDAN1", "CDC16", "CDC45",
      "CDC5L", "CDC73", "CDK11A", "CDK12", "CDK2", "CDK6", "CDK7", "CDK8", "CDK9", "CDKN1B", "CDT1", "CDX1",
      "CDX2", "CDX4", "CDYL", "CDYL2", "CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG", "CEBPZ", "CERS6", "CGGBP1",
      "CHAF1A", "CHAMP1", "CHAT", "CHCHD3", "CHD1", "CHD2", "CHD3", "CHD4", "CHD6", "CHD7", "CHD8", "CHUK",
      "CIC", "CIITA", "CIZ1", "CLOCK", "CNBP", "CNOT3", "COBLL1", "COIL", "COPS2", "CPSF1", "CPSF2", "CPSF3",
      "CPSF4", "CPSF6", "CPSF7", "CREB1", "CREB3", "CREB3L1", "CREB3L2", "CREB5", "CREBBP", "CREBL2", "CREM",
      "CRTC1", "CRTC2", "CRX", "CRY1", "CSDC2", "CSDE1", "CSE1L", "CSNK2A1", "CSRNP3", "CSTF2", "CSTF3", "CTBP1",
      "CTBP2", "CTCF", "CTCFL", "CTNNB1", "CTR9", "CTSG", "CUL3", "CUL4A", "CUL4B", "CUX1", "CXXC1", "CXXC5",
      "DACH1", "DAXX", "DBP", "DCP1A", "DCP2", "DDIT3", "DDX11", "DDX17", "DDX20", "DDX21", "DDX39B", "DDX5",
      "DDX55", "DEAF1", "DEK", "DICER1", "DIDO1", "DKC1", "DLX4", "DLX5", "DLX6", "DMAP1", "DMC1", "DMRT2",
      "DMTF1", "DNA2", "DNAJC2", "DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "DONSON", "DOT1L", "DPF1",
      "DPF2", "DPPA3", "DR1", "DRAP1", "DROSHA", "DTL", "DUX4", "DVL1", "DYRK1A", "DZIP1", "E2F1", "E2F2", "E2F3",
      "E2F4", "E2F5", "E2F6", "E2F7", "E2F8", "E4F1", "EBF1", "EBF3", "EEA1", "EED", "EGFR", "EGR1", "EGR2", "EHF",
      "EHMT1", "EHMT2", "ELANE", "ELAVL1", "ELF1", "ELF2", "ELF3", "ELF4", "ELK1", "ELK3", "ELK4", "ELL", "ELL2",
      "EMSY", "EMX1", "EN1", "EOMES", "EP300", "EP400", "EPAS1", "EPCAM", "EPOP", "ERBB2", "ERCC2", "ERCC3",
      "ERCC6", "ERCC8", "ERF", "ERG", "ERH", "ESR1", "ESR2", "ESRRA", "ESRRB", "ESRRG", "ETS1", "ETS2", "ETV1",
      "ETV2", "ETV3", "ETV4", "ETV5", "ETV6", "EVX1", "EWSR1", "EZH1", "EZH2", "F2RL1", "FANCD2", "FANCM", "FBL",
      "FBP1", "FBXO44", "FERD3L", "FEZF1", "FGFR1", "FIP1L1", "FLI1", "FMR1", "FOS", "FOSB", "FOSL1", "FOSL2",
      "FOXA1", "FOXA2", "FOXA3", "FOXB1", "FOXC1", "FOXC2", "FOXD3", "FOXF1", "FOXF2", "FOXH1", "FOXJ2",
      "FOXJ3", "FOXK1", "FOXK2", "FOXL2", "FOXM1", "FOXN3", "FOXO1", "FOXO3", "FOXO4", "FOXP1", "FOXP2",
      "FOXP3", "FOXP4", "FOXQ1", "FOXR2", "FOXS1", "FUBP1", "FUBP3", "FUS", "FXR1", "FXR2", "G3BP1", "GABPA",
      "GABPB1", "GABPB2", "GATA1", "GATA2", "GATA3", "GATA4", "GATA6", "GATAD1", "GATAD2A", "GATAD2B",
      "GCM1", "GFI1", "GFI1B", "GGCT", "GLI1", "GLI2", "GLI4", "GLIS1", "GLIS2", "GLIS3", "GLMP", "GLYR1", "GMEB1",
      "GMEB2", "GNL3", "GPBP1L1", "GPN1", "GPS2", "GREB1", "GRHL1", "GRHL2", "GRHL3", "GRWD1", "GTF2A2",
      "GTF2B", "GTF2E2", "GTF2F1", "GTF2F2", "GTF2H2", "GTF2H2C", "GTF2I", "GTF3A", "GTF3C1", "GTF3C2",
      "GTF3C5", "GZF1", "HAND1", "HAND2", "HAT1", "HBP1", "HBZ", "HCFC1", "HDAC1", "HDAC2", "HDAC3", "HDAC4",
      "HDAC6", "HDAC8", "HDAC9", "HDGF", "HECTD1", "HES1", "HES2", "HES4", "HES5", "HES7", "HESX1", "HEXIM1",
      "HEY1", "HEY2", "HEYL", "HHEX", "HIC1", "HIC2", "HIF1A", "HIF3A", "HINFP", "HIRA", "HIVEP1", "HIVEP3",
      "HJURP", "HK2", "HLF", "HLTF", "HMBOX1", "HMG20A", "HMG20B", "HMGA1", "HMGA2", "HMGB1", "HMGB2",
      "HMGN1", "HMGN2", "HMGN3", "HMGN4", "HMGXB3", "HMGXB4", "HNF1A", "HNF1B", "HNF4A", "HNF4G",
      "HNRNPA2B1", "HNRNPC", "HNRNPH1", "HNRNPK", "HNRNPL", "HNRNPLL", "HNRNPUL1", "HOMEZ", "HOXA11",
      "HOXA13", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXB1", "HOXB13", "HOXB2", "HOXB4",
      "HOXB5", "HOXB6", "HOXB7", "HOXB8", "HOXC11", "HOXC13", "HOXC5", "HOXC8", "HOXC9", "HOXD1", "HOXD11",
      "HOXD13", "HOXD4", "HOXD9", "HSF1", "HSF2", "HSF4", "HSP90B1", "HSPA1A", "ICE1", "ICE2", "ID1", "ID2", "ID3",
      "ID4", "IFI16", "IKZF1", "IKZF2", "IKZF3", "IKZF4", "IKZF5", "ILF3", "ILK", "ING2", "ING4", "ING5", "INO80", "INSM2",
      "INSR", "INTS10", "INTS12", "INTS13", "INTS3", "INTS5", "INTS6", "IRF1", "IRF2", "IRF2BP2", "IRF3", "IRF4", "IRF5",
      "IRF8", "IRF9", "IRX2", "IRX5", "ISL1", "ISL2", "ISX", "JADE3", "JARID2", "JAZF1", "JDP2", "JMJD1C", "JMJD6", "JRK",
      "JUN", "JUNB", "JUND", "KANSL1", "KANSL3", "KAT2A", "KAT2B", "KAT5", "KAT6A", "KAT6B", "KAT7", "KAT8",
      "KCMF1", "KDM1A", "KDM1B", "KDM2A", "KDM2B", "KDM3A", "KDM3B", "KDM4A", "KDM4B", "KDM4C",
      "KDM5A", "KDM5B", "KDM5C", "KDM5D", "KDM6A", "KDM6B", "KHSRP", "KLF1", "KLF10", "KLF11", "KLF12",
      "KLF13", "KLF14", "KLF15", "KLF16", "KLF17", "KLF3", "KLF4", "KLF5", "KLF6", "KLF7", "KLF8", "KLF9", "KMT2A",
      "KMT2B", "KMT2C", "KMT2D", "L3MBTL2", "L3MBTL3", "LARP7", "LBX2", "LCOR", "LCORL", "LDB1", "LDB2",
      "LEF1", "LEO1", "LHX2", "LHX4", "LHX5", "LHX6", "LIG4", "LIN28B", "LIN54", "LIN9", "LMNA", "LMNB1", "LMO1",
      "LMO2", "LMO3", "LMTK3", "LRIF1", "LRRFIP1", "LRWD1", "LYL1", "MAF", "MAF1", "MAFB", "MAFF", "MAFG",
      "MAFK", "MAML1", "MAML3", "MAP1LC3B", "MAP2K1", "MAPK1", "MAPK14", "MATR3", "MAU2", "MAX", "MAZ",
      "MBD1", "MBD2", "MBD3", "MBD3L2", "MBD4", "MBD5", "MBD6", "MBNL3", "MCM2", "MCM3", "MCM5",
      "MCM7", "MCRS1", "MDM2", "MDN1", "MEAF6", "MECOM", "MECP2", "MED1", "MED12", "MED13", "MED16",
      "MED21", "MED22", "MED25", "MED26", "MED30", "MED31", "MED4", "MED6", "MED7", "MED8", "MED9",
      "MEF2A", "MEF2B", "MEF2C", "MEF2D", "MEIS1", "MEIS2", "MEN1", "MEOX2", "MEPCE", "MESP1", "METTL14",
      "METTL3", "MGA", "MGMT", "MIER1", "MIER2", "MIER3", "MIF", "MITF", "MIXL1", "MKI67", "MLLT1", "MLLT10",
      "MLLT3", "MLLT6", "MLX", "MLXIP", "MNT", "MNX1", "MORC2", "MPHOSPH8", "MSC", "MSH6", "MSX2", "MTA1",
      "MTA2", "MTA3", "MTF1", "MTF2", "MTHFD1", "MTOR", "MXD1", "MXD3", "MXD4", "MXI1", "MYB", "MYBBP1A",
      "MYBL1", "MYBL2", "MYC", "MYCL", "MYCN", "MYF5", "MYF6", "MYH11", "MYNN", "MYOCD", "MYOD1", "MYOG",
      "MYRF", "MZF1", "N6AMT1", "NAB2", "NABP2", "NAIF1", "NANOG", "NAT10", "NBN", "NCAPG", "NCAPG2",
      "NCAPH2", "NCBP1", "NCL", "NCOA1", "NCOA2", "NCOA3", "NCOA4", "NCOA5", "NCOA6", "NCOR1", "NCOR2",
      "NEDD4", "NEIL1", "NELFA", "NELFB", "NELFCD", "NELFE", "NETO2", "NEUROD1", "NEUROD2", "NFAT5", "NFATC1",
      "NFATC2", "NFATC3", "NFATC4", "NFE2", "NFE2L1", "NFE2L2", "NFE2L3", "NFIA", "NFIB", "NFIC", "NFIL3", "NFIX",
      "NFKB1", "NFKB2", "NFKBIA", "NFKBIZ", "NFRKB", "NFXL1", "NFYA", "NFYB", "NFYC", "NHLH1", "NIPBL", "NKRF",
      "NKX2-1", "NKX2-2", "NKX2-3", "NKX2-5", "NKX2-8", "NKX3-1", "NKX6-1", "NLRP2", "NME2", "NONO", "NOS2",
      "NOTCH1", "NOTCH3", "NPAS2", "NPAT", "NPM1", "NR0B1", "NR0B2", "NR1D1", "NR1H2", "NR1H3", "NR1H4",
      "NR1I2", "NR2C1", "NR2C2", "NR2E1", "NR2E3", "NR2F1", "NR2F2", "NR2F6", "NR3C1", "NR3C2", "NR4A1",
      "NR4A2", "NR4A3", "NR5A1", "NR5A2", "NRF1", "NRIP1", "NRL", "NRP2", "NUDT21", "NUFIP1", "NUP153",
      "NUP214", "NUP93", "NUP98", "NUTM1", "NXF1", "OGT", "OLA1", "OLIG2", "ONECUT1", "ONECUT2", "ORC1",
      "ORC2", "OSR1", "OSR2", "OTX1", "OTX2", "OVOL1", "OVOL2", "OVOL3", "PABPC1", "PADI2", "PAF1", "PALB2",
      "PAPOLG", "PARK7", "PARP1", "PATZ1", "PAWR", "PAX2", "PAX5", "PAX6", "PAX7", "PAX8", "PAX9", "PAXIP1",
      "PBRM1", "PBX1", "PBX2", "PBX3", "PBX4", "PBXIP1", "PCBP1", "PCBP2", "PCF11", "PCGF1", "PCGF2", "PCGF3",
      "PCGF5", "PCGF6", "PCNA", "PDHA1", "PDS5A", "PDS5B", "PDX1", "PGBD5", "PGR", "PHB2", "PHC1", "PHF10",
      "PHF12", "PHF13", "PHF19", "PHF2", "PHF20", "PHF20L1", "PHF21A", "PHF5A", "PHF6", "PHF8", "PHIP",
      "PHOX2B", "PHTF2", "PIAS1", "PIAS4", "PITX1", "PITX3", "PKM", "PKNOX1", "PLAG1", "PLRG1", "PLSCR1", "PML",
      "POGZ", "POLA1", "POLH", "POLQ", "POLR2C", "POLR3C", "POU2AF1", "POU2F1", "POU2F2", "POU2F3",
      "POU3F2", "POU5F1", "POU6F1", "PPARA", "PPARD", "PPARG", "PPARGC1A", "PPP1R10", "PPP2CA", "PPP2R1A",
      "PRAME", "PRDM1", "PRDM10", "PRDM11", "PRDM12", "PRDM14", "PRDM15", "PRDM2", "PRDM4", "PRDM6",
      "PREB", "PRKCQ", "PRKDC", "PRMT1", "PRMT3", "PRMT5", "PRMT6", "PROSER1", "PROX1", "PRPF4", "PRRX1",
      "PRRX2", "PRTN3", "PSIP1", "PSMB5", "PSMD1", "PTBP1", "PTBP3", "PTEN", "PTF1A", "PTTG1", "PURB", "PYGO2",
      "QKI", "QRICH1", "QSER1", "RAC1", "RAD21", "RAD21L1", "RAD51", "RAG1", "RAG2", "RARA", "RARB", "RARG",
      "RB1", "RBAK", "RBBP4", "RBBP5", "RBBP7", "RBBP8", "RBFOX2", "RBL1", "RBL2", "RBM14", "RBM15", "RBM17",
      "RBM22", "RBM25", "RBM26", "RBM34", "RBM39", "RBMX", "RBPJ", "RCOR1", "RCOR2", "RCOR3", "REC8", "REL",
      "RELA", "RELB", "REPIN1", "RERE", "REST", "RFX1", "RFX2", "RFX3", "RFX5", "RFX7", "RFXANK", "RFXAP",
      "RHOXF2B", "RING1", "RIOK2", "RLF", "RNASEH1", "RNASEH2A", "RNF2", "RNF20", "RORA", "RORB", "RORC",
      "RPA1", "RPA2", "RREB1", "RTF1", "RUNX1", "RUNX1T1", "RUNX2", "RUNX3", "RUVBL1", "RUVBL2", "RXRA",
      "RXRB", "RYBP", "S100A10", "S100A9", "SAFB", "SAFB2", "SAG", "SALL1", "SALL2", "SALL3", "SALL4", "SAMD1",
      "SAP130", "SAP30", "SATB1", "SATB2", "SCMH1", "SCML2", "SCML4", "SCRT1", "SCRT2", "SERPINE1", "SET",
      "SETD1A", "SETD2", "SETD7", "SETDB1", "SETX", "SF3B1", "SFMBT1", "SFPQ", "SHOX2", "SIN3A", "SIN3B", "SIPA1",
      "SIRT1", "SIRT3", "SIRT6", "SIRT7", "SIX1", "SIX2", "SIX4", "SIX5", "SKI", "SKIL", "SKP2", "SLBP", "SLC30A9", "SLFN5",
      "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD9", "SMARCA1", "SMARCA2", "SMARCA4",
      "SMARCA5", "SMARCAL1", "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1", "SMARCE1", "SMC1A", "SMC1B",
      "SMC2", "SMC3", "SMC4", "SMCHD1", "SMN1", "SMYD3", "SNAI1", "SNAI2", "SNAPC1", "SNAPC2", "SNAPC4",
      "SNAPC5", "SNIP1", "SNRNP70", "SOD1", "SON", "SOX1", "SOX10", "SOX11", "SOX13", "SOX15", "SOX17",
      "SOX18", "SOX2", "SOX21", "SOX3", "SOX4", "SOX5", "SOX6", "SOX8", "SOX9", "SP1", "SP110", "SP140", "SP140L",
      "SP2", "SP3", "SP4", "SP5", "SP7", "SPDEF", "SPEN", "SPI1", "SPIN1", "SRCAP", "SREBF1", "SREBF2", "SRF",
      "SRPK1", "SRPK2", "SRSF1", "SRSF3", "SRSF4", "SRSF7", "SRSF9", "SS18", "SS18L1", "SSRP1", "STAG1", "STAG2",
      "STAG3", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6", "SUMO1", "SUMO2", "SUPT16H",
      "SUPT20H", "SUPT5H", "SUPT6H", "SUV39H1", "SUZ12", "SVIL", "SYNCRIP", "TAF1", "TAF15", "TAF1C", "TAF2",
      "TAF3", "TAF5L", "TAF7", "TAF9B", "TAL1", "TARDBP", "TBL1XR1", "TBL1Y", "TBP", "TBPL1", "TBX18", "TBX2",
      "TBX21", "TBX3", "TBX5", "TCEA2", "TCEAL8", "TCF12", "TCF15", "TCF19", "TCF21", "TCF25", "TCF3", "TCF4",
      "TCF7", "TCF7L1", "TCF7L2", "TCFL5", "TDRD3", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "TEF", "TERF1", "TERF2",
      "TERT", "TET1", "TET2", "TET3", "TEX10", "TFAM", "TFAP2A", "TFAP2B", "TFAP2C", "TFAP4", "TFCP2", "TFCP2L1",
      "TFDP1", "TFDP2", "TFE3", "TFEB", "TGIF2", "TGIF2LX", "THAP1", "THAP11", "THAP12", "THAP4", "THAP7",
      "THAP8", "THAP9", "THOC5", "THRA", "THRAP3", "THRB", "TIGD3", "TIGD6", "TLE3", "TLX1", "TLX3", "TOE1",
      "TOP1", "TOP2A", "TOP2B", "TOX", "TOX2", "TOX4", "TP53", "TP53BP1", "TP63", "TP73", "TRAFD1", "TRERF1",
      "TRIB3", "TRIM22", "TRIM24", "TRIM25", "TRIM28", "TRIM33", "TRIM37", "TRIP13", "TRMT112", "TRMT13",
      "TRPS1", "TRRAP", "TSC1", "TSC22D2", "TSC22D4", "TSHZ1", "TSHZ2", "TTF2", "TUBG1", "TWIST1", "U2AF1",
      "U2AF2", "UBE2A", "UBE2I", "UBE3A", "UBP1", "UBR5", "UBTF", "UPF1", "USF1", "USF2", "USF3", "USP22",
      "USP7", "VCP", "VDR", "VEZF1", "VGLL1", "VPS72", "VRK1", "VSX2", "WAPL", "WAS", "WDHD1", "WDR33",
      "WDR5", "WDR77", "WDR82", "WIZ", "WRN", "WRNIP1", "WT1", "WTAP", "WWTR1", "XBP1", "XPO1", "XRCC1",
      "XRCC3", "XRCC4", "XRCC5", "XRCC6", "XRN2", "YAP1", "YBX1", "YBX3", "YEATS2", "YEATS4", "YTHDC1", "YWHAZ",
      "YY1", "YY2", "ZBED1", "ZBED2", "ZBED4", "ZBED5", "ZBTB1", "ZBTB10", "ZBTB11", "ZBTB12", "ZBTB14", "ZBTB16",
      "ZBTB17", "ZBTB2", "ZBTB20", "ZBTB21", "ZBTB24", "ZBTB25", "ZBTB26", "ZBTB3", "ZBTB33", "ZBTB34", "ZBTB37",
      "ZBTB38", "ZBTB39", "ZBTB4", "ZBTB40", "ZBTB42", "ZBTB43", "ZBTB44", "ZBTB46", "ZBTB48", "ZBTB49", "ZBTB5",
      "ZBTB6", "ZBTB7A", "ZBTB7B", "ZBTB8A", "ZBTB9", "ZC3H10", "ZC3H11A", "ZC3H4", "ZC3H8", "ZEB1", "ZEB2",
      "ZFC3H1", "ZFHX2", "ZFHX4", "ZFP1", "ZFP14", "ZFP2", "ZFP3", "ZFP30", "ZFP36", "ZFP36L2", "ZFP37", "ZFP41",
      "ZFP42", "ZFP62", "ZFP64", "ZFP69B", "ZFP82", "ZFP91", "ZFX", "ZFY", "ZGPAT", "ZHX1", "ZHX2", "ZHX3", "ZIC2",
      "ZIC5", "ZIK1", "ZKSCAN1", "ZKSCAN3", "ZKSCAN5", "ZKSCAN8", "ZMAT3", "ZMAT4", "ZMAT5", "ZMIZ1", "ZMYM2",
      "ZMYM3", "ZMYND11", "ZMYND8", "ZNF10", "ZNF101", "ZNF112", "ZNF114", "ZNF117", "ZNF12", "ZNF121",
      "ZNF124", "ZNF132", "ZNF133", "ZNF134", "ZNF138", "ZNF140", "ZNF142", "ZNF143", "ZNF146", "ZNF148",
      "ZNF155", "ZNF157", "ZNF16", "ZNF165", "ZNF169", "ZNF17", "ZNF174", "ZNF175", "ZNF18", "ZNF180", "ZNF184",
      "ZNF189", "ZNF19", "ZNF195", "ZNF197", "ZNF2", "ZNF20", "ZNF202", "ZNF205", "ZNF207", "ZNF211", "ZNF212",
      "ZNF213", "ZNF214", "ZNF215", "ZNF217", "ZNF219", "ZNF221", "ZNF223", "ZNF224", "ZNF225", "ZNF227",
      "ZNF23", "ZNF230", "ZNF232", "ZNF234", "ZNF239", "ZNF24", "ZNF248", "ZNF25", "ZNF250", "ZNF251", "ZNF253",
      "ZNF256", "ZNF257", "ZNF26", "ZNF260", "ZNF263", "ZNF264", "ZNF266", "ZNF274", "ZNF276", "ZNF277",
      "ZNF280A", "ZNF280B", "ZNF280C", "ZNF280D", "ZNF281", "ZNF282", "ZNF292", "ZNF296", "ZNF3", "ZNF30",
      "ZNF300", "ZNF302", "ZNF311", "ZNF316", "ZNF317", "ZNF318", "ZNF319", "ZNF324", "ZNF326", "ZNF329",
      "ZNF331", "ZNF333", "ZNF334", "ZNF335", "ZNF337", "ZNF33A", "ZNF33B", "ZNF34", "ZNF341", "ZNF343",
      "ZNF346", "ZNF347", "ZNF350", "ZNF354B", "ZNF354C", "ZNF362", "ZNF366", "ZNF367", "ZNF37A", "ZNF382",
      "ZNF383", "ZNF384", "ZNF391", "ZNF394", "ZNF395", "ZNF396", "ZNF397", "ZNF398", "ZNF404", "ZNF407",
      "ZNF408", "ZNF41", "ZNF410", "ZNF414", "ZNF416", "ZNF417", "ZNF423", "ZNF425", "ZNF426", "ZNF430",
      "ZNF431", "ZNF432", "ZNF433", "ZNF436", "ZNF438", "ZNF44", "ZNF441", "ZNF444", "ZNF445", "ZNF446",
      "ZNF449", "ZNF451", "ZNF460", "ZNF462", "ZNF473", "ZNF48", "ZNF483", "ZNF484", "ZNF485", "ZNF486",
      "ZNF488", "ZNF490", "ZNF491", "ZNF493", "ZNF501", "ZNF503", "ZNF507", "ZNF510", "ZNF511", "ZNF512",
      "ZNF512B", "ZNF513", "ZNF514", "ZNF516", "ZNF518A", "ZNF518B", "ZNF521", "ZNF524", "ZNF526", "ZNF529",
      "ZNF530", "ZNF532", "ZNF543", "ZNF544", "ZNF547", "ZNF548", "ZNF549", "ZNF550", "ZNF551", "ZNF552",
      "ZNF555", "ZNF556", "ZNF557", "ZNF558", "ZNF560", "ZNF561", "ZNF565", "ZNF569", "ZNF57", "ZNF570",
      "ZNF571", "ZNF572", "ZNF574", "ZNF576", "ZNF577", "ZNF579", "ZNF580", "ZNF581", "ZNF583", "ZNF584",
      "ZNF585B", "ZNF586", "ZNF587", "ZNF589", "ZNF592", "ZNF593", "ZNF596", "ZNF597", "ZNF598", "ZNF600",
      "ZNF605", "ZNF607", "ZNF608", "ZNF609", "ZNF610", "ZNF611", "ZNF614", "ZNF615", "ZNF616", "ZNF619",
      "ZNF621", "ZNF622", "ZNF623", "ZNF624", "ZNF626", "ZNF629", "ZNF639", "ZNF641", "ZNF644", "ZNF646",
      "ZNF652", "ZNF654", "ZNF655", "ZNF658", "ZNF660", "ZNF662", "ZNF664", "ZNF668", "ZNF670", "ZNF671",
      "ZNF674", "ZNF677", "ZNF678", "ZNF680", "ZNF687", "ZNF691", "ZNF692", "ZNF695", "ZNF696", "ZNF697",
      "ZNF699", "ZNF7", "ZNF70", "ZNF700", "ZNF701", "ZNF703", "ZNF704", "ZNF707", "ZNF709", "ZNF71", "ZNF710",
      "ZNF711", "ZNF713", "ZNF717", "ZNF737", "ZNF740", "ZNF747", "ZNF750", "ZNF75A", "ZNF76", "ZNF761",
      "ZNF764", "ZNF766", "ZNF768", "ZNF77", "ZNF770", "ZNF772", "ZNF773", "ZNF775", "ZNF776", "ZNF777",
      "ZNF778", "ZNF780A", "ZNF781", "ZNF782", "ZNF784", "ZNF785", "ZNF786", "ZNF79", "ZNF791", "ZNF792",
      "ZNF8", "ZNF800", "ZNF816", "ZNF83", "ZNF830", "ZNF837", "ZNF839", "ZNF84", "ZNF843", "ZNF850", "ZNF865",
      "ZNF883", "ZNF891", "ZNF92", "ZNHIT1", "ZSCAN12", "ZSCAN16", "ZSCAN18", "ZSCAN20", "ZSCAN21",
      "ZSCAN22", "ZSCAN23", "ZSCAN26", "ZSCAN29", "ZSCAN30", "ZSCAN31", "ZSCAN32", "ZSCAN4", "ZSCAN5A",
      "ZSCAN5B", "ZSCAN5C", "ZSCAN9", "ZXDB", "ZXDC", "ZZZ3"
    )
  } else if (genome_assembly == "hg19") {
    TFs_in_ChIP_Atlas = c(
      "AATF", "ABI1", "ACSS2", "ACTR5", "ADAR", "ADNP", "ADNP2", "AEBP2", "AFF1", "AFF4", "AGO1", "AGO2",
      "AHCTF1", "AHDC1", "AHR", "AHRR", "AICDA", "AIM2", "AKAP8", "AKNA", "ALKBH1", "ALKBH3", "ALOX5",
      "ALYREF", "ANAPC1", "APC", "APOBEC3B", "AR", "ARHGAP35", "ARID1A", "ARID1B", "ARID2", "ARID3A",
      "ARID3B", "ARID4A", "ARID4B", "ARID5B", "ARNT", "ARNT2", "ARNTL", "ARNTL2", "ARRB1", "ASCL1", "ASCL2",
      "ASF1A", "ASH1L", "ASH2L", "ASPSCR1", "ASXL1", "ASXL2", "ASXL3", "ATAD3A", "ATF1", "ATF2", "ATF3", "ATF4",
      "ATF5", "ATF6", "ATF7", "ATF7IP", "ATM", "ATOH1", "ATOH8", "ATR", "ATRX", "ATXN7L3", "AUTS2", "BACH1",
      "BACH2", "BAHCC1", "BAHD1", "BANF1", "BANP", "BAP1", "BARX1", "BATF", "BATF3", "BBX", "BCL10", "BCL11A",
      "BCL11B", "BCL3", "BCL6", "BCL6B", "BCL9", "BCLAF1", "BCOR", "BCORL1", "BDP1", "BHLHE40", "BHLHE41",
      "BMI1", "BNC2", "BOD1L1", "BPTF", "BRCA1", "BRCA2", "BRD1", "BRD2", "BRD3", "BRD4", "BRD7", "BRD8",
      "BRD9", "BRDT", "BRF1", "BRF2", "BRPF1", "BRPF3", "BRWD1", "BUB3", "CAMTA1", "CAMTA2", "CARM1",
      "CASP3", "CASP8AP2", "CASZ1", "CBFA2T2", "CBFA2T3", "CBFB", "CBX1", "CBX2", "CBX3", "CBX4", "CBX5", "CBX6",
      "CBX7", "CBX8", "CC2D1A", "CCAR2", "CCND2", "CCNE1", "CCNT1", "CCNT2", "CD74", "CDAN1", "CDC16",
      "CDC45", "CDC5L", "CDC73", "CDK11A", "CDK12", "CDK2", "CDK6", "CDK7", "CDK8", "CDK9", "CDKN1B", "CDT1",
      "CDX1", "CDX2", "CDX4", "CDYL", "CDYL2", "CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG", "CEBPZ", "CERS6",
      "CGGBP1", "CHAF1A", "CHAMP1", "CHAT", "CHCHD3", "CHD1", "CHD2", "CHD3", "CHD4", "CHD6", "CHD7",
      "CHD8", "CHUK", "CIC", "CIITA", "CIZ1", "CLOCK", "CNBP", "CNOT3", "COBLL1", "COIL", "COPS2", "CPSF1", "CPSF2",
      "CPSF3", "CPSF4", "CPSF6", "CPSF7", "CREB1", "CREB3", "CREB3L1", "CREB3L2", "CREB5", "CREBBP", "CREBL2",
      "CREM", "CRTC1", "CRTC2", "CRX", "CRY1", "CSDC2", "CSDE1", "CSE1L", "CSNK2A1", "CSRNP3", "CSTF2", "CSTF3",
      "CTBP1", "CTBP2", "CTCF", "CTCFL", "CTNNB1", "CTR9", "CTSG", "CUL3", "CUL4A", "CUL4B", "CUX1", "CXXC1",
      "CXXC5", "DACH1", "DAXX", "DBP", "DCP1A", "DCP2", "DDIT3", "DDX11", "DDX17", "DDX20", "DDX21", "DDX39B",
      "DDX5", "DDX55", "DEAF1", "DEK", "DICER1", "DIDO1", "DKC1", "DLX4", "DLX5", "DLX6", "DMAP1", "DMC1",
      "DMRT2", "DMTF1", "DNA2", "DNAJC2", "DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "DONSON", "DOT1L",
      "DPF1", "DPF2", "DPPA3", "DR1", "DRAP1", "DROSHA", "DTL", "DUX4", "DVL1", "DYRK1A", "DZIP1", "E2F1",
      "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8", "E4F1", "EBF1", "EBF3", "EEA1", "EED", "EGFR", "EGR1",
      "EGR2", "EHF", "EHMT1", "EHMT2", "ELANE", "ELAVL1", "ELF1", "ELF2", "ELF3", "ELF4", "ELK1", "ELK3", "ELK4",
      "ELL", "ELL2", "EMSY", "EMX1", "EN1", "EOMES", "EP300", "EP400", "EPAS1", "EPCAM", "EPOP", "ERBB2",
      "ERCC2", "ERCC3", "ERCC6", "ERCC8", "ERF", "ERG", "ERH", "ESR1", "ESR2", "ESRRA", "ESRRB", "ESRRG", "ETS1",
      "ETS2", "ETV1", "ETV2", "ETV3", "ETV4", "ETV5", "ETV6", "EVX1", "EWSR1", "EZH1", "EZH2", "F2RL1", "FANCD2",
      "FANCM", "FBL", "FBP1", "FBXO44", "FERD3L", "FEZF1", "FGFR1", "FIP1L1", "FLI1", "FMR1", "FOS", "FOSB",
      "FOSL1", "FOSL2", "FOXA1", "FOXA2", "FOXA3", "FOXB1", "FOXC1", "FOXC2", "FOXD3", "FOXF1", "FOXF2",
      "FOXH1", "FOXJ2", "FOXJ3", "FOXK1", "FOXK2", "FOXL2", "FOXM1", "FOXN3", "FOXO1", "FOXO3", "FOXO4",
      "FOXP1", "FOXP2", "FOXP3", "FOXP4", "FOXQ1", "FOXR2", "FOXS1", "FUBP1", "FUBP3", "FUS", "FXR1", "FXR2",
      "G3BP1", "GABPA", "GABPB1", "GABPB2", "GATA1", "GATA2", "GATA3", "GATA4", "GATA6", "GATAD1", "GATAD2A",
      "GATAD2B", "GCM1", "GFI1", "GFI1B", "GGCT", "GLI1", "GLI2", "GLI4", "GLIS1", "GLIS2", "GLIS3", "GLMP", "GLYR1",
      "GMEB1", "GMEB2", "GNL3", "GPBP1L1", "GPN1", "GPS2", "GREB1", "GRHL1", "GRHL2", "GRHL3", "GRWD1",
      "GTF2A2", "GTF2B", "GTF2E2", "GTF2F1", "GTF2F2", "GTF2H2", "GTF2H2C", "GTF2I", "GTF3A", "GTF3C1",
      "GTF3C2", "GTF3C5", "GZF1", "HAND1", "HAND2", "HAT1", "HBP1", "HBZ", "HCFC1", "HDAC1", "HDAC2",
      "HDAC3", "HDAC4", "HDAC6", "HDAC8", "HDAC9", "HDGF", "HECTD1", "HES1", "HES2", "HES4", "HES5", "HES7",
      "HESX1", "HEXIM1", "HEY1", "HEY2", "HEYL", "HHEX", "HIC1", "HIC2", "HIF1A", "HIF3A", "HINFP", "HIRA",
      "HIVEP1", "HIVEP3", "HJURP", "HK2", "HLF", "HLTF", "HMBOX1", "HMG20A", "HMG20B", "HMGA1", "HMGA2",
      "HMGB1", "HMGB2", "HMGN1", "HMGN2", "HMGN3", "HMGN4", "HMGXB3", "HMGXB4", "HNF1A", "HNF1B",
      "HNF4A", "HNF4G", "HNRNPA2B1", "HNRNPC", "HNRNPH1", "HNRNPK", "HNRNPL", "HNRNPLL", "HNRNPUL1",
      "HOMEZ", "HOXA11", "HOXA13", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXB1",
      "HOXB13", "HOXB2", "HOXB4", "HOXB5", "HOXB6", "HOXB7", "HOXB8", "HOXC11", "HOXC13", "HOXC5", "HOXC8",
      "HOXC9", "HOXD1", "HOXD11", "HOXD13", "HOXD4", "HOXD9", "HSF1", "HSF2", "HSF4", "HSP90B1", "HSPA1A",
      "ICE1", "ICE2", "ID1", "ID2", "ID3", "ID4", "IFI16", "IGF1", "IKZF1", "IKZF2", "IKZF3", "IKZF4", "IKZF5", "ILF3", "ILK",
      "ING2", "ING4", "ING5", "INO80", "INSM2", "INSR", "INTS10", "INTS12", "INTS13", "INTS3", "INTS5", "INTS6",
      "IRF1", "IRF2", "IRF2BP2", "IRF3", "IRF4", "IRF5", "IRF8", "IRF9", "IRX2", "IRX5", "ISL1", "ISL2", "ISX", "JADE3",
      "JARID2", "JAZF1", "JDP2", "JMJD1C", "JMJD6", "JRK", "JUN", "JUNB", "JUND", "KANSL1", "KANSL3", "KAT2A",
      "KAT2B", "KAT5", "KAT6A", "KAT6B", "KAT7", "KAT8", "KCMF1", "KDM1A", "KDM1B", "KDM2A", "KDM2B",
      "KDM3A", "KDM3B", "KDM4A", "KDM4B", "KDM4C", "KDM5A", "KDM5B", "KDM5C", "KDM5D", "KDM6A",
      "KDM6B", "KHSRP", "KLF1", "KLF10", "KLF11", "KLF12", "KLF13", "KLF14", "KLF15", "KLF16", "KLF17", "KLF3",
      "KLF4", "KLF5", "KLF6", "KLF7", "KLF8", "KLF9", "KMT2A", "KMT2B", "KMT2C", "KMT2D", "L3MBTL2", "L3MBTL3",
      "LARP7", "LBX2", "LCOR", "LCORL", "LDB1", "LDB2", "LEF1", "LEO1", "LHX2", "LHX4", "LHX5", "LHX6", "LIG4",
      "LIN28B", "LIN54", "LIN9", "LMNA", "LMNB1", "LMO1", "LMO2", "LMO3", "LMTK3", "LRIF1", "LRRFIP1", "LRWD1",
      "LYL1", "MAF", "MAF1", "MAFB", "MAFF", "MAFG", "MAFK", "MAML1", "MAML3", "MAP1LC3B", "MAP2K1",
      "MAPK1", "MAPK14", "MATR3", "MAU2", "MAX", "MAZ", "MBD1", "MBD2", "MBD3", "MBD3L2", "MBD4",
      "MBD5", "MBD6", "MBNL3", "MCM2", "MCM3", "MCM5", "MCM7", "MCRS1", "MDM2", "MDN1", "MEAF6",
      "MECOM", "MECP2", "MED1", "MED12", "MED13", "MED16", "MED21", "MED22", "MED25", "MED26", "MED30",
      "MED31", "MED4", "MED6", "MED7", "MED8", "MED9", "MEF2A", "MEF2B", "MEF2C", "MEF2D", "MEIS1",
      "MEIS2", "MEN1", "MEOX2", "MEPCE", "MESP1", "METTL14", "METTL3", "MGA", "MGMT", "MIER1", "MIER2",
      "MIER3", "MIF", "MITF", "MIXL1", "MKI67", "MLLT1", "MLLT10", "MLLT3", "MLLT6", "MLX", "MLXIP", "MNT",
      "MNX1", "MORC2", "MPHOSPH8", "MSC", "MSH6", "MSX2", "MTA1", "MTA2", "MTA3", "MTF1", "MTF2",
      "MTHFD1", "MTOR", "MXD1", "MXD3", "MXD4", "MXI1", "MYB", "MYBBP1A", "MYBL1", "MYBL2", "MYC",
      "MYCL", "MYCN", "MYF5", "MYF6", "MYH11", "MYNN", "MYOCD", "MYOD1", "MYOG", "MYRF", "MZF1",
      "N6AMT1", "NAB2", "NABP2", "NAIF1", "NANOG", "NAT10", "NBN", "NCAPG", "NCAPG2", "NCAPH2", "NCBP1",
      "NCL", "NCOA1", "NCOA2", "NCOA3", "NCOA4", "NCOA5", "NCOA6", "NCOR1", "NCOR2", "NEDD4", "NEIL1",
      "NELFA", "NELFB", "NELFCD", "NELFE", "NETO2", "NEUROD1", "NEUROD2", "NFAT5", "NFATC1", "NFATC2",
      "NFATC3", "NFATC4", "NFE2", "NFE2L1", "NFE2L2", "NFE2L3", "NFIA", "NFIB", "NFIC", "NFIL3", "NFIX", "NFKB1",
      "NFKB2", "NFKBIA", "NFKBIZ", "NFRKB", "NFXL1", "NFYA", "NFYB", "NFYC", "NHLH1", "NIPBL", "NKRF", "NKX2-1",
      "NKX2-2", "NKX2-3", "NKX2-5", "NKX2-8", "NKX3-1", "NKX6-1", "NLRP2", "NME2", "NONO", "NOS2", "NOTCH1",
      "NOTCH3", "NPAS2", "NPAT", "NPM1", "NR0B1", "NR0B2", "NR1D1", "NR1H2", "NR1H3", "NR1H4", "NR1I2",
      "NR2C1", "NR2C2", "NR2E3", "NR2F1", "NR2F2", "NR2F6", "NR3C1", "NR3C2", "NR4A1", "NR4A2", "NR4A3",
      "NR5A1", "NR5A2", "NRF1", "NRIP1", "NRL", "NRP2", "NUDT21", "NUFIP1", "NUP153", "NUP214", "NUP93",
      "NUP98", "NUTM1", "NXF1", "OGT", "OLA1", "OLIG2", "ONECUT1", "ONECUT2", "ORC1", "ORC2", "OSR1", "OSR2",
      "OTX1", "OTX2", "OVOL1", "OVOL2", "OVOL3", "PABPC1", "PADI2", "PAF1", "PALB2", "PAPOLG", "PARK7", "PARP1",
      "PATZ1", "PAWR", "PAX2", "PAX5", "PAX6", "PAX7", "PAX8", "PAX9", "PAXIP1", "PBRM1", "PBX1", "PBX2", "PBX3",
      "PBX4", "PBXIP1", "PCBP1", "PCBP2", "PCF11", "PCGF1", "PCGF2", "PCGF3", "PCGF5", "PCGF6", "PCNA", "PDHA1",
      "PDS5A", "PDS5B", "PDX1", "PGBD5", "PGR", "PHB2", "PHC1", "PHF10", "PHF12", "PHF13", "PHF19", "PHF2",
      "PHF20", "PHF20L1", "PHF21A", "PHF5A", "PHF6", "PHF8", "PHIP", "PHOX2B", "PHTF2", "PIAS1", "PIAS4", "PITX1",
      "PITX3", "PKM", "PKNOX1", "PLAG1", "PLRG1", "PLSCR1", "PML", "POGZ", "POLA1", "POLH", "POLQ", "POLR2C",
      "POLR3C", "POU2AF1", "POU2F1", "POU2F2", "POU2F3", "POU3F2", "POU5F1", "POU6F1", "PPARA", "PPARD",
      "PPARG", "PPARGC1A", "PPP1R10", "PPP2CA", "PPP2R1A", "PRAME", "PRDM1", "PRDM10", "PRDM11",
      "PRDM12", "PRDM14", "PRDM15", "PRDM2", "PRDM4", "PRDM6", "PREB", "PRKCQ", "PRKDC", "PRMT1",
      "PRMT3", "PRMT5", "PRMT6", "PROSER1", "PROX1", "PRPF4", "PRRX1", "PRRX2", "PRTN3", "PSIP1", "PSMB5",
      "PSMD1", "PTBP1", "PTBP3", "PTEN", "PTF1A", "PTTG1", "PURB", "PYGO2", "QKI", "QRICH1", "QSER1", "RAC1",
      "RAD21", "RAD21L1", "RAD51", "RAG1", "RAG2", "RARA", "RARB", "RARG", "RB1", "RBAK", "RBBP4", "RBBP5",
      "RBBP7", "RBBP8", "RBFOX2", "RBL1", "RBL2", "RBM14", "RBM15", "RBM17", "RBM22", "RBM25", "RBM26",
      "RBM34", "RBM39", "RBMX", "RBPJ", "RCOR1", "RCOR2", "RCOR3", "REC8", "REL", "RELA", "RELB", "REPIN1",
      "RERE", "REST", "RFX1", "RFX2", "RFX3", "RFX5", "RFX7", "RFXANK", "RFXAP", "RHOXF2B", "RING1", "RIOK2",
      "RLF", "RNASEH1", "RNASEH2A", "RNF2", "RNF20", "RORA", "RORB", "RORC", "RPA1", "RPA2", "RREB1", "RTF1",
      "RUNX1", "RUNX1T1", "RUNX2", "RUNX3", "RUVBL1", "RUVBL2", "RXRA", "RXRB", "RYBP", "S100A10", "S100A9",
      "SAFB", "SAFB2", "SAG", "SALL1", "SALL2", "SALL3", "SALL4", "SAMD1", "SAP130", "SAP30", "SATB1", "SATB2",
      "SCMH1", "SCML2", "SCML4", "SCRT1", "SCRT2", "SERPINE1", "SET", "SETD1A", "SETD2", "SETD7", "SETDB1",
      "SETX", "SF3B1", "SFMBT1", "SFPQ", "SHOX2", "SIN3A", "SIN3B", "SIPA1", "SIRT1", "SIRT3", "SIRT6", "SIRT7",
      "SIX1", "SIX2", "SIX4", "SIX5", "SKI", "SKIL", "SKP2", "SLBP", "SLC30A9", "SLFN5", "SMAD1", "SMAD2", "SMAD3",
      "SMAD4", "SMAD5", "SMAD9", "SMARCA1", "SMARCA2", "SMARCA4", "SMARCA5", "SMARCAL1", "SMARCB1",
      "SMARCC1", "SMARCC2", "SMARCD1", "SMARCE1", "SMC1A", "SMC1B", "SMC2", "SMC3", "SMC4", "SMCHD1",
      "SMN1", "SMYD3", "SNAI1", "SNAI2", "SNAPC1", "SNAPC2", "SNAPC4", "SNAPC5", "SNIP1", "SNRNP70", "SOD1",
      "SON", "SOX1", "SOX10", "SOX11", "SOX13", "SOX15", "SOX18", "SOX2", "SOX21", "SOX3", "SOX4", "SOX5",
      "SOX6", "SOX8", "SOX9", "SP1", "SP110", "SP140", "SP140L", "SP2", "SP3", "SP4", "SP5", "SP7", "SPDEF", "SPEN",
      "SPI1", "SPIN1", "SRCAP", "SREBF1", "SREBF2", "SRF", "SRPK1", "SRPK2", "SRSF1", "SRSF3", "SRSF4", "SRSF7",
      "SRSF9", "SS18", "SS18L1", "SSRP1", "STAG1", "STAG2", "STAG3", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A",
      "STAT5B", "STAT6", "SUMO1", "SUMO2", "SUPT16H", "SUPT20H", "SUPT5H", "SUPT6H", "SUV39H1", "SUZ12",
      "SVIL", "SYNCRIP", "TAF1", "TAF15", "TAF1C", "TAF2", "TAF3", "TAF5L", "TAF7", "TAF9B", "TAL1", "TARDBP",
      "TBL1XR1", "TBL1Y", "TBP", "TBPL1", "TBX18", "TBX2", "TBX21", "TBX3", "TBX5", "TCEA2", "TCEAL8", "TCF12",
      "TCF15", "TCF19", "TCF21", "TCF25", "TCF3", "TCF4", "TCF7", "TCF7L1", "TCF7L2", "TCFL5", "TDRD3", "TEAD1",
      "TEAD2", "TEAD3", "TEAD4", "TEF", "TERF1", "TERF2", "TERT", "TET1", "TET2", "TET3", "TEX10", "TFAM",
      "TFAP2A", "TFAP2B", "TFAP2C", "TFAP4", "TFCP2", "TFCP2L1", "TFDP1", "TFDP2", "TFE3", "TFEB", "TGIF2",
      "TGIF2LX", "THAP1", "THAP11", "THAP12", "THAP4", "THAP7", "THAP8", "THAP9", "THOC5", "THRA", "THRAP3",
      "THRB", "TIGD3", "TIGD6", "TLE3", "TLX1", "TLX3", "TOE1", "TOP1", "TOP2A", "TOP2B", "TOX", "TOX2", "TOX4",
      "TP53", "TP53BP1", "TP63", "TP73", "TRAFD1", "TRERF1", "TRIB3", "TRIM22", "TRIM24", "TRIM25", "TRIM28",
      "TRIM33", "TRIM37", "TRIP13", "TRMT112", "TRMT13", "TRPS1", "TRRAP", "TSC1", "TSC22D2", "TSC22D4",
      "TSHZ1", "TSHZ2", "TTF2", "TUBG1", "TWIST1", "U2AF1", "U2AF2", "UBE2A", "UBE2I", "UBE3A", "UBP1", "UBR5",
      "UBTF", "UPF1", "USF1", "USF2", "USF3", "USP22", "USP7", "VCP", "VDR", "VEZF1", "VGLL1", "VPS72", "VRK1",
      "VSX2", "WAPL", "WAS", "WDHD1", "WDR33", "WDR5", "WDR77", "WDR82", "WIZ", "WRN", "WRNIP1", "WT1",
      "WTAP", "WWTR1", "XBP1", "XPO1", "XRCC1", "XRCC3", "XRCC4", "XRCC5", "XRCC6", "XRN2", "YAP1", "YBX1",
      "YBX3", "YEATS2", "YEATS4", "YTHDC1", "YWHAZ", "YY1", "YY2", "ZBED1", "ZBED2", "ZBED4", "ZBED5", "ZBTB1",
      "ZBTB10", "ZBTB11", "ZBTB12", "ZBTB14", "ZBTB16", "ZBTB17", "ZBTB2", "ZBTB20", "ZBTB21", "ZBTB24",
      "ZBTB25", "ZBTB26", "ZBTB3", "ZBTB33", "ZBTB34", "ZBTB37", "ZBTB38", "ZBTB39", "ZBTB4", "ZBTB40", "ZBTB42",
      "ZBTB43", "ZBTB44", "ZBTB46", "ZBTB48", "ZBTB49", "ZBTB5", "ZBTB6", "ZBTB7A", "ZBTB7B", "ZBTB8A", "ZBTB9",
      "ZC3H10", "ZC3H11A", "ZC3H4", "ZC3H8", "ZEB1", "ZEB2", "ZFC3H1", "ZFHX2", "ZFHX4", "ZFP1", "ZFP14", "ZFP2",
      "ZFP3", "ZFP30", "ZFP36", "ZFP36L2", "ZFP37", "ZFP41", "ZFP42", "ZFP62", "ZFP64", "ZFP69B", "ZFP82", "ZFP91",
      "ZFX", "ZFY", "ZGPAT", "ZHX1", "ZHX2", "ZHX3", "ZIC2", "ZIC5", "ZIK1", "ZKSCAN1", "ZKSCAN3", "ZKSCAN5",
      "ZKSCAN8", "ZMAT3", "ZMAT4", "ZMAT5", "ZMIZ1", "ZMYM2", "ZMYM3", "ZMYND11", "ZMYND8", "ZNF10",
      "ZNF101", "ZNF112", "ZNF114", "ZNF117", "ZNF12", "ZNF121", "ZNF124", "ZNF132", "ZNF133", "ZNF134",
      "ZNF138", "ZNF140", "ZNF142", "ZNF143", "ZNF146", "ZNF148", "ZNF155", "ZNF157", "ZNF16", "ZNF165",
      "ZNF169", "ZNF17", "ZNF174", "ZNF175", "ZNF18", "ZNF180", "ZNF184", "ZNF189", "ZNF19", "ZNF195", "ZNF197",
      "ZNF2", "ZNF20", "ZNF202", "ZNF205", "ZNF207", "ZNF211", "ZNF212", "ZNF213", "ZNF214", "ZNF215", "ZNF217",
      "ZNF219", "ZNF221", "ZNF223", "ZNF224", "ZNF225", "ZNF227", "ZNF23", "ZNF230", "ZNF232", "ZNF234",
      "ZNF239", "ZNF24", "ZNF248", "ZNF25", "ZNF250", "ZNF251", "ZNF253", "ZNF256", "ZNF257", "ZNF26", "ZNF260",
      "ZNF263", "ZNF264", "ZNF266", "ZNF274", "ZNF276", "ZNF277", "ZNF280A", "ZNF280B", "ZNF280C", "ZNF280D",
      "ZNF281", "ZNF282", "ZNF292", "ZNF296", "ZNF3", "ZNF30", "ZNF300", "ZNF302", "ZNF311", "ZNF316", "ZNF317",
      "ZNF318", "ZNF319", "ZNF324", "ZNF326", "ZNF329", "ZNF331", "ZNF333", "ZNF334", "ZNF335", "ZNF337",
      "ZNF33A", "ZNF33B", "ZNF34", "ZNF341", "ZNF343", "ZNF346", "ZNF347", "ZNF350", "ZNF354B", "ZNF354C",
      "ZNF362", "ZNF366", "ZNF367", "ZNF37A", "ZNF382", "ZNF383", "ZNF384", "ZNF391", "ZNF394", "ZNF395",
      "ZNF396", "ZNF397", "ZNF398", "ZNF404", "ZNF407", "ZNF408", "ZNF41", "ZNF410", "ZNF414", "ZNF416",
      "ZNF417", "ZNF423", "ZNF425", "ZNF426", "ZNF430", "ZNF431", "ZNF432", "ZNF433", "ZNF436", "ZNF438",
      "ZNF44", "ZNF441", "ZNF444", "ZNF445", "ZNF446", "ZNF449", "ZNF451", "ZNF460", "ZNF462", "ZNF473",
      "ZNF48", "ZNF483", "ZNF484", "ZNF485", "ZNF486", "ZNF488", "ZNF490", "ZNF491", "ZNF493", "ZNF501",
      "ZNF503", "ZNF507", "ZNF510", "ZNF511", "ZNF512", "ZNF512B", "ZNF513", "ZNF514", "ZNF516", "ZNF518A",
      "ZNF518B", "ZNF521", "ZNF524", "ZNF526", "ZNF529", "ZNF530", "ZNF532", "ZNF543", "ZNF544", "ZNF547",
      "ZNF548", "ZNF549", "ZNF550", "ZNF551", "ZNF552", "ZNF555", "ZNF556", "ZNF557", "ZNF558", "ZNF560",
      "ZNF561", "ZNF565", "ZNF569", "ZNF57", "ZNF570", "ZNF571", "ZNF572", "ZNF574", "ZNF576", "ZNF577",
      "ZNF579", "ZNF580", "ZNF581", "ZNF583", "ZNF584", "ZNF585B", "ZNF586", "ZNF587", "ZNF589", "ZNF592",
      "ZNF593", "ZNF596", "ZNF597", "ZNF598", "ZNF600", "ZNF605", "ZNF607", "ZNF608", "ZNF609", "ZNF610",
      "ZNF611", "ZNF614", "ZNF615", "ZNF616", "ZNF619", "ZNF621", "ZNF622", "ZNF623", "ZNF624", "ZNF626",
      "ZNF629", "ZNF639", "ZNF641", "ZNF644", "ZNF646", "ZNF652", "ZNF654", "ZNF655", "ZNF658", "ZNF660",
      "ZNF662", "ZNF664", "ZNF668", "ZNF670", "ZNF671", "ZNF674", "ZNF677", "ZNF678", "ZNF680", "ZNF687",
      "ZNF691", "ZNF692", "ZNF695", "ZNF696", "ZNF697", "ZNF699", "ZNF7", "ZNF70", "ZNF700", "ZNF701", "ZNF703",
      "ZNF704", "ZNF707", "ZNF709", "ZNF71", "ZNF710", "ZNF711", "ZNF713", "ZNF717", "ZNF737", "ZNF740",
      "ZNF747", "ZNF750", "ZNF75A", "ZNF76", "ZNF761", "ZNF764", "ZNF766", "ZNF768", "ZNF77", "ZNF770",
      "ZNF772", "ZNF773", "ZNF775", "ZNF776", "ZNF777", "ZNF778", "ZNF780A", "ZNF781", "ZNF782", "ZNF784",
      "ZNF785", "ZNF786", "ZNF79", "ZNF791", "ZNF792", "ZNF8", "ZNF800", "ZNF816", "ZNF83", "ZNF830", "ZNF837",
      "ZNF839", "ZNF84", "ZNF843", "ZNF850", "ZNF865", "ZNF883", "ZNF891", "ZNF92", "ZNHIT1", "ZSCAN12",
      "ZSCAN16", "ZSCAN18", "ZSCAN20", "ZSCAN21", "ZSCAN22", "ZSCAN23", "ZSCAN26", "ZSCAN29", "ZSCAN30",
      "ZSCAN31", "ZSCAN32", "ZSCAN4", "ZSCAN5A", "ZSCAN5B", "ZSCAN5C", "ZSCAN9", "ZXDB", "ZXDC", "ZZZ3"
    )
  } else if (genome_assembly == "mm10") {
    TFs_in_ChIP_Atlas = c(
      "Acaa2", "Acss2", "Actb", "Adnp", "Aebp2", "Aff3", "Aff4", "Ago2", "Ahcy", "Ahr", "Aicda", "Aire",
      "Alkbh5", "Anp32e", "Ar", "Arid1a", "Arid1b", "Arid2", "Arid3a", "Arid4b", "Arid5a", "Arnt",
      "Arntl", "Art1", "Arx", "Ascl1", "Ascl2", "Ash1l", "Ash2l", "Asxl1", "Asxl3", "Atf2", "Atf3", "Atf4",
      "Atf7", "Atf7ip", "Atm", "Atoh1", "Atrx", "Auts2", "Bach1", "Bach2", "Banp", "Bap1", "Batf",
      "Batf3", "Baz1a", "Baz1b", "Baz2a", "Bbx", "Bcl11a", "Bcl11b", "Bcl6", "Bcor", "Bend3", "Bhlhe40",
      "Bhlhe41", "Bmi1", "Bnc1", "Brca1", "Brca2", "Brd2", "Brd3", "Brd4", "Brd7", "Brd9", "Brdt",
      "Brwd1", "Btaf1", "Carm1", "Cbfa2t2", "Cbfa2t3", "Cbfb", "Cbx1", "Cbx2", "Cbx3", "Cbx4", "Cbx5",
      "Cbx7", "Cbx8", "Ccar1", "Ccnk", "Cd38", "Cdc73", "Cdk7", "Cdk8", "Cdk9", "Cdx2", "Cdyl",
      "Cebpa", "Cebpb", "Cebpd", "Cebpe", "Cebpg", "Chaf1b", "Chd1", "Chd2", "Chd4", "Chd7",
      "Chd8", "Chuk", "Cic", "Clock", "Creb1", "Creb3l2", "Crebbp", "Crtc2", "Crtc3", "Crx", "Cry1",
      "Ctbp2", "Ctcf", "Ctcfl", "Ctf1", "Ctnnb1", "Ctr9", "Cux2", "Cxxc1", "Cxxc5", "Daxx", "Dbp", "Ddb1",
      "Ddit3", "Ddx1", "Ddx18", "Ddx21", "Ddx5", "Dgcr8", "Dhx9", "Dicer1", "Dlx1", "Dlx3", "Dmap1",
      "Dmc1", "Dmrt1", "Dmrtb1", "Dnajc2", "Dnmt1", "Dnmt3a", "Dnmt3b", "Dnmt3l", "Dnttip1",
      "Dot1l", "Dpf2", "Dpf3", "Dppa2", "Dppa3", "Dppa4", "Dpy30", "Dr1", "Drosha", "Duxbl1",
      "Dyrk1a", "E2f1", "E2f3", "E2f4", "E2f6", "E4f1", "Ebf1", "Ebf2", "Ebf4", "Eed", "Egr1", "Egr2",
      "Ehf", "Ehmt2", "Elf1", "Elf4", "Elf5", "Elk1", "Elk3", "Elk4", "Ell", "Ell2", "Ell3", "Eloa", "Elob",
      "Emx2", "Eomes", "Ep300", "Ep400", "Epas1", "Epop", "Ercc3", "Erf", "Erg", "Esr1", "Esr2", "Esrra",
      "Esrrb", "Esrrg", "Ets1", "Ets2", "Etv2", "Etv5", "Etv6", "Ezh1", "Ezh2", "Fbxo21", "Fev", "Fgfr1",
      "Fli1", "Fos", "Fosb", "Fosl1", "Fosl2", "Foxa1", "Foxa2", "Foxa3", "Foxc1", "Foxc2", "Foxf1",
      "Foxg1", "Foxh1", "Foxj2", "Foxk1", "Foxk2", "Foxl2", "Foxo1", "Foxo3", "Foxp1", "Foxp2", "Foxp3",
      "Gabpa", "Gad1", "Gadd45a", "Gapdh", "Gata1", "Gata2", "Gata3", "Gata4", "Gata6", "Gckr",
      "Gcm1", "Gfi1", "Gfi1b", "Gli1", "Gli2", "Gli3", "Glis3", "Glyr1", "Gmnn", "Gps2", "Grhl2", "Grhl3",
      "Grip1", "Gsk3b", "Gtf2a2", "Gtf2b", "Gtf2e2", "Gtf2f1", "Gtf3c1", "Hand1", "Hand2", "Hcfc1",
      "Hdac1", "Hdac2", "Hdac3", "Hdac4", "Hdac5", "Hdac6", "Hdac7", "Hells", "Hes1", "Hexim1",
      "Hey1", "Heyl", "Hic2", "Hif1a", "Hira", "Hmga2", "Hmgn2", "Hnf1a", "Hnf1b", "Hnf4a", "Hnf4g",
      "Hnrnpk", "Hnrnpl", "Hnrnpu", "Hnrnpul1", "Hopx", "Hormad1", "Hoxa11", "Hoxa13", "Hoxa5",
      "Hoxa9", "Hoxb1", "Hoxb13", "Hoxb4", "Hoxc10", "Hoxc13", "Hoxc6", "Hoxc8", "Hoxc9",
      "Hoxd13", "Hoxd9", "Hsf1", "Hsf2", "Htatsf1", "Htt", "Id2", "Iho1", "Ikzf1", "Ikzf2", "Ikzf3", "Ing1",
      "Ing5", "Ino80", "Insm1", "Ints3", "Ints5", "Irf1", "Irf2", "Irf2bp2", "Irf3", "Irf4", "Irf5", "Irf7", "Irf8",
      "Irf9", "Isl1", "Jarid2", "Jmjd1c", "Jmjd6", "Jun", "Junb", "Jund", "Kansl3", "Kat2a", "Kat2b", "Kat5",
      "Kat7", "Kat8", "Kdm1a", "Kdm2a", "Kdm2b", "Kdm3a", "Kdm3b", "Kdm4a", "Kdm4b", "Kdm4c",
      "Kdm4d", "Kdm5a", "Kdm5b", "Kdm5c", "Kdm6a", "Kdm6b", "Klf1", "Klf15", "Klf2", "Klf3", "Klf4",
      "Klf5", "Klf6", "Klf9", "Kmt2a", "Kmt2b", "Kmt2c", "Kmt2d", "L3mbtl2", "L3mbtl3", "Ldb1", "Lef1",
      "Leo1", "Lhx2", "Lhx3", "Lhx5", "Lhx6", "Lin9", "Lmna", "Lmnb1", "Lmo2", "Lmx1b", "Lrrfip1",
      "Lyl1", "Macf1", "Maf", "Mafa", "Mafb", "Maff", "Mafg", "Mafk", "Mapk3", "Mapk8", "Matr3",
      "Mau2", "Max", "Maz", "Mbd1", "Mbd2", "Mbd3", "Mcrs1", "Mecp2", "Med1", "Med12",
      "Med14", "Med15", "Med23", "Med24", "Med26", "Med4", "Mef2a", "Mef2c", "Mef2d", "Mei4",
      "Meiosin", "Meis1", "Men1", "Meox2", "Mesp1", "Mettl14", "Mettl3", "Mitf", "Mllt10", "Mllt3",
      "Mlxipl", "Mnt", "Mnx1", "Morc3", "Mphosph8", "Mpnd", "Mre11a", "Msh2", "Msh6", "Msl1",
      "Msl2", "Mta1", "Mta2", "Mtf2", "Mtor", "Mxi1", "Myb", "Mybl1", "Myc", "Mycn", "Myf5",
      "Myo1c", "Myod1", "Myog", "Mysm1", "Myt1l", "N6amt1", "Naa10", "Naca", "Nacc1", "Nanog",
      "Nbn", "Ncapd3", "Ncapg", "Ncaph2", "Ncl", "Ncoa1", "Ncoa2", "Ncoa3", "Ncor1", "Ncor2",
      "Nelfa", "Nelfb", "Nelfcd", "Nelfe", "Neurod1", "Neurod2", "Neurog2", "Nfat5", "Nfatc1",
      "Nfatc2", "Nfe2", "Nfe2l1", "Nfe2l2", "Nfia", "Nfib", "Nfic", "Nfil3", "Nfix", "Nfkb1", "Nfkbia",
      "Nfya", "Nfyb", "Nfyc", "Nipbl", "Nkx2-1", "Nkx2-2", "Nkx2-5", "Nkx3-1", "Nkx3-2", "Nkx6-1",
      "Nlrc5", "Nlrp3", "Nono", "Notch1", "Npas4", "Npm1", "Nr0b1", "Nr0b2", "Nr1d2", "Nr1h2",
      "Nr1h3", "Nr1h4", "Nr1i2", "Nr2f1", "Nr2f2", "Nr3c1", "Nr4a1", "Nr5a1", "Nr5a2", "Nrf1", "Nsl1",
      "Nucks1", "Nup153", "Nup98", "Obox1", "Ogdh", "Ogt", "Olig2", "Onecut1", "Onecut2", "Otx2",
      "Ovol1", "Ovol2", "Paf1", "Pagr1a", "Park7", "Parp1", "Patz1", "Pawr", "Pax3", "Pax5", "Pax6",
      "Pax7", "Pax8", "Paxip1", "Pbrm1", "Pbx1", "Pcgf1", "Pcgf2", "Pcgf3", "Pcgf5", "Pcgf6", "Pcna",
      "Pds5a", "Pds5b", "Pdx1", "Peg3", "Pgr", "Phc1", "Phc2", "Phf10", "Phf19", "Phf2", "Phf20",
      "Phf5a", "Phf6", "Phf7", "Phf8", "Phrf1", "Pitx1", "Pknox1", "Plagl1", "Pml", "Polr1a", "Pou2af1",
      "Pou2f1", "Pou2f2", "Pou2f3", "Pou3f1", "Pou3f2", "Pou5f1", "Ppara", "Ppard", "Pparg",
      "Ppargc1b", "Prdm1", "Prdm10", "Prdm13", "Prdm14", "Prdm15", "Prdm16", "Prdm5", "Prdm9",
      "Prm1", "Prmt5", "Prox1", "Prrx1", "Psip1", "Pspc1", "Ptbp1", "Ptf1a", "Pth1r", "Pwp1", "Pygo2",
      "Rad21", "Rad23b", "Rad51", "Raf1", "Rag1", "Rag2", "Rai1", "Raly", "Rara", "Rarb", "Rarg", "Rb1",
      "Rbbp4", "Rbbp5", "Rbbp7", "Rbfox2", "Rbp2", "Rbpj", "Rbpjl", "Rcor2", "Rec8", "Rel", "Rela",
      "Relb", "Rest", "Rfx1", "Rfx2", "Rfx3", "Rif1", "Rnf2", "Rnf20", "Rora", "Rorc", "Rpa1", "Rpa2",
      "Rreb1", "Rtf1", "Runx1", "Runx2", "Runx3", "Ruvbl1", "Ruvbl2", "Rxra", "Rybp", "Safb", "Safb2",
      "Sall4", "Samd7", "Sap130", "Sap18", "Satb1", "Satb2", "Scml2", "Setd1a", "Setd2", "Setd5",
      "Setdb1", "Sfmbt1", "Sfmbt2", "Shox2", "Sin3a", "Sin3b", "Sirt1", "Sirt6", "Sirt7", "Six1", "Six2",
      "Six4", "Smad1", "Smad2", "Smad3", "Smad4", "Smad6", "Smad7", "Smarca2", "Smarca4",
      "Smarca5", "Smarcad1", "Smarcb1", "Smarcc1", "Smc1a", "Smc3", "Smc5", "Smchd1", "Smyd2",
      "Snai1", "Snai2", "Socs1", "Sox10", "Sox11", "Sox13", "Sox15", "Sox17", "Sox2", "Sox3", "Sox4",
      "Sox5", "Sox6", "Sox9", "Sp1", "Sp110", "Sp3", "Sp6", "Sp7", "Sp9", "Spen", "Spi1", "Spib", "Spin1",
      "Srcap", "Srebf1", "Srebf2", "Srf", "Srsf1", "Srsf2", "Ss18", "Ssrp1", "Ssu72", "Stab1", "Stag1",
      "Stag2", "Stat1", "Stat2", "Stat3", "Stat4", "Stat5a", "Stat5b", "Stat6", "Stra8", "Sumo1", "Sumo2",
      "Supt16", "Supt4a", "Supt5", "Supt6", "Suv39h2", "Suz12", "Sycp3", "T", "Taf1", "Taf10", "Taf12",
      "Taf15", "Taf2", "Taf3", "Taf7l", "Taf9b", "Tal1", "Tbl1x", "Tbp", "Tbpl1", "Tbr1", "Tbx19", "Tbx20",
      "Tbx21", "Tbx3", "Tbx4", "Tbx5", "Tcf12", "Tcf3", "Tcf4", "Tcf7", "Tcf7l1", "Tcf7l2", "Tcfl5", "Tdg",
      "Tead1", "Tead2", "Tead4", "Terf1", "Terf2ip", "Tet1", "Tet2", "Tet3", "Tfam", "Tfap2a", "Tfap2c",
      "Tfap4", "Tfcp2l1", "Tfe3", "Tfeb", "Tgif1", "Tgif2", "Thap1", "Thap11", "Thoc5", "Thra", "Thrap3",
      "Thrb", "Timeless", "Tirap", "Tle1", "Tle3", "Tlx1", "Tmem131", "Tmpo", "Tnni2", "Top1", "Top2a",
      "Top2b", "Tox", "Tox3", "Trim24", "Trim28", "Trim33", "Trmt112", "Trp53", "Trp63", "Trps1",
      "Tsc22d3", "Ttf1", "Twist1", "Ubtf", "Uhrf1", "Uhrf2", "Usf1", "Usf2", "Usp16", "Usp22", "Usp7",
      "Utf1", "Utp3", "Utp6", "Vdr", "Vsx2", "Wapl", "Was", "Wdr11", "Wdr5", "Wdr82", "Wiz", "Wt1",
      "Wtap", "Xbp1", "Xpo1", "Xrn2", "Yaf2", "Yap1", "Ybx1", "Yeats4", "Ythdc1", "Yy1", "Yy2",
      "Zbtb11", "Zbtb16", "Zbtb17", "Zbtb2", "Zbtb24", "Zbtb7b", "Zc3h10", "Zc3h11a", "Zc3h4",
      "Zeb1", "Zfp143", "Zfp217", "Zfp281", "Zfp322a", "Zfp335", "Zfp384", "Zfp407", "Zfp410",
      "Zfp42", "Zfp541", "Zfp57", "Zfp64", "Zfp708", "Zfp961", "Zfpm1", "Zfpm2", "Zic1", "Zic2", "Zic3",
      "Zkscan1", "Zmiz1", "Zmynd11", "Zmynd8", "Zscan4b", "Zscan4c", "Zscan5b"
    )
  } else if (genome_assembly == "mm9") {
    TFs_in_ChIP_Atlas = c(
      "Acaa2", "Acss2", "Actb", "Adnp", "Aebp2", "Aff3", "Aff4", "Ago2", "Ahcy", "Ahr", "Aicda", "Aire",
      "Alkbh5", "Anp32e", "Ar", "Arid1a", "Arid1b", "Arid2", "Arid3a", "Arid4b", "Arid5a", "Arnt",
      "Arntl", "Art1", "Arx", "Ascl1", "Ascl2", "Ash1l", "Ash2l", "Asxl1", "Asxl3", "Atf2", "Atf3", "Atf4",
      "Atf7", "Atf7ip", "Atm", "Atoh1", "Atrx", "Auts2", "Bach1", "Bach2", "Banp", "Bap1", "Batf",
      "Batf3", "Baz1a", "Baz1b", "Baz2a", "Bbx", "Bcl11a", "Bcl11b", "Bcl6", "Bcor", "Bend3", "Bhlhe40",
      "Bhlhe41", "Bmi1", "Bnc1", "Brca1", "Brca2", "Brd2", "Brd3", "Brd4", "Brd7", "Brd9", "Brdt",
      "Brwd1", "Btaf1", "Carm1", "Cbfa2t2", "Cbfa2t3", "Cbfb", "Cbx1", "Cbx2", "Cbx3", "Cbx4", "Cbx5",
      "Cbx7", "Cbx8", "Ccar1", "Ccnk", "Cd38", "Cdc73", "Cdk7", "Cdk8", "Cdk9", "Cdx2", "Cdyl",
      "Cebpa", "Cebpb", "Cebpd", "Cebpe", "Cebpg", "Chaf1b", "Chd1", "Chd2", "Chd4", "Chd7",
      "Chd8", "Chuk", "Cic", "Clock", "Creb1", "Creb3l2", "Crebbp", "Crtc2", "Crtc3", "Crx", "Cry1",
      "Ctbp2", "Ctcf", "Ctcfl", "Ctf1", "Ctnnb1", "Ctr9", "Cux2", "Cxxc1", "Cxxc5", "Daxx", "Dbp", "Ddb1",
      "Ddit3", "Ddx1", "Ddx18", "Ddx21", "Ddx5", "Dgcr8", "Dhx9", "Dicer1", "Dlx1", "Dlx3", "Dmap1",
      "Dmc1", "Dmrt1", "Dmrtb1", "Dnajc2", "Dnmt1", "Dnmt3a", "Dnmt3b", "Dnmt3l", "Dnttip1",
      "Dot1l", "Dpf2", "Dpf3", "Dppa2", "Dppa3", "Dppa4", "Dpy30", "Dr1", "Drosha", "Duxbl1",
      "Dyrk1a", "E2f1", "E2f3", "E2f4", "E2f6", "E4f1", "Ebf1", "Ebf2", "Ebf4", "Eed", "Egr1", "Egr2",
      "Ehf", "Ehmt2", "Elf1", "Elf4", "Elf5", "Elk1", "Elk3", "Elk4", "Ell", "Ell2", "Ell3", "Eloa", "Elob",
      "Emx2", "Eomes", "Ep300", "Ep400", "Epas1", "Epop", "Ercc3", "Erf", "Erg", "Esr1", "Esr2", "Esrra",
      "Esrrb", "Esrrg", "Ets1", "Ets2", "Etv2", "Etv5", "Etv6", "Ezh1", "Ezh2", "Fbxo21", "Fev", "Fgfr1",
      "Fli1", "Fos", "Fosb", "Fosl1", "Fosl2", "Foxa1", "Foxa2", "Foxa3", "Foxc1", "Foxc2", "Foxf1",
      "Foxg1", "Foxh1", "Foxj2", "Foxk1", "Foxk2", "Foxl2", "Foxo1", "Foxo3", "Foxp1", "Foxp2", "Foxp3",
      "Gabpa", "Gad1", "Gadd45a", "Gapdh", "Gata1", "Gata2", "Gata3", "Gata4", "Gata6", "Gckr",
      "Gcm1", "Gfi1", "Gfi1b", "Gli1", "Gli2", "Gli3", "Glis3", "Glyr1", "Gmnn", "Gps2", "Grhl2", "Grhl3",
      "Grip1", "Gsk3b", "Gtf2a2", "Gtf2b", "Gtf2e2", "Gtf2f1", "Gtf3c1", "Hand1", "Hand2", "Hcfc1",
      "Hdac1", "Hdac2", "Hdac3", "Hdac4", "Hdac5", "Hdac6", "Hdac7", "Hells", "Hes1", "Hexim1",
      "Hey1", "Heyl", "Hic2", "Hif1a", "Hira", "Hmga2", "Hmgn2", "Hnf1a", "Hnf1b", "Hnf4a", "Hnf4g",
      "Hnrnpk", "Hnrnpl", "Hnrnpu", "Hnrnpul1", "Hopx", "Hormad1", "Hoxa11", "Hoxa13", "Hoxa5",
      "Hoxa9", "Hoxb1", "Hoxb13", "Hoxb4", "Hoxc10", "Hoxc13", "Hoxc6", "Hoxc8", "Hoxc9",
      "Hoxd13", "Hoxd9", "Hsf1", "Hsf2", "Htatsf1", "Htt", "Id2", "Iho1", "Ikzf1", "Ikzf2", "Ikzf3", "Ing1",
      "Ing5", "Ino80", "Insm1", "Ints3", "Ints5", "Irf1", "Irf2", "Irf2bp2", "Irf3", "Irf4", "Irf5", "Irf7", "Irf8",
      "Irf9", "Isl1", "Jarid2", "Jmjd1c", "Jmjd6", "Jun", "Junb", "Jund", "Kansl3", "Kat2a", "Kat2b", "Kat5",
      "Kat7", "Kat8", "Kdm1a", "Kdm2a", "Kdm2b", "Kdm3a", "Kdm3b", "Kdm4a", "Kdm4b", "Kdm4c",
      "Kdm4d", "Kdm5a", "Kdm5b", "Kdm5c", "Kdm6a", "Kdm6b", "Klf1", "Klf15", "Klf2", "Klf3", "Klf4",
      "Klf5", "Klf6", "Klf9", "Kmt2a", "Kmt2b", "Kmt2c", "Kmt2d", "L3mbtl2", "L3mbtl3", "Ldb1", "Lef1",
      "Leo1", "Lhx2", "Lhx3", "Lhx5", "Lhx6", "Lin9", "Lmna", "Lmnb1", "Lmo2", "Lmx1b", "Lrrfip1",
      "Lyl1", "Macf1", "Maf", "Mafa", "Mafb", "Maff", "Mafg", "Mafk", "Mapk3", "Mapk8", "Matr3",
      "Mau2", "Max", "Maz", "Mbd1", "Mbd2", "Mbd3", "Mcrs1", "Mecp2", "Med1", "Med12",
      "Med14", "Med15", "Med23", "Med24", "Med26", "Med4", "Mef2a", "Mef2c", "Mef2d", "Mei4",
      "Meiosin", "Meis1", "Men1", "Meox2", "Mesp1", "Mettl14", "Mettl3", "Mitf", "Mllt10", "Mllt3",
      "Mlxipl", "Mnt", "Mnx1", "Morc3", "Mphosph8", "Mpnd", "Mre11a", "Msh2", "Msh6", "Msl1",
      "Msl2", "Mta1", "Mta2", "Mtf2", "Mtor", "Mxi1", "Myb", "Mybl1", "Myc", "Mycn", "Myf5",
      "Myo1c", "Myod1", "Myog", "Mysm1", "Myt1l", "N6amt1", "Naa10", "Naca", "Nacc1", "Nanog",
      "Nbn", "Ncapd3", "Ncapg", "Ncaph2", "Ncl", "Ncoa1", "Ncoa2", "Ncoa3", "Ncor1", "Ncor2",
      "Nelfa", "Nelfb", "Nelfcd", "Nelfe", "Neurod1", "Neurod2", "Neurog2", "Nfat5", "Nfatc1",
      "Nfatc2", "Nfe2", "Nfe2l1", "Nfe2l2", "Nfia", "Nfib", "Nfic", "Nfil3", "Nfix", "Nfkb1", "Nfkbia",
      "Nfya", "Nfyb", "Nfyc", "Nipbl", "Nkx2-1", "Nkx2-2", "Nkx2-5", "Nkx3-1", "Nkx3-2", "Nkx6-1",
      "Nlrc5", "Nlrp3", "Nono", "Notch1", "Npas3", "Npas4", "Npm1", "Nr0b1", "Nr0b2", "Nr1d1",
      "Nr1d2", "Nr1h2", "Nr1h3", "Nr1h4", "Nr1i2", "Nr2f1", "Nr2f2", "Nr3c1", "Nr4a1", "Nr5a1",
      "Nr5a2", "Nrf1", "Nsl1", "Nucks1", "Nup153", "Nup98", "Obox1", "Ogdh", "Ogt", "Olig2",
      "Onecut1", "Onecut2", "Otx2", "Ovol1", "Ovol2", "Paf1", "Pagr1a", "Park7", "Parp1", "Patz1",
      "Pawr", "Pax3", "Pax5", "Pax6", "Pax7", "Pax8", "Paxip1", "Pbrm1", "Pbx1", "Pcgf1", "Pcgf2",
      "Pcgf3", "Pcgf5", "Pcgf6", "Pcna", "Pds5a", "Pds5b", "Pdx1", "Peg3", "Pgr", "Phc1", "Phc2",
      "Phf10", "Phf19", "Phf2", "Phf20", "Phf5a", "Phf6", "Phf7", "Phf8", "Phrf1", "Pitx1", "Pknox1",
      "Plagl1", "Pml", "Polr1a", "Pou2af1", "Pou2f1", "Pou2f2", "Pou2f3", "Pou3f1", "Pou3f2", "Pou5f1",
      "Ppara", "Ppard", "Pparg", "Ppargc1b", "Prdm1", "Prdm10", "Prdm13", "Prdm14", "Prdm15",
      "Prdm16", "Prdm5", "Prdm9", "Prm1", "Prmt5", "Prox1", "Prrx1", "Psip1", "Pspc1", "Ptbp1",
      "Ptf1a", "Pth1r", "Pwp1", "Pygo2", "Rad21", "Rad23b", "Rad51", "Raf1", "Rag1", "Rag2", "Rai1",
      "Raly", "Rara", "Rarb", "Rarg", "Rb1", "Rbbp4", "Rbbp5", "Rbbp7", "Rbfox2", "Rbp2", "Rbpj",
      "Rbpjl", "Rcor2", "Rec8", "Rel", "Rela", "Relb", "Rest", "Rfx1", "Rfx2", "Rfx3", "Rif1", "Rnf2",
      "Rnf20", "Rora", "Rorc", "Rpa1", "Rpa2", "Rreb1", "Rtf1", "Runx1", "Runx2", "Runx3", "Ruvbl1",
      "Ruvbl2", "Rxra", "Rybp", "Safb", "Safb2", "Sall1", "Sall4", "Samd7", "Sap130", "Sap18", "Satb1",
      "Satb2", "Scml2", "Setd1a", "Setd2", "Setd5", "Setdb1", "Sfmbt1", "Sfmbt2", "Shox2", "Sin3a",
      "Sin3b", "Sirt1", "Sirt6", "Sirt7", "Six1", "Six2", "Six4", "Smad1", "Smad2", "Smad3", "Smad4",
      "Smad6", "Smad7", "Smarca2", "Smarca4", "Smarca5", "Smarcad1", "Smarcb1", "Smarcc1",
      "Smc1a", "Smc3", "Smc5", "Smchd1", "Smyd2", "Snai1", "Snai2", "Socs1", "Sox10", "Sox11",
      "Sox13", "Sox15", "Sox17", "Sox2", "Sox3", "Sox30", "Sox4", "Sox5", "Sox6", "Sox9", "Sp1",
      "Sp110", "Sp3", "Sp6", "Sp7", "Sp9", "Spen", "Spi1", "Spib", "Spin1", "Srcap", "Srebf1", "Srebf2",
      "Srf", "Srsf1", "Srsf2", "Ss18", "Ssrp1", "Ssu72", "Stab1", "Stag1", "Stag2", "Stat1", "Stat2", "Stat3",
      "Stat4", "Stat5a", "Stat5b", "Stat6", "Stra8", "Sumo1", "Sumo2", "Supt16", "Supt4a", "Supt5",
      "Supt6", "Suv39h2", "Suz12", "Sycp3", "T", "Taf1", "Taf10", "Taf12", "Taf15", "Taf2", "Taf3", "Taf7l",
      "Taf9b", "Tal1", "Tbl1x", "Tbp", "Tbpl1", "Tbr1", "Tbx19", "Tbx20", "Tbx21", "Tbx3", "Tbx4",
      "Tbx5", "Tcf12", "Tcf3", "Tcf4", "Tcf7", "Tcf7l1", "Tcf7l2", "Tcfl5", "Tdg", "Tead1", "Tead2", "Tead4",
      "Terf1", "Terf2ip", "Tet1", "Tet2", "Tet3", "Tfam", "Tfap2a", "Tfap2c", "Tfap4", "Tfcp2l1", "Tfe3",
      "Tfeb", "Tgif1", "Tgif2", "Thap1", "Thap11", "Thoc5", "Thra", "Thrap3", "Thrb", "Timeless", "Tirap",
      "Tle1", "Tle3", "Tlx1", "Tmem131", "Tmpo", "Tnni2", "Top1", "Top2a", "Top2b", "Tox", "Tox3",
      "Trim24", "Trim28", "Trim33", "Trmt112", "Trp53", "Trp63", "Trps1", "Tsc22d3", "Ttf1", "Twist1",
      "Ubtf", "Uhrf1", "Uhrf2", "Usf1", "Usf2", "Usp16", "Usp22", "Usp7", "Utf1", "Utp3", "Utp6", "Vdr",
      "Vsx2", "Wapl", "Was", "Wdr11", "Wdr5", "Wdr82", "Wiz", "Wt1", "Wtap", "Xbp1", "Xpo1",
      "Xrn2", "Yaf2", "Yap1", "Ybx1", "Yeats4", "Ythdc1", "Yy1", "Yy2", "Zbtb11", "Zbtb16", "Zbtb17",
      "Zbtb2", "Zbtb24", "Zbtb7b", "Zc3h10", "Zc3h11a", "Zc3h4", "Zeb1", "Zfp143", "Zfp217",
      "Zfp281", "Zfp322a", "Zfp335", "Zfp384", "Zfp407", "Zfp410", "Zfp42", "Zfp541", "Zfp57", "Zfp64",
      "Zfp708", "Zfp961", "Zfpm1", "Zfpm2", "Zic1", "Zic2", "Zic3", "Zkscan1", "Zmiz1", "Zmynd11",
      "Zmynd8", "Zscan4b", "Zscan4c", "Zscan5b"
    )
  } else if (genome_assembly == "rn6") {
    TFs_in_ChIP_Atlas = c(
      "Ahr", "Ar", "Arnt", "Arnt2", "Brd4", "Chd7", "Creb1", "Creb3l1", "Crebbp", "Crtc2", "Ctcf",
      "Ctnnb1", "Dbp", "Dlx4", "Egr1", "Egr2", "Ep300", "Esr1", "Foxe1", "Foxo1", "Hdac2", "Hdac3",
      "Hils1", "Hmgb2", "Hnf1a", "Hnf4a", "Ino80", "Jun", "Kdm1a", "Lmnb1", "Matr3", "Mbd3",
      "Med1", "Mef2a", "Mkx", "Mlxipl", "Myt1", "Neurod1", "Nfia", "Nkx2-1", "Npas4", "Nr1h4",
      "Nr3c1", "Nr3c2", "Olig2", "Parp1", "Pax8", "Pdx1", "Pgr", "Pou1f1", "Rest", "Runx1", "Satb1",
      "Seh1l", "Sox10", "Sox6", "Sox9", "Sp1", "Srebf1", "Stat1", "Suz12", "Taf3", "Tcf4", "Tcf7l2",
      "Xrcc1", "Ybx1", "Zhx2"
    )
  } else if (genome_assembly == "dm6") {
    TFs_in_ChIP_Atlas = c(
      "AGO2", "Abd-B", "Acf", "Act5C", "Ada2b", "Adf1", "Aef1", "Art1", "BEAF-32", "BRWD3", "BigH1",
      "Br140", "CBP", "CG10543", "CG10949", "CG11902", "CG11906", "CG42741", "CG7099",
      "CG7691", "CG8478", "CG9609", "CTCF", "Cap-D3", "Cap-H2", "Cdc6", "Cdk12", "Cdk8", "Cdk9",
      "Cf2", "Chd1", "Chro", "Clamp", "Clk", "Cp190", "D", "Dfd", "Dll", "DnaJ-1", "Dp1", "Dref", "Dsp1",
      "E(bx)", "E(z)", "EcR", "Elba2", "Elba3", "Etl1", "Fer3", "FoxK", "Gcn5", "Gug", "HDAC4", "HDAC6",
      "HGTX", "HP1b", "HP1c", "HP4", "Hmr", "Hnf4", "Hr39", "Hr4", "Hr78", "Hrb27C", "Hrb87F",
      "Hsp83", "Ibf1", "Ibf2", "Ice1", "Irbp", "Iswi", "JIL-1", "Jarid2", "Kdm2", "Kdm4A", "Kr", "Kr-h1",
      "Lam", "Lgr1", "Lint-O", "Lpt", "LysE", "M1BP", "MBD-R2", "MRG15", "Mad", "Max", "Mcm2",
      "Med", "Mef2", "Mes2", "Mi-2", "Mlf", "Myc", "N", "NELF-B", "Nelf-A", "Nelf-E", "Nup107",
      "Nup93-1", "Nup98-96", "Oga", "Orc2", "PCID2", "Pc", "Pcl", "Pof", "Psc", "Ptip", "Rbf", "Rbp1",
      "Rrp40", "Rrp6", "Rsf1", "RunxB", "SF1", "SMC3", "Saf6", "Sce", "Scm", "Set1", "Sfmbt", "Sgf11",
      "Sidpn", "Sin3A", "Smox", "SoxN", "Spps", "Spt3", "Spt5", "Ssrp", "Stat92E", "Su(H)", "Su(Tpl)",
      "Su(var)2-HP2", "Su(var)205", "Su(var)3-3", "Su(var)3-7", "Su(var)3-9", "Su(z)12", "SuUR", "Tbp",
      "Tet", "TfAP-2", "TfIIB", "TfIIFalpha", "Trf2", "Trl", "Ubqn", "Ubx", "Upf1", "Utx", "ZIPIC", "Zaf1",
      "ab", "abd-A", "aop", "apt", "arm", "ash1", "ash2", "atms", "bap", "bcd", "bin", "brk", "brm", "bs",
      "cad", "cbt", "cg", "cic", "cnc", "cuff", "disco", "dl", "dre4", "drm", "dve", "dwg", "elav", "en",
      "enok", "erm", "esg", "exd", "ey", "fd19B", "feo", "foxo", "fs(1)h", "gfzf", "grh", "gro", "gt", "hb",
      "hkb", "hng3", "hth", "insv", "kis", "kn", "l(1)sc", "l(3)mbt", "lilli", "lmd", "lola", "lov", "mago",
      "mid", "mle", "mod(mdg4)", "mof", "msl-1", "msl-2", "msl-3", "mxc", "not", "nsl1", "nub", "opa",
      "osa", "pan", "pdm2", "ph-d", "ph-p", "pho", "pita", "polybromo", "prd", "psq", "pzg", "rhi",
      "row", "rump", "sc", "scro", "sd", "shep", "shn", "slp1", "sna", "so", "stil", "su(Hw)", "term", "tll",
      "trem", "trr", "trx", "tsu", "twi", "upSET", "ush", "usp", "velo", "vtd", "wda", "wds", "woc", "yki",
      "z", "zen"
    )
  } else if (genome_assembly == "dm3") {
    TFs_in_ChIP_Atlas = c(
      "AGO2", "Abd-B", "Acf", "Act5C", "Ada2b", "Adf1", "Aef1", "Art1", "BEAF-32", "BRWD3", "BigH1",
      "Br140", "CBP", "CG10543", "CG10949", "CG11902", "CG11906", "CG42741", "CG7099",
      "CG7691", "CG8478", "CG9609", "CTCF", "Cap-D3", "Cap-H2", "Cdc6", "Cdk12", "Cdk8", "Cdk9",
      "Cf2", "Chd1", "Chro", "Clamp", "Clk", "Cp190", "D", "Dfd", "Dll", "DnaJ-1", "Dp1", "Dref", "Dsp1",
      "E(bx)", "E(z)", "EcR", "Elba2", "Elba3", "Etl1", "FoxK", "Gcn5", "Gug", "HDAC4", "HDAC6", "HGTX",
      "HP1b", "HP1c", "HP4", "Hmr", "Hr39", "Hrb27C", "Hrb87F", "Hsp83", "Ibf1", "Ibf2", "Ice1",
      "Irbp", "Iswi", "JIL-1", "Jarid2", "Kdm2", "Kdm4A", "Kr", "Kr-h1", "Lam", "Lgr1", "Lint-O", "Lpt",
      "LysE", "M1BP", "MBD-R2", "MRG15", "Mad", "Max", "Mcm2", "Med", "Mef2", "Mes2", "Mi-2",
      "Mlf", "Myc", "N", "NELF-B", "Nelf-A", "Nelf-E", "Nup107", "Nup93-1", "Nup98-96", "Oga", "Orc2",
      "PCID2", "Pc", "Pcl", "Pof", "Psc", "Rbf", "Rbp1", "Rrp40", "Rrp6", "Rsf1", "RunxB", "SF1", "SMC3",
      "Saf6", "Sce", "Scm", "Set1", "Sfmbt", "Sgf11", "Sidpn", "Sin3A", "Smox", "SoxN", "Spps", "Spt3",
      "Spt5", "Ssrp", "Stat92E", "Su(H)", "Su(Tpl)", "Su(var)2-HP2", "Su(var)205", "Su(var)3-3",
      "Su(var)3-7", "Su(var)3-9", "Su(z)12", "SuUR", "Tbp", "Tet", "TfAP-2", "TfIIB", "TfIIFalpha", "Trf2",
      "Trl", "Ubqn", "Ubx", "Upf1", "Utx", "ZIPIC", "Zaf1", "ab", "abd-A", "aop", "apt", "arm", "ash1",
      "ash2", "atms", "bap", "bcd", "bin", "brk", "brm", "bs", "cad", "cbt", "cg", "cic", "cuff", "disco",
      "dl", "dre4", "drm", "dve", "dwg", "elav", "en", "enok", "erm", "esg", "eve", "exd", "ey", "fd19B",
      "feo", "fs(1)h", "gfzf", "grh", "gro", "gt", "hb", "hng3", "hth", "insv", "kis", "kni", "l(1)sc", "l(3)mbt",
      "lilli", "lmd", "lola", "mago", "mid", "mle", "mod(mdg4)", "mof", "msl-1", "msl-2", "msl-3", "mxc",
      "not", "nub", "opa", "osa", "pan", "ph-d", "ph-p", "pho", "pita", "polybromo", "psq", "pzg", "rhi",
      "row", "rump", "sc", "sd", "shep", "shn", "slp1", "sna", "so", "stil", "su(Hw)", "trem", "trr", "trx",
      "tsu", "twi", "upSET", "ush", "usp", "velo", "vtd", "wda", "wds", "woc", "yki", "z", "zen"
    )
  } else if (genome_assembly == "ce11") {
    TFs_in_ChIP_Atlas = c(
      "aha-1", "alr-1", "aly-2", "ama-1", "ash-2", "atfs-1", "blmp-1", "capg-1", "capg-2", "cbp-1", "cec-3",
      "cec-5", "cec-7", "ceh-16", "ceh-28", "ceh-30", "ceh-37", "ceh-39", "ceh-51", "ceh-82", "ces-1",
      "chd-3", "coh-1", "coh-3", "csr-1", "daf-12", "daf-16", "dpl-1", "dpy-26", "dpy-27", "dpy-28",
      "dpy-30", "dve-1", "efl-1", "egl-27", "egl-5", "elt-1", "elt-2", "elt-3", "eor-1", "epc-1", "fkh-10",
      "fkh-2", "fkh-7", "fos-1", "gmeb-1", "ham-1", "hcp-3", "hcp-4", "hcp-6", "hda-1", "hif-1", "him-1",
      "him-17", "him-3", "him-5", "him-8", "hlh-1", "hlh-30", "hpl-2", "htp-3", "imb-1", "jmjd-2",
      "jun-1", "kle-2", "lag-1", "lem-2", "let-418", "let-607", "lin-11", "lin-13", "lin-15B", "lin-35",
      "lin-36", "lin-37", "lin-39", "lin-42", "lin-52", "lin-53", "lin-54", "lin-61", "lin-9", "lsy-2", "mab-5",
      "mau-2", "mdl-1", "mdt-15", "med-1", "mef-2", "mep-1", "met-2", "mis-12", "mix-1", "mml-1",
      "mre-11", "mrg-1", "msh-5", "nfya-1", "nhr-10", "nhr-11", "nhr-116", "nhr-12", "nhr-129",
      "nhr-2", "nhr-21", "nhr-23", "nhr-237", "nhr-25", "nhr-28", "nhr-6", "nhr-67", "nhr-76", "nhr-77",
      "npp-13", "npp-3", "nurf-1", "pax-1", "peb-1", "pes-1", "pha-4", "pop-1", "pqm-1", "pqn-85",
      "rec-8", "rpc-1", "sax-3", "scc-1", "sdc-1", "sdc-2", "sdc-3", "sea-2", "sem-4", "sin-3", "skn-1",
      "sma-3", "sma-4", "sma-9", "smc-3", "smc-4", "smc-6", "sptf-1", "taf-1", "tbp-1", "tlp-1", "top-2",
      "tra-1", "unc-119", "unc-130", "unc-39", "unc-55", "unc-62", "wdr-5", "zag-1", "zfp-1", "zhp-3",
      "zim-1", "zim-3", "zip-2", "zip-8", "ztf-11", "ztf-4", "ztf-7"
    )
  } else if (genome_assembly == "ce10") {
    TFs_in_ChIP_Atlas = c(
      "aha-1", "alr-1", "aly-2", "ama-1", "ash-2", "atfs-1", "blmp-1", "capg-1", "capg-2", "cbp-1", "cec-3",
      "cec-5", "cec-7", "ceh-16", "ceh-28", "ceh-30", "ceh-37", "ceh-39", "ceh-51", "ceh-82", "ces-1",
      "chd-3", "coh-1", "coh-3", "csr-1", "daf-12", "daf-16", "dpl-1", "dpy-26", "dpy-27", "dpy-28",
      "dpy-30", "dve-1", "efl-1", "egl-27", "egl-5", "elt-1", "elt-2", "elt-3", "eor-1", "epc-1", "fkh-10",
      "fkh-2", "fkh-7", "fos-1", "gmeb-1", "ham-1", "hcp-3", "hcp-4", "hcp-6", "hda-1", "hif-1", "him-1",
      "him-17", "him-3", "him-5", "him-8", "hlh-1", "hlh-30", "hpl-2", "htp-3", "imb-1", "jmjd-2",
      "jun-1", "kle-2", "lag-1", "lem-2", "let-418", "let-607", "lin-11", "lin-13", "lin-15B", "lin-35",
      "lin-36", "lin-37", "lin-39", "lin-42", "lin-52", "lin-53", "lin-54", "lin-61", "lin-9", "lsy-2", "mab-5",
      "mau-2", "mdl-1", "mdt-15", "med-1", "mef-2", "mep-1", "met-2", "mis-12", "mix-1", "mml-1",
      "mre-11", "mrg-1", "msh-5", "nfya-1", "nhr-10", "nhr-11", "nhr-116", "nhr-12", "nhr-129",
      "nhr-2", "nhr-21", "nhr-23", "nhr-237", "nhr-25", "nhr-28", "nhr-6", "nhr-67", "nhr-76", "nhr-77",
      "npp-13", "npp-3", "nurf-1", "pax-1", "peb-1", "pes-1", "pha-4", "pop-1", "pqm-1", "pqn-85",
      "rec-8", "rpc-1", "sax-3", "scc-1", "sdc-1", "sdc-2", "sdc-3", "sea-2", "sem-4", "sin-3", "skn-1",
      "sma-3", "sma-4", "sma-9", "smc-3", "smc-4", "smc-6", "sptf-1", "taf-1", "tbp-1", "tlp-1", "top-2",
      "tra-1", "unc-119", "unc-130", "unc-39", "unc-55", "unc-62", "wdr-5", "zag-1", "zfp-1", "zhp-3",
      "zim-1", "zim-3", "zip-2", "zip-8", "ztf-11", "ztf-4", "ztf-7"
    )
  } else if (genome_assembly == "sacCer3") {
    TFs_in_ChIP_Atlas = c(
      "ABF1", "ASF1", "CHD1", "CSE4", "CST6", "CTF4", "DPB3", "ESA1", "FHL1", "GAL11", "GCN4",
      "HMO1", "HSF1", "HST2", "HST3", "HST4", "IFH1", "INO80", "ISW1", "LEU3", "MCD1", "MCM1",
      "MCM2", "MED7", "NCB2", "ORC1", "PCH2", "PGD1", "RAD51", "RAD52", "RAP1", "REB1", "REC8",
      "RED1", "RFA1", "RGR1", "RGT1", "RLF2", "ROX3", "RPL25", "RRP6", "RTT106", "SET1", "SFP1",
      "SIR3", "SLD2", "SMC3", "SMC4", "SMC6", "SMT3", "SNF2", "SPT15", "SPT16", "SRB4", "SSN8",
      "SUA7", "TAF1", "TAF10", "TAF12", "TAF2", "TAF4", "TAF5", "TAF8", "TFG1", "TOA2", "TOP2",
      "UME6"
    )
  } else {
    stop("Parameter genome_assembly input error.")
  }

  ### Get candidate TFs from the built-in ChIP-Atlas list.
  message("Select ChIP-Atlas TFs present in scRNA-seq data.")
  TFs_from_priors = base::sort(base::unique(base::intersect(
    TFs_in_ChIP_Atlas,
    base::rownames(scRNAseq)
  )))

  # Keep only ChIP-Atlas TFs that are detected in the scRNA-seq data as the candidate pool.
  message("Extract raw gene-by-cell count matrix for expression quantification.")
  gene_counts_all = Seurat::GetAssayData(scRNAseq, layer = "counts")

  ### Filter candidate TFs by expression quality.
  message("Compute per-gene expression statistics (non-zero frequency, mean in non-zero cells, overall mean).")
  n_cells = base::ncol(gene_counts_all)
  row_sums = Matrix::rowSums(gene_counts_all)
  n_nonzero = Matrix::rowSums(gene_counts_all != 0)
  TF_features_all = base::data.frame(
    nonzero_freq = n_nonzero / n_cells,
    # Set the conditional mean to 0 for genes detected in no cell (avoid 0/0 = NaN).
    nonzero_mean = base::ifelse(n_nonzero == 0, 0, row_sums / n_nonzero),
    overall_mean = row_sums / n_cells,
    row.names = base::rownames(gene_counts_all),
    check.names = FALSE
  )
  base::rm(gene_counts_all, row_sums, n_nonzero, n_cells)

  message("Add percentile ranks for each feature (0 = lowest, 1 = highest).")
  TF_features_all[base::paste0(base::names(TF_features_all), "_quantile")] =
    base::lapply(TF_features_all, function(x) base::rank(x) / base::length(x))
  # Restrict the genome-wide feature table to candidate TFs, then flag TF-side columns.
  TF_features_all = TF_features_all[base::intersect(TFs_from_priors, base::rownames(TF_features_all)), ]
  base::colnames(TF_features_all) = base::paste0("TF_", base::colnames(TF_features_all))

  message("Apply expression quality filters.")
  TF_features_all$keep = base::with(
    TF_features_all,
    TF_nonzero_freq >= min_TF_nonzero_freq &
      TF_nonzero_mean >= min_TF_nonzero_mean &
      TF_overall_mean >= min_TF_overall_mean &
      TF_nonzero_freq_quantile >= min_TF_nonzero_freq_quantile &
      TF_nonzero_mean_quantile >= min_TF_nonzero_mean_quantile &
      TF_overall_mean_quantile >= min_TF_overall_mean_quantile
  )

  message("Extract TFs passing all quality thresholds.")
  TFs_filtered_by_features = base::sort(base::rownames(TF_features_all[TF_features_all$keep, ]))

  ### Download ChIP-Atlas target-gene priors with local caching.
  # Initialize the per-TF download status tracker.
  message("Initialize download status tracker.")
  nTF = base::length(TFs_filtered_by_features)
  download_information = base::data.frame(base::matrix(
    "", nrow = nTF, ncol = 4,
    dimnames = base::list(TFs_filtered_by_features, base::c("url", "destination", "status", "time"))
  ))

  # Create the local directory that stores downloaded TSV files.
  message("Set up ChIP-Atlas prior download directory.")
  if (!base::dir.exists("./ChIP_Atlas_Priors")) {
    base::dir.create("./ChIP_Atlas_Priors", recursive = TRUE)
  }

  # Download one target-gene table per TF, reusing cached files when they already exist.
  message("Download ChIP-Atlas target gene data for each TF.")
  original_timeout = base::getOption("timeout")
  base::options(timeout = timeout)
  time0 = base::Sys.time()
  for (i in base::seq_along(TFs_filtered_by_features)) {
    time1 = base::Sys.time()
    tf = TFs_filtered_by_features[i]

    # Build the request URL and the local destination path.
    download_information[tf, "url"] = base::sprintf(
      "https://chip-atlas.dbcls.jp/data/%s/target/%s.%s.tsv",
      genome_assembly, tf, distance_from_TSS
    )
    download_information[tf, "destination"] = base::file.path(
      base::getwd(), "ChIP_Atlas_Priors", base::paste0(tf, ".tsv")
    )

    # Reuse an existing local file instead of downloading it again.
    if (base::file.exists(download_information[tf, "destination"])) {
      download_information[tf, "status"] = "already downloaded"
      download_information[tf, "time"] = base::format(base::Sys.time() - time1)
      base::Sys.sleep(delay_time)
      message(base::sprintf(
        "[%d/%d] [skip] %s - already downloaded - %s",
        i, nTF, tf, base::format(base::Sys.time() - time0)
      ))
      next
    }

    # Attempt the download; treat any error or warning as a failed attempt.
    download_successful = tryCatch({
      utils::download.file(
        download_information[tf, "url"],
        download_information[tf, "destination"],
        mode = "wb", quiet = TRUE
      )
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)

    # Accept only a non-empty file; otherwise flag the TF as failed.
    if (download_successful && base::file.size(download_information[tf, "destination"]) > 0) {
      download_information[tf, "status"] = "successful"
      download_information[tf, "time"] = base::format(base::Sys.time() - time1)
      message(base::sprintf(
        "[%d/%d] [ok] %s - successful - %s",
        i, nTF, tf, base::format(base::Sys.time() - time0)
      ))
    } else {
      # Delete any corrupted or empty file left behind by the failed attempt.
      if (base::file.exists(download_information[tf, "destination"])) {
        base::file.remove(download_information[tf, "destination"])
      }
      download_information[tf, "status"] = "URL not found or download error"
      download_information[tf, "time"] = base::format(base::Sys.time() - time1)
      message(base::sprintf(
        "[%d/%d] [fail] %s - URL not found or download error - %s",
        i, nTF, tf, base::format(base::Sys.time() - time0)
      ))
    }

    base::Sys.sleep(delay_time)
  }
  base::options(timeout = original_timeout)

  message("Save download information.")
  base::saveRDS(download_information, base::paste0("./ChIP_Atlas_Priors/download_information.rds"))

  ### Process priors and validate their quality.
  message("Process raw priors: remove STRING column, filter cell types, compute mean intensity, normalize & transform, filter low-quality data.")
  priors = base::list()
  for (tf in TFs_filtered_by_features) {
    # Skip TFs whose download failed.
    if (download_information[tf, "status"] == "URL not found or download error") {
      priors[[tf]] = NULL
      next
    }

    # Read the raw target-gene table (row names are target gene symbols).
    priors[[tf]] = utils::read.delim(
      download_information[tf, "destination"],
      row.names = 1, stringsAsFactors = FALSE
    )

    # Discard empty tables.
    if (base::nrow(priors[[tf]]) == 0) {
      priors[[tf]] = NULL
      next
    }

    # Drop the STRING association column; only ChIP-seq sample columns are used as prior.
    priors[[tf]][["STRING"]] = NULL

    # Remove sample columns whose name matches any excluded cell-type keyword.
    exclude_col_mask = base::grepl(
      base::paste(excluded_cell_types_in_priors, collapse = "|"),
      base::colnames(priors[[tf]]),
      ignore.case = ignore_case_when_filtering_priors
    )
    priors[[tf]] = priors[[tf]][, !exclude_col_mask, drop = FALSE]

    # Require at least 2 columns: column 1 becomes the mean, columns 2+ are averaged.
    if (base::ncol(priors[[tf]]) <= 1) {
      priors[[tf]] = NULL
      next
    }

    # Average intensity across the retained samples, storing the result in column 1.
    priors[[tf]][, 1] = base::rowMeans(priors[[tf]][, -1, drop = FALSE])

    # Remove target genes with zero mean intensity.
    priors[[tf]] = priors[[tf]][priors[[tf]][, 1] > 0, ]

    # Discard priors with too few valid target genes.
    if (base::nrow(priors[[tf]]) < prior_min_TG_number) {
      priors[[tf]] = NULL
      next
    }

    # Rank target genes by decreasing mean intensity.
    priors[[tf]] = priors[[tf]][base::order(-priors[[tf]][, 1]), , drop = FALSE]

    # Min-max normalize the mean (with a small floor/ceiling), then power-transform it into 'intensity'.
    min_average = base::min(priors[[tf]][, 1]) * 0.2
    max_average = base::max(priors[[tf]][, 1]) + 0.0001
    priors[[tf]][, "intensity"] = ((priors[[tf]][, 1] - min_average) / (max_average - min_average)) ** intensity_scale_factor
  }
  message("Save processed priors.")
  base::saveRDS(priors, base::paste0("./ChIP_Atlas_Priors/processed_priors.rds"))

  ### Simplify priors and extract the final TFs.
  message("Simplify priors into TF-TG intensity.")
  priors_simplify = base::list()
  for (tf in TFs_filtered_by_features) {
    if (!base::is.null(priors[[tf]])) {
      # Retain only the normalized intensity column.
      priors_simplify[[tf]] = priors[[tf]][, "intensity", drop = FALSE]

      # Keep top targets under the dual absolute/proportional cap, always retaining at least one.
      priors_simplify[[tf]] = priors_simplify[[tf]][
        1:base::max(
          1, base::round(base::min(
            prior_top_n_keep,
            base::nrow(priors[[tf]]) * prior_top_q_keep
          ))
        ), , drop = FALSE
      ]
    }
  }
  message("Save simplified processed priors.")
  base::saveRDS(priors_simplify, base::paste0("./ChIP_Atlas_Priors/simplified_processed_priors.rds"))

  # TFs with a valid simplified prior constitute the final screened TF set.
  message("Extract final TFs for downstream analysis.")
  TFs = base::names(priors_simplify)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Retrieve TFs from prior database. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(TFs)
}
