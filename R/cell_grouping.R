#' @title Adaptive Cell Grouping Along Pseudotime Trajectories
#'
#' @description Performs adaptive cell grouping along pseudotime trajectories for each cell type.
#' It first fits a logarithmic model to determine the minimum number of cells per group from the total cell count,
#' then iteratively enlarges groups so that every retained gene has sufficient non-zero expression values in each group.
#' It returns aggregated gene expression profiles and pseudotime endpoints for downstream dynamic modeling.
#'
#' @param interest_cell_type_genes_pseudotime_info List. Output from \code{\link{get_genes_pseudotime_info}}.
#' @param points_x_for_fitting_nc_and_nmin Integer vector of length 3, strictly increasing positive values.
#' Total cell counts (x-axis) used to fit the minimum-cells-per-group logarithmic model.
#' @param points_y_for_fitting_nc_and_nmin Integer vector of length 3, strictly increasing positive values.
#' Corresponding minimum cells per group (y-axis) for the logarithmic model.
#' Defaults to \code{base::c(5, 10, 15)}.
#' @param ncores Positive integer. Number of cores for parallel computation via \code{\link[parallel]{mclapply}}.
#' Defaults to \code{1}. Fork-based parallelism is supported only on Unix-like systems;
#' set to \code{1} when running on Windows.
#'
#' @return Nested list named by cell type; each entry contains:
#' \itemize{
#'   \item \code{Branches_n_group}: List of integers. Number of valid groups per branch (\code{0} if grouping fails for that branch).
#'   \item \code{Branches_n_k}: List of integer vectors. Number of cells in each group per branch.
#'   \item \code{n_min}: Numeric vector. Fitted minimum cell threshold for each branch.
#'   \item \code{Branches_T_k}: List of numeric vectors. Pseudotime value at the endpoint of each group.
#'   \item \code{Branches_Tslot_k}: List of numeric vectors. Pseudotime interval length covered by each group.
#'   \item \code{Branches_TFE_T}: List of data frames (genes as rows, groups as columns).
#'     Mean of non-zero TF expression values within each group; \code{NA} when a TF is zero throughout a group.
#'   \item \code{Branches_TGE_T}: List of data frames (genes as rows, groups as columns).
#'     Mean of non-zero TG expression values within each group; \code{NA} when a TG is zero throughout a group.
#' }
#'
#' @details For each cell type and branch, the workflow executes as follows:
#' \enumerate{
#'   \item Validates the three fitting points so that they are strictly increasing positive integers.
#'   \item Fits \eqn{y = c \cdot \ln(ax + b)} once from the three calibration points
#'         (x = total cells in a branch, y = minimum cells per group) and rounds the prediction up to obtain \code{n_min}.
#'   \item Groups pseudotime-ordered cells adaptively: each group starts with \code{n_min} cells;
#'         if any gene in the TFE or TGE matrix has fewer than \code{n_min} non-zero values in the
#'         current group, one cell is added at a time until the constraint is met or cells are exhausted.
#'   \item Handles edge cases: fewer than \code{n_min} remaining cells are merged into the previous group;
#'         if the first group cannot meet the constraint, that branch returns \code{0} groups.
#'   \item Computes the endpoint pseudotime and interval length of each group.
#'   \item Averages the non-zero expression values of each gene within each group and returns the aggregated TFE/TGE matrices.
#' }
#'
#' @note
#' \itemize{
#'   \item The three fitting points determine \code{n_min} and hence the number of groups.
#'         A smaller \code{n_min} yields more groups but increases sparsity risk.
#'   \item A branch that cannot satisfy the constraint returns \code{0} groups together with empty
#'         aggregated objects; other branches of the same cell type are still processed.
#'         To recover such a branch, lower the y calibration values to reduce \code{n_min},
#'         or tune upstream parameters in \code{\link{get_genes_pseudotime_info}} to reduce gene sparsity.
#'   \item For reliable downstream modeling, 8-15 groups per branch is recommended: too few groups
#'         reduce temporal resolution, while too many groups weaken per-group statistical power.
#'   \item The function creates the \code{2.2 Data Processing - Cell Grouping} directory, temporarily
#'         switches the working directory (restored on exit), and saves one logarithmic-fit plot per cell type.
#'   \item Parallel computation relies on \code{\link[parallel]{mclapply}} and is not supported on Windows.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_genes_pseudotime_info = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds"
#' )
#' points_x_for_fitting_nc_and_nmin = base::c(280, 450, 650)
#' points_y_for_fitting_nc_and_nmin = base::c(5, 10, 15)
#' ncores = parallel::detectCores() - 1 # In Linux.
#' # ncores = 1 # In Windows.
#' interest_cell_type_group = cell_grouping(
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   points_x_for_fitting_nc_and_nmin = points_x_for_fitting_nc_and_nmin,
#'   points_y_for_fitting_nc_and_nmin = points_y_for_fitting_nc_and_nmin,
#'   ncores = ncores
#' )
#' base::saveRDS(
#'   interest_cell_type_group,
#'   file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds"
#' )
#' }
cell_grouping = function(
    interest_cell_type_genes_pseudotime_info,
    points_x_for_fitting_nc_and_nmin,
    points_y_for_fitting_nc_and_nmin = base::c(5, 10, 15),
    ncores = 1
) {
  ### Validate parameters.
  # The calibration points must be strictly increasing positive integers of length 3.
  validate_points = function(x, name = "x") {
    if (!base::is.numeric(x)) {
      base::stop("Parameter '", name, "' must be a numeric vector.")
    }
    if (base::length(x) != 3) {
      base::stop("Parameter '", name, "' must have exactly 3 elements.")
    }
    if (base::any(x < 0)) {
      base::stop("Parameter '", name, "' must contain only positive values.")
    }
    if (base::any(x %% 1 != 0)) {
      base::stop("Parameter '", name, "' must contain only integers.")
    }
    if (!base::all(base::diff(x) > 0)) {
      base::stop("Parameter '", name, "' must be strictly increasing.")
    }
    message("Validation passed for parameter '", name, "'.")
  }

  message("Validating points_x_for_fitting_nc_and_nmin.")
  validate_points(points_x_for_fitting_nc_and_nmin, name = "points_x_for_fitting_nc_and_nmin")

  message("Validating points_y_for_fitting_nc_and_nmin.")
  validate_points(points_y_for_fitting_nc_and_nmin, name = "points_y_for_fitting_nc_and_nmin")

  if (!(base::is.numeric(ncores) && base::length(ncores) == 1 && ncores >= 1)) {
    base::stop("ncores must be a single positive number.")
  }

  # Fork-based parallelism is unavailable on Windows.
  if (.Platform$OS.type == "windows" && ncores > 1) {
    message("Parallel computing is not supported on Windows; ncores forced to 1.")
    ncores = 1
  }

  # Every cell type must carry the outputs of get_genes_pseudotime_info().
  required_slots = base::c("Branches_n_cell", "Branches_TFE", "Branches_TGE", "Branches")
  for (cell_type in base::names(interest_cell_type_genes_pseudotime_info)) {
    if (!base::all(required_slots %in% base::names(interest_cell_type_genes_pseudotime_info[[cell_type]]))) {
      base::stop("interest_cell_type_genes_pseudotime_info[['", cell_type,
                 "']] is malformed; use the output of get_genes_pseudotime_info().")
    }
  }

  ### Start.
  t_start = base::Sys.time()
  message("Run: Grouping cells based on pseudotime points and valid (non-zero) counts of expression and activity values. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "2.2 Data Processing - Cell Grouping"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Fit n_min from total cell count: y = c*ln(a*x+b).
  # The three calibration points are identical for every branch,
  # so the curve is fitted once per cell type and then evaluated at each branch's cell count;
  # n_min is the rounded-up prediction.
  calculate_n_min = function(Branches_n_cell,
                             x = base::c(200, 1200, 2000),
                             y = base::c(5, 10, 15),
                             xmax = 5000,
                             cell_type) {
    fit_data = base::data.frame(x = x, y = y)

    message("[", cell_type, "] Fitting the n_min model y = c*ln(a*x+b).")
    fit = minpack.lm::nlsLM(
      y ~ coef_c * base::log(coef_a * x + coef_b),
      data = fit_data,
      start = base::list(coef_a = 1, coef_b = 1, coef_c = stats::median(fit_data$y)),
      control = minpack.lm::nls.lm.control(maxiter = 1000)
    )
    coef_a = stats::coef(fit)["coef_a"]
    coef_b = stats::coef(fit)["coef_b"]
    coef_c = stats::coef(fit)["coef_c"]
    message("[", cell_type, "] Fitted parameters: a = ", base::round(coef_a, 9),
            ", b = ", base::round(coef_b, 9), ", c = ", base::round(coef_c, 9), ".")

    # Save one calibration-curve plot per cell type.
    x_grid = base::seq(0, xmax, length.out = 400)
    curve_df = base::data.frame(
      x = x_grid,
      y = base::ceiling(coef_c * base::log(coef_a * x_grid + coef_b))
    )
    x_annotate = 200
    step = 300
    while ((next_val = x_annotate[base::length(x_annotate)] + step) <= xmax) {
      x_annotate = base::c(x_annotate, next_val)
      step = step + 100
    }
    annotate_df = base::data.frame(
      x = x_annotate,
      y = base::ceiling(coef_c * base::log(coef_a * x_annotate + coef_b)),
      label = base::sprintf("(%.0f, %.0f)", base::ceiling(x_annotate),
                            base::ceiling(coef_c * base::log(coef_a * x_annotate + coef_b)))
    )
    p = ggplot2::ggplot(curve_df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_line(color = "blue", linewidth = 1) +
      ggplot2::geom_hline(yintercept = 0, color = "gray50") +
      ggplot2::geom_vline(xintercept = 0, color = "gray50") +
      ggplot2::geom_point(data = annotate_df, ggplot2::aes(x = x, y = y),
                          color = "red", fill = NA, shape = 21, size = 3, stroke = 1.5) +
      ggrepel::geom_text_repel(data = annotate_df, ggplot2::aes(label = label),
                               nudge_x = 50, nudge_y = 0.2 * base::max(curve_df$y),
                               size = 3.5, segment.color = "gray60",
                               box.padding = 0.5, max.overlaps = Inf) +
      ggplot2::labs(
        title = "y = ceiling(c*ln(a*x+b))",
        subtitle = base::paste0("a = ", base::format(base::round(coef_a, 4), nsmall = 4),
                                ", b = ", base::format(base::round(coef_b, 4), nsmall = 4),
                                ", c = ", base::format(base::round(coef_c, 4), nsmall = 4)),
        x = "Number of cells in the branch",
        y = "Minimum number of cells per group in the branch") +
      ggplot2::scale_x_continuous(labels = function(x) base::format(x, big.mark = ",")) +
      ggplot2::scale_y_continuous(labels = function(y) base::format(base::round(y, 1), nsmall = 1)) +
      ggplot2::theme_minimal(base_family = "sans") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "gray40"),
        axis.text = ggplot2::element_text(size = 10),
        panel.grid.minor = ggplot2::element_blank()
      )
    ggplot2::ggsave(base::paste0(cell_type, " --- y = ceiling(c*ln(ax+b)).png"),
                    plot = p, width = 6, height = 5)

    # Predict n_min per branch and clamp it to at least one cell.
    base::vapply(Branches_n_cell, function(n_i) {
      base::max(1L, base::as.integer(base::ceiling(coef_c * base::log(coef_a * n_i + coef_b) + 0.01)))
    }, base::integer(1))
  }

  ### Defining a function for individual cell type processing.
  cell_grouping_process = function(
    cell_type,
    Branches_n_cell,
    Branches_TFE, Branches_TGE,
    Branches,
    calculate_n_min,
    points_x_for_fitting_nc_and_nmin, points_y_for_fitting_nc_and_nmin
  ) {
    # n_min for every branch.
    n_min = calculate_n_min(
      Branches_n_cell = Branches_n_cell,
      x = points_x_for_fitting_nc_and_nmin,
      y = points_y_for_fitting_nc_and_nmin,
      xmax = base::ceiling(base::max(base::unlist(Branches_n_cell) * 1.5)),
      cell_type = cell_type
    )
    message("[", cell_type, "] n_min per branch: ", base::paste(n_min, collapse = ", "), ".")

    # Preallocate the maximum possible number of groups (each valid group holds >= n_min cells).
    Branches_n_k = base::list() # Cell count of each group per branch.
    for (i in base::seq_along(Branches)) {
      Branches_n_k[[i]] = base::numeric(base::ceiling(Branches_n_cell[[i]] / n_min[i]))
    }
    Branches_n_group = base::list() # Final number of groups.
    Branches_k = base::list()       # Current group index k.
    Branches_n_r = base::list()     # Remaining ungrouped cell count.

    for (i in base::seq_along(Branches)) {
      message("[", cell_type, "] Branch ", i, ": grouping cells (n_cell = ",
              Branches_n_cell[[i]], ", n_min = ", n_min[i], ").")

      # n_min above the total makes grouping impossible: return 0 groups for this branch.
      if (n_min[i] > Branches_n_cell[[i]]) {
        message("[", cell_type, "] Branch ", i, ": n_min exceeds total cells; 0 groups.")
        Branches_n_group[[i]] = 0
        Branches_n_r[[i]] = Branches_n_cell[[i]]
        Branches_k[[i]] = 0
        Branches_n_k[[i]] = 0
        next
      }

      # State machine: 1 init, 2 assign, 3 check, 4 prepare next group.
      current_step = 1
      repeat {
        if (current_step == 1) {
          # Step 1: start group k = 1 with all n_cell cells ungrouped (n_r = n_cell).
          Branches_k[[i]] = 1 # Preparing to partition group k=1.
          Branches_n_r[[i]] = Branches_n_cell[[i]] # Currently there are n_r = n_cell ungrouped cells (all cells remain ungrouped).
          current_step = 2 # Proceeding to Step 2.

        } else if (current_step == 2) {
          # Step 2: pre-assign n_min cells to group k.
          Branches_n_k[[i]][Branches_k[[i]]] = n_min[i] # Assigning n_k = n_min cells to group k
          Branches_n_r[[i]] = Branches_n_r[[i]] - Branches_n_k[[i]][Branches_k[[i]]]  # Currently there are n_r = n_r - n_k ungrouped cells remaining.
          current_step = 3 # Proceeding to Step 3.

        } else if (current_step == 3) {
          # Step 3: every gene in TFE and TGE must have >= n_min non-zero values within the current group's pseudotime-ordered cell window.
          # If not, enlarge the group;
          # if cells are exhausted, merge into the previous group, or fail with 0 groups for the first group.
          group_start = base::sum(Branches_n_k[[i]]) - Branches_n_k[[i]][Branches_k[[i]]] + 1
          group_end = base::sum(Branches_n_k[[i]])
          tfe_ok = base::all(base::rowSums(Branches_TFE[[i]][, group_start:group_end, drop = FALSE] != 0) >= n_min[i])
          tge_ok = base::all(base::rowSums(Branches_TGE[[i]][, group_start:group_end, drop = FALSE] != 0) >= n_min[i])
          if (!tfe_ok || !tge_ok) {

            # No cell left to enlarge the current group.
            if (Branches_n_r[[i]] == 0) {
              if (Branches_k[[i]] == 1) {
                # The first group already fails: grouping impossible, return 0 groups.
                Branches_n_group[[i]] = 0
                Branches_n_r[[i]] = Branches_n_cell[[i]]
                Branches_k[[i]] = 0
                Branches_n_k[[i]] = 0
                break

              } else {
                # Merge the unfinished group into the previous group, yielding k-1 groups.
                 Branches_n_k[[i]][Branches_k[[i]] - 1] =
                  Branches_n_k[[i]][Branches_k[[i]] - 1] + Branches_n_k[[i]][Branches_k[[i]]]
                Branches_k[[i]] = Branches_k[[i]] - 1
                Branches_n_r[[i]] = 0
                while (base::length(Branches_n_k[[i]]) > 0 &&
                       Branches_n_k[[i]][base::length(Branches_n_k[[i]])] == 0) {
                  Branches_n_k[[i]] = Branches_n_k[[i]][-base::length(Branches_n_k[[i]])] # Drop preallocated trailing zeros.
                }
                Branches_n_group[[i]] = Branches_k[[i]]
                Branches_n_k[[i]] = Branches_n_k[[i]][-base::length(Branches_n_k[[i]])] # Drop the merged group's slot.
                break
              }
            }

            # Enlarge the current group by one cell and re-evaluate.
            if (Branches_n_r[[i]] > 0) {
              Branches_n_k[[i]][Branches_k[[i]]] = Branches_n_k[[i]][Branches_k[[i]]] + 1
              Branches_n_r[[i]] = Branches_n_r[[i]] - 1
              current_step = 3 # Proceeding to Step 3.
            }

            # Defensive guard; should never be reached under normal flow.
            if (Branches_n_r[[i]] < 0) {break}

          } else {
            # The current group meets the threshold.
            if (Branches_n_r[[i]] == 0) {
              # Cells exhausted exactly when this group completed.
              while (base::length(Branches_n_k[[i]]) > 0 &&
                     Branches_n_k[[i]][base::length(Branches_n_k[[i]])] == 0) {
                Branches_n_k[[i]] = Branches_n_k[[i]][-base::length(Branches_n_k[[i]])] # Drop preallocated trailing zeros.
              }
              Branches_n_group[[i]] = Branches_k[[i]]
              Branches_n_r[[i]] = 0
              break

            } else {
              # Start the next group.
              Branches_k[[i]] = Branches_k[[i]] + 1
              current_step = 4 # Proceeding to Step 4.
            }
          }

        } else if (current_step == 4) {
          # Step 4: fewer than n_min cells cannot open a new group; merge them into the previous group.
         if (Branches_n_r[[i]] < n_min[i]) {
           # If n_r < n_min, it indicates insufficient remaining cells to form a new group.
           # Merge these cells into group k-1, resulting in k-1 total groups, and terminate grouping.
            Branches_n_k[[i]][Branches_k[[i]]] = Branches_n_r[[i]]
            Branches_n_k[[i]][Branches_k[[i]] - 1] =
              Branches_n_k[[i]][Branches_k[[i]] - 1] + Branches_n_k[[i]][Branches_k[[i]]]
            Branches_k[[i]] = Branches_k[[i]] - 1
            Branches_n_r[[i]] = 0
            while (base::length(Branches_n_k[[i]]) > 0 &&
                   Branches_n_k[[i]][base::length(Branches_n_k[[i]])] == 0) {
              Branches_n_k[[i]] = Branches_n_k[[i]][-base::length(Branches_n_k[[i]])] # Drop preallocated trailing zeros.
            }
            Branches_n_k[[i]] = Branches_n_k[[i]][-base::length(Branches_n_k[[i]])] # Drop the unfinished new-group slot.
            Branches_n_group[[i]] = Branches_k[[i]] # Recording final group count.
            break

          } else {
            # Sufficient remaining cells are available to form a new group.
            current_step = 2 # Proceeding to Step 2.
          }

        }

      }
      message("[", cell_type, "] Branch ", i, ": ", Branches_n_group[[i]], " group(s).")
    }

    # Mean of the non-zero values of one gene over the segments defined by n_k (NA if all zero).
    calculate_segment_means = function(gene_value_vector, n_k_vector) {
      gene_value_vector = base::as.numeric(gene_value_vector)
      n_k_vector = base::as.numeric(n_k_vector)
      m = base::length(n_k_vector)
      starts = base::cumsum(base::c(1, n_k_vector[-m]))
      ends = base::cumsum(n_k_vector)
      base::vapply(base::seq_len(m), function(k) {
        segment = gene_value_vector[starts[k]:ends[k]]
        segment = segment[segment != 0]
        if (base::length(segment) == 0) {return(NA_real_)}
        base::mean(segment)
      }, base::numeric(1))
    }

    # Apply the segmented non-zero mean to every gene (row) of a branch matrix.
    apply_segment_means = function(X, n_k_vector) {
      X = base::as.data.frame(X)
      n_k_vector = base::as.numeric(n_k_vector)
      m = base::length(n_k_vector)
      if (m == 0 || base::all(n_k_vector == 0)) {
        return(X[, base::integer(0), drop = FALSE])
      }
      Y = base::as.data.frame(base::t(base::apply(X, 1, function(row) {
        calculate_segment_means(base::as.numeric(row), n_k_vector)
      })))
      base::colnames(Y) = base::as.character(base::seq_len(base::ncol(Y)))
      return(Y)
    }

    ### Build pseudotime endpoints and aggregated TFE/TGE profiles.
    Branches_T_k = base::list()
    Branches_Tslot_k = base::list()
    Branches_TFE_T = base::list()
    Branches_TGE_T = base::list()
    for (i in base::seq_along(Branches)) {
      # Failed branch: return empty aggregates but keep the branch entry.
      if (Branches_n_group[[i]] == 0) {
        Branches_T_k[[i]] = base::numeric(0)
        Branches_Tslot_k[[i]] = base::numeric(0)
        Branches_TFE_T[[i]] = Branches_TFE[[i]][, base::integer(0), drop = FALSE]
        Branches_TGE_T[[i]] = Branches_TGE[[i]][, base::integer(0), drop = FALSE]
        next
      }

      # Integrity checks for successfully grouped branches.
      if (base::length(Branches_n_k[[i]]) != Branches_n_group[[i]]) {
        base::stop("[", cell_type, "] Branch ", i, ": invalid group count record.")
      }
      if (base::sum(Branches_n_k[[i]]) != Branches_n_cell[[i]]) {
        base::stop("[", cell_type, "] Branch ", i, ": grouped cells do not sum to the total cell count.")
      }
      if (Branches_n_r[[i]] != 0) {
        base::stop("[", cell_type, "] Branch ", i, ": ungrouped cells remain.")
      }

      # Pseudotime at the last cell of each group, and each group's interval length.
      end_labels = base::colnames(Branches_TFE[[i]])[base::cumsum(Branches_n_k[[i]])]
      Branches_T_k[[i]] = Branches[[i]]$Pseudotime[
        base::match(end_labels, Branches[[i]]$`k~`)
      ]
      Branches_Tslot_k[[i]] = base::diff(base::c(0, base::as.numeric(Branches_T_k[[i]])))

      # Aggregated non-zero mean expression per group.
      Branches_TFE_T[[i]] = apply_segment_means(Branches_TFE[[i]], Branches_n_k[[i]])
      Branches_TGE_T[[i]] = apply_segment_means(Branches_TGE[[i]], Branches_n_k[[i]])
    }

    return(base::list("Branches_n_group" = Branches_n_group,
                      "Branches_n_k" = Branches_n_k,
                      "n_min" = n_min,
                      "Branches_T_k" = Branches_T_k,
                      "Branches_Tslot_k" = Branches_Tslot_k,
                      "Branches_TFE_T" = Branches_TFE_T,
                      "Branches_TGE_T" = Branches_TGE_T))
  }

  ### Invoking cell grouping function in parallel.
  interest_cell_type_group = parallel::mclapply(
    X = base::names(interest_cell_type_genes_pseudotime_info),
    FUN = function(cell_type) {
      message("Preparing for cell grouping across all branches in cell type ", cell_type, ".")
      cell_grouping_process(
        cell_type = cell_type,
        Branches_n_cell = interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_n_cell"]],
        Branches_TFE = interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]],
        Branches_TGE = interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGE"]],
        Branches = interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches"]],
        calculate_n_min = calculate_n_min,
        points_x_for_fitting_nc_and_nmin = points_x_for_fitting_nc_and_nmin,
        points_y_for_fitting_nc_and_nmin = points_y_for_fitting_nc_and_nmin
      )
    },
    mc.cores = ncores
  )
  base::names(interest_cell_type_group) = base::names(interest_cell_type_genes_pseudotime_info)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Grouping cells based on pseudotime points and valid (non-zero) counts of expression and activity values. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_group)
}
