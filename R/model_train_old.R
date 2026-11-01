#' @title Train scDARI Gene Regulatory Network Models
#'
#' @description Trains a model for every TF-TG pair across cell types and branches
#' with a hybrid strategy combining genetic-algorithm global search and gradient-ascent local refinement.
#' Multiple early stopping criteria are supported, and regulons can be trained in parallel.
#'
#' @param interest_cell_type_branch_training_set List.
#' Output of \code{\link{select_training_set}},
#' containing branch-level training data with expression and activity features.
#' @param initial_params Data frame. Output of \code{\link{initialize_params}},
#' containing the initialization, upper bounds and lower bounds of the model parameters.
#' @param interest_cell_type_group List. Output of \code{\link{cell_grouping}},
#' containing pseudotime interval lengths for each branch.
#' @param interest_cell_type_tau Named numeric vector.
#' Output of \code{\link{get_cell_type_tau}},
#' the cell-type-specific binarization thresholds.
#' @param interest_cell_type_genes_pseudotime_info List.
#' Output of \code{\link{get_genes_pseudotime_info}},
#' containing branch cell counts used to weight the integrated results.
#' @param max_epochs Positive integer.
#' Maximum number of training epochs.
#' Defaults to 1000.
#' @param early_stop Positive integer.
#' Number of consecutive epochs without loss improvement that triggers early stopping.
#' Defaults to 10.
#' @param eps_theta Tiny positive numeric.
#' Threshold at which the predicted regulatory strength is considered converged to the observed strength.
#' Defaults to 1e-3.
#' @param eps_loss Tiny positive numeric.
#' Threshold for negligible loss improvement.
#' Defaults to 1e-5.
#' @param popSize_ga Integer. Population size of the genetic algorithm.
#' Defaults to 50.
#' See \code{popSize} in \code{\link[GA]{ga}}.
#' @param maxiter_ga Positive integer.
#' Maximum generations per genetic algorithm run.
#' Defaults to 150.
#' See \code{maxiter} in \code{\link[GA]{ga}}.
#' @param pcrossover_ga Numeric between 0 and 1. Crossover probability.
#' Defaults to 0.9.
#' See \code{pcrossover} in \code{\link[GA]{ga}}.
#' @param pmutation_ga Numeric between 0 and 1. Mutation probability.
#' Defaults to 0.5.
#' See \code{pmutation} in \code{\link[GA]{ga}}.
#' @param parallel_ga Logical. Whether the genetic algorithm runs in parallel.
#' Defaults to FALSE.
#' See \code{parallel} in \code{\link[GA]{ga}}.
#' @param seed_ga Integer.
#' Random seed for reproducible genetic algorithm results.
#' Defaults to 123.
#' See \code{seed} in \code{\link[GA]{ga}}.
#' @param iterations_grad Positive integer.
#' Number of gradient ascent steps per epoch.
#' Defaults to 30.
#' See \code{iterations} in \code{\link{gradient_ascent_train_old}}
#' @param ncores_grad Positive integer.
#' Cores for parallel gradient evaluation within gradient ascent.
#' Defaults to 1.
#' See \code{ncores} in \code{\link{gradient_ascent_train_old}}
#' @param ncores Positive integer.
#' Cores for regulon-level parallel evaluation via \code{\link[parallel]{mclapply}}.
#' Defaults to 1; set to 1 on Windows.
#'
#' @return Nested list organized by cell type, each entry containing:
#' \itemize{
#'   \item \code{best_pred_result}: integrated data frame across branches, with
#'     predicted regulatory strength, binarized result, fitting loss and
#'     training metadata.
#'   \item \code{branch*}: detailed results for every regulon in the branch,
#'     including training data, pseudotime intervals, threshold, objective
#'     values, parameter trajectories, loss history, best epoch, best
#'     parameters, best prediction, stop cause, training time and total epochs.
#' }
#'
#' @details For each cell type and branch, every TF-TG pair is trained
#' independently:
#' \enumerate{
#'   \item The genetic algorithm population is initialized from the shared parameter template.
#'   \item Each epoch runs two stages:
#'     \itemize{
#'       \item Stage 1: genetic algorithm for global exploration,
#'         initialized with the previous epoch's best solution.
#'       \item Stage 2: gradient ascent for local refinement from the best GA solution.
#'     }
#'   \item After each epoch, stopping conditions are checked in priority order:
#'     maximum epochs reached, no loss improvement for consecutive epochs,
#'     predicted strength converged to observed strength, negligible loss improvement.
#'   \item Results are integrated across branches: binary relationships use union logic,
#'     and continuous values are weighted by branch cell counts.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates and switches to the \code{4 Build Prediction Model} directory,
#'     and restores the original working directory on exit.
#'   \item Regulons are trained independently,
#'     enabling regulon-level parallel computing through \code{ncores}.
#'   \item When \code{parallel_ga = TRUE},
#'     an additional serial GA iteration is run after the parallel GA to reduce result instability.
#'   \item Larger population,
#'     iteration count or mutation probability can help escape local optima at the cost of runtime.
#' }
#'
#' @seealso \code{\link{f_train_cal_old}}, \code{\link{loss_cal_old}}, \code{\link{theta_p_cal_old}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_training_set = base::readRDS(
#'   "./4 Build Prediction Model/interest_cell_type_branch_training_set.rds"
#' )
#' initial_params = base::readRDS(
#'   "./4 Build Prediction Model/initial_params.rds"
#' )
#' interest_cell_type_group = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds"
#' )
#' interest_cell_type_tau = base::readRDS("./3 get iGRN/interest_cell_type_tau.rds")
#' interest_cell_type_genes_pseudotime_info = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds"
#' )
#'
#' max_epochs = 50
#' early_stop = 5
#'
#' popSize_ga = 50
#' maxiter_ga = 150
#' pcrossover_ga = 0.9
#' pmutation_ga = 0.5
#' parallel_ga = FALSE
#' seed_ga = 123 # Reproducible
#' # seed_ga = stats::runif(1) # Not reproducible
#'
#' iterations_grad = 30
#' ncores_grad = 1
#'
#' ncores = parallel::detectCores() - 1 # Linux
#' # ncores = 1 # Windows
#'
#' interest_cell_type_branch_model_train = model_train_old(
#'   interest_cell_type_branch_training_set = interest_cell_type_branch_training_set,
#'   initial_params = initial_params,
#'   interest_cell_type_group = interest_cell_type_group,
#'   interest_cell_type_tau = interest_cell_type_tau,
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   max_epochs = max_epochs,
#'   early_stop = early_stop,
#'   popSize_ga = popSize_ga,
#'   maxiter_ga = maxiter_ga,
#'   pcrossover_ga = pcrossover_ga,
#'   pmutation_ga = pmutation_ga,
#'   parallel_ga = parallel_ga,
#'   seed_ga = seed_ga,
#'   iterations_grad = iterations_grad,
#'   ncores_grad = ncores_grad,
#'   ncores = ncores
#' )
#' base::saveRDS(
#'   interest_cell_type_branch_model_train,
#'   file = "./4 Build Prediction Model/interest_cell_type_branch_model_train.rds"
#' )
#' }
model_train_old = function(
    interest_cell_type_branch_training_set,
    initial_params,
    interest_cell_type_group,
    interest_cell_type_tau,
    interest_cell_type_genes_pseudotime_info,
    max_epochs = 1000, early_stop = 10,
    eps_theta = 1e-3, eps_loss = 1e-5,
    popSize_ga = 50, maxiter_ga = 150, pcrossover_ga = 0.9, pmutation_ga = 0.5,
    parallel_ga = FALSE, seed_ga = 123,
    iterations_grad = 30, ncores_grad = 1,
    ncores = 1
) {
  ### Start.
  t_start = base::Sys.time()
  message("Run: Training models for each branch. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "4 Build Prediction Model"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {
    message("Folder already exists: ", new_folder, ".")
  }
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Process.
  interest_cell_type_branch_model_train = base::list()
  for (cell_type in base::names(interest_cell_type_branch_training_set)) {
    message("Starting to train prediction models for all branches in cell type: ", cell_type, ".")
    interest_cell_type_branch_model_train[[cell_type]] = base::list()

    # Exclude the "all" element when counting branches.
    n_branch = base::length(interest_cell_type_branch_training_set[[cell_type]]) - 1
    best_pred_result = interest_cell_type_branch_training_set[[cell_type]][["all"]]
    best_pred_result["theta_p"] = 0
    best_pred_result[base::paste0("theta_p_branch", 1:n_branch)] = 0
    best_pred_result$theta_i_bin = base::as.integer(best_pred_result$theta_i > 0)
    best_pred_result["theta_p_bin"] = 0
    best_pred_result[base::paste0("theta_p_bin_branch", 1:n_branch)] = 0
    best_pred_result["fit_loss"] = Inf
    best_pred_result[base::paste0("fit_loss_branch", 1:n_branch)] = Inf
    interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]] = best_pred_result

    for (n in 1:n_branch) {
      message("Starting training for branch ", n, " of cell type ", cell_type, ".")
      branch_train_data = interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]]
      message("There are ", base::nrow(branch_train_data), " regulons for training.")
      interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]] = base::list()

      # Train one regulon.
      model_train_process = function(re,
                                     n,
                                     n_re = base::match(re, base::rownames(branch_train_data)),
                                     N_re = base::nrow(branch_train_data),
                                     N_branch = n_branch,
                                     initial_params,
                                     train_data = base::as.matrix(branch_train_data[re, , drop = FALSE]),
                                     Tslot_K = interest_cell_type_group[[cell_type]][["Branches_Tslot_k"]][[n]],
                                     tau = base::as.numeric(interest_cell_type_tau[cell_type]),
                                     popSize_ga,
                                     max_epochs,
                                     pcrossover_ga,
                                     pmutation_ga,
                                     maxiter_ga,
                                     parallel_ga,
                                     seed_ga,
                                     best_pred_result_re = best_pred_result[re, , drop = FALSE],
                                     early_stop) {
        result = base::list()
        result[["re"]] = re
        signal = base::paste0(n_re, "/", N_re)
        message("Progress ", signal, "...")
        Signal = base::paste0("[Branch ", n, "/", N_branch, " - Train Regulon ", signal, "] ")
        message(Signal, "Extracting initial parameter values, upper bounds, lower bounds and names.")

        params = initial_params$initialization
        params_init = initial_params$initialization
        params_upper = initial_params$upper
        params_lower = initial_params$lower
        params_names = base::rownames(initial_params)
        base::names(params) = params_names
        base::names(params_init) = params_names
        base::names(params_upper) = params_names
        base::names(params_lower) = params_names
        message(Signal, "Extracting the training set.")

        result[["train_data"]] = train_data
        message(Signal, "Extracting pseudotime interval lengths for all groups.")

        result[["Tslot_K"]] = Tslot_K
        message(Signal, "Extracting the binarization threshold for regulatory strength.")

        result[["tau"]] = tau
        message(Signal, "Initializing the genetic algorithm population.")

        prev_population = base::matrix(
          data = base::rep(base::as.numeric(params), popSize_ga),
          nrow = popSize_ga, ncol = base::length(params), byrow = TRUE
        )
        message(Signal, "Initializing parameter storage for each epoch.")

        para_data_frame = base::as.data.frame(base::matrix(
          data = 0, nrow = max_epochs, ncol = base::length(params),
          dimnames = base::list(1:max_epochs, params_names)
        ))
        message(Signal, "Initializing loss storage for each epoch.")

        loss_train = base::rep(0, max_epochs)
        base::names(loss_train) = 1:max_epochs
        f_value_train_vec = base::rep(0, max_epochs)
        base::names(f_value_train_vec) = 1:max_epochs

        t1_re = base::Sys.time()
        for (epoch in 1:max_epochs) {
          # Each epoch: GA global search, then gradient ascent refinement,
          # and the refined parameters seed the next epoch's GA.
          SIGNAL = base::paste0(base::gsub("] ", "", base::paste0(Signal, " - Epoch ", epoch), fixed = TRUE), "] ")
          if (epoch %% 10 == 1) message(SIGNAL, "Epoch ", epoch, "/", max_epochs, " in progress.")

          # Genetic algorithm global search.
          result_ga = GA::ga(
            type = "real-valued", fitness = scDARI::f_train_cal_old, lower = params_lower, upper = params_upper,
            popSize = popSize_ga, pcrossover = pcrossover_ga, pmutation = pmutation_ga,
            elitism = base::ceiling(base::min(0.2 * maxiter_ga, popSize_ga / 5)),
            maxiter = maxiter_ga, suggestions = prev_population, parallel = parallel_ga,
            seed = seed_ga, monitor = FALSE,
            params_names = params_names, train_data = train_data, tau = tau, Tslot_K = Tslot_K,
            params_init = params_init, params_upper = params_upper, params_lower = params_lower
          )
          # Parallel GA may be unstable; rebuild the population from its best
          # solution and run one serial GA iteration for stabilization.
          if (parallel_ga) {
            prev_population = base::matrix(
              data = base::rep(base::as.numeric(result_ga@solution[1, ]), popSize_ga),
              nrow = popSize_ga, ncol = base::length(params), byrow = TRUE
            )
            result_ga = GA::ga(
              type = "real-valued", fitness = scDARI::f_train_cal_old, lower = params_lower, upper = params_upper,
              popSize = popSize_ga, pcrossover = pcrossover_ga, pmutation = pmutation_ga,
              elitism = base::ceiling(base::min(0.2 * maxiter_ga, popSize_ga / 5)),
              maxiter = 1, suggestions = prev_population, parallel = FALSE,
              seed = seed_ga, monitor = FALSE,
              params_names = params_names, train_data = train_data, tau = tau, Tslot_K = Tslot_K,
              params_init = params_init, params_upper = params_upper, params_lower = params_lower
            )
            base::gc()
          }

          # Gradient ascent local refinement.
          result_gradient = scDARI::gradient_ascent_train_old(
            iterations = iterations_grad, alpha_lower = 1e-5, alpha_upper = 1e-3, alpha_guess = 1e-4,
            init_theta = base::as.numeric(result_ga@solution[1, ]),
            h_max = 1e-9, h_min_percent = 0.01, ncores = ncores_grad,
            params_lower = params_lower, params_upper = params_upper, monitor = FALSE,
            params_names = params_names, train_data = train_data, tau = tau, Tslot_K = Tslot_K,
            params_init = params_init
          )
          f_value_train_vec[epoch] = -result_gradient$objective_value
          result[["f_value_train"]] = f_value_train_vec
          prev_population = base::matrix(
            data = base::rep(result_gradient$parameters, popSize_ga),
            nrow = popSize_ga, ncol = base::length(params), byrow = TRUE
          )
          para_data_frame[epoch, ] = result_gradient$parameters
          result[["para_data_frame"]] = para_data_frame

          # Pure prediction loss at the refined point.
          loss_train[epoch] = scDARI::loss_cal_old(
            data = train_data,
            params = result_gradient$parameters, params_names = params_names,
            tau = tau, Tslot_K = Tslot_K
          )

          # Stopping conditions, checked in priority order.
          if (epoch >= max_epochs) {
            message(SIGNAL, "Reaching the maximum number of training epochs.")
            loss_vec = loss_train[1:base::ifelse(
              base::length(base::which(loss_train != 0)) > 0, base::max(base::which(loss_train != 0)), 0
            )]
            result[["loss_train"]] = loss_vec
            best_epoch = base::as.numeric(base::names(base::which.min(loss_vec)))
            result[["best_epoch"]] = best_epoch
            best_params = base::as.numeric(para_data_frame[best_epoch, , drop = FALSE])
            base::names(best_params) = params_names
            result[["best_params"]] = best_params
            best_pred_result_re[re, base::paste0("fit_loss_branch", n)] = loss_vec[best_epoch]
            best_pred_result_re[re, base::paste0("theta_p_branch", n)] = scDARI::theta_p_cal_old(
              data = train_data,
              best_params = best_params, params_names = params_names,
              tau = tau, Tslot_K = Tslot_K
            )
            best_pred_result_re[re, base::paste0("theta_p_bin_branch", n)] =
              base::ifelse(best_pred_result_re[re, base::paste0("theta_p_branch", n)] >= tau, 1, 0)
            result[["best_pred_result_re"]] = best_pred_result_re
            stop_cause = "maximum training epochs reached"
            message(SIGNAL, "Meeting the early stopping condition: ", stop_cause, ".")
            result[["stop_cause"]] = stop_cause
            break
          } else if (epoch >= early_stop) {
            loss_vec = loss_train[1:base::ifelse(
              base::length(base::which(loss_train != 0)) > 0, base::max(base::which(loss_train != 0)), 0
            )]
            result[["loss_train"]] = loss_vec
            best_epoch = base::as.numeric(base::names(base::which.min(loss_vec)))
            result[["best_epoch"]] = best_epoch
            best_params = base::as.numeric(para_data_frame[best_epoch, , drop = FALSE])
            base::names(best_params) = params_names
            result[["best_params"]] = best_params
            best_pred_result_re[re, base::paste0("fit_loss_branch", n)] = loss_vec[best_epoch]
            best_pred_result_re[re, base::paste0("theta_p_branch", n)] = scDARI::theta_p_cal_old(
              data = train_data,
              best_params = best_params, params_names = params_names,
              tau = tau, Tslot_K = Tslot_K
            )
            best_pred_result_re[re, base::paste0("theta_p_bin_branch", n)] =
              base::ifelse(best_pred_result_re[re, base::paste0("theta_p_branch", n)] >= tau, 1, 0)
            result[["best_pred_result_re"]] = best_pred_result_re
            if (best_epoch <= base::length(loss_vec) - early_stop + 1) {
              stop_cause = "no loss improvement for consecutive early_stop epochs"
              message(SIGNAL, "Meeting the early stopping condition: ", stop_cause, ".")
              result[["stop_cause"]] = stop_cause
              break
            } else if (base::abs(best_pred_result_re[re, base::paste0("theta_p_branch", n)] -
                                 best_pred_result_re[re, "theta_i"]) < eps_theta) {
              stop_cause = "predicted theta_p converges to observed theta_i"
              message(SIGNAL, "Meeting the early stopping condition: ", stop_cause, ".")
              result[["stop_cause"]] = stop_cause
              break
            } else if (base::max(loss_vec[(base::length(loss_vec) - early_stop + 1):base::length(loss_vec)]) -
                       loss_vec[best_epoch] < eps_loss) {
              stop_cause = "loss improvement is below the threshold"
              message(SIGNAL, "Meeting the early stopping condition: ", stop_cause, ".")
              result[["stop_cause"]] = stop_cause
              break
            }
          }
        }
        t2_re = base::Sys.time()
        message(Signal, "Training completed. Duration: ", base::format(t2_re - t1_re), ".")

        result[["train_time"]] = base::sprintf("%.2f %s", base::as.numeric(t2_re - t1_re),
                                               base::units(t2_re - t1_re))
        result[["train_epoch"]] = epoch
        return(result)
      }

      interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]] = parallel::mclapply(
        X = base::rownames(branch_train_data),
        FUN = function(re) {
          model_train_process(
            re,
            n = n,
            n_re = base::match(re, base::rownames(branch_train_data)),
            N_re = base::nrow(branch_train_data),
            N_branch = n_branch,
            initial_params = initial_params,
            train_data = base::as.matrix(branch_train_data[re, , drop = FALSE]),
            Tslot_K = interest_cell_type_group[[cell_type]][["Branches_Tslot_k"]][[n]],
            tau = base::as.numeric(interest_cell_type_tau[cell_type]),
            popSize_ga = popSize_ga,
            max_epochs = max_epochs,
            pcrossover_ga = pcrossover_ga,
            pmutation_ga = pmutation_ga,
            maxiter_ga = maxiter_ga,
            parallel_ga = parallel_ga,
            seed_ga = seed_ga,
            best_pred_result_re = best_pred_result[re, , drop = FALSE],
            early_stop = early_stop
          )
        },
        mc.cores = ncores
      )
      base::names(interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]]) =
        base::rownames(branch_train_data)

      # Copy per-regulon branch results into the integrated table.
      for (re in base::rownames(branch_train_data)) {
        re_result = interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["best_pred_result_re"]]
        best_pred_result[re, base::paste0("theta_p_branch", n)] = re_result[re, base::paste0("theta_p_branch", n)]
        best_pred_result[re, base::paste0("theta_p_bin_branch", n)] = re_result[re, base::paste0("theta_p_bin_branch", n)]
        best_pred_result[re, base::paste0("fit_loss_branch", n)] = re_result[re, base::paste0("fit_loss_branch", n)]
      }
    }
    message("All branch trainings have been completed.")

    # Integrate branches: binary relationships use union logic, and continuous
    # values are weighted by branch cell counts.
    message("Integrating the results from all branches.")
    if (n_branch > 1) {
      best_pred_result$theta_p_bin = base::as.integer(base::rowSums(
        best_pred_result[, base::grep("^theta_p_bin_branch\\d+$", base::names(best_pred_result))]) > 0
      )
      weights = base::unlist(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_n_cell"]])
      best_pred_result$theta_p = base::apply(
        best_pred_result[, base::grep("^theta_p_branch([0-9]{1,2})$", base::names(best_pred_result), value = TRUE)],
        1, function(x) {
          valid = base::which(x != 0)
          if (base::length(valid) == 0) 0
          else base::sum(x[valid] * weights[valid]) / base::sum(weights[valid])
        })
      best_pred_result$fit_loss = base::apply(
        best_pred_result[, base::grep("^fit_loss_branch([0-9]{1,2})$", base::names(best_pred_result), value = TRUE)],
        1, function(x) {
          valid = base::which(x != Inf)
          if (base::length(valid) == 0) Inf
          else base::sum(x[valid] * weights[valid]) / base::sum(weights[valid])
        })
    } else if (n_branch == 1) {
      best_pred_result$theta_p_bin = best_pred_result$theta_p_bin_branch1
      best_pred_result$theta_p = best_pred_result$theta_p_branch1
      best_pred_result$fit_loss = best_pred_result$fit_loss_branch1
    } else {
      base::stop("The number of branches is less than 1.")
    }

    # Insert the integrated absolute difference after the observed strength.
    best_pred_result = base::cbind(
      best_pred_result[, 1:4],
      base::abs(best_pred_result$theta_i - best_pred_result$theta_p),
      best_pred_result[, 5:base::ncol(best_pred_result)]
    )
    base::colnames(best_pred_result)[5] = "theta_difference"
    best_pred_result[base::paste0("theta_difference_branch", 1:n_branch)] = Inf
    best_pred_result[base::paste0("train_time_branch", 1:n_branch)] = 0
    best_pred_result[base::paste0("train_epoch_branch", 1:n_branch)] = 0
    best_pred_result[base::paste0("best_epoch_branch", 1:n_branch)] = 0
    best_pred_result[base::paste0("stop_cause_branch", 1:n_branch)] = NA
    for (n in 1:n_branch) {
      branch_train_data = interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]]
      for (re in base::rownames(branch_train_data)) {
        re_detail = interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]]
        best_pred_result[re, base::paste0("theta_difference_branch", n)] =
          base::abs(best_pred_result[re, "theta_i"] - best_pred_result[re, base::paste0("theta_p_branch", n)])
        best_pred_result[re, base::paste0("train_time_branch", n)] = re_detail[["train_time"]]
        best_pred_result[re, base::paste0("train_epoch_branch", n)] = re_detail[["train_epoch"]]
        best_pred_result[re, base::paste0("best_epoch_branch", n)] = re_detail[["best_epoch"]]
        best_pred_result[re, base::paste0("stop_cause_branch", n)] = re_detail[["stop_cause"]]
      }
    }
    interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]] = best_pred_result
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Training models for each branch. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_branch_model_train)
}
