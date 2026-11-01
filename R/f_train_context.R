#' @title Build the Fixed Training Context for One Regulon
#'
#' @description Precomputes every quantity that stays constant while one TF-TG
#' pair (regulon) is trained, and returns them in a single context list.
#' Building the context once and reusing it across the many genetic-algorithm
#' and gradient-ascent evaluations avoids repeated column matching, sorting,
#' lag-1 indexing and memory allocation.
#'
#' @param params_names Character vector. Names of the model parameters, in the
#' same order as the elements of a candidate parameter vector. It is matched
#' against the canonical parameter order to build stable positional indices,
#' so every later evaluation can read a parameter by name robustly regardless
#' of vector ordering. Closely related to every candidate \code{params} vector
#' passed to \code{\link{f_train_core}}.
#' @param train_data Data frame or one-row matrix. Training data for one
#' regulon. The first three columns are \code{TF}, \code{TG} and the prior
#' regulatory strength \code{theta_i}; from column 4 onward the columns
#' alternate between TF expression \code{TFE_k} and target expression
#' \code{TGE_k}, one pair per group. This fixed layout is used to extract and
#' lag-1 shift the two expression series. Closely related to \code{Tslot_K},
#' which gives the interval lengths for the same groups.
#' @param tau Numeric strictly between 0 and 1. Cell-type-specific
#' binarization threshold for regulatory strength. It is stored unchanged in
#' the context and later passed to \code{\link{R_cal}}. Closely related to the
#' threshold used in \code{\link{model_train}}.
#' @param Tslot_K Positive numeric vector with one entry per group.
#' Pseudotime interval lengths between consecutive groups. Its first
#' \code{n_group - 1} entries become the lag-1 step lengths \code{Tslot}, and
#' its type-7 75th percentile \code{q75} is the reference pseudotime step used
#' by the discrete-stability penalty. Closely related to the
#' \code{TFE_k}/\code{TGE_k} columns in \code{train_data}.
#'
#' @return A list holding the fixed context:
#' \itemize{
#'   \item \code{pos}: integer vector mapping the canonical parameter order to
#'     positions in \code{params_names}.
#'   \item \code{Length}: number of groups, equal to the length of each
#'     expression series.
#'   \item \code{TFE_cur}, \code{TFE_next}, \code{TGE_cur}, \code{TGE_next}:
#'     lag-1 shifted expression series, each of length \code{n_group - 1}.
#'   \item \code{Tslot}: lag-1 interval lengths of length \code{n_group - 1}.
#'   \item \code{q75}: type-7 75th percentile of \code{Tslot_K}.
#'   \item \code{theta}: prior regulatory strength \code{theta_i}.
#'   \item \code{tau}: the binarization threshold.
#' }
#'
#' @details The context is assembled in a single pass:
#' \enumerate{
#'   \item Matches the canonical parameter order against \code{params_names}
#'     to obtain stable positional indices.
#'   \item Extracts the TFE and TGE group series from the fixed alternating
#'     column layout and pre-shifts them into current and next lag-1 vectors.
#'   \item Subsets the lag-1 interval lengths and computes the 75th-percentile
#'     pseudotime step with a nested type-7 routine that is numerically
#'     identical to \code{stats::quantile(x, probs = 0.75)}.
#'   \item Stores the prior strength and threshold, which are also constant.
#' }
#'
#' @note
#' \itemize{
#'   \item This function performs no optimization; it only organizes constant
#'     data, so rebuilding it more than once per regulon is wasted work.
#'   \item The percentile is computed by the nested \code{q75_fast} helper
#'     using the same type-7 interpolation as \code{\link[stats]{quantile}}.
#'   \item The fixed column layout is assumed to match the output of
#'     \code{\link{select_training_set}}.
#' }
#'
#' @seealso Called by \code{\link{f_train_cal}}, \code{\link{f_train_cal_gr}},
#'   \code{\link{gradient_ascent_train}} and \code{\link{model_train}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ctx = f_train_context(
#'   params_names = params_names,
#'   train_data = train_data,
#'   tau = tau,
#'   Tslot_K = Tslot_K
#' )
#' }
f_train_context = function(params_names, train_data, tau, Tslot_K) {
  # Fixed parameter order used to map the named vector to positional indices.
  .scdari_canon = c("r2", "alpha1", "alpha2", "beta1", "beta2", "s",
                    "u1i", "u1j", "u2i", "u2j", "v_tilde")

  pos = base::match(.scdari_canon, params_names)

  # Package layout: TF, TG, theta_i, then alternating TFE_k / TGE_k columns.
  nc = base::ncol(train_data)
  tfe_cols = base::seq(4L, nc, 2L)
  tge_cols = base::seq(5L, nc, 2L)
  TFE = base::as.numeric(train_data[tfe_cols])
  TGE = base::as.numeric(train_data[tge_cols])

  Length = base::length(TFE)
  K_1 = base::seq_len(Length - 1L)

  # Pre-shifted series so the hot loop never repeats lag-1 indexing.
  TFE_cur = TFE[K_1]; TFE_next = TFE[K_1 + 1L]
  TGE_cur = TGE[K_1]; TGE_next = TGE[K_1 + 1L]

  # Internal fast type-7 75th percentile, equivalent to
  # stats::quantile(x, probs = 0.75) with the default settings.
  q75_fast = function(x) {
    xs = base::sort(x)
    n = base::length(xs)
    h = (n - 1) * 0.75 + 1
    il = base::floor(h)
    if (il == h) return(xs[il])
    xs[il] + (h - il) * (xs[il + 1L] - xs[il])
  }

  list(
    pos = pos,
    Length = Length,
    TFE_cur = TFE_cur, TFE_next = TFE_next,
    TGE_cur = TGE_cur, TGE_next = TGE_next,
    Tslot = Tslot_K[K_1],
    q75 = q75_fast(Tslot_K),
    theta = base::as.numeric(train_data[, "theta_i"]),
    tau = tau
  )
}
