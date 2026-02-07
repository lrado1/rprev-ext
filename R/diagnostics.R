# This function compares the actual variance of the yearly incidence rates with
# rates simulated from a Poisson process with overall rate equal to the overall
# mean rate.
#
# param inc Vector of absolute incidence values for each included year.
# param N_sim Number of simulations to perform.
# return Vector of p-values for over- and under-dispersion based on the
#   position of the observed sequence variance in the distribution.
test_dispersion <- function(inc, N_sim = 1e5) {
    var_sim <- vapply(seq(N_sim), function(i) var(rpois(length(inc), mean(inc))), numeric(1))
    c(sum(var_sim > var(inc))/N_sim, sum(var_sim <= var(inc))/N_sim)
}

#' Test simulated prevalence fit.
#'
#' Calculates a Chi squared test between predicted yearly contributions to
#' prevalence, and the observed values obtained from the registry, indicating
#' whether the simulated prevalence values are accurate.
#'
#' @param object A \code{prevalence} object.
#' @return P-value from a chi-squared test of difference between prevalence
#'   prediction and counted prevalence at the index date.
#' @examples
#' data(prevsim)
#'
#' \dontrun{
#'
#' obj <- prevalence(Surv(time, status) ~ age(age) + sex(sex) + entry(entrydate) + event(eventdate),
#'                   data=prevsim, num_years_to_estimate = c(5, 10), population_size=1e6,
#'                   start = "2005-09-01",
#'                   num_reg_years = 8, cure = 5)
#'
#' test_prevalence_fit(obj)
#' }
#' @import data.table
#'
#' @export
#' @family prevalence functions
test_prevalence_fit <- function(object) {
    # Needed for CRAN
    incident_date <- NULL
    sim <- NULL

    sim_dt <- object$simulated
    idx_dates <- object$index_dates
    if (is.null(idx_dates) || length(idx_dates) == 0) {
        stop("Error: prevalence object must contain non-empty 'index_dates'.")
    }
    if (is.null(sim_dt) || is.null(object$counted)) {
        return(NA_real_)
    }

    # Per-index p-value(s), for both single- and multi-index objects
    pvals <- vapply(seq_along(idx_dates), function(k) {
        alive_col <- sprintf("alive_k%03d", k)
        if (!alive_col %in% names(sim_dt)) {
            stop("Error: simulated prevalence data must contain column '", alive_col, "'.")
        }
        contribs <- sim_dt[, sum((incident_date >= object$registry_start) & (get(alive_col) == 1)), by=sim][[2]]
        predicted <- round(mean(contribs))
        counted_k <- if (!is.null(names(object$counted))) object$counted[[as.character(idx_dates[k])]] else object$counted[k]
        poisson.test(c(counted_k, predicted))$p.value
    }, numeric(1))
    names(pvals) <- as.character(idx_dates)
    if (length(pvals) == 1) unname(pvals) else pvals
}
