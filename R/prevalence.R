INCIDENCE_MARGIN <- 1.5
MIN_INCIDENCE <- 10

#' Estimate point prevalence at an index date.
#'
#' Point prevalence at a specific index date is estimated using contributions to
#' prevalence from both available registry data, and from Monte Carlo
#' simulations of the incidence and survival process, as outlined by Crouch et
#' al (2004) (see References).
#'
#' The most important parameter is \code{num_years_to_estimate}, which governs
#' the number of previous years of data to use when estimating the prevalence at
#' the index date. If this parameter is greater than the number of years of
#' known incident cases available in the supplied registry data (specified with
#' argument \code{num_registry_years}), then the remaining
#' \code{num_years_to_estimate - num_registry_years} years of incident data will
#' be simulated using Monte Carlo simulation.
#'
#' The larger \code{num_years_to_estimate}, the more accurate the prevalence
#' estimate will be, provided an adequate survival model can be fitted to the
#' registry data. It is therefore important to provide as much clean registry
#' data as possible.
#'
#' Prevalence arises from two stochastic processes: incidence and survival.
#' This is reflected in the function arguments by multiple options for
#' each of these processes.
#'
#' The incidence process is specified by an object
#' that has an associated \code{draw_incident_population} method, which produces the new
#' incident population. The default implementation is a homogeneous Poisson process,
#' whereby interarrival times are distributed according to an exponential distribution.
#' The \code{inc_formula} argument specifies the nature of this process, see the
#' description for more details. See the vignette for guidance on providing a custom incidence
#' object.
#'
#' The survival process is characterised by a method \code{predict_survival_probability},
#' that estimates the probability of a given individual being alive at the index date.
#' The default object is a parametric distribution with the functional form being specified
#' in \code{surv_formula} and distribution given in \code{dist}. See the vignette for guidance
#' on providing a custom survival model.
#'
#' @param index Deprecated. Use \code{index_dates} instead.
#' @param index_dates The date(s) at which to estimate point prevalence as a string in the format
#' YYYY-MM-DD.
#' @param num_years_to_estimate Number of years of data to consider when
#'   estimating point prevalence; multiple values can be specified in a vector.
#'   If any values are greater than the number of years of registry data
#'   available before the requested \code{index_dates}, incident cases
#'   for the difference will be simulated.
#' @param data A data frame with the corresponding column names provided in
#'   \code{form}.
#' @param inc_formula A formula specifying the columns used in the incidence process.
#'     The LHS should be the name of the column holding the incident dates,
#'     with the RHS specifying any variables that should be stratified by, or 1 if no
#'     stratification. For example, with the supplied \code{prevsim} data set, it could
#'     be used as follows:
#'
#'     \code{entrydate ~ 1} for a non-stratified process.
#'     \code{entrydate ~ sex} for a process that will stratify incidence by sex.
#'
#' @param inc_model An object that has a \code{draw_incident_population}
#'     method. See the vignette for further guidance.
#' @param surv_formula A formula used to specify a survival model, where the
#' LHS a Surv object, as used by \code{flexsurvreg}.
#' @param dist The distribution used by the default parametric survival model.
#' @param surv_model An object that has a \code{predict_survival_probability}
#'     method. See the vignette for further guidance.
#' @param registry_start_date The starting date of the registry. If not supplied
#'   then defaults to the earliest incidence date in the supplied data set.
#' @param death_column A string providing the name of the column which holds the death
#'     date information. If not provided then prevalence cannot be counted and estimates
#'     will be solely derived from simulation.
#' @param incident_column A string providing the name of the column which holds the diagnosis
#'     date. If not provided either in this argument or in \code{inc_formula},
#'     then prevalence cannot be counted and estimates will be solely derived from simulation.
#' @param age_column A string providing the name of the column that holds patient age. If provided
#'     then patients alive at \code{age_dead} are set to die. This helps combat 'immortal' patients.
#' @param age_dead The age at which patients are set to be dead if they are still alive, to prevent
#'     'immortal' patients. Used in conjunction with \code{age_column}.
#' @param status_column A string providing the name of the column that holds patient event status at
#'     the event time. If not provided in \code{surv_formula} or in this argument then prevalence
#'     cannot be counted.
#' @param N_boot Number of bootstrapped calculations to perform.
#' @param population_size Integer corresponding to the size of the population at
#'   risk.
#' @param proportion The population ratio to estimate prevalence for.
#' @param level Double representing the desired confidence interval width.
#' @param precision Integer representing the number of decimal places required.
#'
#' @return A \code{prevalence} object containing the following attributes:
#'   \item{estimates}{Prevalence estimates at the specified years as both absolute and rates.}
#'   \item{simulated}{A \code{data.table} containing simulated incident cases from each bootstrap iteration
#'     Each row corresponds to a simulated incident case with their simulated attributes and survival status.
#'     Binary flags are provided beginning \code{prev_}, which indicate whether that person contributed
#'     to the prevalence for the specified time-period. The \code{prev_registry} flag indicates whether that
#'     person was incident during the registry time-span and alive at the index. These cases are used to
#'     assess the model fit, as the numbers can be simply compared to the known registry prevalence.}
#'   \item{counted}{The number of incident cases present in the registry data set.}
#'   \item{full_surv_model}{The survival model built on the complete registry data set.}
#'   \item{full_inc_model}{The incidence model built on the complete registry data set.}
#'   \item{surv_models}{A list of the survival models fitted to each bootstrap iteration.}
#'   \item{inc_models}{A list of the incidence models fitted to each bootstrap iteration.}
#'   \item{index_dates}{The index date(s).}
#'   \item{est_years}{The years at which prevalence is estimated at.}
#'   \item{counted_incidence_rate}{The overall incidence rate in the registry data set.}
#'   \item{registry_start}{The date the registry was identified at starting at.}
#'   \item{proportion}{The denominator to use for estimating prevalence rates.}
#'   \item{status_col}{The column in the registry data containing the survival status.}
#'   \item{N_boot}{The number of bootstrap iterations that were run.}
#'   \item{means}{Covariate means, used when plotting Kaplan-Meier estimators using \code{survfit}.}
#'   \item{max_event_time}{The maximum time-to-event in the registry data. Again, used in
#'     \code{survfit} to scale the time-axis.}
#'   \item{pval}{The p-value resulting from a hypothesis test on the difference between the
#'   simulated and counted prevalence on the time-span covered by the registry. Tests the
#'   prevalence fit; if a significant result is found then further diagnostics are required.}
#'
#' @references Crouch, Simon, et al. "Determining disease prevalence from
#'   incidence and survival using simulation techniques." Cancer epidemiology
#'   38.2 (2014): 193-199.
#' @examples
#' data(prevsim)
#'
#' \dontrun{
#' data(prevsim)
#'
#' prevalence(index_dates='2013-01-30',
#'            num_years_to_estimate=c(3, 5, 10, 20),
#'            data=prevsim,
#'            inc_formula = entrydate ~ sex,
#'            surv_formula = Surv(time, status) ~ age + sex,
#'            dist='weibull',
#'            population_size = 1e6,
#'            death_column = 'eventdate')
#' }
#'
#' @family prevalence functions
#' @import data.table
#' @export
prevalence <- function(index=NULL, index_dates=NULL, num_years_to_estimate,
                       data,
                       inc_formula=NULL,
                       inc_model=NULL,
                       surv_formula=NULL,
                       surv_model=NULL,
                       registry_start_date=NULL,
                       death_column=NULL,
                       incident_column=NULL,
                       age_column='age',
                       age_dead=100,
                       status_column='status',
                       N_boot=1000,
                       population_size=NULL, proportion=1e5,
                       level=0.95,
                       dist=c('exponential', 'weibull', 'lognormal'),
                       precision=2) {

    # Needed for CRAN check
    alive_at_index <- NULL
    incident_date <- NULL

    if (is.null(incident_column)) {
        if (!is.null(inc_formula)) {
            incident_column <- all.vars(update(inc_formula, .~0))
        }
    }

    dist <- match.arg(dist)

    # Is it right for surv_formula to have precedence over surv_formula?
    if (!is.null(surv_formula)) {
        surv_LHS <- all.vars(update(surv_formula, .~0))
        status_column <- surv_LHS[length(surv_LHS)]
    }

    # Form formula for counted prevalence
    # extract entry column from incidence formula
    if (is.null(death_column)) {
        message("death_column not provided so prevalence cannot be counted over the registry. Estimates will be solely from simulation.")
        counted_formula <- NULL
    } else if (is.null(incident_column)) {
        message("incident_column not provided so prevalence cannot be counted over the registry. Estimates will be solely from simulation.")
        counted_formula <- NULL
    } else {
        counted_formula <- as.formula(paste(death_column, incident_column, sep='~'))
    }

    if (!is.null(incident_column)) {
        if (!incident_column %in% colnames(data)) {
            stop("Error: Cannot find incident column '", incident_column, "' in supplied data set.")
        }
        data[[incident_column]] <- lubridate::ymd(data[[incident_column]])
    }
    if (!is.null(death_column)) {
        if (!death_column %in% colnames(data)) {
            warning("Death column '", death_column, "' not found in supplied data set so estimates will be solely from simulation.")
            death_column <- NULL
        } else {
            data[[death_column]] <- lubridate::ymd(data[[death_column]])
        }
    }

    # This argument allows the user to specify when their registry started. I.e. it could have
    # started a month before received first incident case, in which case would want that taking into account when
    # estimating incidence rate and prevalence
    if (is.null(registry_start_date)) {
        if (is.null(incident_column)) {
            stop("Error: Unknown registry starting date. Please either provide 'registry_start_date' or the incident column in 'inc_formula'.")
        }
        registry_start_date <- min(data[[incident_column]])
    }

    index_provided <- !missing(index) && !is.null(index)
    index_dates_provided <- !missing(index_dates) && !is.null(index_dates)
    if (!index_dates_provided) {
        if (!index_provided) {
            stop("Error: Please provide 'index_dates'.")
        }
        .Deprecated("index", new="index_dates", package="rprev")
        index_dates <- index
        index_dates_provided <- TRUE
    }

    raw_index_dates <- index_dates
    index_dates <- suppressWarnings(lubridate::ymd(index_dates))
    if (anyNA(index_dates)) {
        bad_inputs <- raw_index_dates[is.na(index_dates)]
        stop(
            "Error: Index date(s) '",
            paste(bad_inputs, collapse = ", "),
            "' cannot be parsed as a date. Please enter it as a string in %Y%m%d or %Y-%m-%d format."
        )
    }
    index_dates <- sort(unique(index_dates))
    if (index_dates_provided && index_provided) {
        parsed_index <- suppressWarnings(lubridate::ymd(index))
        if (anyNA(parsed_index)) {
            bad_inputs <- index[is.na(parsed_index)]
            stop(
                "Error: Index date(s) '",
                paste(bad_inputs, collapse = ", "),
                "' cannot be parsed as a date. Please enter it as a string in %Y%m%d or %Y-%m-%d format."
            )
        }
        parsed_index <- sort(unique(parsed_index))
        if (!identical(parsed_index, index_dates)) {
            stop("Error: Both 'index' (deprecated) and 'index_dates' supplied but do not match.")
        }
        warning("Both 'index' (deprecated) and 'index_dates' supplied; using 'index_dates'.", call.=FALSE)
    }
    K <- length(index_dates)
    index_max <- index_dates[length(index_dates)]
    registry_start_date <- lubridate::ymd(registry_start_date)
    sim_start_date <- min(index_dates) - lubridate::years(max(num_years_to_estimate))

    # NEED SIMULATION if registry window is too short OR counted prevalence cannot be computed    
    #   - have N years > R registry years available
    #   - haven't provided date of death for registry data (not always available)
    can_count_registry <- (sim_start_date >= registry_start_date) && !is.null(counted_formula)
    if (!can_count_registry) {

        # Incidence models
        if (is.null(inc_model) & is.null(inc_formula)) {
            stop("Error: Please provide one of inc_model and inc_formula.")
        }
        if (!is.null(inc_model) & !(is.null(inc_formula))) {
            stop("Error: Please provide only one of inc_model and inc_formula.")
        }
        if (is.null(inc_model)) {
            inc_model <- fit_exponential_incidence(inc_formula, data)
        }

        # Survival models
        if (!is.null(surv_model) & !(is.null(surv_formula))) {
            warning("warning: both surv_model and surv_formula provided, survival model will be built using surv_model and surv_formula ignored.")
        }

        if (missing(surv_model) & missing (surv_formula)) {
            stop("Error: Please provide one of surv_model or surv_formula.")
        }

        if (!missing(surv_formula)) {
            surv_model <- build_survreg(surv_formula, data, dist)
        }

        prev_sim <- sim_prevalence(data=data,
                                   index_dates=index_dates,
                                   starting_date=sim_start_date,
                                   inc_model=inc_model,
                                   surv_model=surv_model,
                                   age_column=age_column,
                                   age_dead=age_dead,
                                   N_boot=N_boot)

        # Create column indicating whether contributed to prevalence for each year of interest
        if (K == 1) {
            for (year in num_years_to_estimate) {
                # Determine starting incident date
                starting_incident_date <- index_max - lubridate::years(year)

                # We'll create a new column to hold a binary indicator of whether that observation contributes to prevalence
                col_name <- paste0("prev_", year, "yr")

                # Determine prevalence as incident date is in range and alive at index date
                prev_sim$results[, (col_name) := as.numeric((incident_date > starting_incident_date & incident_date < index_max) & alive_at_index)]
            }
        } else {
            if (!all(c("k_start", "k_end", "incident_date") %in% colnames(prev_sim$results))) {
                stop("Error: multi-index simulation requires k_start, k_end, and incident_date in simulated results.")
            }
        }

    } else {
        prev_sim <- NULL
    }

    # Determine point estimates of prevalence by combining simulated and counted values
    names <- sapply(num_years_to_estimate, function(x) paste('y', x, sep=''))
    if (length(index_dates) == 1) {
        estimates <- lapply(setNames(num_years_to_estimate, names),
                            new_point_estimate,  # Function
                            prev_sim$results, index_dates, data,
                            counted_formula,
                            registry_start_date,
                            status_column,
                            population_size, proportion, level, precision)
    } else {
        sim_agg <- NULL
        if (!is.null(prev_sim)) {
            sim_prev_counts <- build_prev_counts_multiindex(prev_sim$results,
                                                            index_dates,
                                                            num_years_to_estimate)
            sim_prev_counts_pre_registry <- build_prev_counts_multiindex(prev_sim$results[incident_date < registry_start_date],
                                                                         index_dates,
                                                                         num_years_to_estimate)

            total <- as.data.frame(sim_prev_counts)
            total$index_dates <- index_dates[total$k]
            total$contrib_total <- total$prev_count
            total <- total[, c("sim", "index_dates", "year", "contrib_total")]

            pre <- as.data.frame(sim_prev_counts_pre_registry)
            pre$index_dates <- index_dates[pre$k]
            pre$contrib_pre_registry <- pre$prev_count
            pre <- pre[, c("sim", "index_dates", "year", "contrib_pre_registry")]
            sim_agg <- merge(total, pre, by=c("sim", "index_dates", "year"), all=TRUE)
        }

        estimates <- lapply(setNames(num_years_to_estimate, names),
                            new_point_estimate_multiindex,  # Function
                            sim_agg=sim_agg,
                            index_dates=index_dates,
                            registry_data=data,
                            prev_formula=counted_formula,
                            registry_start_date=registry_start_date,
                            status_col=status_column,
                            population_size=population_size,
                            proportion=proportion,
                            level=level,
                            precision=precision,
                            N_boot=N_boot)
    }

    full_surv_model <- if (!is.null(prev_sim)) prev_sim$full_surv_model else NULL
    full_inc_model <- if (!is.null(prev_sim)) prev_sim$full_inc_model else NULL
    surv_models <- if (!is.null(prev_sim)) prev_sim$surv_models else NULL
    inc_models <- if (!is.null(prev_sim)) prev_sim$inc_models else NULL

    if (!is.null(counted_formula)) {
        if (!status_column %in% colnames(data)) {
            stop("Error: cannot find status column '", status_column, "' in data frame.")
        }
        if (K > 1) {
            counted_vec <- vapply(index_dates,
                                  function(tk) counted_prevalence(counted_formula, tk, data, registry_start_date, status_column),
                                  numeric(1))
            counted_prev <- data.frame(index_dates=index_dates,
                                       counted=unname(counted_vec))
        } else {
            counted_prev <- counted_prevalence(counted_formula, index_dates, data, registry_start_date, status_column)
        }
    } else {
        counted_prev <- NULL
    }
    t_ref <- max(index_dates)
    counted_legacy <- counted_prev
    if (K > 1 && is.data.frame(counted_prev)) {
        idx <- match(t_ref, counted_prev$index_dates)
        counted_legacy <- if (!is.na(idx)) counted_prev$counted[idx] else NA_real_
    }
    counted_incidence_rate <- nrow(data) / as.numeric(difftime(t_ref,
                                                               registry_start_date,
                                                               units='days'))

    object <- list(estimates=estimates,
                   simulated=prev_sim$results,
                   counted_by_index=counted_prev,
                   counted=counted_legacy,
                   full_surv_model=full_surv_model,
                   full_inc_model=full_inc_model,
                   surv_models=surv_models,
                   inc_models=inc_models,
                   index_dates=index_dates,
                   est_years=num_years_to_estimate,
                   counted_incidence_rate=counted_incidence_rate,
                   registry_start=registry_start_date,
                   proportion=proportion,
                   status_col=status_column,
                   N_boot=N_boot)

    # Calculate covariate averages for survfit later on
    if (!is.null(prev_sim)) {
        model_for_means <- full_surv_model
        covars <- extract_covars(model_for_means)
        missing_covars <- setdiff(covars, colnames(data))
        if (length(missing_covars) > 0) {
            stop("Error: cannot find covariate(s) in data: ", paste(missing_covars, collapse = ", "))
        }
        # Obtain if continuous or categorical
        is_factor <- sapply(covars, function(x) is.factor(data[[x]]) || is.character(data[[x]]))
        fact_cols <- covars[is_factor]
        cont_cols <- covars[!is_factor]
        cont_means <- sapply(cont_cols, function(x) mean(data[[x]], na.rm=TRUE))
        cat_modes <- sapply(fact_cols, function(x) {
            vals <- stats::na.omit(data[[x]])
            if (length(vals) == 0) {
                return(NA)
            }
            names(which.max(table(vals)))
        })

        # Save into data frame
        means <- data.frame()
        for (var in cont_cols) {
            means[1, var] <- cont_means[var]
        }
        for (var in fact_cols) {
            means[1, var] <- cat_modes[var]
            levels_vec <- if (is.factor(data[[var]])) {
                levels(data[[var]])
            } else {
                sort(unique(stats::na.omit(data[[var]])))
            }
            means[[var]] <- factor(means[[var]], levels=levels_vec)
        }
        object$means <- means
    } else {
        object$means <- NULL
    }

    # Add max time if possible
    if (!is.null(incident_column) & !is.null(death_column)) {
        time_vals <- as.numeric(difftime(data[[death_column]],
                                         data[[incident_column]],
                                         units='days'))
        time_vals <- time_vals[!is.na(time_vals)]
        if (length(time_vals) == 0) {
            object$max_event_time <- NULL
        } else {
            if (any(time_vals < 0)) {
                warning("Negative event times detected in registry data.")
            }
            object$max_event_time <- max(time_vals)
        }
    }


    if (!is.null(prev_sim)) {
        if (K > 1) {
            pval_ref <- test_prevalence_fit(object)
            pval_by_index <- setNames(rep(NA_real_, K), as.character(index_dates))
            pval_by_index[as.character(t_ref)] <- pval_ref
            object$pval <- pval_ref
            object$pval_by_index <- pval_by_index
        } else {
            object$pval <- test_prevalence_fit(object)
        }
    }

    attr(object, 'class') <- 'prevalence'
    object

}

# Build prevalence counts across multiple index dates from interval contributions.
build_prev_counts_multiindex <- function(results, index_dates, years) {
    # Needed for CRAN check
    delta <- incident_date <- k <- k_end <- k_start <- prev_count <- sim <- year <- NULL

    required <- c("sim", "incident_date", "k_start", "k_end")
    missing_cols <- setdiff(required, colnames(results))
    if (length(missing_cols) > 0) {
        stop("Error: missing columns in results: ", paste(missing_cols, collapse = ", "))
    }

    if (length(index_dates) == 0 || length(years) == 0 || nrow(results) == 0) {
        return(data.table::data.table(sim=integer(),
                                      year=integer(),
                                      k=integer(),
                                      prev_count=integer()))
    }

    idx_dates <- as.Date(index_dates)
    K <- length(idx_dates)

    dt <- data.table::as.data.table(results)[, ..required]
    dt[, incident_date := as.Date(incident_date)]

    # Expand by year (sim, case) x year.
    dt_exp <- dt[rep(seq_len(nrow(dt)), each = length(years))]
    dt_exp[, year := rep(years, times = nrow(dt))]

    # Strict window: incident_date < t_k and incident_date > t_k - years(year).
    dt_exp[, k_inc_start := findInterval(incident_date, idx_dates) + 1L]
    dt_exp[, k_inc_end := findInterval(incident_date + lubridate::years(year) - 1, idx_dates)]

    dt_exp[, k_contrib_start := pmax(k_start, k_inc_start)]
    dt_exp[, k_contrib_end := pmin(k_end, k_inc_end)]

    valid <- !is.na(dt_exp$k_contrib_start) & !is.na(dt_exp$k_contrib_end) &
        dt_exp$k_contrib_start <= dt_exp$k_contrib_end &
        dt_exp$k_contrib_start <= K & dt_exp$k_contrib_end >= 1L

    if (!any(valid)) {
        return(data.table::data.table(sim=integer(),
                                      year=integer(),
                                      k=integer(),
                                      prev_count=integer()))
    }

    events <- data.table::rbindlist(list(
        dt_exp[valid, .(sim, year, k = k_contrib_start, delta = 1L)],
        dt_exp[valid & k_contrib_end < K, .(sim, year, k = k_contrib_end + 1L, delta = -1L)]
    ), use.names = TRUE)

    events <- events[, .(delta = sum(delta)), by = .(sim, year, k)]

    grid <- data.table::CJ(sim = unique(dt$sim), year = years, k = seq_len(K))
    grid <- events[grid, on = .(sim, year, k)]
    grid[is.na(delta), delta := 0L]
    data.table::setorder(grid, sim, year, k)
    grid[, prev_count := cumsum(delta), by = .(sim, year)]
    grid[, delta := NULL]
    grid
}

#' Count prevalence from registry data.
#'
#' Counts contribution to prevalence at a specific index from each year of a
#' registry. A person is included as contributing to disease prevalence if they
#' are incident within the specified time-span, and are either alive or censored
#' at the index date. The rationale for including censored cases in prevalence
#' estimation is that such cases have typically been lost to follow-up, and are
#' often more likely to have been alive at the index date than not.
#'
#' @inheritParams prevalence
#' @param formula A formula of the form <event date column> ~ <entry date column>.
#' @param start_date The initial date to start counting prevalence from as a \code{Date} object.
#'     Typically the index date - (Nyears * 365.25). Allows for non-whole year prevalence estimations.
#' @param status_col The name of the column holding a binary indicator variable
#'     of whether the individual experienced an event at their event time or was
#'     censored.
#'
#' @return The number of prevalent cases at the specified index date(s).
#'   Returns a single integer for length-1 \code{index_dates} and a numeric
#'   vector for multiple dates.
counted_prevalence <- function(formula, index_dates, data, start_date, status_col) {
    death_col <- all.vars(update(formula, .~0))
    entry_col <- all.vars(update(formula, 0~.))

    if (length(death_col) != 1 || length(entry_col) != 1) {
        stop("Error: formula must contain exactly one event column and one entry column.")
    }

    index_dates <- as.Date(index_dates)
    start_date <- as.Date(start_date)
    K <- length(index_dates)

    if (length(start_date) == 1) {
        start_date <- rep(start_date, K)
    }
    if (length(start_date) != K) {
        stop("Error: start_date must have length 1 or match length(index_dates).")
    }

    entry <- data[[entry_col]]
    death <- data[[death_col]]
    status <- data[[status_col]]

    if (!inherits(entry, "Date")) {
        entry <- as.Date(entry)
    }
    if (!inherits(death, "Date")) {
        death <- as.Date(death)
    }

    count_one <- function(idx_date, start_date_k) {
        incident <- !is.na(entry) & entry >= start_date_k & entry < idx_date
        dead_at_index <- (!is.na(status) & status == 1) & !is.na(death) & (death <= idx_date)
        sum(incident & !dead_at_index, na.rm=TRUE)
    }

    if (K == 1) {
        count_one(index_dates[1], start_date[1])
    } else {
        out <- vapply(seq_len(K),
                      function(i) count_one(index_dates[i], start_date[i]),
                      numeric(1))
        names(out) <- as.character(index_dates)
        out
    }
}

#' Estimate prevalence using Monte Carlo simulation.
#'
#' Estimates prevalent cases at a specific index date by use of Monte Carlo
#' simulation. Simulated cases are marked with age and sex to enable agreement
#' with population survival data where a cure model is used, and calculation of
#' the posterior distributions of each.
#' @inheritParams prevalence
#' @param starting_date The initial date to start simulating prevalence from as a \code{Date} object.
#'     Typically the index date - (Nyears * 365.25). Allows for non-whole year prevalence estimations.
#'
#' @return A list with the following attributes:
#'   \item{results}{A data.table containing the simulated incident populations from each
#'   simulation along with their covariates and survival status at the index.}
#'   \item{full_surv_model}{The survival model built on the full registry data set.}
#'   \item{full_inc_model}{The incidence model built on the full registry data set.}
#'   \item{surv_models}{A list containing survival models built on each bootstrap sample.}
#'   \item{inc_models}{A list containing incidence models built on each bootstrap sample.}
#' @importFrom survival survfit
sim_prevalence <- function(data, index_dates, starting_date,
                           inc_model, surv_model,
                           age_column='age',
                           N_boot=1000,
                           age_dead=100)
                           {

    # Needed for CRAN check
    alive_at_index <- NULL
    time_to_index <- NULL
    k_start <- NULL
    k_end <- NULL
    xi <- NULL

    raw_index_dates <- index_dates
    index_dates <- suppressWarnings(lubridate::ymd(index_dates))
    if (anyNA(index_dates)) {
        bad_inputs <- raw_index_dates[is.na(index_dates)]
        stop("Error: Index date(s) '", paste(bad_inputs, collapse=", "),
             "' cannot be parsed as a date. Please enter it as a string in %Y%m%d or %Y-%m-%d format.")
    }
    index_dates <- sort(unique(index_dates))
    K <- length(index_dates)
    if (K == 0) {
        stop("Error: No valid index dates provided.")
    }

    starting_date <- lubridate::ymd(starting_date)
    if (anyNA(starting_date)) {
        stop("Error: starting_date cannot be parsed as a date.")
    }
    index_max <- max(index_dates)

    data <- data[complete.cases(data), ]
    full_data <- data

    covars <- extract_covars(surv_model)
    number_incident_days <- as.numeric(difftime(index_max, starting_date, units='days'))
    if (number_incident_days <= 0) {
        stop("Error: starting_date must be earlier than index date(s).")
    }

    run_sample <- function() {
        # bootstrap dataset
        data <- full_data[sample(seq(nrow(full_data)), replace=T), ]

        # fit incidence and survival models.
        bs_inc <- eval(inc_model$call)
        bs_surv <- eval(surv_model$call)

        # Draw the incident population using the fitted model and predict their death times
        incident_population <- draw_incident_population(bs_inc, full_data, number_incident_days, extract_covars(bs_surv))
        data.table::setDT(incident_population)

        # For each individual determine the time between incidence and each index
        incident_date <- as.Date(starting_date + incident_population[[1]])
        time_to_index_mat <- t(outer(index_dates,
                                     incident_date,
                                     FUN=function(tk, dj) as.numeric(difftime(tk, dj, units='days'))))
        mask <- time_to_index_mat >= 0
        times_mat <- pmax(time_to_index_mat, 0)

        if (K == 1) {
            time_to_index <- time_to_index_mat[, 1]
            surv_prob <- predict_survival_probability(bs_surv, incident_population[, -1], times_mat[, 1])
            alive_at_index <- rbinom(length(surv_prob), size=1, prob=surv_prob)
            alive_mat <- matrix(alive_at_index, ncol=1)
        } else {
            p_mat <- matrix(NA_real_, nrow=nrow(time_to_index_mat), ncol=K)
            for (k in seq_len(K)) {
                p_mat[, k] <- predict_survival_probability(bs_surv, incident_population[, -1], times_mat[, k])
            }
            xi <- runif(nrow(p_mat))
            alive_mat <- (xi <= p_mat) & mask
            alive_at_index <- as.integer(alive_mat[, K])
        }

        any_alive <- rowSums(alive_mat) > 0
        k_start <- ifelse(any_alive, max.col(alive_mat, ties.method="first"), NA_integer_)
        idx_rev <- max.col(alive_mat[, K:1, drop=FALSE], ties.method="first")
        k_end <- K - idx_rev + 1
        k_end[!any_alive] <- NA_integer_

        if (!is.null(age_column) & age_column %in% colnames(incident_population)) {
            cap_days <- (age_dead - incident_population[[age_column]]) * DAYS_IN_YEAR
            cap_days[is.na(cap_days)] <- Inf
            k_age_end <- ifelse(is.infinite(cap_days),
                                K,
                                findInterval(incident_date + cap_days, index_dates))
            k_age_end[cap_days < 0] <- 0
            k_end <- pmin(k_end, k_age_end)
            invalid <- is.na(k_start) | is.na(k_end) | k_end < k_start
            k_start[invalid] <- NA_integer_
            k_end[invalid] <- NA_integer_
        }

        alive_at_index <- as.integer(!is.na(k_start) & !is.na(k_end) & k_start <= K & k_end >= K)

        incident_population[, 'incident_date' := incident_date]
        incident_population[, 'k_start' := k_start]
        incident_population[, 'k_end' := k_end]
        incident_population[, 'alive_at_index' := alive_at_index]
        if (K == 1) {
            incident_population[, 'time_to_index' := time_to_index]
        }
        list(pop=incident_population,
             surv=bs_surv,
             inc=bs_inc)
    }

    # TODO Implement this and turn into user facing argument
    n_cores <- 1
    if (n_cores > 1) {
        message("Multi-core functionality not currently implemented, defaulting to single-core.")
    }
    all_results <- replicate(N_boot, run_sample(), simplify=FALSE)
    names(all_results) <- seq_len(N_boot)

    # Combine incident population into single table
    pops <- lapply(all_results, function(x) {
        pop <- x$pop
        data.table::setDT(pop)
        pop
    })

    if (length(pops) > 0) {
        ref_cols <- names(pops[[1]])
        for (i in seq_along(pops)) {
            cols_i <- names(pops[[i]])
            if (!identical(cols_i, ref_cols)) {
                stop("Error: inconsistent columns in bootstrap population at sim ",
                     i, ". Expected: ", paste(ref_cols, collapse=", "),
                     "; got: ", paste(cols_i, collapse=", "), ".")
            }
            required_cols <- c("incident_date", "k_start", "k_end", "alive_at_index")
            missing_cols <- setdiff(required_cols, cols_i)
            if (length(missing_cols) > 0) {
                stop("Error: missing required columns in bootstrap population at sim ",
                     i, ": ", paste(missing_cols, collapse=", "), ".")
            }
            if (!inherits(pops[[i]]$incident_date, "Date")) {
                stop("Error: incident_date must be a Date in bootstrap population at sim ", i, ".")
            }
            if (!is.numeric(pops[[i]]$k_start) || !is.numeric(pops[[i]]$k_end)) {
                stop("Error: k_start/k_end must be numeric in bootstrap population at sim ", i, ".")
            }
            if (!is.numeric(pops[[i]]$alive_at_index)) {
                stop("Error: alive_at_index must be numeric in bootstrap population at sim ", i, ".")
            }
        }
    }

    results <- data.table::rbindlist(pops, idcol='sim', use.names=TRUE, fill=FALSE)

    # Force death at 100 if possible (handled in run_sample); warn once if age unavailable
    if (is.null(age_column) | !age_column %in% colnames(results)) {
        message("No column found for age in ", age_column, ", so cannot assume death at 100 years of age. Be careful of 'infinite' survival times.")
    }

    # These intermediary columns aren't useful for the user and would just clutter up the output
    cols_drop <- intersect(c('time_to_index', 'time_to_entry'), colnames(results))
    if (length(cols_drop) > 0) {
        results[, (cols_drop) := NULL]
    }

    list(results=results,
         full_surv_model=surv_model,
         full_inc_model=inc_model,
         surv_models=lapply(all_results, function(x) x$surv),
         inc_models=lapply(all_results, function(x) x$inc))
}

#' @export
print.prevalence <- function(x, ...) {
    if (is.null(x$index_dates)) {
        stop("Error: prevalence object missing index_dates.")
    }
    index_dates <- x$index_dates
    years <- x$est_years

    get_rate_name <- function(est_names) {
        rate_names <- grep("^per[^.]*$", est_names, value=TRUE)
        if (length(rate_names) == 0) {
            return(NULL)
        }
        if (length(rate_names) > 1) {
            stop("Error: multiple rate columns found in estimates: ", paste(rate_names, collapse=", "))
        }
        rate_names
    }

    if (length(index_dates) <= 1) {
        cat(paste0("Estimated prevalence at ", index_dates, ":\n"))
        lapply(paste0("y", years), function(item) {
            est <- x$estimates[[item]]
            abs_prev_est <- est[["absolute.prevalence"]]
            rate_name <- get_rate_name(names(est))
            if (!is.null(rate_name)) {
                rel_prev <- est[[rate_name]]
                rel_prev_est <- paste0("(", rel_prev, " per ", x$proportion, ")")
            } else {
                rel_prev_est <- NULL
            }
            year <- gsub("^y", "", item)
            cat(paste(year, "years:", abs_prev_est, rel_prev_est, '\n'))
        })
    } else {
        abs_table <- data.frame(index_dates=index_dates)
        for (year in years) {
            item <- paste0("y", year)
            est <- x$estimates[[item]]
            abs_vec <- est[["absolute.prevalence"]]
            if (length(abs_vec) != length(index_dates)) {
                stop("Error: length mismatch for ", item, " absolute.prevalence.")
            }
            abs_table[[paste0(year, "y")]] <- abs_vec
        }

        cat("Estimated prevalence:\n")
        print(abs_table, row.names=FALSE)

        rate_name <- get_rate_name(names(x$estimates[[paste0("y", years[1])]]))
        if (!is.null(rate_name)) {
            rate_table <- data.frame(index_dates=index_dates)
            for (year in years) {
                item <- paste0("y", year)
                est <- x$estimates[[item]]
                rate_vec <- est[[rate_name]]
                if (length(rate_vec) != length(index_dates)) {
                    stop("Error: length mismatch for ", item, " rate column.")
                }
                rate_table[[paste0(year, "y")]] <- rate_vec
            }
            cat(paste0("Prevalence rate (per ", x$proportion, "):\n"))
            print(rate_table, row.names=FALSE)
        }
    }
    invisible(x)
}

#' @export
summary.prevalence <- function(object, ...) {

    # Required to pass R CMD CHECK
    incident_date <- NULL
    sim <- NULL
    V1 <- NULL

    if (is.null(object$index_dates)) {
        stop("Error: prevalence object missing index_dates.")
    }
    index_dates <- object$index_dates
    K <- length(index_dates)

    cat("Prevalence \n~~~~~~~~~~\n")
    print(object)
    cat("\n")

    cat("Registry Data\n~~~~~~~~~~~~~\n")
    if (K > 1) {
        cat("Index dates:", paste0(min(index_dates), " to ", max(index_dates), " (K=", K, ")"), "\n")
    } else {
        cat("Index date:", as.character(index_dates), "\n")
    }
    cat("Start date:", as.character(object$registry_start), "\n")
    cat("Overall incidence rate:", round(object$counted_incidence_rate, 3), " (per day, up to ", max(index_dates), ")\n")

    cat("Counted prevalent cases:\n")
    if (is.null(object$counted_by_index)) {
        cat(object$counted, "\n")
    } else {
        counted_tbl <- object$counted_by_index
        if (is.numeric(counted_tbl)) {
            counted_tbl <- data.frame(index_dates=index_dates, counted=counted_tbl)
        }
        print(counted_tbl, row.names=FALSE)
    }

    if (!is.null(object$simulated) && !all(is.na(object$simulated))) {
        cat("\nSimulation\n~~~~~~~~~~\n")
        cat("Iterations:", object$N_boot, "\n")
        cat("Average incidence rate:",
            round(object$simulated[, length(incident_date), by=sim][, mean(V1)] / (max(object$est_years)*DAYS_IN_YEAR), 3),
            "\n")
        if (!is.null(object$pval_by_index)) {
            cat("P-values:\n")
            p_tbl <- data.frame(index_dates=index_dates, p_value=as.numeric(object$pval_by_index))
            print(p_tbl, row.names=FALSE)
        } else if (!is.null(object$pval)) {
            cat("P-value:", object$pval)
        }
    }
}
