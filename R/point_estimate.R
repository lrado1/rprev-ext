new_point_estimate <- function(year, sim_results, index_dates, registry_data, prev_formula, registry_start_date, status_col,
                               population_size=NULL, proportion=1e5,
                               level=0.95, precision=2) {
    if (length(index_dates) > 1) {
        return(new_point_estimate_multiindex(year=year,
                                             sim_agg=sim_results,
                                             index_dates=index_dates,
                                             registry_data=registry_data,
                                             prev_formula=prev_formula,
                                             registry_start_date=registry_start_date,
                                             status_col=status_col,
                                             population_size=population_size,
                                             proportion=proportion,
                                             level=level,
                                             precision=precision))
    }

    if (year <= 0) {
        warning("Cannot estimate prevalence for a non-positive value of num_year_to_estimate.")
        return(list(absolute.prevalence=0))
    }

    # CRAN check
    incident_date <- NULL
    sim <- NULL

    # See if need simulation if have less registry data than required
    index_ref <- index_dates[1]
    initial_date <- index_ref - lubridate::years(year)
    need_simulation <- initial_date < registry_start_date

    # Only count prevalence if formula isn't null
    if (!is.null(prev_formula)) {
        count_prev <- counted_prevalence(prev_formula, index_ref, registry_data, max(initial_date, registry_start_date), status_col)

        # See if appending prevalence to simulation data or it's entirely counted
        if (initial_date < registry_start_date) {
            stopifnot(!is.null(sim_results))

            col_name <- paste0("prev_", year, "yr")
            sim_contributions <- sim_results[incident_date < registry_start_date][, sum(get(col_name)), by=sim][[2]]  # Return results column
            the_estimate <- count_prev + mean(sim_contributions)

            # Closure to calculate combined standard error
            se_func <- build_se_func(counted_contribs=count_prev, sim_contribs=sim_contributions)

        } else {
            the_estimate <- count_prev
            # Closure to calculate standard error of counted data
            se_func <- build_se_func(counted_contribs=count_prev)
        }
    } else {
        # If don't have counted data then prevalence estimates are entirely simulated
        col_name <- paste0("prev_", year, "yr")
        sim_contributions <- sim_results[, sum(get(col_name)), by=sim][[2]]  # Return results column
        the_estimate <- mean(sim_contributions)

        # Closure to calculate standard error of simulated data
        se_func <- build_se_func(sim_contribs=sim_contributions)
    }

    result <- list(absolute.prevalence=the_estimate)

    if (!is.null(population_size)) {
        the_proportion <- (the_estimate / population_size) * proportion
        se <- se_func(population_size)

        z_level <- qnorm((1+level)/2)
        CI <- z_level * se * proportion

        # Setup labels for proportion list outputs
        est_lab <- paste0('per', proportion_label(proportion))
        upper_lab <- paste(est_lab, '.upper', sep='')
        lower_lab <- paste(est_lab, '.lower', sep='')
        result[[est_lab]] <- the_proportion
        result[[upper_lab]] <- the_proportion + CI
        result[[lower_lab]] <- the_proportion - CI
    }

    lapply(result, round, precision)
}

new_point_estimate_multiindex <- function(year,
                                          sim_agg,
                                          index_dates,
                                          registry_data,
                                          prev_formula,
                                          registry_start_date,
                                          status_col,
                                          population_size=NULL,
                                          proportion=1e5,
                                          level=0.95,
                                          precision=2,
                                          N_boot=NULL) {
    if (year <= 0) {
        warning("Cannot estimate prevalence for a non-positive value of num_year_to_estimate.")
        return(data.frame(index_dates=index_dates, absolute.prevalence=0))
    }

    if (length(index_dates) < 1) {
        stop("Error: index_dates must have length >= 1.")
    }
    if (!inherits(index_dates, "Date")) {
        stop("Error: index_dates must be a Date vector.")
    }

    if (!is.null(sim_agg)) {
        sim_agg <- as.data.frame(sim_agg)
        required <- c("sim", "index_dates", "year", "contrib_total")
        missing_cols <- setdiff(required, names(sim_agg))
        if (length(missing_cols) > 0) {
            stop("Error: sim_agg missing required columns: ", paste(missing_cols, collapse=", "))
        }
        if (!inherits(sim_agg$index_dates, "Date")) {
            sim_agg$index_dates <- as.Date(sim_agg$index_dates)
        }
        if (is.null(N_boot) && nrow(sim_agg) > 0) {
            N_boot <- max(sim_agg$sim)
        }
    }

    if (!is.null(prev_formula) && !status_col %in% colnames(registry_data)) {
        stop("Error: cannot find status column '", status_col, "' in registry data.")
    }

    initial_dates <- index_dates - lubridate::years(year)
    need_sim <- initial_dates < registry_start_date

    if (!is.null(prev_formula)) {
        count_prev_vec <- vapply(seq_along(index_dates), function(i) {
            counted_prevalence(prev_formula,
                               index_dates[i],
                               registry_data,
                               max(initial_dates[i], registry_start_date),
                               status_col)
        }, numeric(1))
    } else {
        count_prev_vec <- rep(0, length(index_dates))
    }

    need_sim_any <- is.null(prev_formula) || any(need_sim)
    if (need_sim_any) {
        if (is.null(sim_agg)) {
            stop("Error: sim_agg is required for simulation contributions.")
        }
        if (is.null(N_boot)) {
            stop("Error: N_boot is required to build simulation contribution vectors.")
        }
        if (!"contrib_pre_registry" %in% names(sim_agg) && any(need_sim)) {
            stop("Error: sim_agg missing required column: contrib_pre_registry.")
        }
    }

    get_sim_contribs <- function(idx_date, col_name) {
        if (is.null(sim_agg)) {
            stop("Error: sim_agg is required for simulation contributions.")
        }
        sub <- sim_agg[sim_agg$year == year & sim_agg$index_dates == idx_date, c("sim", col_name)]
        contribs <- rep(0, N_boot)
        if (nrow(sub) > 0) {
            contribs[sub$sim] <- sub[[col_name]]
        }
        contribs
    }

    rows <- vector("list", length(index_dates))
    for (i in seq_along(index_dates)) {
        index_k <- index_dates[i]
        if (is.null(prev_formula)) {
            sim_contribs <- get_sim_contribs(index_k, "contrib_total")
            the_estimate <- mean(sim_contribs)
            se_func <- build_se_func(sim_contribs=sim_contribs)
        } else if (need_sim[i]) {
            sim_contribs <- get_sim_contribs(index_k, "contrib_pre_registry")
            the_estimate <- count_prev_vec[i] + mean(sim_contribs)
            se_func <- build_se_func(counted_contribs=count_prev_vec[i], sim_contribs=sim_contribs)
        } else {
            the_estimate <- count_prev_vec[i]
            se_func <- build_se_func(counted_contribs=count_prev_vec[i])
        }

        row <- data.frame(index_dates=index_k, absolute.prevalence=the_estimate)
        if (!is.null(population_size)) {
            the_proportion <- (the_estimate / population_size) * proportion
            se <- se_func(population_size)
            z_level <- qnorm((1+level)/2)
            CI <- z_level * se * proportion

            est_lab <- paste0('per', proportion_label(proportion))
            upper_lab <- paste(est_lab, '.upper', sep='')
            lower_lab <- paste(est_lab, '.lower', sep='')
            row[[est_lab]] <- the_proportion
            row[[upper_lab]] <- the_proportion + CI
            row[[lower_lab]] <- the_proportion - CI
        }

        rows[[i]] <- row
    }

    result <- do.call(rbind, rows)
    result[] <- lapply(result, function(x) {
        if (is.numeric(x) && !inherits(x, "Date")) round(x, precision) else x
    })
    result
}

build_se_func <- function(counted_contribs=NULL, sim_contribs=NULL) {
    # Pure simulated
    if (is.null(counted_contribs)) {
        function(pop_size) {
            calculate_se_sim(pop_size, sim_contribs)
        }
    }  else if (is.null(sim_contribs)) {
        # Pure counted
        function(pop_size) {
            calculate_se_counted(pop_size, counted_contribs)
        }
    } else {
        # Combination
        function(pop_size) {
            calculate_se_combined(pop_size, counted_contribs, sim_contribs)
        }
    }
}

calculate_se_combined <- function(population_size, counted_contribs, sim_contribs) {
    calculate_se_sim(population_size, sim_contribs) +
        calculate_se_counted(population_size, counted_contribs)
}

calculate_se_sim <- function(population_size, sim_contribs) {
    sd(sim_contribs) / population_size
}
calculate_se_counted <- function(population_size, counted_contribs) {
    raw_proportion <- counted_contribs / population_size
    sqrt((raw_proportion * (1 - raw_proportion)) / population_size)
}
