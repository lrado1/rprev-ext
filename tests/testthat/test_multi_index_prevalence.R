library(rprev)

data(prevsim)

make_counted_prevalence <- function(index = c("2010-01-01", "2011-01-01", "2012-01-01"),
                                    years = c(1, 2),
                                    population_size = NULL) {
    prevalence(
        index = index,
        num_years_to_estimate = years,
        data = prevsim,
        incident_column = "entrydate",
        death_column = "eventdate",
        population_size = population_size
    )
}

test_that("a single index date retains single-date behaviour", {
    scalar <- make_counted_prevalence(index = "2012-01-01", years = 2)
    length_one_vector <- make_counted_prevalence(index = c("2012-01-01"), years = 2)

    expect_s3_class(scalar, "prevalence")
    expect_equal(scalar, length_one_vector)
    expect_equal(scalar$index_dates, as.Date("2012-01-01"))
    expect_equal(nrow(scalar$estimates$y2), 1)
    expect_null(scalar$simulated)
})

test_that("multi-index estimates have one row per date and period", {
    object <- make_counted_prevalence(population_size = 1e6)
    expected_dates <- as.Date(c("2010-01-01", "2011-01-01", "2012-01-01"))

    expect_equal(names(object$estimates), c("y1", "y2"))
    expect_true(all(vapply(object$estimates, nrow, integer(1)) == 3))
    expect_true(all(vapply(object$estimates, function(x) identical(x$index_date, expected_dates), logical(1))))
    expect_true(all(vapply(
        object$estimates,
        function(x) identical(names(x), c("index_date", "absolute.prevalence", "per100K", "per100K.upper", "per100K.lower")),
        logical(1)
    )))
    expect_true(all(vapply(object$estimates, function(x) {
        all(x$per100K.lower <= x$per100K & x$per100K <= x$per100K.upper)
    }, logical(1))))
})

test_that("index dates are sorted before estimation", {
    unsorted <- make_counted_prevalence(
        index = c("2012-01-01", "2010-01-01", "2011-01-01"),
        years = 2
    )
    sorted <- make_counted_prevalence(
        index = c("2010-01-01", "2011-01-01", "2012-01-01"),
        years = 2
    )

    expect_equal(unsorted, sorted)
    expect_equal(unsorted$index_dates, sort(unsorted$index_dates))
})

test_that("duplicate index dates are removed before estimation", {
    duplicated <- make_counted_prevalence(
        index = c("2011-01-01", "2010-01-01", "2011-01-01", "2012-01-01"),
        years = 2
    )
    unique_dates <- make_counted_prevalence(
        index = c("2010-01-01", "2011-01-01", "2012-01-01"),
        years = 2
    )

    expect_equal(duplicated, unique_dates)
    expect_equal(length(duplicated$index_dates), 3)
})

test_that("Date vectors are accepted as index input", {
    date_input <- make_counted_prevalence(
        index = as.Date(c("2010-01-01", "2011-01-01")),
        years = 2
    )
    character_input <- make_counted_prevalence(
        index = c("2010-01-01", "2011-01-01"),
        years = 2
    )
    compact_character_input <- make_counted_prevalence(
        index = c("20100101", "20110101"),
        years = 2
    )

    expect_equal(date_input, character_input)
    expect_equal(compact_character_input, character_input)
})

test_that("an invalid value in an index vector reports the offending date", {
    expect_error(
        make_counted_prevalence(
            index = c("2010-01-01", "not-a-date", "2012-01-01"),
            years = 2
        ),
        "not-a-date.*cannot be parsed"
    )
    expect_error(
        make_counted_prevalence(
            index = c("2010-01-01", NA_character_, "2012-01-01"),
            years = 2
        ),
        "NA.*cannot be parsed"
    )
})

test_that("the counted-only branch returns per-date registry estimates", {
    object <- make_counted_prevalence(
        index = c("2011-01-01", "2012-01-01"),
        years = 2
    )
    registry <- prevsim
    registry$entrydate <- as.Date(registry$entrydate)
    registry$eventdate <- as.Date(registry$eventdate)
    expected <- vapply(object$index_dates, function(index_date) {
        rprev:::counted_prevalence(
            eventdate ~ entrydate,
            index_date,
            registry,
            index_date - lubridate::years(2),
            "status"
        )
    }, numeric(1))

    expect_null(object$simulated)
    expect_null(object$inc_models)
    expect_null(object$surv_models)
    expect_equal(names(object$counted), as.character(object$index_dates))
    expect_equal(object$estimates$y2$absolute.prevalence, expected)
})

test_that("counted prevalence includes the period start and excludes the index date", {
    registry <- data.frame(
        entrydate = as.Date(c("2018-01-01", "2019-12-31", "2020-01-01")),
        eventdate = as.Date(rep("2021-01-01", 3)),
        status = rep(1, 3)
    )

    counted <- rprev:::counted_prevalence(
        eventdate ~ entrydate,
        as.Date("2020-01-01"),
        registry,
        as.Date("2018-01-01"),
        "status"
    )

    expect_equal(counted, 2)
})

make_simulated_prevalence <- function(index = c("2010-01-01", "2011-01-01", "2012-01-01"),
                                      years = c(5, 8),
                                      seed = 2026,
                                      n_boot = 2,
                                      death_column = "eventdate",
                                      age_dead = 100) {
    set.seed(seed)
    suppressMessages(prevalence(
        index = index,
        num_years_to_estimate = years,
        data = prevsim,
        inc_formula = entrydate ~ 1,
        surv_formula = Surv(time, status) ~ age,
        dist = "weibull",
        death_column = death_column,
        age_dead = age_dead,
        N_boot = n_boot
    ))
}

multi_index_cache <- new.env(parent = emptyenv())

get_mixed_prevalence <- function() {
    if (!exists("mixed", envir = multi_index_cache, inherits = FALSE)) {
        multi_index_cache$mixed <- make_simulated_prevalence()
    }
    multi_index_cache$mixed
}

test_that("the simulation gate uses the earliest date and longest period", {
    counted_only <- make_counted_prevalence(
        index = c("2010-01-01", "2012-01-01"),
        years = 5
    )
    requires_simulation <- get_mixed_prevalence()

    expect_null(counted_only$simulated)
    expect_false(is.null(requires_simulation$simulated))
    expect_equal(sort(unique(requires_simulation$simulated$sim)), 1:2)
})

test_that("one simulated population is shared across index dates", {
    multi <- make_simulated_prevalence(
        index = c("2010-01-01", "2012-01-01"),
        years = 8,
        seed = 501
    )
    equivalent_single <- make_simulated_prevalence(
        index = "2012-01-01",
        years = 10,
        seed = 501
    )
    population_columns <- c("sim", "incident_date", "age")

    expect_equal(
        as.data.frame(multi$simulated[, population_columns, with = FALSE]),
        as.data.frame(equivalent_single$simulated[, population_columns, with = FALSE])
    )
})

test_that("survival states are masked before incidence and consistent over time", {
    object <- get_mixed_prevalence()
    alive_columns <- sprintf("alive_k%03d", seq_along(object$index_dates))
    alive <- as.matrix(object$simulated[, alive_columns, with = FALSE])
    eligible <- outer(
        as.numeric(object$simulated$incident_date),
        as.numeric(object$index_dates),
        "<="
    )

    expect_true(all(alive[!eligible] == 0))
    monotone_after_incidence <- vapply(seq_len(nrow(alive)), function(i) {
        states <- alive[i, eligible[i, ]]
        length(states) < 2 || all(diff(states) <= 0)
    }, logical(1))
    expect_true(all(monotone_after_incidence))
})

test_that("prevalence contribution columns match every date-period pair", {
    object <- get_mixed_prevalence()
    contribution_columns <- grep("^prev_[0-9]+yr_k[0-9]+$", names(object$simulated), value = TRUE)
    expected_columns <- as.vector(outer(
        object$est_years,
        seq_along(object$index_dates),
        function(year, k) sprintf("prev_%dyr_k%03d", year, k)
    ))

    expect_equal(sort(contribution_columns), sort(expected_columns))
    for (year in object$est_years) {
        for (k in seq_along(object$index_dates)) {
            expected <- as.numeric(
                object$simulated$incident_date > object$index_dates[k] - lubridate::years(year) &
                    object$simulated$incident_date < object$index_dates[k] &
                    object$simulated[[sprintf("alive_k%03d", k)]] == 1
            )
            expect_equal(object$simulated[[sprintf("prev_%dyr_k%03d", year, k)]], expected)
        }
    }
})

test_that("mixed estimates combine counted and pre-registry simulated cases", {
    object <- get_mixed_prevalence()
    index_date <- object$index_dates[1]
    registry <- prevsim
    registry$entrydate <- as.Date(registry$entrydate)
    registry$eventdate <- as.Date(registry$eventdate)
    counted <- rprev:::counted_prevalence(
        eventdate ~ entrydate,
        index_date,
        registry,
        object$registry_start,
        "status"
    )
    simulated_by_bootstrap <- vapply(
        split(object$simulated, object$simulated$sim),
        function(sample) {
            sum(sample$prev_8yr_k001[sample$incident_date < object$registry_start])
        },
        numeric(1)
    )

    expect_equal(
        object$estimates$y8$absolute.prevalence[1],
        round(counted + mean(simulated_by_bootstrap), 2)
    )
})

test_that("simulation-only estimates use all simulated contributions", {
    object <- make_simulated_prevalence(
        index = c("2010-01-01", "2012-01-01"),
        years = 5,
        death_column = NULL,
        seed = 711
    )

    expect_null(object$counted)
    for (k in seq_along(object$index_dates)) {
        column <- sprintf("prev_5yr_k%03d", k)
        expected <- mean(vapply(
            split(object$simulated, object$simulated$sim),
            function(sample) sum(sample[[column]]),
            numeric(1)
        ))
        expect_equal(object$estimates$y5$absolute.prevalence[k], round(expected, 2))
    }
})

test_that("age_dead is applied at every index date", {
    object <- make_simulated_prevalence(
        index = c("2010-01-01", "2011-01-01", "2012-01-01"),
        years = 5,
        death_column = NULL,
        age_dead = 40,
        seed = 912
    )

    for (k in seq_along(object$index_dates)) {
        age_at_index <- object$simulated$age + as.numeric(
            object$index_dates[k] - object$simulated$incident_date
        ) / 365.25
        too_old <- age_at_index > 40
        expect_true(any(too_old))
        expect_true(all(object$simulated[[sprintf("alive_k%03d", k)]][too_old] == 0))
    }
})

test_that("multi-index goodness-of-fit values are keyed by index date", {
    object <- get_mixed_prevalence()

    expect_equal(length(object$pval), length(object$index_dates))
    expect_equal(names(object$pval), as.character(object$index_dates))
    expect_equal(object$pval, test_prevalence_fit(object))
})

test_that("fixed seeds reproduce multi-index simulation results", {
    first <- make_simulated_prevalence(seed = 31415)
    second <- make_simulated_prevalence(seed = 31415)

    expect_equal(first$index_dates, second$index_dates)
    expect_equal(first$estimates, second$estimates)
    expect_equal(first$simulated, second$simulated)
    expect_equal(first$pval, second$pval)
})

test_that("pre-fitted incidence and survival models support multiple dates", {
    data <- prevsim
    incidence_model <- rprev:::fit_exponential_incidence(entrydate ~ 1, data)
    survival_model <- rprev:::build_survreg(Surv(time, status) ~ age, data, "weibull")
    set.seed(804)
    object <- prevalence(
        index = c("2010-01-01", "2012-01-01"),
        num_years_to_estimate = 8,
        data = data,
        inc_model = incidence_model,
        surv_model = survival_model,
        registry_start_date = min(as.Date(data$entrydate)),
        incident_column = "entrydate",
        death_column = "eventdate",
        N_boot = 2
    )

    expect_s3_class(object, "prevalence")
    expect_equal(object$index_dates, as.Date(c("2010-01-01", "2012-01-01")))
    expect_equal(nrow(object$estimates$y8), 2)
    expect_false(is.null(object$simulated))
})

test_that("print displays every index date and estimation period", {
    object <- get_mixed_prevalence()
    output <- paste(capture.output(print(object)), collapse = "\n")

    expect_match(output, "Estimated prevalence")
    expect_match(output, "2010-01-01")
    expect_match(output, "2011-01-01")
    expect_match(output, "2012-01-01")
    expect_match(output, "2010-01-01\\s+5")
    expect_match(output, "2012-01-01\\s+8")
})

test_that("summary displays multi-index registry and simulation results", {
    object <- get_mixed_prevalence()
    output <- paste(capture.output(summary(object)), collapse = "\n")

    expect_match(output, "Index date\\(s\\): 2010-01-01, 2011-01-01, 2012-01-01")
    expect_match(output, "Counted prevalent cases")
    expect_match(output, "Simulation")
    expect_match(output, "Iterations: 2")
    expect_match(output, "P-values")
})

test_that("survfit supports prevalence objects with multiple index dates", {
    object <- get_mixed_prevalence()
    fitted <- survival::survfit(object)

    expect_s3_class(fitted, "survfit.prev")
    expect_equal(nrow(fitted$surv), object$N_boot)
    expect_equal(ncol(fitted$surv), length(fitted$time))
    expect_equal(length(fitted$fullsurv), length(fitted$time))
    expect_true(all(diff(as.numeric(fitted$fullsurv)) <= 0))
})
