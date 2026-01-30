library(rprev)
context("multiindex interface")

counted_df <- data.frame(
    incident = as.Date(c("2019-06-01", "2019-07-01", "2019-08-01")),
    death = as.Date(c(NA, "2020-01-15", NA)),
    status = c(0, 1, 0)
)

counted_args <- list(
    num_years_to_estimate = 1,
    data = counted_df,
    incident_column = "incident",
    death_column = "death",
    status_column = "status",
    registry_start_date = as.Date("2019-01-01")
)

test_that("index is deprecated alias and index_dates is canonical", {
    withr::with_seed(1, {
        expect_warning(
            res <- do.call(prevalence, c(list(index = "2020-03-01"), counted_args)),
            regexp = "deprecated|index_dates"
        )
        expect_equal(res$index_dates, as.Date("2020-03-01"))
        expect_false("index_date" %in% names(res))

        expect_warning(
            res2 <- do.call(prevalence, c(list(index = "2020-03-01", index_dates = "2020-03-01"), counted_args)),
            regexp = "index_dates"
        )
        expect_equal(res2$index_dates, as.Date("2020-03-01"))

        expect_error(
            do.call(prevalence, c(list(index = "2020-03-01", index_dates = "2020-04-01"), counted_args)),
            "do not match"
        )
    })
})

test_that("index_dates is sorted and unique", {
    withr::with_seed(1, {
        res <- do.call(prevalence, c(list(index_dates = c("2020-06-01", "2020-01-01", "2020-01-01")), counted_args))
        expect_equal(res$index_dates, as.Date(c("2020-01-01", "2020-06-01")))
        expect_false("index_date" %in% names(res))
    })
})

test_that("index and single index_dates give same result", {
    withr::with_seed(2, {
        res1 <- suppressWarnings(do.call(prevalence, c(list(index = "2020-03-01"), counted_args)))
        res2 <- do.call(prevalence, c(list(index_dates = "2020-03-01"), counted_args))
        expect_equal(res1$estimates, res2$estimates)
        expect_equal(res1$index_dates, res2$index_dates)
    })
})

test_that("counted-only path returns expected estimates without simulation", {
    withr::with_seed(3, {
        res <- do.call(prevalence, c(list(index_dates = "2020-03-01"), counted_args))
        expect_true(is.null(res$simulated))

        expected_registry <- rprev:::counted_prevalence(death ~ incident,
                                                        index_dates = as.Date("2020-03-01"),
                                                        data = counted_df,
                                                        start_date = as.Date("2019-01-01"),
                                                        status_col = "status")
        expected_year <- rprev:::counted_prevalence(death ~ incident,
                                                    index_dates = as.Date("2020-03-01"),
                                                    data = counted_df,
                                                    start_date = as.Date("2019-03-01"),
                                                    status_col = "status")
        expect_equal(res$counted, expected_registry)
        expect_equal(res$estimates$y1$absolute.prevalence, expected_year)
    })
})

test_that("multi-index simulation returns valid k_start/k_end", {
    data(prevsim)
    withr::with_seed(4, {
        res <- suppressWarnings(suppressMessages(prevalence(index_dates = c("2013-01-01", "2013-06-01", "2014-01-01"),
                                                            num_years_to_estimate = 5,
                                                            data = prevsim,
                                                            inc_formula = entrydate ~ sex,
                                                            surv_formula = Surv(time, status) ~ sex + age,
                                                            dist = "weibull",
                                                            death_column = NULL,
                                                            N_boot = 3)))
        sim <- res$simulated
        K <- length(res$index_dates)
        expect_true(all(c("k_start", "k_end") %in% names(sim)))
        expect_true(all(is.na(sim$k_start) | (sim$k_start >= 1 & sim$k_start <= K)))
        expect_true(all(is.na(sim$k_end) | (sim$k_end >= 1 & sim$k_end <= K)))
        ok <- is.na(sim$k_start) | is.na(sim$k_end) | (sim$k_start <= sim$k_end)
        expect_true(all(ok))
    })
})
