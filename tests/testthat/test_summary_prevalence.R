library(rprev)

test_that("summary.prevalence preserves K=1 output", {
    obj <- list(
        index_date=as.Date("2020-01-01"),
        est_years=c(5),
        estimates=list(y5=list(absolute.prevalence=10, per100000=2)),
        proportion=100000,
        registry_start=as.Date("2019-01-01"),
        counted=10,
        counted_incidence_rate=0.1,
        simulated=NA,
        N_boot=10
    )
    class(obj) <- "prevalence"

    out <- capture.output(summary(obj))
    expect_true(any(grepl("Registry Data", out)))
    expect_true(any(grepl("Index date:", out)))
})

test_that("summary.prevalence prints multi-index tables", {
    idx <- as.Date(c("2020-01-01", "2020-02-01"))
    obj <- list(
        index_dates=idx,
        index_date=idx[2],
        est_years=c(3, 5),
        estimates=list(
            y3=data.frame(index_date=idx,
                          absolute.prevalence=c(1, 2),
                          per100000=c(10, 20)),
            y5=data.frame(index_date=idx,
                          absolute.prevalence=c(3, 4),
                          per100000=c(30, 40))
        ),
        proportion=100000,
        registry_start=as.Date("2019-01-01"),
        counted_by_index=data.frame(index_date=idx, counted=c(1, 2)),
        counted_incidence_rate=0.2,
        simulated=NA,
        N_boot=10,
        pval_by_index=setNames(c(0.1, 0.2), as.character(idx))
    )
    class(obj) <- "prevalence"

    out <- capture.output(summary(obj))
    expect_true(any(grepl("Index dates:", out)))
    expect_true(any(grepl("Counted prevalent cases:", out)))
    expect_true(all(vapply(as.character(idx),
                           function(d) any(grepl(d, out)),
                           logical(1))))
})
