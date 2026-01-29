library(rprev)

test_that("print.prevalence preserves K=1 format", {
    obj <- list(
        index_date=as.Date("2020-01-01"),
        est_years=c(5),
        estimates=list(y5=list(absolute.prevalence=10, per100000=2)),
        proportion=100000
    )
    class(obj) <- "prevalence"

    out <- capture.output(print(obj))
    expect_true(any(grepl("Estimated prevalence at 2020-01-01", out)))
    expect_true(any(grepl("5 years:", out)))
})

test_that("print.prevalence prints multi-index tables", {
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
        proportion=100000
    )
    class(obj) <- "prevalence"

    out <- capture.output(print(obj))
    expect_true(any(grepl("Estimated prevalence:", out)))
    expect_true(any(grepl("Prevalence rate", out)))
    expect_true(all(vapply(as.character(idx),
                           function(d) any(grepl(d, out)),
                           logical(1))))
})
