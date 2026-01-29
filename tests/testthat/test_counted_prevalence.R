library(rprev)

test_that("counted_prevalence returns scalar for single index", {
    df <- data.frame(
        entry = as.Date(c("2020-01-01", "2020-01-15", "2020-02-01", NA)),
        death = as.Date(c("2020-03-01", NA, "2020-02-15", "2020-02-01")),
        status = c(1, 0, 1, 0)
    )
    res <- rprev:::counted_prevalence(death ~ entry,
                                      index=as.Date("2020-02-01"),
                                      data=df,
                                      start_date=as.Date("2020-01-01"),
                                      status_col="status")
    expect_equal(res, 2)
})

test_that("counted_prevalence handles multi-index and start_date vectors", {
    df <- data.frame(
        entry = as.Date(c("2020-01-01", "2020-01-15", "2020-02-01", NA)),
        death = as.Date(c("2020-03-01", NA, "2020-02-15", "2020-02-01")),
        status = c(1, 0, 1, 0)
    )
    idx <- as.Date(c("2020-02-01", "2020-03-01"))
    res <- rprev:::counted_prevalence(death ~ entry,
                                      index=idx,
                                      data=df,
                                      start_date=as.Date("2020-01-01"),
                                      status_col="status")
    expect_equal(as.numeric(res), c(2, 1))
    expect_equal(names(res), as.character(idx))

    res2 <- rprev:::counted_prevalence(death ~ entry,
                                       index=idx,
                                       data=df,
                                       start_date=as.Date(c("2020-01-01", "2020-02-01")),
                                       status_col="status")
    expect_equal(as.numeric(res2), c(2, 0))
})

test_that("counted_prevalence is NA-robust", {
    df <- data.frame(
        entry = as.Date(c("2020-01-01", NA)),
        death = as.Date(c(NA, NA)),
        status = c(0, 0)
    )
    res <- rprev:::counted_prevalence(death ~ entry,
                                      index=as.Date(c("2020-02-01", "2020-03-01")),
                                      data=df,
                                      start_date=as.Date("2020-01-01"),
                                      status_col="status")
    expect_false(anyNA(res))
})
