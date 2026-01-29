library(rprev)

test_that("new_point_estimate_multiindex handles sim-only estimates", {
    idx <- as.Date(c("2020-01-10", "2020-02-10"))
    sim_agg <- data.frame(
        sim=c(1, 2, 1, 2),
        index_date=rep(idx, each=2),
        year=1,
        contrib_total=c(4, 6, 2, 4),
        contrib_pre_registry=c(1, 1, 1, 1)
    )

    res <- rprev:::new_point_estimate_multiindex(year=1,
                                                 sim_agg=sim_agg,
                                                 index_dates=idx,
                                                 registry_data=data.frame(),
                                                 prev_formula=NULL,
                                                 registry_start_date=as.Date("2019-01-01"),
                                                 status_col="status")
    expect_equal(as.numeric(res$absolute.prevalence), c(5, 3))
})

test_that("new_point_estimate_multiindex combines counted and pre-registry sim", {
    idx <- as.Date(c("2020-01-10", "2020-02-10"))
    sim_agg <- data.frame(
        sim=c(1, 2, 1, 2),
        index_date=rep(idx, each=2),
        year=1,
        contrib_total=c(4, 6, 2, 4),
        contrib_pre_registry=c(2, 4, 0, 0)
    )
    registry_data <- data.frame(
        entry=as.Date(c("2019-02-10", "2019-03-01")),
        death=as.Date(c(NA, "2019-03-15")),
        status=c(0, 1)
    )

    res <- rprev:::new_point_estimate_multiindex(year=1,
                                                 sim_agg=sim_agg,
                                                 index_dates=idx,
                                                 registry_data=registry_data,
                                                 prev_formula=death ~ entry,
                                                 registry_start_date=as.Date("2019-02-01"),
                                                 status_col="status")
    expect_equal(as.numeric(res$absolute.prevalence), c(4, 1))
})

test_that("new_point_estimate_multiindex validates sim_agg columns", {
    idx <- as.Date(c("2020-01-10", "2020-02-10"))
    sim_bad <- data.frame(sim=1, index_date=idx[1], year=1, contrib_total=1)
    registry_data <- data.frame(
        entry=as.Date("2019-02-10"),
        death=as.Date(NA),
        status=0
    )

    expect_error(rprev:::new_point_estimate_multiindex(year=1,
                                                       sim_agg=sim_bad,
                                                       index_dates=idx,
                                                       registry_data=registry_data,
                                                       prev_formula=death ~ entry,
                                                       registry_start_date=as.Date("2019-02-01"),
                                                       status_col="status"))
})
