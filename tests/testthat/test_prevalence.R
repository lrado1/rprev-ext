library(rprev)
context('Prevalence')
data(prevsim)

test_that("prevalence outputs the same values as before", {
    # This function tests that the estimates haven't changed, not that the estimates or code are correct.
    # It's useful when large changes have been added that could affect how the prevalence estimates are produced.

    # Build models
    suppressWarnings(RNGversion("3.5.0"))
    set.seed(17)
    counted_nonstrat <- prevalence(index_dates="2013-01-01",
                                   num_years_to_estimate = c(5, 8),
                                   data=prevsim,
                                   inc_formula=entrydate ~ 1,
                                   surv_formula=Surv(time, status) ~ 1,
                                   dist='weibull',
                                   death_column='eventdate')
    wei_nonstrat <- prevalence(index_dates="2013-01-01",
                               num_years_to_estimate = 20,
                               data=prevsim,
                               inc_formula=entrydate ~ 1,
                               surv_formula=Surv(time, status) ~ 1,
                               dist='weibull',
                               death_column='eventdate',
                               N_boot = 10)
    wei_strat <- prevalence(index_dates="2013-01-01",
                            num_years_to_estimate = 20,
                            data=prevsim,
                            inc_formula=entrydate ~ sex,
                            surv_formula=Surv(time, status) ~ 1,
                            dist='weibull',
                            death_column='eventdate',
                            N_boot = 10)
    wei_agesex <- prevalence(index_dates="2013-01-01",
                             num_years_to_estimate = 20,
                             data=prevsim,
                             inc_formula=entrydate ~ sex,
                             surv_formula=Surv(time, status) ~ sex + age,
                             dist='weibull',
                             death_column='eventdate',
                             N_boot = 10)
    lnorm_age <- prevalence(index_dates="2013-01-01",
                            num_years_to_estimate = 17,
                            data=prevsim,
                            inc_formula=entrydate ~ sex,
                            surv_formula=Surv(time, status) ~ age,
                            dist='lognormal',
                            death_column='eventdate',
                            N_boot = 10)
    lnorm_sex <- prevalence(index_dates="2013-01-01",
                            num_years_to_estimate = 13,
                            data=prevsim,
                            inc_formula=entrydate ~ sex,
                            surv_formula=Surv(time, status) ~ sex,
                            dist='lognormal',
                            death_column='eventdate',
                            N_boot = 5)
    exp_full <- prevalence(index_dates="2013-01-01",
                           num_years_to_estimate = 13,
                           data=prevsim,
                           inc_formula=entrydate ~ sex,
                           surv_formula=Surv(time, status) ~ sex + age,
                           dist='exponential',
                           death_column='eventdate',
                           N_boot = 20)
    exp_sex <- prevalence(index_dates="2013-01-01",
                           num_years_to_estimate = 13,
                           data=prevsim,
                           inc_formula=entrydate ~ sex,
                           surv_formula=Surv(time, status) ~ sex,
                           dist='exponential',
                           death_column='eventdate',
                           N_boot = 15)


    # Compare estimates (absolute prevalence for K=1)
    expect_equal(c(counted_nonstrat$estimates$y5$absolute.prevalence,
                   counted_nonstrat$estimates$y8$absolute.prevalence),
                 c(305, 435))
    expect_equal(wei_nonstrat$estimates$y20$absolute.prevalence, c(755.5), tolerance = 0.5)
    expect_equal(wei_strat$estimates$y20$absolute.prevalence, c(769.9))
    expect_equal(wei_agesex$estimates$y20$absolute.prevalence, c(768.2))
    expect_equal(lnorm_age$estimates$y17$absolute.prevalence, c(742.2))
    expect_equal(lnorm_sex$estimates$y13$absolute.prevalence, c(630.4))
    expect_equal(exp_full$estimates$y13$absolute.prevalence, c(566.85), tolerance = 0.5)
    expect_equal(exp_sex$estimates$y13$absolute.prevalence, c(566.87))

    # Compare # simulated individuals
    expect_equal(nrow(wei_nonstrat$simulated), 20019)
    expect_equal(nrow(wei_strat$simulated), 19767)
    expect_equal(nrow(wei_agesex$simulated), 19889)
    expect_equal(nrow(lnorm_age$simulated), 17066)
    expect_equal(nrow(lnorm_sex$simulated), 6306)
    expect_equal(nrow(exp_full$simulated), 25928)
    expect_equal(nrow(exp_sex$simulated), 19333)

    # Compare p-vals
    expect_equal(round(wei_nonstrat$pval, 3), 0.431)
    expect_equal(round(wei_strat$pval, 3), 0.282)
    expect_equal(round(wei_agesex$pval, 3), 0.312)
    expect_equal(round(lnorm_age$pval, 3), 0.594)
    expect_equal(round(lnorm_sex$pval, 3), 0.216)
    expect_equal(round(exp_full$pval, 3), 0.469)
    expect_equal(round(exp_sex$pval, 3), 0.254)
})

test_that("prevalence function handles incorrectly specified inputs", {
    # Index date isn't formatted as date
    expect_error(prevalence(index_dates="201311",
                            num_years_to_estimate = 13,
                            data=prevsim,
                            inc_formula=entrydate ~ sex,
                            surv_formula=Surv(time, status) ~ sex,
                            dist='exponential',
                            death_column='eventdate',
                            N_boot = 15))

    # Incorrect specification of num_years_to_estimate
    expect_warning(prevalence(index_dates="20130101",
                              num_years_to_estimate = -10,
                              data=prevsim,
                              inc_formula=entrydate ~ sex,
                              surv_formula=Surv(time, status) ~ sex,
                              dist='exponential',
                              death_column='eventdate',
                              N_boot = 15))

    # Incidence formula containing variables that aren't in the data set
    expect_error(prevalence(index_dates="20130101",
                            num_years_to_estimate = 15,
                            data=prevsim,
                            inc_formula=entry ~ 1,
                            surv_formula=Surv(time, status) ~ sex,
                            dist='exponential',
                            death_column='eventdate',
                            N_boot = 15))
    expect_error(prevalence(index_dates="20130101",
                            num_years_to_estimate = 15,
                            data=prevsim,
                            inc_formula=entrydate ~ Sex,
                            surv_formula=Surv(time, status) ~ sex,
                            dist='exponential',
                            death_column='eventdate',
                            N_boot = 15))

    # Survival formula containing variables that aren't in the data set
    expect_error(prevalence(index_dates="20130101",
                            num_years_to_estimate = 15,
                            data=prevsim,
                            inc_formula=entrydate ~ 1,
                            surv_formula=Surv(stime, status) ~ 1,
                            dist='exponential',
                            death_column='eventdate',
                            N_boot = 15))
    expect_error(prevalence(index_dates="20130101",
                            num_years_to_estimate = 15,
                            data=prevsim,
                            inc_formula=entrydate ~ sex,
                            surv_formula=Surv(time, stat) ~ 1,
                            dist='exponential',
                            death_column='eventdate',
                            N_boot = 15))
    expect_error(prevalence(index_dates="20130101",
                            num_years_to_estimate = 15,
                            data=prevsim,
                            inc_formula=entrydate ~ sex,
                            surv_formula=Surv(time, status) ~ Sex,
                            dist='exponential',
                            death_column='eventdate',
                            N_boot = 15))

    # Distribution isn't in accepted list
    expect_error(prevalence(index_dates="20130101",
                            num_years_to_estimate = 15,
                            data=prevsim,
                            inc_formula=entrydate ~ sex,
                            surv_formula=Surv(time, status) ~ sex,
                            dist='exponental',
                            death_column='eventdate',
                            N_boot = 15))
    expect_error(prevalence(index_dates="20130101",
                            num_years_to_estimate = 15,
                            data=prevsim,
                            inc_formula=entrydate ~ sex,
                            surv_formula=Surv(time, status) ~ sex,
                            dist='gompertz',
                            death_column='eventdate',
                            N_boot = 15))
    expect_error(prevalence(index_dates="20130101",
                            num_years_to_estimate = 15,
                            data=prevsim,
                            inc_formula=entrydate ~ sex,
                            surv_formula=Surv(time, status) ~ sex,
                            dist='lnorm',
                            death_column='eventdate',
                            N_boot = 15))

    # Death column incorrectly specified
    expect_warning(prevalence(index_dates="20130101",
                              num_years_to_estimate = 15,
                              data=prevsim,
                              inc_formula=entrydate ~ sex,
                              surv_formula=Surv(time, status) ~ sex,
                              dist='weibull',
                              death_column='event',
                              N_boot = 15))
})

test_that("build_prev_counts_multiindex respects strict incident window", {
    index_dates <- as.Date(c("2020-01-01", "2020-02-01", "2020-03-01"))
    results <- data.table::data.table(sim = 1L,
                                      incident_date = as.Date("2020-01-15"),
                                      k_start = 1L,
                                      k_end = 3L)

    counts <- rprev:::build_prev_counts_multiindex(results, index_dates, years = 1)
    prev_vec <- counts[sim == 1 & year == 1][order(k)]$prev_count

    expect_equal(prev_vec, c(0, 1, 1))
})
