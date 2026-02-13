# Synthetic registry data generation helpers (yearly incidence specification).
#
# This script implements the updated synthetic-data logic where incidence is
# generated year-by-year using:
# - deterministic trend: lambda_bar_y = lambda0 * (1 + gamma)^(y - Y_min)
# - multiplicative log-normal uncertainty
# - Poisson yearly counts
#
# It also includes a prevalence counting helper for a selected index date.
# No example run is included in this file; usage is intended from a notebook.

options(stringsAsFactors = FALSE)

# Convert arbitrary date-like input to Date.
synth_as_date <- function(x) {
  as.Date(x)
}

# Convert Date(s) to integer day indices from origin date.
synth_date_to_day_index <- function(date_value, t_orig) {
  as.integer(synth_as_date(date_value) - synth_as_date(t_orig))
}

# Convert integer day indices back to Date values from origin date.
synth_day_index_to_date <- function(day_index, t_orig) {
  as.Date(as.integer(day_index), origin = synth_as_date(t_orig))
}

# Build validated time-axis configuration.
synth_prepare_time_scale <- function(t_orig, t_min, t_max, t_end) {
  t_orig <- synth_as_date(t_orig)
  t_min <- synth_as_date(t_min)
  t_max <- synth_as_date(t_max)
  t_end <- synth_as_date(t_end)

  if (t_min >= t_max) stop("t_min must be earlier than t_max.")
  if (t_end < t_max) stop("t_end must be on or after t_max.")

  list(
    t_orig = t_orig,
    t_min = t_min,
    t_max = t_max,
    t_end = t_end,
    t_min_idx = as.integer(t_min - t_orig),
    t_max_idx = as.integer(t_max - t_orig),
    t_end_idx = as.integer(t_end - t_orig)
  )
}

# Build calendar-year windows for [t_min, t_max).
#
# strict_calendar_bounds = TRUE enforces exact year boundaries:
# - t_min must be YYYY-01-01
# - t_max must be YYYY-01-01
synth_build_year_windows <- function(time_cfg, strict_calendar_bounds = TRUE) {
  if (isTRUE(strict_calendar_bounds)) {
    if (format(time_cfg$t_min, "%m-%d") != "01-01") {
      stop("With strict_calendar_bounds=TRUE, t_min must be January 1.")
    }
    if (format(time_cfg$t_max, "%m-%d") != "01-01") {
      stop("With strict_calendar_bounds=TRUE, t_max must be January 1.")
    }
  }

  t_max_minus_one <- time_cfg$t_max - 1L
  y_min <- as.integer(format(time_cfg$t_min, "%Y"))
  y_max <- as.integer(format(t_max_minus_one, "%Y"))
  years <- seq.int(y_min, y_max)

  rows <- vector("list", length(years))
  k <- 0L

  for (y in years) {
    year_start <- as.Date(sprintf("%04d-01-01", y))
    next_start <- as.Date(sprintf("%04d-01-01", y + 1L))

    a_date <- max(year_start, time_cfg$t_min)
    b_date <- min(next_start, time_cfg$t_max)

    if (a_date >= b_date) next

    if (isTRUE(strict_calendar_bounds) && (a_date != year_start || b_date != next_start)) {
      stop("Year windows are not full calendar years under strict mode.")
    }

    k <- k + 1L
    rows[[k]] <- data.frame(
      year = y,
      a_idx = synth_date_to_day_index(a_date, time_cfg$t_orig),
      b_idx = synth_date_to_day_index(b_date, time_cfg$t_orig),
      a_date = a_date,
      b_date = b_date
    )
  }

  if (k == 0L) {
    return(data.frame(
      year = integer(0),
      a_idx = integer(0),
      b_idx = integer(0),
      a_date = as.Date(character(0)),
      b_date = as.Date(character(0))
    ))
  }

  do.call(rbind, rows[seq_len(k)])
}

# Sample yearly incidence and diagnosis-day indices under the updated model.
synth_sample_yearly_incidence <- function(year_windows, inc_lam0, inc_g, inc_sdlog) {
  if (inc_lam0 <= 0) stop("inc_lam0 must be > 0.")
  if (inc_sdlog < 0) stop("inc_sdlog must be >= 0.")
  if (1 + inc_g <= 0) stop("inc_g must satisfy 1 + inc_g > 0.")

  if (nrow(year_windows) == 0L) {
    return(list(
      diagnosis_day = integer(0),
      yearly_table = data.frame(
        year = integer(0),
        lambda_bar = numeric(0),
        u_year = numeric(0),
        lambda_year = numeric(0),
        n_year = integer(0)
      )
    ))
  }

  y0 <- min(year_windows$year)
  year_offset <- year_windows$year - y0

  lambda_bar <- inc_lam0 * (1 + inc_g) ^ year_offset
  u_year <- exp(rnorm(
    n = nrow(year_windows),
    mean = -0.5 * inc_sdlog ^ 2,
    sd = inc_sdlog
  ))
  lambda_year <- lambda_bar * u_year
  n_year <- rpois(n = nrow(year_windows), lambda = lambda_year)

  total_n <- sum(n_year)
  diagnosis_day <- integer(total_n)

  if (total_n > 0L) {
    pos <- 1L
    for (i in seq_len(nrow(year_windows))) {
      ni <- n_year[i]
      if (ni == 0L) next

      support <- year_windows$a_idx[i]:(year_windows$b_idx[i] - 1L)
      draws <- sample(support, size = ni, replace = TRUE)
      diagnosis_day[pos:(pos + ni - 1L)] <- draws
      pos <- pos + ni
    }

    diagnosis_day <- sample(diagnosis_day, size = length(diagnosis_day), replace = FALSE)
  }

  yearly_table <- data.frame(
    year = year_windows$year,
    lambda_bar = lambda_bar,
    u_year = u_year,
    lambda_year = lambda_year,
    n_year = n_year
  )

  list(
    diagnosis_day = diagnosis_day,
    yearly_table = yearly_table
  )
}

# Sample covariates X1, X2, X3 according to the specification.
synth_sample_covariates <- function(
  n,
  p_x1,
  beta_a_x1 = c(2.0, 5.0),
  beta_b_x1 = c(5.0, 2.0),
  x2_range = c(20, 85),
  x3_mu = c(0.0, 0.6, 0.03),
  x3_sd = 1.0
) {
  if (n < 0L) stop("n must be >= 0.")
  if (p_x1 <= 0 || p_x1 >= 1) stop("p_x1 must be in (0,1).")
  if (length(beta_a_x1) != 2 || length(beta_b_x1) != 2) {
    stop("beta_a_x1 and beta_b_x1 must have length 2 (for x1=0,1).")
  }
  if (x2_range[1] >= x2_range[2]) stop("x2_range must satisfy min < max.")
  if (length(x3_mu) != 3) stop("x3_mu must have length 3: (eta0, eta1, eta2).")
  if (x3_sd <= 0) stop("x3_sd must be > 0.")

  x1 <- rbinom(n, size = 1, prob = p_x1)
  u <- numeric(n)

  idx0 <- which(x1 == 0L)
  idx1 <- which(x1 == 1L)
  if (length(idx0) > 0L) u[idx0] <- rbeta(length(idx0), beta_a_x1[1], beta_b_x1[1])
  if (length(idx1) > 0L) u[idx1] <- rbeta(length(idx1), beta_a_x1[2], beta_b_x1[2])

  x2 <- x2_range[1] + (x2_range[2] - x2_range[1]) * u
  mean_x3 <- x3_mu[1] + x3_mu[2] * x1 + x3_mu[3] * x2
  x3 <- rnorm(n, mean = mean_x3, sd = x3_sd)

  data.frame(
    x1 = x1,
    x2 = x2,
    x3 = x3
  )
}

# Draw from Prentice-style generalized gamma parameterization.
synth_rgengamma_prentice <- function(n, mu, sigma, gg_q) {
  if (sigma <= 0) stop("sigma must be > 0.")
  if (n <= 0L) return(numeric(0))
  mu <- rep(mu, length.out = n)

  if (abs(gg_q) < sqrt(.Machine$double.eps)) {
    return(exp(rnorm(n, mean = mu, sd = sigma)))
  }

  k <- 1 / (gg_q ^ 2)
  y <- rgamma(n, shape = k, rate = k)
  exp(mu + (sigma / gg_q) * log(y))
}

# Sample latent event times T (continuous and day-level).
synth_sample_latent_event_time <- function(covariates, gg_lam, gg_k, gg_b, gg_q) {
  if (gg_lam <= 0) stop("gg_lam must be > 0.")
  if (gg_k <= 0) stop("gg_k must be > 0.")
  if (length(gg_b) != 3) stop("gg_b must have length 3 for x1,x2,x3.")

  n <- nrow(covariates)
  x_mat <- as.matrix(covariates[, c("x1", "x2", "x3")])
  linpred <- as.vector(x_mat %*% gg_b)

  mu <- log(gg_lam) - linpred / gg_k
  sigma <- 1 / gg_k
  t_cont <- synth_rgengamma_prentice(n = n, mu = mu, sigma = sigma, gg_q = gg_q)
  t_day <- as.integer(floor(pmax(t_cont, 0)))

  list(
    T_cont = t_cont,
    T_day = t_day,
    mu = mu,
    sigma = sigma
  )
}

# Sample LTFU censoring time uniformly before min(T, C_adm).
synth_sample_pre_event_ltfu_time <- function(T_day, C_adm_day) {
  upper <- pmin(T_day, C_adm_day)
  out <- rep(Inf, length(upper))

  eligible <- which(upper > 0L)
  if (length(eligible) > 0L) {
    out[eligible] <- vapply(
      upper[eligible],
      function(u) sample.int(as.integer(u), size = 1L) - 1L,
      FUN.VALUE = integer(1)
    )
  }

  out
}

# Apply administrative and random-LTFU censoring.
synth_apply_censoring <- function(D_day, T_day, t_end_idx, p_ltfu) {
  if (p_ltfu < 0 || p_ltfu > 1) stop("p_ltfu must be in [0,1].")

  n <- length(D_day)
  C_adm_day <- pmax(as.integer(t_end_idx - D_day), 0L)

  L <- rbinom(n, size = 1, prob = p_ltfu)
  C_ltfu_day <- rep(Inf, n)
  idx_ltfu <- which(L == 1L)
  if (length(idx_ltfu) > 0L) {
    C_ltfu_day[idx_ltfu] <- synth_sample_pre_event_ltfu_time(
      T_day = T_day[idx_ltfu],
      C_adm_day = C_adm_day[idx_ltfu]
    )
  }

  C_day <- pmin(C_adm_day, C_ltfu_day)
  Y_day <- pmin(T_day, C_day)
  delta <- as.integer(T_day <= C_day)
  F_day <- as.integer(D_day + Y_day)

  list(
    C_adm_day = C_adm_day,
    C_ltfu_day = C_ltfu_day,
    C_day = C_day,
    Y_day = as.integer(Y_day),
    delta = delta,
    F_day = F_day,
    L = L
  )
}

# Main synthetic registry generator (yearly incidence trend version).
generate_synthetic_registry_yearly <- function(
  t_orig = "2000-01-01",
  t_min = "2010-01-01",
  t_max = "2015-01-01",
  t_end = "2016-12-31",
  inc_lam0 = 900,
  inc_g = 0.03,
  inc_sdlog = 0.20,
  strict_calendar_bounds = TRUE,
  p_x1 = 0.45,
  beta_a_x1 = c(2.0, 5.0),
  beta_b_x1 = c(5.0, 2.0),
  x2_range = c(20, 85),
  x3_mu = c(0.0, 0.6, 0.03),
  x3_sd = 1.0,
  gg_lam = 1800,
  gg_k = 1.2,
  gg_b = c(-0.5, 0.015, 0.6),
  gg_q = 1.0,
  p_ltfu = 0.25,
  seed = 20260213
) {
  set.seed(seed)

  time_cfg <- synth_prepare_time_scale(
    t_orig = t_orig,
    t_min = t_min,
    t_max = t_max,
    t_end = t_end
  )

  year_windows <- synth_build_year_windows(
    time_cfg = time_cfg,
    strict_calendar_bounds = strict_calendar_bounds
  )

  incidence <- synth_sample_yearly_incidence(
    year_windows = year_windows,
    inc_lam0 = inc_lam0,
    inc_g = inc_g,
    inc_sdlog = inc_sdlog
  )

  D_day <- incidence$diagnosis_day
  N <- length(D_day)

  covars <- synth_sample_covariates(
    n = N,
    p_x1 = p_x1,
    beta_a_x1 = beta_a_x1,
    beta_b_x1 = beta_b_x1,
    x2_range = x2_range,
    x3_mu = x3_mu,
    x3_sd = x3_sd
  )

  latent <- synth_sample_latent_event_time(
    covariates = covars,
    gg_lam = gg_lam,
    gg_k = gg_k,
    gg_b = gg_b,
    gg_q = gg_q
  )

  T_day <- latent$T_day
  F_true_day <- D_day + T_day
  delta_true <- rep(1L, N)

  observed <- synth_apply_censoring(
    D_day = D_day,
    T_day = T_day,
    t_end_idx = time_cfg$t_end_idx,
    p_ltfu = p_ltfu
  )

  out <- data.frame(
    case_id_id = seq_len(N),
    diagnosis_date_D = synth_day_index_to_date(D_day, time_cfg$t_orig),
    covariate_binary_x1 = covars$x1,
    covariate_scaled_x2 = covars$x2,
    covariate_hidden_x3 = covars$x3,
    last_followup_date_F = synth_day_index_to_date(observed$F_day, time_cfg$t_orig),
    event_observed_delta = observed$delta,
    true_event_date_F_true = synth_day_index_to_date(F_true_day, time_cfg$t_orig),
    true_event_indicator_delta_true = delta_true
  )

  attr(out, "time_config") <- time_cfg
  attr(out, "year_windows") <- year_windows
  attr(out, "yearly_incidence") <- incidence$yearly_table
  attr(out, "seed") <- seed

  out
}

# Count prevalent cases at one or more selected index dates:
# - truth-based count from F_true
# - observed count from F
count_prevalence_at_index <- function(df, index_date, population_size = NULL, include_index_day = FALSE) {
  required_cols <- c(
    "diagnosis_date_D",
    "last_followup_date_F",
    "true_event_date_F_true"
  )
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  idx_vec <- synth_as_date(index_date)
  if (length(idx_vec) == 0L) {
    stop("index_date must contain at least one date.")
  }
  if (any(is.na(idx_vec))) {
    stop("index_date contains invalid date value(s).")
  }

  pop_vec <- rep(NA_real_, length(idx_vec))
  if (!is.null(population_size)) {
    if (!is.numeric(population_size) || any(!is.finite(population_size))) {
      stop("population_size must be a finite numeric value.")
    }
    if (any(population_size <= 0)) {
      stop("population_size must be > 0.")
    }

    if (length(population_size) == 1L) {
      pop_vec <- rep(as.numeric(population_size), length(idx_vec))
    } else if (length(population_size) == length(idx_vec)) {
      pop_vec <- as.numeric(population_size)
    } else {
      stop("population_size must be length 1 or the same length as index_date.")
    }
  }

  cmp <- if (isTRUE(include_index_day)) {
    function(x, y) x >= y
  } else {
    function(x, y) x > y
  }

  summarize_one_index <- function(idx, pop) {
    diagnosed <- df$diagnosis_date_D <= idx
    n_diagnosed <- sum(diagnosed)

    prevalent_true <- diagnosed & cmp(df$true_event_date_F_true, idx)
    prevalent_observed <- diagnosed & cmp(df$last_followup_date_F, idx)

    n_true <- sum(prevalent_true)
    n_observed <- sum(prevalent_observed)

    prev_rate <- if (is.na(pop)) {
      rep(NA_real_, 2)
    } else {
      c(n_true / pop, n_observed / pop)
    }

    data.frame(
      index_date = rep(idx, 2),
      followup_basis = c("true_F", "observed_F"),
      prevalent_cases = c(n_true, n_observed),
      diagnosed_cases = rep(n_diagnosed, 2),
      prevalence_rate = prev_rate
    )
  }

  result_list <- Map(summarize_one_index, idx_vec, pop_vec)
  result <- do.call(rbind, result_list)
  rownames(result) <- NULL
  result
}
