
# Get minimum of a vector, return Inf if vector is null or all elements are NA
min_na_inf <- function(x){
  if(all(is.na(x)) | is.null(x)){
  # if(is.null(x)){
    return(Inf)
  } else {
    return(min(x, na.rm = TRUE))
  }
}


# Extract first time where status == k by time t
# output Inf if no such time exists
# Argument: df = tibble(time, status)
extact_time_with_status <- function(time, status, k, t = Inf){
  time[status == k & time <= t] |> min_na_inf()
}


# df1 <- df_y |>
#   filter(id == uid[1])
# df2 <- df_y |>
#   filter(id == uid[2])
#
# t <- 220
#
# wp(df1, df2, deathcode = 1, nfcode = 2, t = Inf)

# Transform outcomes data by standard win function (Pocock's setting)
poc_win <- function(df_y, n){

  df_y_wide <- df_y |>
    group_by(id) |>
    reframe(
      results = long_to_wide_death_nf(time, status)
    ) |>
    mutate(
      names = rep(c("death_time", "nonfatal_time", "fu_time"), n)
    ) |>
    pivot_wider(
      names_from = names,
      values_from = results
    )

  # Generate all pairwise combinations
  pairwise_tibble <- crossing(i = df_y_wide, j = df_y_wide) %>%
    filter(i |> row_number() < j |> row_number()) |>
    unnest(c(i, j), names_sep = "_")


  # Compute pairwise wins and losses
  pw_win <- pairwise_tibble |>
    mutate(
      win_death = case_when(
        j_death_time < pmin(i_death_time, i_fu_time) ~ 1,
        i_death_time < pmin(j_death_time, j_fu_time) ~ -1,
        TRUE ~ 0
      ),
      win_nonfatal = if_else(win_death == 0, case_when(
        j_nonfatal_time < pmin(i_nonfatal_time, i_fu_time) ~ 1,
        i_nonfatal_time < pmin(j_nonfatal_time, j_fu_time) ~ -1,
        TRUE ~ 0
      ), 0),
      win = win_death + win_nonfatal,
      component = case_when(
        win_death != 0 ~ "death",
        win_nonfatal != 0 ~ "nonfatal",
        TRUE ~ "tie"
      )
    ) |>
    select(
      id_i = i_id,
      id_j = j_id,
      win,
      component
    ) |>
    filter(win != 0) # Remove pairs with win == 0

  # Remove pairs with win == 0 in pairwise tibble as well
  # by left joining with
  pairwise_tibble <- pw_win |> select(id_i, id_j) |>
    left_join(pairwise_tibble, by = c("id_i" = "i_id", "id_j" = "j_id"))

  list(pairwise_tibble = pairwise_tibble, pw_win = pw_win)
}

# wtest <- function(df1, df2, deathcode = 1, nfcode = 2, t = Inf) {
#   list(win = 1, component = "tie")
# }

# Reorganize long outcomes (death, nonfatal, censor) to wide

long_to_wide_death_nf <- function(time, status, deathcode = 1, nfcode = 2, t = Inf){

 fu_time <- min(max(time), t)
 death_time <- extact_time_with_status(time, status, deathcode, fu_time)
 nonfatal_time <- extact_time_with_status(time, status, nfcode, fu_time)

 obj <- list(
   death_time = death_time,
   nonfatal_time = nonfatal_time,
   fu_time = fu_time
 )
  # return(obj)
 return(c(death_time, nonfatal_time, fu_time))

}

# Compute the derivative of score function
pw_score_fun <- function(beta, delta_pm, Zd, N){

  # Exp(beta^T Zd)
  exp_beta_Zd <- exp(Zd %*% beta)
  mu <- as.vector(exp_beta_Zd / (1 + exp_beta_Zd))
  # Compute score function
  M <- pmax(delta_pm, 0) - abs(delta_pm) * mu
  # M <- (delta_pm == 1) - abs(delta_pm) * mu
  score_matrix <- Zd * M
  score <- colSums(score_matrix) / N
  # Compute its derivative A
  A <- - t(Zd * mu) %*% (Zd * (1 - mu)) / N

  return(list(score = score, A = A, mu = mu, M = M, score_matrix = score_matrix))
}

newton_pw <- function(pw_df, N, p, beta = NULL, maxiter = 50, eps = 1e-5){

  # Signed win vector
  delta_pm <- pw_df |>
    select(win) |>
    pull()

  # Matrix of Zd (N x p)
  Zd <- pw_df |> select(5: (4 + p)) |> as.matrix()


  # Initialize beta
  if(is.null(beta)){
    beta <- rep(0, ncol(Zd))
  }
  # Initialize score function
  obj <- pw_score_fun(beta, delta_pm, Zd, N)
  score <- obj$score
  A <- obj$A
  # Initialize iteration counter
  iter <- 0
  # Initialize convergence indicator
  conv <- FALSE
  # Newton-Raphson iteration
  while(!conv & iter < maxiter){
    # Update beta
    beta <- beta - solve(A) %*% score
    # Update score function
    obj <- pw_score_fun(beta, delta_pm, Zd, N)
    score <- obj$score
    A <- obj$A
    # Update iteration counter
    iter <- iter + 1
    # Check convergence
    conv <- max(abs(score)) < eps
  }
  # Return the estimated beta
  beta <- beta[, 1]
  # If conv == FALSE, print a warning message
  if(!conv){
    warning("Newton-Raphson algorithm did not converge.")
  }

  return(list(beta = beta, conv = conv, score = score, A = A, iter = iter,
              mu = obj$mu, M = obj$M, score_matrix = obj$score_matrix))
}


# Function to get individualized quantities
# delta, R, Lambda, M, and scores by projection

proj_scores <- function(pw_df, score_matrix, mu, M, n, p){

  # Combine score matrix with win loss id's
  score_tibble <- tibble(
    pw_df |> select(id_i, id_j, win), mu, M,
    as_tibble(score_matrix)
  )

  # ncol(score_tibble)
  # Function to sum and divide by (n - 1)
  proj_n <- function(x){
    sum(x) / (n - 1)
  }

  # Combine score matrix with win loss id's
  score_tibble <- tibble(
    pw_df |> select(id_i, id_j, win), mu, M,
    as_tibble(score_matrix)
  )

  # ncol(score_tibble)
  # Function to sum and divide by (n - 1)
  proj_n <- function(x){
    sum(x) / (n - 1)
  }


  # Score is symmetric w.r.t id_i and id_j
  ## Summarize by id_i
  score_avg_i <- score_tibble |>
    mutate(
      delta = pmax(win, 0),
      R = abs(win), # should all be 1
      Lambda = R * mu,
      M = M,
      .after = id_i,
    ) |>
    select(-c(id_j, win, mu)) |> # include win, scores
    group_by(id_i) |>
    summarize(
      across(everything(), proj_n)
    ) |>
    rename(id = id_i)

  ## Summarize by id_j
  score_avg_j <- score_tibble |>
    mutate(
      delta = pmax(-win, 0),
      R = abs(win), # should all be 1
      Lambda = R * (1 - mu),
      M = - M,
      .after = id_j,
    ) |>
    select(-c(id_i, win, mu)) |> # include win, scores
    group_by(id_j) |>
    summarize(
      across(everything(), proj_n)
    ) |>
    rename(id = id_j)

  # Combine score_sum_i and score_sum_j
  score_avg_i |>
    bind_rows(score_avg_j) |>
    group_by(id) |>
    summarize(
      across(everything(), sum)
    ) |>
    arrange(id)

}





# A test version of pwreg (for v1.1)

#' A test version of pwreg
#'
#' This function is a test version of the \code{pwreg} function.
#' It is intended to simplify the code and make computations more efficient.
#'
#' @param id A vector of subject identifiers.
#' @param time A vector of event times.
#' @param status A vector of event indicators.
#' @param Z A matrix or data frame (tibble) of covariates.
#' @param wfun A function for the win function (Default is standard Pocock's win function).
#' @param Y An optional matrix or data frame (tibble) of continuous response variables.
#' @param strata A vector of strata indicators.
#' @param fixedL A logical value indicating whether the frailty variance is fixed.
#' @param eps A numeric value for the convergence criterion.
#' @param maxiter A numeric value for the maximum number of iterations.
#' @export
#' @import dplyr tidyr stringr
pwreg1 <- function(id, time, status, Z, wfun = NULL, Y = NULL, strata = NULL, fixedL = TRUE, eps = 1e-6,
                   maxiter = 50){


# Check input format and organize them into tibble ------------------------

 # Convert Z to tibble
  Z <- as_tibble(Z)
 # Check the input data
    if (is.null(Y)){
      if(length(id) != length(time) | length(id) != length(status) | length(id) != nrow(Z)){
        stop("The lengths of id, time, status, and Z must be the same.")
        }
        # number of Y columns
        q <- 0
      } else{
        # Convert Y to a tibble and get column number
        Y <- as_tibble(Y)
        # q is dimension of continuous response
        q <- ncol(Y)
      if(length(id) != length(time) | length(id) != length(status) | length(id) != nrow(Y) | length(id) != nrow(Z)){
        stop("The lengths of id, time, status, Y, and Z must be the same.")
        }
      }

  # Organize data into a tibble
    df <- tibble(
      id, time, status, Y, Z
     ) |>
      arrange(
        id, time
      ) |>
      drop_na()
  # Get attributes of the data
    # Unique patient id
    uid <- df |> distinct(id) |> pull(id)
    # Number of patients
    n <- length(uid)
    # Number of pairs
    N <- n * (n - 1) / 2
  # De-duplicate covariate matrix
    # Zn is tibble with id column and p covariate columns
    Zn <- df |>
      group_by(id) |>
      slice(1) |>
      select(-(2: (3 + q))) |>
      ungroup()
    # Number of covariates
    p <- ncol(Zn) - 1
  # Outcomes data (long format)
    df_y <- df[, 1 : (3 + q)]


# Flatten outcome data and self-pair --------------------------------------
# Specific to win function wfun
# df_y -> pairwise_tibble, pw_win
   if (is.null(wfun)){
     poc_win_obj <- poc_win(df_y, n)
    } else {
      stop("wfun is not implemented yet.")
    }
  # Extract pairwise data
    pairwise_tibble <- poc_win_obj$pairwise_tibble
    pw_win <- poc_win_obj$pw_win
  # Self-pair covariates
    pw_Z <- crossing(i = Zn, j = Zn) %>%
      filter(i |> row_number() < j |> row_number()) |>
      unnest(c(i, j), names_sep = "_") |>
      # Right joint with pw_win to remove pairs with win ==0
      right_join(pw_win |> select(id_i, id_j),
                 by = c("i_id" = "id_i", "j_id" = "id_j"))
  # Compute difference in covariates
    pw_Zd <- pw_Z |>
      select(id_i = i_id, id_j = j_id) |>
      bind_cols(
        # Z_i - Z_j
        pw_Z |> select(starts_with("i_")) |> rename_with(~ str_remove(., "i_")) |> select(-id) -
          pw_Z |> select(starts_with("j_")) |> rename_with(~ str_remove(., "j_")) |> select(-id)
      )

    # Merge Zd with pw_win by id_i and id_j
    pw_df <- pw_win |>
      left_join(pw_Zd, by = c("id_i", "id_j"))




    # Newton-Raphson iteration to estimate beta
    newton_obj <- newton_pw(pw_df, N, p, maxiter = maxiter, eps = eps)
    # Extract elements
    beta <- newton_obj$beta # fitted bete
    conv <- newton_obj$conv # convergence indicator
    iter <- newton_obj$iter # number of iterations
    A <- newton_obj$A # estimated A matrix
    score_matrix <- newton_obj$score_matrix # score matrix
    mu <- newton_obj$mu # estimated mu
    M <- newton_obj$M # estimated win residual

# Variance calculation ----------------------------------------------------

    ## check if M = delta - Lambda
    # score_sum |>
    #   mutate(M_check = delta - Lambda) |>
    #   select(id, M, M_check) |>
    #   summarize(error = sum(abs(M - M_check)))

    score_sum <- proj_scores(pw_df, score_matrix, mu, M, n, p)

    # Kappa matrix (n x p) with influence function of the score
    kappa_matrix <- score_sum |>
      select(- c(id, delta,   R, Lambda, M)) |>
      as.matrix()
    # Compute A inverse
    A_inv <- solve(A)
    # Influence function of beta
    psi_matrix <- - 2 * A_inv %*% t(kappa_matrix) # p x n
    Var <- psi_matrix %*% t(psi_matrix) / (n * (n - p - 1))

    ## Output id-based influence function of beta
    psi_df_wl <- score_sum |> select(id) |> bind_cols(as_tibble(t(psi_matrix)))

    psi_df <- tibble(
      id = uid
      ) |>
      left_join(
        psi_df_wl,
        by = c("id" = "id")
      )|>
      mutate(
        across(-id, ~ replace_na(., 0))
      )

    psi_matrix <- psi_df |> select(-id) |> as.matrix() |> t()
    # Take projected, individualized quantities: delta, R, Lambda,  M, score residuals
    resids_n <- score_sum

    obj <- list(beta = beta, Var = Var, conv = conv, iter = iter, n = n, N = N, # basic model fit
                resids_n = resids_n, Zn = Zn,  # individualized quantities (w. id column)
                A = A, kappa_matrix = kappa_matrix, psi_matrix = psi_matrix, # influence-related quantities
                pairwise_tibble = pairwise_tibble,  pw_win = pw_win,# pairwise data for possible re-calculation of win-loss (by t)
                df = df, # original data
                call = match.call(), wfun = wfun, fixedL = fixedL # call and control parameters
                )

  class(obj) <- "pwreg1"

  return(obj)
}


#' @export
residuals.pwreg1 <- function(x, ...) {

  call <- x$call
  # Contains delta, R, Lambda, M (win residuals)
  resids_n <- x$resids_n
  # Contains Zn
  Zn <- x$Zn
  p <- ncol(Zn) - 1
  n <- nrow(Zn)
  resids_n <- x$resids_n |>
    right_join(Zn |> select(id), by = "id") |>
    mutate(
      across(everything(), ~ replace_na(., 0))
    )
# Compute hat matrix and Cook's distance ----------------------------------
  # Compute hat matrix
  Zbar <- resids_n |> select(id, Lambda) |>
    left_join(Zn, by = c("id" = "id")) |>
    summarize(
      across(3 : (2+p), \(x) weighted.mean(x, Lambda))
    )


  # Compute Z - Zbar
  Zn_cen <- Zn |> select(-id) |> as.matrix() - matrix(rep(Zbar |> as.matrix(), n), nrow = n, byrow = TRUE)
  # Compute hat matrix
  H <- Zn_cen %*% solve(t(Zn_cen) %*% Zn_cen) %*% t(Zn_cen)
  hii <- diag(H)

  # Compute Cook's distance
  M <- resids_n$M
  sigma2 <- sum(M^2) / (n - p - 1)
  cook_d <- {M^2 / (sigma2 * (p + 1))} * {hii / (1 - hii)^2}
# max(cook_d)

  df_resids <- resids_n |>
    mutate(
      r = M / R
    ) |>
    bind_cols(
      hii = hii,
      cook_d = cook_d
    ) |>
    relocate(
      id, delta, R, Lambda, M, r, hii, cook_d
    )

  return(df_resids)

  }


#' @export
scores.pwreg1 <- function(x, ...) {

  call <- x$call
  # Contains delta, R, Lambda, M (win residuals)
  resids_n <- x$resids_n
  # Contains Zn
  Zn <- x$Zn
  p <- ncol(Zn) - 1
  n <- nrow(Zn)
  resids_n <- x$resids_n |>
    right_join(Zn |> select(id), by = "id") |>
    mutate(
      across(everything(), ~ replace_na(., 0))
    )
  # Compute hat matrix and Cook's distance ----------------------------------
  # Compute hat matrix
  Zbar <- resids_n |> select(id, Lambda) |>
    left_join(Zn, by = c("id" = "id")) |>
    summarize(
      across(3 : (2+p), \(x) weighted.mean(x, Lambda))
    )


  # Compute Z - Zbar
  Zn_cen <- Zn |> select(-id) |> as.matrix() - matrix(rep(Zbar |> as.matrix(), n), nrow = n, byrow = TRUE)
  # Compute hat matrix
  H <- Zn_cen %*% solve(t(Zn_cen) %*% Zn_cen) %*% t(Zn_cen)
  hii <- diag(H)

  # Compute Cook's distance
  M <- resids_n$M
  sigma2 <- sum(M^2) / (n - p - 1)
  cook_d <- {M^2 / (sigma2 * (p + 1))} * {hii / (1 - hii)^2}
  # max(cook_d)

  df_resids <- resids_n |>
    mutate(
      r = M / R
    ) |>
    bind_cols(
      hii = hii,
      cook_d = cook_d
    ) |>
    relocate(
      id, delta, R, Lambda, M, r, hii, cook_d
    )

  return(df_resids)

}

#' @export
predict.pwreg1 <- function(x, z1, z2, alpha = 0.05, contrast = FALSE, ...) {

  # z1 <- non_ischemic[1, 4:ncol(non_ischemic)]
  # z2 <- z1
  # z1$trt_ab = 1

  # Check if x$Y is not NULL
  if(!is.null(x$Y)){
    stop("Not applicable if there is continuous response.")
  }
  call <- x$call
  n <- x$n
  # Extract PW model data -------------------------------------------------------
  df <- x$df
  df_n <- df |>
    group_by(id) |>
    arrange(
      time
    ) |>
    slice(1) |>
    ungroup()
  # Compute beta
  beta <- x$beta
  expz1beta <- exp(sum(z1 * beta))
  expz2beta  <- exp(sum(z2 * beta))
  mu <- expz1beta  / (expz1beta  + expz2beta)
  # Influence function of beta
  psi_matrix <- x$psi_matrix

  # Fit Cox model to first event --------------------------------------------
  # Extract first event data
  # Fit Cox model
  cox_fit <- coxph(Surv(time, status > 0) ~ ., data = df_n |>
                     select(-id))
  # Get model matrix
  # Z <- model.matrix(cox_fit)
  Zn <- df_n |> select(-c(time, status))
  # Extract coefficients
  gamma <- coef(cox_fit)
  p <- length(gamma)
  # Extract efficient information
  v <- vcov(cox_fit)
  # sqrt(diag(v))
  invI <- v * n
  # Extract baseline cumulative hazard
  cox_base <- basehaz(cox_fit, centered = FALSE) |> as_tibble()
  # ?basehaz
  # Unique times
  ts <- cox_base$time
  m <- length(ts)
  # Extract baseline hazard
  Lambda <- cox_base$hazard
  # Take increments
  dLambda <- c(Lambda[1], diff(Lambda))

  # Point estimates of win-loss probs ---------------------------------------
  # Survival probabilities
  expz1 <- exp(sum(gamma * z1))
  expz2 <- exp(sum(gamma * z2))

  Sz1 <- exp(- expz1 * Lambda)
  Sz2 <- exp(- expz2 * Lambda)
  # Comparability probabilities
  comp_probs <- 1 - Sz1 * Sz2
  # WL probs
  win_prob <- comp_probs * mu
  loss_prob <- comp_probs - win_prob
  # Compute standard errors -------------------------------------------------
  # Function to compute counting, at-risk, martingale process

  counting_at_risk <- function(X, status){
    # Counting process and increments
    Nt <- (X <= ts) * status
    dNt <- c(Nt[1], diff(Nt))
    # At-risk process
    Yt <- (X >= ts) + 0
    return(list(dNt = dNt, Yt = Yt, ts = ts, dLambda = dLambda))
  }

  df_mart_stats <- df_n |>
    select(
      id,
      X = time,
      status
    ) |>
    mutate(
      status =  (status > 0) + 0
    ) |> bind_cols(
      nu = exp(Zn |> select(-id) |> as.matrix() %*% gamma |> as.vector())
    ) |>
    mutate(
      martingales = map2(X, status, counting_at_risk)
    ) |>
    unnest_wider(martingales) |>
    unnest_longer(c(dNt, Yt, dLambda, ts)) |>
    mutate(
      dMt = dNt - Yt * nu * dLambda,
    )
  # Compue s0t and s1t
  df_s0t <- df_mart_stats |>
    group_by(ts) |>
    summarize(
      s0t = mean(Yt * nu),
    )

  df_s1t <-
    df_mart_stats |>
    select(id, ts, Yt, nu) |>
    left_join(Zn, by = c("id" = "id")) |>
    group_by(ts) |>
    summarize(
      across(-c(id, Yt, nu), \(x) mean(x * Yt * nu)) # weighted average by Y(t) * nu
    )

  # E-function of t
  df_Et <- df_s0t |>
    left_join(df_s1t, by = "ts") |>
    mutate(
      across(-c(ts, s0t), \(x) if_else(s0t > 0, x / s0t, 0))
    )


  df_all_stats <- df_mart_stats |>
    select(id, ts, dLambda, dMt) |>
    left_join(df_Et, by = "ts")


  Et_dM_n <-  df_all_stats |>
    group_by(id) |>
    summarize(
      across(-c(ts, dLambda, dMt, s0t), \(x)  - sum(x * dMt))
    )

  Z_dM_n <- df_mart_stats|>
    select(id, ts, dMt) |>
    left_join(Zn, by = "id") |>
    group_by(id) |>
    summarize(
      across(-c(ts, dMt), \(x) sum(x * dMt))
    )


  cox_score_psi <- Z_dM_n |>
    bind_rows(Et_dM_n) |>
    group_by(id) |>
    summarize(
      across(everything(), sum)
    )


  cox_score_psi_mat <- cox_score_psi |>
    select(-id) |>
    as.matrix() %*% invI


  # Calculate H function
  # Joint with baseline hazards by ts
  Et_mat <-  df_Et |> select(-s0t) |>
    left_join(cox_base, by = c("ts" = "time"))
  # Get Lambda and dLambda
  Lambda <- Et_mat$hazard
  dLambda <- c(Lambda[1], diff(Lambda))
  Et_mat <- Et_mat |>
    bind_cols(
      dLambda = dLambda
    )
  # Extract Et
  Et <- Et_mat |> select(-c(ts, dLambda, hazard)) |> as.matrix()
  # Hz1
  z1 <- as.numeric(z1)
  H1mat <- Et_mat |>
    select(ts, dLambda) |>
    bind_cols(
      matrix(rep(z1, m), byrow = TRUE, nrow = m) - Et
    ) |>
    arrange(ts) |>
    mutate(
      across(-c(ts, dLambda),  \(x) expz1 * cumsum(x * dLambda))
    ) |>
    select(-c(ts, dLambda)) |>
    as.matrix()
  # Hz2
  z2 <- as.numeric(z2)
  H2mat <- Et_mat |>
    select(ts, dLambda) |>
    bind_cols(
      matrix(rep(z2, m), byrow = TRUE, nrow = m) - Et
    ) |>
    arrange(ts) |>
    mutate(
      across(-c(ts, dLambda),  \(x) expz2 * cumsum(x * dLambda))
    ) |>
    select(-c(ts, dLambda)) |>
    as.matrix()

  # Second part of Omega psi
  omega_mat2 <- df_all_stats |>
    select(id, ts, s0t, dMt) |>
    group_by(id) |>
    arrange(ts) |>
    mutate(
      s0t_inv = if_else(s0t > 0, 1 / s0t, 0),
      omega_psi_p2 = (expz1 + expz2) * cumsum(s0t_inv * dMt)
    ) |>
    pivot_wider(
      id_cols = ts,
      names_from = id,
      values_from = omega_psi_p2
    ) |>
    select(-ts) |>
    as.matrix()


  # Adding two parts to obtain omega psi
  # m x n matrix
  omega_psi_mat1 <- (H1mat + H2mat) %*% t(cox_score_psi_mat)
  omega_psi_mat <- Sz1 * Sz2 * (omega_psi_mat1 + omega_mat2)

  # rowMeans(omega_psi_mat)

  # PW part of psu
  # Compute psu
  # Win part
  win_if_beta <- (1 - mu) * win_prob %*% t(z1 - z2) %*% psi_matrix
  # loss_if_beta <- - win_if_beta
  win_if_mat <- mu * omega_psi_mat + win_if_beta
  loss_if_mat <- (1 - mu) * omega_psi_mat - win_if_beta
  # Pointwise variance
  win_prob_se <- sqrt(rowMeans(win_if_mat^2) / (n - p - 1))
  loss_prob_se <- sqrt(rowMeans(loss_if_mat^2) / (n - p - 1))


# Estimate and make inference on contrasts --------------------------------
if (contrast){
  # Variance matrix for beta
  Var <- x$Var
  # log-WR
  zd <- z1 - z2
  log_wr <- sum(zd * beta)
  log_wr_se <- sqrt(zd %*% Var %*% zd) |> as.numeric()






}


  # Return results
  result <- tibble(
    time = ts,
    win_prob = win_prob,
    win_prob_se = win_prob_se,
    loss_prob = loss_prob,
    loss_prob_se = loss_prob_se
  )
  return(result)
}



#' Print model information for a proportional win-fractions regression
#' @description Print model information for a proportional win-fractions regression.
#' @param x an object of class \code{pwreg1}.
#' @param ... further arguments passed to or from other methods.
#' @return Print the results of \code{pwreg1} object
#' @seealso \code{\link{pwreg1}}
#' @import dplyr
#' @export
#' @keywords pwreg1
#' @examples
#' # see the example for pwreg1
print.pwreg1 <- function(x,...){

  cat("Call:\n")
  print(x$call)
  cat("\n")
  strata <- x$strata
  beta <- x$beta
  # t <- x$t
  n <- x$n
  N <- x$N
  win_comp_tab <- x$pw_win |>
    count(component) |>
    mutate(prop = n / N)


  if(is.null(strata)){
    cat("Proportional win-fractions regression analysis:\n")
    tn <- n*(n-1)/2
    cat("Total number of pairs:", N,"\n")
    print(win_comp_tab)
    }else{
    cat("Stratified proportional win-fractions regression analysis:\n\n")
    cat("Total number of strata:", length(levels(strata)), "\n")
    }

  print(beta)

}






# Test another version of pw score ----------------------------------
# Tested
#
#
# pw_score_fun_alt <- function(pw_win, pw_Zd, Zn, N){
# ## mu
# delta_pm_f <- pw_win |> select(win) |> pull()
# Zd_f <- pw_Zd |> select(-id_i, -id_j) |> as.matrix()
# exp_beta_Zd_f <- exp(Zd_f %*% beta)
# mu_f <- as.vector(exp_beta_Zd_f / (1 + exp_beta_Zd_f))
#
# Z_i_f <- pw_win |> select(id_i) |> left_join(Zn, by = c("id_i" = "id")) |> select(-id_i) |> as.matrix()
# Z_j_f <- pw_win |> select(id_j) |> left_join(Zn, by = c("id_j" = "id")) |> select(-id_j) |> as.matrix()
#
# Rbar <- sum(abs(delta_pm_f))
#
# Zbar_f <- rep(0, p)
#
# for (k in 1:N){
#   if (abs(delta_pm_f[k]) == 1){
#     Zbar_f  <- Zbar_f + (mu_f[k] * Z_i_f[k, ]  + (1 - mu_f[k]) * Z_j_f[k, ]) / Rbar
#
#   }
# }
#
#
#
# delta_n_w <- pw_win |> group_by(id_i) |> summarize(delta_n = sum(win == 1) / (n - 1))
# delta_n_l <- pw_win |> group_by(id_j) |> summarize(delta_n = sum(win == -1) / (n - 1))
#
# delta_n <- delta_n_w |> full_join(delta_n_l, by = c("id_i" = "id_j")) |>
#   replace_na(list(delta_n.x = 0, delta_n.y = 0)) |> mutate(delta_n = delta_n.x + delta_n.y) |> pull(delta_n)
#
# score_f <- colMeans((Zn |> select(-id) |> as.matrix() - matrix(rep(Zbar_f, n), nrow = n, byrow = TRUE)) * delta_n)
#   score_f * 2
#
#   # score - score_f * 2
#   # trt_ab            age            sex Black.vs.White Other.vs.White            bmi       bipllvef       hyperten
#   # -6.422293e-14  -1.473488e-12  -4.957701e-13  -9.142687e-14  -4.748806e-17   3.168132e-12   7.837342e-13  -1.700168e-13
#   # COPD       diabetes           acei          betab      smokecurr
#   # -1.604619e-16   4.309053e-15  -3.320001e-13  -5.959018e-13   7.114101e-15
#
# }

# # test Zbar (verified)
#
# Zbar <- score_sum |> select(id, Lambda) |>
#   left_join(Zn, by = c("id" = "id")) |>
#   summarize(
#     across(3 : (2+p), \(x) weighted.mean(x, Lambda))
#   )
