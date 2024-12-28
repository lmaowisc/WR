
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
    Var <- psi_matrix %*% t(psi_matrix) / n^2

    # Take projected, individualized quantities: delta, R, Lambda,  M, score residuals
    resids_n <- score_sum

    obj <- list(beta = beta, Var = Var, conv = conv, iter = iter, # basic model fit
                resids_n = resids_n, Zn = Zn,  # individualized quantities (w. id row)
                A = A, kappa_matrix = kappa_matrix, psi_matrix = psi_matrix, # influence-related quantities
                pairwise_tibble = pairwise_tibble, # pairwise data for possible re-calculation of win-loss (by t)
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
