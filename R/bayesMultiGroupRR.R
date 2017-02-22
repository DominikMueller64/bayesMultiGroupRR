#' Bayesian hierarchical multi-group random regression model for genomic prediction.
#'
#' @return A \code{\link[rstan]{stanfit-class}} object.
#'
#' @param data A list. The data for fitting the model:
#'
##' \tabular{ll}{
##'    \strong{name}   \tab \strong{structure}                    \cr
##'    group           \tab index vector with group memberships   \cr
##'    X               \tab marker matrix with genotypes in rows  \cr
##'    y               \tab vector with phenotypic data
##' }
#'
#'
#' @details This function is a minimal wrapper around \code{\link[rstan]{stan}} and it only
#' supplies the model file argument and the data. All other argument must be taken from
#' \code{\link[rstan]{stan}}.
#'
#' This implementation can only handle a very limited number of markers (<= 1000).
#'
#' @seealso See \code{\link[rstan]{stan}}, for which this function is a minimal wrapper.
#'
#' @examples
#' ### The model file can be found here:
#' stanfile <- system.file('extdat', 'bayesMultiGroupRR.stan', package = 'bayesMultiGroupRR')
#'
#' ### This is a simple example, using simulated data.
#'
#' ## Load packages and set options. ------------------------------------------------------------
#'
#' library('rstan')
#' library('magrittr')
#' rstan_options(auto_write = TRUE)
#' options(mc.cores = parallel::detectCores())
#' rm(list = ls())
#'
#' ## copied from https://github.com/DominikMueller64/dmisc/blob/master/R/equal_split.R.
#' equal_split <- function(x, n, random = TRUE, beginning = FALSE) {
#'   len <- length(x)
#'   base_size <- len %/% n
#'   rem <- len %% n
#'   sizes <- rep(base_size, times = n)
#'   if (rem > 0L) {
#'     ix <- if (beginning) seq.int(from = 1L, to = rem) else sample.int(n, size = rem)
#'     sizes[ix] <- sizes[ix] + 1L
#'   }
#'
#'   s <- seq_len(length(x))
#'   ret <- purrr::map(sizes, function(size) {
#'     if (random) {
#'       smp <- sample(x = s, size = size)
#'     } else {
#'       smp <- s[seq.int(from = 1L, to = size)]
#'     }
#'     s <<- setdiff(s, smp)
#'     x[smp]
#'   })
#'   if (rem > 0L)
#'     attr(ret, which = 'indices') <- ix
#'   ret
#' }
#'
#' ## Load data ---------------------------------------------------------------------------------
#'
#' data('wheat', package = 'BGLR')
#' X <- 2 * wheat.X
#' rm(list = grep(pattern = 'wheat*', x = ls(), value = TRUE))
#'
#' ## Set parameters ----------------------------------------------------------------------------
#'
#' h2 <- 0.5  ## heritability
#' n_loci <- min(ncol(X), 50L)  ## number of (known) loci
#' n_groups <- 3L  ## number of groups (sub-populations)
#' base <- 0.75
#'   ## The true correlation matrix of the locus effects.
#' (true_cor_eff <- base^abs(outer(seq_len(n_groups), seq_len(n_groups), FUN = `-`)))
#' true_sigma_eff <- seq_len(n_groups)
#' (true_cov_eff <- diag(true_sigma_eff) %*% true_cor_eff %*% diag(true_sigma_eff))
#'   ## Sampling true locus effects.
#' true_eff <- mvtnorm::rmvnorm(n = n_loci, sigma = true_cov_eff) %>% split(., f = col(.))
#'   ## Empirical covariances.
#' (emp_cov_eff <- cov(do.call(what = cbind, true_eff)))
#'
#' ## Prepare data ------------------------------------------------------------------------------
#'
#'   ## Subset data to the desired number of loci.
#' X <- X[, sort(sample.int(n = ncol(X), size = n_loci)), drop = FALSE]
#'   ## Split the genotypes into (randomly sampled!) groups.
#' X <- dmisc::equal_split(x = split(x = X, f = row(X)), n = n_groups, random = TRUE)
#'   ## Compute genetic values.
#' g <- purrr::map2(X, true_eff, function(x, e) purrr::map_dbl(x, ~.x %*% e))
#'   ## Vector indication group-membership of individual observations.
#' group_size <- X %>% purrr::map_dbl(length)
#' group <- rep(seq_along(group_size), times = group_size)
#' n <- sum(group_size)  ## Total number of observations.
#'   ## Flatten genotypes to a list of matrices.
#' X <- X %>% purrr::map(~do.call(what = 'rbind', args = .x))
#'   ## Standardize genotypes within groups. Remove first homozygous loci!
#' ## X <- X %>% purrr::map(scale)
#'   ## Flatten genotypes to matrix.
#' X <- X %>% do.call(what = 'rbind')
#'   ## Flatten genetic values to vector
#' g <- purrr::flatten_dbl(g)
#'   ## Add some noise.
#' y <- g + rnorm(n = n, sd = sqrt(var(g) * (1 - h2) / h2))
#'
#' dat <- list(group = group,
#'             y = y,
#'             X = X)
#' \dontrun{
#' ## Fit model.
#' fit <- bayesMultiGroupRR(data = dat, iter = 300L, warmup = 100L,
#'                          chains = 4L, verbose = TRUE)
#'
#' ## Inspect traceplots of chains.
#' rstan::traceplot(fit, 'mu') # overall mean
#' rstan::traceplot(fit, 'beta')  # group-specific means
#' rstan::traceplot(fit, 'alpha_mu')  # overall locus effects
#' ## rstan::traceplot(fit, 'alpha') # group-specific locus effects (too much)
#' rstan::traceplot(fit, 'sigma_e') # residual standard deviation
#' rstan::traceplot(fit, 'sigma_alpha_mu') # overall locus effect standard deviation
#' rstan::traceplot(fit, 'sigma_alpha') # group-specific locus effect standard deviation
#' rstan::traceplot(fit, 'Omega') # correlations
#' rstan::traceplot(fit, 'h2') # heritability (locus - level!)
#'
#' ## Compute posterior covariance matrix, correlation matrix and linear predictors (GEBVs).
#' tmp <- rstan::get_posterior_mean(fit, 'Sigma')[, 'mean-all chains']
#' (post_mean_cov_eff <- matrix(data = tmp, ncol = n_groups))
#' cov2cor(post_mean_cov_eff)
#' tmp <- rstan::get_posterior_mean(fit, 'Omega')[, 'mean-all chains']
#' (post_mean_cor_eff <- matrix(data = tmp, ncol = n_groups))
#' post_theta <- rstan::get_posterior_mean(fit, 'theta')[, 'mean-all chains']
#'
#' ## Compute the prediction accuracy on per-group basis.
#' purrr::map2(split(g, f = group), split(post_theta, f = group), ~cor(.x, .y))
#'
#' }
#' @export
bayesMultiGroupRR <- function(data, ...) {
  stanfile <- system.file('extdat', 'bayesMultiGroupRR.stan', package = 'bayesMultiGroupRR')


  req_fields <- c('group', 'X', 'y')

  sdf <- setdiff(x = names(data), y = req_fields)

  if (length(sdf) > 0L)
    stop(sprintf(fmt = "The following data are missing: %s.",
                 paste(sdf, collapse = ', ')))

  n <- length(data$group)
  n_groups <- length(ug <- unique(data$group))
  n_loci <- ncol(data$X)
  group <- match(x = data$group, table = ug)  # 1-1 map of group to 1, ..., n_groups

  if (length(unique(n, nrow(data$X), length(data$y))) > 1L)
    stop("Dimensions of 'group', 'X' and 'y' are not conformable.")

  data <- list(n = n, n_groups = n_groups, n_loci = n_loci, group = group,
               y = data$y, X = data$X)
  rstan::stan(file = stanfile, data = data, ...)
}
