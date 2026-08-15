#' Wrapper for inference: binomial response with Gneiting-class kernel
#' 
#' (Most hyperparameters are hardcoded with the exception of number of samples/
#' length of warmup/number of chains. Option to provide folds for cross valiation)
#'
#' @param mut_data data.frame from `setup_mut_data()`
#' @param covariates SpatRaster
#' @param pfpr_years vector of numerics (years represented in `covariates`)
#' @param out_dir character: directory in "output/" to write to ("marker/model/")
#' @param fold numeric: indexes `folds`
#' @param folds list of folds from `caret::createFolds()`
#' @param warmup numeric
#' @param nsamples numeric
#' @param nchains numeric
#'
#' @returns
#' @export
#'
#' @examples
#' 
#' # without CV folds:
#' fit_binom(mut_data = mut_data,
#' covariates = covariates,
#' pfpr_years = pfpr_years,
#' out_dir = out_dir)
#' 
#' # with CV folds:
#' fit_binom(mut_data = mut_data,
#' covariates = covariates,
#' pfpr_years = pfpr_years,
#' out_dir = out_dir,
#' fold = 1,
#' folds = folds)
#' 
fit_binom <- function(mut_data,
                     covariates,
                     pfpr_years,
                     out_dir,
                     fold = NULL, # numeric !
                     folds = NULL,
                     warmup = 100, 
                     nsamples = 100, 
                     nchains = 6,
                     buffer = 0){
  
  message(nrow(mut_data))
  
  out <- build_design_matrix(covariates,
                             coords = mut_data,
                             scale = FALSE,
                             temporal_var = TRUE,
                             temporal_covt_range = pfpr_years,
                             degs_to_rads = TRUE,
                             buffer = buffer)
  if(!is.null(fold)){
    train <- unlist(folds[-c(fold)])
    mut_data = mut_data[train,]
    fold = paste0("_", fold)  # faff
  } else {
    fold = ""
  }

  message("here")
  
  X_obs <- out$df %>% 
    # not previously here but feels sensible:
    filter(!is.na(pfpr))
  
  scaled_years <- out$scaled_years
  
  coord_cols <- c("x_rd", "y_rd", "year_scaled")
  
  design_cols <- c("intercept", "year_scaled", "pfpr")
  
  # hyperparameters
  gneiting_len <- normal(0, 3, truncation = c(0, Inf))
  gneiting_tim <- normal(0, 3, truncation = c(0, Inf))
  gneiting_sd <- normal(0, 2, truncation = c(0, Inf))
  white_sd <- normal(0, 3, truncation = c(0, Inf))
  
  # kernel & GP
  kernel <- gneiting(lengthscale = gneiting_len, 
                     timescale = gneiting_tim,
                     variance = gneiting_sd ** 2, 
                     columns = 1:3) + 
    white(white_sd ** 2)
  
  kmn <- kmeans(X_obs[,coord_cols], centers = 40)
  random_field <- gp(x = X_obs[,coord_cols], 
                     kernel = kernel,
                     inducing = kmn$centers)
  
  beta <- normal(0, 1, dim = 3)
  gp_mean_obs <- X_obs[,design_cols] %*% beta + random_field
  X_prob_obs <- ilogit(gp_mean_obs)
  
  # likelihood
  distribution(X_obs$present) <- binomial(X_obs$tested, X_prob_obs)
  
  # fit the model by Hamiltonian Monte Carlo
  m <- model(gneiting_len, gneiting_tim, gneiting_sd, white_sd, beta)
  
  # hardcoding hmc hyperparams - have had a play with these earlier
  {start <- Sys.time()
    set.seed(seed)
    draws <- mcmc(m,
                  sampler = hmc(Lmin = 10, Lmax = 15),
                  chains = nchains,
                  warmup = warmup,
                  n_samples = nsamples,
                  initial_values = initials(gneiting_len = 1,
                                            gneiting_tim = 3,
                                            gneiting_sd = 13,
                                            white_sd = 0.5,
                                            beta = rep(0, 3)),
                  verbose = FALSE)
    end <- Sys.time()
    end - start}
  
  r_hats <- coda::gelman.diag(draws,
                              autoburnin = FALSE,
                              multivariate = FALSE)
  message(summary(r_hats$psrf))
  
  parameters <- list(gneiting_len, gneiting_tim, gneiting_sd, white_sd, beta)
  names(parameters) <- c("gneiting_len", "gneiting_tim", "gneiting_sd",
                         "white_sd", "beta")

  message(paste0("Writing", out_dir, fold))
  
  # save everything and do the prediction separately
  write_rds(parameters, paste0("output/", out_dir, "parameters", fold, ".rds"))
  write_rds(kernel, paste0("output/", out_dir, "kernel", fold, ".rds"))
  write_rds(random_field, paste0("output/", out_dir, "random_field", fold, ".rds"))
  write_rds(m, paste0("output/", out_dir, "m", fold, ".rds"))
  write_rds(draws, paste0("output/", out_dir, "draws", fold, ".rds"))
  
  return(TRUE)
}



#' Wrapper for inference: beta-binomial response with Gneiting-class kernel
#' 
#' (Most hyperparameters are hardcoded with the exception of number of samples/
#' length of warmup/number of chains. Option to provide folds for cross valiation)
#'
#' @param mut_data data.frame from `setup_mut_data()`
#' @param covariates SpatRaster
#' @param pfpr_years vector of numerics (years represented in `covariates`)
#' @param out_dir character: directory in "output/" to write to ("marker/model/")
#' @param fold numeric: indexes `folds`
#' @param folds list of folds from `caret::createFolds()`
#' @param warmup numeric
#' @param nsamples numeric
#' @param nchains numeric
#'
#' @returns
#' @export
#'
#' @examples
#' 
#' # without CV folds:
#' fit_betabinom(mut_data = mut_data,
#' covariates = covariates,
#' pfpr_years = pfpr_years,
#' out_dir = out_dir)
#' 
#' # with CV folds:
#' fit_betabinom(mut_data = mut_data,
#' covariates = covariates,
#' pfpr_years = pfpr_years,
#' out_dir = out_dir,
#' fold = 1,
#' folds = folds)
#' 
fit_betabinom <- function(mut_data,
                      covariates,
                      pfpr_years,
                      out_dir,
                      fold = NULL,
                      folds = NULL,
                      lfo = FALSE, # flag for leave-future-out CV
                      warmup = 100, 
                      nsamples = 100, 
                      nchains = 6,
                      buffer = 0){

  message("here!")
  
  if(!is.null(fold)){
    # remove testing data from training set
    if(!lfo){
      train <- unlist(folds[-c(fold)])
      # response-stratified CV/spatially-stratified CV: 
      # remove test set from full set
    } else {
      message("lfo!")
      train <- unlist(folds[fold])  
      # LFO: fold is everything up to current year
      # (setdiff on next fold for test set)
    }
    
    mut_data = mut_data[train,]
    fold = paste0("_", fold)  # faff
    
  } else {
    fold = ""
  }

  message(nrow(mut_data))
  
  out <- build_design_matrix(covariates,
                             coords = mut_data,
                             scale = FALSE,
                             temporal_var = TRUE,
                             temporal_covt_range = pfpr_years,
                             degs_to_rads = TRUE,
                             buffer = buffer)
  X_obs <- out$df%>% 
    # not previously here but feels sensible:
    filter(!is.na(pfpr))
  
  scaled_years <- out$scaled_years
  
  coord_cols <- c("x_rd", "y_rd", "year_scaled")
  
  design_cols <- c("intercept", "year_scaled", "pfpr")
  
  # hyperparameters
  gneiting_len <- normal(0, 3, truncation = c(0, Inf))
  gneiting_tim <- normal(0, 3, truncation = c(0, Inf))
  gneiting_sd <- normal(0, 2, truncation = c(0, Inf))
  white_sd <- normal(0, 3, truncation = c(0, Inf)) # Median :1.041  Mean   :1.115
  rho <- greta::lognormal(0, 1)
  
  # kernel & GP
  # could potentially give stricter priors to variance/nugget here - trouble with IDability?
  kernel <- gneiting(lengthscale = gneiting_len, 
                     timescale = gneiting_tim,
                     variance = gneiting_sd ** 2, 
                     columns = 1:3) + 
    white(white_sd ** 2)
  
  kmn <- kmeans(X_obs[,coord_cols], centers = 40)
  random_field <- greta.gp::gp(x = X_obs[,coord_cols], 
                               kernel = kernel,
                               inducing = kmn$centers)
  
  beta <- normal(0, 1, dim = 3)
  gp_mean_obs <- X_obs[,design_cols] %*% beta + random_field
  X_prob_obs <- ilogit(gp_mean_obs)
  
  # likelihood
  distribution(X_obs$present) <- betabinomial_p_rho(X_obs$tested, X_prob_obs, rho)
  
  # fit the model by Hamiltonian Monte Carlo
  m <- model(gneiting_len, gneiting_tim, gneiting_sd, white_sd, beta, rho)
  
  # hardcoding hmc hyperparams - have had a play with these earlier
  {start <- Sys.time()
    set.seed(seed)
    draws <- mcmc(m,
                  sampler = hmc(Lmin = 10, Lmax = 15),
                  chains = nchains,
                  warmup = warmup,
                  n_samples = nsamples,
                  initial_values = initials(gneiting_len = 1,
                                            gneiting_tim = 3,
                                            gneiting_sd = 5,
                                            white_sd = 0.5,
                                            beta = rep(0, 3),
                                            rho = 0.5),
                  verbose = FALSE)
    end <- Sys.time()
    end - start}
  
  r_hats <- coda::gelman.diag(draws,
                              autoburnin = FALSE,
                              multivariate = FALSE)
  message(summary(r_hats$psrf))
  
  parameters <- list(gneiting_len, gneiting_tim, gneiting_sd, white_sd, beta, rho)
  names(parameters) <- c("gneiting_len", "gneiting_tim", "gneiting_sd",
                         "white_sd", "beta", "rho")
  
  # would be nice to put a traceplot back in here ...
  message(paste0("Writing", out_dir, fold))
  
  # save everything and do the prediction separately
  write_rds(parameters, paste0("output/", out_dir, "parameters", fold, ".rds"))
  write_rds(kernel, paste0("output/", out_dir, "kernel", fold, ".rds"))
  write_rds(random_field, paste0("output/", out_dir, "random_field", fold, ".rds"))
  write_rds(m, paste0("output/", out_dir, "m", fold, ".rds"))
  write_rds(draws, paste0("output/", out_dir, "draws", fold, ".rds"))
  
  return(TRUE)
}