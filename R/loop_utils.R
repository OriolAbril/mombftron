#' Create data.table objects
#'
#' @param n_vec Array of numbers of observations
#' @param n_rep Number of repetitions to done of the simulation
#' @param prior_vec Array of prior names to compare
#' @param covariates Array with covariate names
#' @export
create_data.table_n_rep_prior <- function(n_vec, n_rep, prior_vec, covariates) {
  n_col <- rep(n_vec, each=length(prior_vec)*n_rep)
  prior_col <- rep(prior_vec, each=n_rep, times=length(n_vec))
  rep_col <- rep(seq(n_rep), times=length(n_vec)*length(prior_vec))
  model_dt <- data.table::data.table(
    n=n_col,
    prior=prior_col,
    repetition=rep_col,
    model_top=0,
    model_prob=0
  )

  marg_df_aux <- data.frame(matrix(0, nrow=length(prior_col), ncol=length(covariates)))
  colnames(marg_df_aux) <- covariates
  marg_dt <- data.table::data.table(
    cbind(
        data.frame(
            n=n_col,
            prior=prior_col,
            repetition=rep_col
        ),
        marg_df_aux
    )
  )
  return(list(model_dt=model_dt, marg_dt=marg_dt))
}

#' Store fit results in data.tables
#'
#' @param model_dt Data table to store model probability and model frequency
#' @param marg_dt Data table to store marginal inclusion probabilities
#' @param row Row of the data.table to modify
#' @param fit Result of [mombf::modelSelection()]
#' @param model_idxs String containing the indexs of the model whose details have to be stored
#' @param covariates Array with covariate names
#' @export
update_data.table <- function(model_dt, marg_dt, row, fit, model_idxs, covariates) {
  if (missing(covariates)) { covariates <- colnames(fit$xstd) }
  marg_dt[row, (covariates) := subset(data.frame(t(fit$margpp)), select=covariates)]
  pprobs <- mombf::postProb(fit)
  modelids <- pprobs$modelid
  model_chosen <- model_idxs == as.character(modelids[1])
  model_dt[row, "model_top" := model_chosen]
  mask <- model_idxs == modelids
  if (!sum(mask)) {
    model_dt[row, "model_prob" := 0]
  } else {
    model_dt[row, "model_prob" := pprobs$pp[mask]]
  }
}

#' Convert prior code to mobmf priorspec object
#'
#' Automatically creates priorspec objects from a string. For mom `tau=0.34`;
#' for normid and groupmom `tau=1` and for zellner and groupzellner `tau=n`
#'
#' @param priorname One of the following strings: mom, normid, zell, gzell, gmom
#' @param n Number of observations in data where prior is to be used
#' @export
get_prior <- function(priorname, n) {
  if (priorname == "mom") {
    prior <- mombf::momprior(tau=.34)
  } else if (priorname == "normid") {
    prior <- mombf::normalidprior(tau=1)
  } else if (priorname == "zell") {
    prior <- mombf::zellnerprior(tau=n)
  } else if (priorname == "gzell") {
    prior <- mombf::groupzellnerprior(tau=n)
  } else if (priorname == "gmom") {
    prior <- mombf::groupmomprior(tau=1)
  }
  return(prior)
}

#' Convert prior pair code to priorspec objects
#'
#' Convert a `_` separated pair of prior codes (as defined in [get_prior()])
#'
#' @param pair_name `_` separated prior codes.
#' @param n Number of observations
#' @seealso [get_prior()]
#' @export
get_prior_pair <- function(pair_name, n) {
  priorname <- strsplit(pair_name, "_")[[1]]
  return(list(pCoef=get_prior(priorname[1],  n), pGroup=get_prior(priorname[2], n)))
}