#~******************************************************************************
#~* OBJET: Estimation bayesienne d'une distribution à partir d'un échantillon 
#~*        mixte pouvant contenir des valeurs parfaitement connues ou censurées, 
#~*        avec une prise en compte des erreurs systématiques dues aux courbes
#~*        de tarage. Utile pour l'intégration de crues historiques.
#~******************************************************************************
#~* PROGRAMMEUR: Benjamin Renard, INRAE Aix-en-Provence
#~******************************************************************************
#~* CREE/MODIFIE: 15/04/2024
#~******************************************************************************
#~* PRINCIPALES FONCTIONS
#~*    1. GetEstimate_HBay
#~******************************************************************************
#~* REF.: Neppel, L., Renard, B., Lang, M., Ayral, P. A., Coeur, D., 
#~*       Gaume, E., … Vinet, F. (2010). Flood frequency analysis using 
#~*       historical data: accounting for random and systematic errors. 
#~*       Hydrological Sciences Journal, 55(2), 192–208.
#~*       https://doi.org/10.1080/02626660903546092
#~******************************************************************************
#~* A FAIRE: -
#~******************************************************************************
#~* COMMENTAIRES: Ré-implémentation en R de l'éxécutable HBay, voir
#~*               https://github.com/benRenard/BMSL/tree/main/cli/HBay
#~******************************************************************************

#****************************
# Fonctions principales ----
#****************************

#' Bayesian estimation of a distribution using historical data
#'
#' XXX 
#' Returns MCMC samples from the posterior distribution.
#'
#' @param y numeric 2-column matrix, data. 
#'     The first column gives the lower bound, the second column gives the upper bound.
#'     Where y[i,1]==y[i,2], the value is assumed perfectly known (up to systematic errors, see below).
#'     Where y[i,1]<y[i,2], the value is assumed to be in the interval [y[i,1];y[i,2]]
#'     -Inf and +Inf are allowed for data being only right- or left-censored 
#'     (i.e. values known to be smaller than or larger than some threshold).
#' @param dist character, distribution name. Only distributions 'GEV', 'Gumbel' and 'LogNormal' are supported.
#' @param prior list of lists, prior distributions. For each parameter to be estimated, the prior
#'     is a list of the form pr=list(dist=..., par=...). See example below.
#' @param SystErrorIndex integer vector, length NROW(y). Index of systematic errors.
#'     Rows where SystErrorIndex==k are all affected by the same multiplicative error gamma_k, 
#'     typically induced by the kth rating curve. SystErrorIndex==0 means no systematic error.
#'     Should only contain integer values between 0 and N_\{systematic errors\}.
#' @param SystErrorPrior list of lists,  prior distribution for each systematic error.
#'     For instance for a systematic error in the range +/- 20\%, you may use a Uniform
#'     between 0.8 and 1.2, or a triangular distribution with the same bounds and peak at 1.
#' @param par0 numeric vector, initial parameter guess.
#' @param SystError0 numeric vector, initial guess for systematic errors. Typically a vector of 1.
#' @param mult numeric, initial jump standard deviations are set to mult * abs(par0)
#' @param eps numeric, where par0 is zero, initial jump standard deviations are set to eps (to avoid jumps of size zero)
#' @param batch.length integer, MCMC parameter: length of each non-adaptive batch 
#' @param batch.n integer, MCMC parameter: number of batches (= adaptation period). Total number of simulations is nsim=batch.n*batch.length
#' @param moverate.min numeric in (0;1), MCMC parameter: lower bound for the desired move rate interval
#' @param moverate.max numeric in (0;1), MCMC parameter: upper bound for the desired move rate interval
#' @param mult.down numeric in (0;1), MCMC parameter: multiplication factor used to decrease jump size when move rate is too low.
#' @param mult.up numeric (>1, avoid 1/mult.down), MCMC parameter: multiplication factor used to increase jump size when move rate is too high.
#' @return A list with the following components:
#'     \item{x}{numeric matrix nsim * (length(par0)+length(SystError0)), MCMC simulations}
#'     \item{fx}{numeric vector, corresponding values f(x)}
#' @examples
#' @export
GetEstimate_HBay<-function(y,dist,prior,
                           SystErrorIndex=rep(0,NROW(y)),SystErrorPrior=list(),
                           par0=GetEstimate_ROUGH(0.5*(y[,1]+y[,2]),dist)$par,SystError0=rep(1,length(SystErrorPrior)),
                           mult=0.1,eps=0.1,
                           batch.length=100,batch.n=100,
                           moverate.min=0.1,moverate.max=0.5,
                           mult.down=0.9, mult.up=1.1){

  # sdjump=mult*par0;sdjump[sdjump==0]=eps
  # out=Metropolis_OAAT_adaptive(GetLogPost,par0,sdjump,
  #                              y,dist,prior,
  #                              batch.length=batch.length,batch.n=batch.n,
  #                              moverate.min=moverate.min,moverate.max=moverate.max,
  #                              mult.down=mult.down, mult.up=mult.up)
  # return(out)
}
