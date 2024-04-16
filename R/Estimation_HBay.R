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

#' Bayesian estimation using historical data
#'
#' Bayesian estimation of a GEV or Gumbel distribution based on a mixed sample containing 
#' point (i.e. perfectly known) or interval (i.e. known to be within bounds) data.
#' Systematic errors induced by rating curve errors can also be accounted for. 
#' Returns MCMC samples from the posterior distribution.
#'
#' @param y numeric 2-column matrix, data. 
#'     The first column gives the lower bound, the second column gives the upper bound.
#'     Where y[i,1]==y[i,2], the value is assumed perfectly known (up to systematic errors, see below).
#'     Where y[i,1]<y[i,2], the value is assumed to be in the interval [y[i,1];y[i,2]]
#'     -Inf and +Inf are allowed for data being only right- or left-censored 
#'     (i.e. values known to be smaller than or larger than some threshold).
#' @param dist character, distribution name. Only distributions 'GEV' and 'Gumbel' are supported.
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
#' set.seed(98765)
#' y0=Generate('GEV',c(100,50,-0.2),100)
#' y=cbind(y0,y0)
#' # Mimics censoring between 0 and 300 for first 70 years
#' y[1:70,1][y0[1:70]<300]=0
#' y[1:70,2][y0[1:70]<300]=300
#' plot(y[,1]);points(y[,2])
#' # Systematoc errors
#' SystErrorIndex=c(rep(1,70),rep(2,30))
#' SystErrorPrior=list(list(dist="Triangle",par=c(1,0.7,1.3)),
#'                     list(dist="Triangle",par=c(1,0.95,1.05)))
#' # Priors on GEV parameters
#' prior=list(list(dist="FlatPrior",par=NULL),
#'            list(dist="FlatPrior",par=NULL),
#'            list(dist="Normal",par=c(0,0.25)))
#' # Go!
#' mcmc=GetEstimate_HBay(y=y,dist='GEV',prior=prior,
#'                       SystErrorIndex=SystErrorIndex,
#'                       SystErrorPrior=SystErrorPrior,
#'                       batch.length=25,batch.n=40)
#' par(mfrow=c(2,3));for(i in 1:5){hist(mcmc$x[,i])}
#' @export
GetEstimate_HBay<-function(y,dist,prior,
                           SystErrorIndex=rep(0,NROW(y)),SystErrorPrior=list(),
                           par0=GetEstimate_ROUGH(0.5*(y[,1]+y[,2]),dist)$par,SystError0=rep(1,length(SystErrorPrior)),
                           mult=0.1,eps=0.1,
                           batch.length=100,batch.n=100,
                           moverate.min=0.1,moverate.max=0.5,
                           mult.down=0.9, mult.up=1.1){
  
  # Preliminary checks ----
  if( !(dist %in% c('GEV','Gumbel')) ){
    stop(paste0('Unsupported distribution: "',dist,'". Only "GEV" and "Gumbel" are allowed'))
  }
  if( any(y[,1]>y[,2]) ){
    stop(paste0('Bounds inversion y[i,1]>y[i,2]: lower bound is larger than upper bound for at least one row in y'))
  }
  if( any(SystErrorIndex <0 ) ){
    stop(paste0('Invalid systematic error index: SystErrorIndex <0'))
  } 
  if( max(SystErrorIndex)>length(SystErrorPrior) ){
    stop(paste0('The number of prior distributions provided in SystErrorPrior ',
                'is too small compared with the number of systematic errors',
                ': max(SystErrorIndex)>length(SystErrorPrior)'))
  } 
  if(length(prior)>0){
    if( any(sapply(prior,class)!='list') ){
      stop(paste0(
        'Badly formed prior. A list of lists is expected, e.g. for 2 parameters:\n ',
        'prior=list( list(dist="FlatPrior",par=NULL) , list(dist="Normal",par=c(0,0.25)) )')
        )
    }
  }
  if(length(SystErrorPrior)>0){
    if( any(sapply(SystErrorPrior,class)!='list') ){
      stop(paste0(
        'Badly formed SystErrorPrior A list of lists is expected, e.g. for 2 systematic errors:\n ',
        'SystErrorPrior=list( list(dist="Triangle",par=c(1,0.7,1.3)) , list(dist="Triangle",par=c(1,0.9,1.1)) )')
      )
    }
  }
  
  # MCMC sampling ----
  sdjump=mult*c(par0,SystError0)
  sdjump[sdjump==0]=eps
  out=Metropolis_OAAT_adaptive(GetHBayLogPost,c(par0,SystError0),sdjump,
                               y=y,dist=dist,prior=prior,
                               SystErrorIndex=SystErrorIndex,SystErrorPrior=SystErrorPrior,
                               batch.length=batch.length,batch.n=batch.n,
                               moverate.min=moverate.min,moverate.max=moverate.max,
                               mult.down=mult.down, mult.up=mult.up)
  return(out)
}

GetHBayLogPost<-function(par,y,dist,prior,SystErrorIndex,SystErrorPrior){
  # Number of systematic errors gamma
  Nk=length(SystErrorPrior)
  # Interpret par
  theta=par[1:(length(par)-Nk)] # distribution parameters
  if(Nk>0){
    gamma=par[(length(par)-Nk+1):(length(par))] # systematic errors
    if(any(gamma<=0)){return(NA)}
  } else {
    gamma=NULL
  }
  # Priors------
  prior=GetLogPrior(theta,prior)
  if(Nk>0){prior=prior+GetLogPrior(gamma,SystErrorPrior)}
  
  # Likelihood-------
  # Get data type
  dif=y[,2]-y[,1]
  if(any(dif<0)){return(NA)}
  dataType=as.integer(dif>0) # 0=perfectly known, 1=interval
  # Start computations
  lkh=0
  for(k in 0:Nk){
    indx=which(SystErrorIndex==k)
    if(length(indx)>0){
      z=y[indx,]
      dt=dataType[indx]
      if(k==0) {g=1} else {g=gamma[k]}
      distPar=c(theta[1]/g,theta[2]/g)
      if(dist=='GEV'){distPar=c(distPar,theta[3])}
      # Perfectly-known data
      w=z[dt==0,]
      if(NROW(w)>0){
        lkh=lkh+GetLogLikelihood(distPar,w[,1],dist)
      }
      # Interval data
      w=z[dt==1,]
      if(NROW(w)>0){
        upper=sapply(w[,2],GetCdf,dist,distPar)
        lower=sapply(w[,1],GetCdf,dist,distPar)
        lkh=lkh+sum(log(upper-lower))
      }
    }
  }
  
  # Compute post and return ------
  post=prior+lkh
  return(post)
}
