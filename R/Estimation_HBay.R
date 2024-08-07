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
#~*    1. GetEstimate_HBay: Bayesian+MCMC estimation, returns MCMC samples
#~*    2. Hydro3_HBay: wraps the function above to return an Hydro3 object
#~*    3. Import_HBayConfig: Imports a configuration folder used by HBay executable
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
#' n=50;n_censored=30
#' y0=Generate('GEV',c(100,50,-0.2),n)
#' y=cbind(y0,y0)
#' # Mimics censoring between 0 and 300 for first n_censored years
#' y[1:n_censored,1][y0[1:n_censored]<300]=0
#' y[1:n_censored,2][y0[1:n_censored]<300]=300
#' plot(y[,1]);points(y[,2])
#' # Systematic errors
#' SystErrorIndex=c(rep(1,n_censored),rep(2,n-n_censored))
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
#'                       # The values below aim at making this example fast to run.
#'                       # In practice, it is recommended to use the default values
#'                       # (batch.length=100,batch.n=100) or larger.
#'                       batch.length=25,batch.n=20) 
#' graphicalpar=par(mfrow=c(2,3))
#' for(i in 1:5){hist(mcmc$x[,i])}
#' par(graphicalpar)
#' @export
GetEstimate_HBay<-function(y,dist,prior,
                           SystErrorIndex=rep(0,NROW(y)),SystErrorPrior=list(),
                           par0=GetEstimate_ROUGH(0.5*(y[,1]+y[,2])[is.finite(y[,1]+y[,2])],dist)$par,
                           SystError0=rep(1,length(SystErrorPrior)),
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


#' Bayesian estimation using historical data
#'
#' Bayesian estimation of a GEV or Gumbel distribution based on a mixed sample containing 
#' point (i.e. perfectly known) or interval (i.e. known to be within bounds) data.
#' Systematic errors induced by rating curve errors can also be accounted for. 
#' Returns an Hydro3 object
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
#' @param options list, options, see details below.
#' @param mcmcoptions list, MCMC options, see details below.
#' @param do.KS,do.MK,do.Pettitt logical, perform KS/MK/Pettitt tests?
#' @return A list with the following components:
#'     \item{dist}{character, estimated distribution.}
#'     \item{ok}{logical, did estimation succeed?}
#'     \item{err}{integer, error code (0 if ok).}
#'     \item{message}{error message.}
#'     \item{empirical}{data frame, sorted data and empirical estimates 
#'         (nonexceedance frequency, return period and reduced variate). 
#'         NOTE: interval data are replaced by a value randomly sampled from
#'         a GEV constrained in this interval. See ?GenerateWithinBounds.}
#'     \item{pcdf}{data frame, estimated pdf and cdf}
#'     \item{quantile}{data frame, estimated quantiles and uncertainty intervals}
#'     \item{par}{data frame, estimated parameters and uncertainty intervals}
#'     \item{KS}{list, result of the Kolmogorov-Smirnov test, see ?KS.
#'         NOTE: interval data are replaced by a value randomly sampled from
#'         a GEV constrained in this interval. See ?HBay_simGEV.}
#'     \item{MK}{list, result of the Mann-Kendall test, see ?MK. Same note as KS test.}
#'     \item{Pettitt}{list, result of the Pettitt test, see ?Pettitt. Same note as KS test.}
#'     \item{u}{list, parameter uncertainty in the form of a covariance matrix ($cov)
#'         and simulated parameter replicates ($sim). Also contains error-handling flags 
#'         $ok, $err and $message.}
#' @details The argument 'options' allows controlling various properties of the analysis and results.
#'     It is a list with the following components:
#'     \itemize{
#'     \item{FreqFormula, character, formula for computing nonexceedance frequency, see ?GetEmpFreq.}
#'     \item{pgrid, numeric vector, probabilities defining the x values where pdf f(x) and cdf F(x) 
#'         are computed. These x values are quantiles from the estimated distribution 
#'         with probabilities pgrid.}
#'     \item{Tgrid, numeric vector, return periods where quantile function q(T) is computed.}
#'     \item{IClevel, numeric, level of uncertainty interval.}
#'     \item{p2T, numeric, conversion factor between nonexceedance probability p and return period T.
#'         p=1-1/(p2T*T). Here p2T=1 in general since GEV/Gumbel are applied to annual maxima in general.}
#'     \item{invertT, logical, when invertT=TRUE, LARGE return periods correspond to SMALL data values.
#'         This is typically used for low-flow statistics. Unused here.}
#'     \item{splitZeros, logical, when splitZeros=TRUE zero and negative values are removed from the data y before 
#'         estimating the distribution,and are used to estimate the probability of zeros p0. This is 
#'         typically used for low-flow statistics to estimate the probability of zero streamflow. Unused here.}
#'     \item{lang, chanracter, language ('fr' or 'en').}
#'     \item{nsim, integer, number of replicated parameters representing uncertainty. Unused here (derives from mcmc options)}
#'     }
#'     The argument 'mcmcoptions' is a list controlling MCMC properties: 
#'     \itemize{
#'     \item{mult, numeric, see ?Metropolis_OAAT_adaptive}
#'     \item{eps, numeric, see ?Metropolis_OAAT_adaptive}
#'     \item{batch.length, integer, see ?Metropolis_OAAT_adaptive}
#'     \item{batch.n, integer, see ?Metropolis_OAAT_adaptive}
#'     \item{moverate.min, numeric, see ?Metropolis_OAAT_adaptive}
#'     \item{moverate.max, numeric, see ?Metropolis_OAAT_adaptive}
#'     \item{mult.down, numeric, see ?Metropolis_OAAT_adaptive}
#'     \item{mult.up, numeric, see ?Metropolis_OAAT_adaptive}
#'     \item{burn, numeric, burn-in factor, e.g. if burn=0.2 the first 20 percents of MCMC samples are discarded}
#'     \item{slim, integer, sliming factor, e.g. if slim=5 only one MCMC sample every 5 is kept (after burn-in)}
#'     }
#' @examples
#' set.seed(98765)
#' n=50;n_censored=30
#' y0=Generate('GEV',c(100,50,-0.2),n)
#' y=cbind(y0,y0)
#' # Mimics censoring between 0 and 300 for first n_censored years
#' y[1:n_censored,1][y0[1:n_censored]<300]=0
#' y[1:n_censored,2][y0[1:n_censored]<300]=300
#' plot(y[,1]);points(y[,2])
#' # Systematic errors
#' SystErrorIndex=c(rep(1,n_censored),rep(2,n-n_censored))
#' SystErrorPrior=list(list(dist="Triangle",par=c(1,0.7,1.3)),
#'                     list(dist="Triangle",par=c(1,0.95,1.05)))
#' # Priors on GEV parameters
#' prior=list(list(dist="FlatPrior",par=NULL),
#'            list(dist="FlatPrior",par=NULL),
#'            list(dist="Normal",par=c(0,0.25)))
#' # Handle MCMC options
#' # The values below aim at making this example fast to run.
#' # In practice, it is recommended to use the default values
#' # (batch.length=100,batch.n=100) or larger.
#' mcmcoptions=mcmcoptions_def
#' mcmcoptions$batch.length=25
#' mcmcoptions$batch.n=20
#' # Go!
#' H3=Hydro3_HBay(y=y,dist='GEV',prior=prior,
#'                SystErrorIndex=SystErrorIndex,
#'                SystErrorPrior=SystErrorPrior,
#'                mcmcoptions=mcmcoptions) 
#' Hydro3_Plot(H3)
#' @export
Hydro3_HBay <- function(y,dist,prior=GetDefaultPrior(GetParNumber(dist)),
                        SystErrorIndex=rep(0,NROW(y)),SystErrorPrior=list(),
                        options=options_def,mcmcoptions=mcmcoptions_def,
                        do.KS=TRUE,do.MK=TRUE,do.Pettitt=TRUE){
  # get MCMC samples
  mcmc=GetEstimate_HBay(y=y,dist=dist,prior=prior,
                        SystErrorIndex=SystErrorIndex,SystErrorPrior=SystErrorPrior,
                        mult=mcmcoptions$mult,eps=mcmcoptions$eps,
                        batch.length=mcmcoptions$batch.length,batch.n=mcmcoptions$batch.n,
                        moverate.min=mcmcoptions$moverate.min,moverate.max=mcmcoptions$moverate.max,
                        mult.down=mcmcoptions$mult.down, mult.up=mcmcoptions$mult.up)
  
  # build H3 object ----
  par.ncol=8
  empirical.ncol=4
  pcdf.ncol=3
  quantile.ncol=6
  out=H3_def_success
  out$dist=dist
  w=GetMode(mcmc)
  if(w$ok==FALSE) {out=H3_def_fail;out$message="estimation:echec";return(out)}
  # fill in $par
  nDpar=GetParNumber(dist)
  npar=NCOL(mcmc$x)
  out$par=data.frame(matrix(NA,npar,par.ncol))
  names(out$par)<-names(H3_def_success$par)
  out$par$index=1:npar
  if(nDpar==npar){
    out$par$name=GetParName(dist,options$lang)
  } else {
    out$par$name=c(GetParName(dist,options$lang),
                  paste0('SystError_',1:(npar-nDpar)))
  }
  out$par$estimate=w$par
  # fill in $empirical
  ny=NROW(y)
  out$empirical=data.frame(matrix(NA,ny,empirical.ncol))
  names(out$empirical)<-names(H3_def_success$empirical)
  # Generate one possible replication for interval data
  y_oneRep=y[,1]
  ny=NROW(y)
  for (i in 1:ny){
    if(y[i,2]>y[i,1]){
      y_oneRep[i]=GenerateWithinBounds(dist=dist,par=out$par$estimate,n=1,
                                       lowerBound=y[i,1],higherBound=y[i,2])
    }
  }
  foo=sort.int(y_oneRep,index.return=TRUE)
  ixEmpirical=foo$ix
  out$empirical$y=foo$x
  out$empirical$freq=sapply(1:ny,GetEmpFreq,ny,options$FreqFormula)
  out$empirical$T=sapply(out$empirical$freq,p2T,options$p2T,options$invertT)
  out$empirical$u=sapply(out$empirical$freq,GetReducedVariate,dist)
  # fill in $pcdf
  xgrid=sapply(options$pgrid,GetQuantile,dist,w$par)
  nx=length(xgrid)
  out$pcdf=data.frame(matrix(NA,nx,pcdf.ncol))
  names(out$pcdf)<-names(H3_def_success$pcdf)
  out$pcdf$x=xgrid
  out$pcdf$pdf=sapply(xgrid,GetPdf,dist,w$par)
  out$pcdf$cdf=sapply(xgrid,GetCdf,dist,w$par)
  # fill in $quant
  pgrid=sapply(options$Tgrid,T2p,options$p2T,options$invertT)
  np=length(pgrid)
  out$quantile=data.frame(matrix(NA,np,quantile.ncol))
  names(out$quantile)<-names(H3_def_success$quantile)
  out$quantile$T=options$Tgrid
  out$quantile$p=pgrid
  out$quantile$u=sapply(pgrid,GetReducedVariate,dist)
  out$quantile$q=sapply(pgrid,GetQuantile,dist,w$par)
  # Tests
  if(do.KS){out$KS=KS(y=y_oneRep,dist=dist,par=w$par)}
  if(do.MK){out$MK=MK(y=y_oneRep)}
  if(do.Pettitt){out$Pettitt=Pettitt(y=y_oneRep)}
  
  # Uncertainties ----
  u=GetUncertainty(mcmc,burn=mcmcoptions$burn,slim=mcmcoptions$slim)
  if(u$ok==TRUE){
    # actually not ok if simulated pars encompass infinity
    if( any(is.infinite(u$sim)) ) {u$ok=FALSE}
  }
  if(u$ok==FALSE) {
    out$err=u$err
    out$message=paste("incertitude:echec:",u$message,sep="")
    out$u=Uncertainty_fail
    return(out)
  }
  # fill in $par
  out$u=u
  out$par$mean=apply(u$sim,2,mean,na.rm=TRUE)
  out$par$median=apply(u$sim,2,stats::median,na.rm=TRUE)
  out$par$sdev=apply(u$sim,2,stats::sd,na.rm=TRUE)
  out$par$IC.low=apply(u$sim,2,stats::quantile,probs=0.5*(1-options$IClevel),na.rm=TRUE)
  out$par$IC.high=apply(u$sim,2,stats::quantile,probs=1-0.5*(1-options$IClevel),na.rm=TRUE)
  # fill in $quantile
  Q=matrix(NA,NROW(u$sim),length(pgrid))
  for (i in 1:NROW(u$sim)) {
    Q[i,]=sapply(pgrid,GetQuantile,dist,u$sim[i,])
  }
  out$quantile$IC.low=apply(Q,2,stats::quantile,probs=0.5*(1-options$IClevel),na.rm=TRUE)
  out$quantile$IC.high=apply(Q,2,stats::quantile,probs=1-0.5*(1-options$IClevel),na.rm=TRUE)
  
  # Additional fields, not in regular H3 object
  out$y=y
  out$SystErrorIndex=SystErrorIndex
  out$ixEmpirical=ixEmpirical
  
  return(out)
}

#' Import HBay Configuration folder
#'
#' Imports configuration data as specified with HBay executable.
#' Returns NULL if configuration folder is not found
#'
#' @param path character, path to configuration folder. 
#' @return A list with the following components (see ?Hydro3_HBay for details):
#'     \item{y}{numeric matrix, data.}
#'     \item{dist}{character, distribution name.}
#'     \item{prior}{list of lists, prior distributions.}
#'     \item{SystErrorIndex}{integer vector, index of systematic errors.}
#'     \item{SystErrorPrior}{list of lists, prior distribution for each systematic error.}
#'     \item{options}{list, inference options.}
#'     \item{mcmcoptions}{list, MCMC options.}
#'     \item{year}{numeric vector, years.}
#' @examples
#' config=Import_HBayConfig('path/to/config')
#' if(!is.null(config)){
#'   H3=Hydro3_HBay(y=config$y,dist=config$dist,prior=config$prior,
#'                SystErrorIndex=config$SystErrorIndex,
#'                SystErrorPrior=config$SystErrorPrior,
#'                options=config$options,
#'                mcmcoptions=config$mcmcoptions)
#'   Hydro3_Plot(H3)
#' }
#' @importFrom utils read.table
#' @export
Import_HBayConfig <- function(path){
  if(!dir.exists(path)){
    message(paste0('folder not found: ',path))
    return(NULL)
  }
  folder=strsplit(path,.Platform$file.sep)[[1]]
  folder=folder[length(folder)]
  
  # Read Config_data
  foo=readLines(file.path(path,'Config_Data.txt'),n=1)
  foo=strsplit(foo,'!')[[1]][1] # remove end-of-line comment
  foosplit=strsplit(foo,'[/\\]')[[1]] # split at windows or unix separators
  if(length(foosplit)==1){ # implicit path is used
    fname=file.path(path,foosplit)
  } else if(foosplit[1]==folder){ # implicit path is used
    fname=file.path(path,do.call(file.path,as.list(foosplit[2:length(foosplit)])))
  } else { # full path is used
    fname=foo
  }

  # Read data
  D=utils::read.table(fname,header=TRUE)
  year=D[,1]
  SystErrorIndex=D[,9]
  y=matrix(NA,NROW(D),2)
  mask=D[,2]==0;y[mask,1]=D[mask,3];y[mask,2]=D[mask,3] # "equal-to" data
  mask=D[,2]==1;y[mask,1]=-Inf;     y[mask,2]=D[mask,4] # "smaller-than" data
  mask=D[,2]==2;y[mask,1]=D[mask,5];y[mask,2]=Inf       # "larger-than" data
  mask=D[,2]==3;y[mask,1]=D[mask,6];y[mask,2]=D[mask,7] # "between" data
  y[y==-9999]=NA # replace NA codes with actual NA's
  # remove NA's
  mask=stats::complete.cases(y)
  year=year[mask]
  SystErrorIndex=SystErrorIndex[mask]
  y=y[mask,]
  
  # read Config_MCMC
  D=utils::read.table(file.path(path,'Config_MCMC.txt'),header=FALSE,comment.char='!')
  mcmcoptions=mcmcoptions_def
  mcmcoptions$batch.length=D[1,1]
  mcmcoptions$batch.n=D[2,1]
  mcmcoptions$burn=D[3,1]
  mcmcoptions$slim=D[4,1]
  mcmcoptions$moverate.min=D[6,1]
  mcmcoptions$moverate.max=D[7,1]
  mcmcoptions$mult.down=D[8,1]
  mcmcoptions$mult.up=D[9,1]
  mcmcoptions$mult=D[10,1]
  mcmcoptions$eps=D[11,1]
  
  # read Config_Inference
  D=utils::read.table(file.path(path,'Config_Inference.txt'),header=FALSE,comment.char='!',blank.lines.skip=FALSE)
  dist=as.character(D[1,1])
  prior=list()
  prior[[1]]=readOnePrior(D[2:4,1])
  prior[[2]]=readOnePrior(D[5:7,1])
  if(dist=='GEV'){
    prior[[3]]=readOnePrior(D[8:10,1])
  }
  
  # read Config_SystematicErrors
  D=utils::read.table(file.path(path,'Config_SystematicErrors.txt'),header=FALSE,comment.char='!',blank.lines.skip=FALSE)
  nK=max(SystErrorIndex)
  SystErrorPrior=list()
  if(nK>0){
    for(i in 1:nK){
      SystErrorPrior[[i]]=readOnePrior(D[3*(i-1)+(1:3),1])
    }
  }
  
  # read Config_ResultOptions
  D=utils::read.table(file.path(path,'Config_ResultOptions.txt'),header=FALSE,comment.char='!',blank.lines.skip=FALSE)
  options=options_def
  pmin=as.numeric(D[1,1])
  pmax=as.numeric(D[2,1])
  pn=as.numeric(D[3,1])
  options$pgrid=seq(pmin,pmax,length.out=pn)
  options$FreqFormula=as.character(D[4,1])
  options$IClevel=as.numeric(D[9,1])
  options$p2T=as.numeric(D[10,1])
  
  return(list(y=y,dist=dist,prior=prior,
              SystErrorIndex=SystErrorIndex,SystErrorPrior=SystErrorPrior,
              options=options,mcmcoptions=mcmcoptions,
              year=year))
}
#' HBay plot
#'
#' Plot summarizing the results of Hydro3_HBay()
#'
#' @param H3 list, resulting from a call to Hydro3_HBay()
#' @param curve_color color, color used for quantile curve
#' @return nothing (just creates a plot)
#' @examples
#' set.seed(98765)
#' n=50;n_censored=30
#' y0=Generate('GEV',c(100,50,-0.2),n)
#' y=cbind(y0,y0)
#' # Mimics censoring between 0 and 300 for first n_censored years
#' y[1:n_censored,1][y0[1:n_censored]<300]=0
#' y[1:n_censored,2][y0[1:n_censored]<300]=300
#' plot(y[,1]);points(y[,2])
#' # Systematic errors
#' SystErrorIndex=c(rep(1,n_censored),rep(2,n-n_censored))
#' SystErrorPrior=list(list(dist="Triangle",par=c(1,0.7,1.3)),
#'                     list(dist="Triangle",par=c(1,0.95,1.05)))
#' # Priors on GEV parameters
#' prior=list(list(dist="FlatPrior",par=NULL),
#'            list(dist="FlatPrior",par=NULL),
#'            list(dist="Normal",par=c(0,0.25)))
#' # Handle MCMC options
#' # The values below aim at making this example fast to run.
#' # In practice, it is recommended to use the default values
#' # (batch.length=100,batch.n=100) or larger.
#' mcmcoptions=mcmcoptions_def
#' mcmcoptions$batch.length=25
#' mcmcoptions$batch.n=20
#' # Go!
#' H3=Hydro3_HBay(y=y,dist='GEV',prior=prior,
#'                SystErrorIndex=SystErrorIndex,
#'                SystErrorPrior=SystErrorPrior,
#'                mcmcoptions=mcmcoptions) 
#' # HBay plot
#' HBay_Plot(H3)
#' @importFrom graphics plot lines points legend
#' @export
HBay_Plot <- function(H3,curve_color='black'){
  mini=min(H3$quantile$IC.low)
  maxi=max(H3$quantile$IC.high)
  graphics::plot(H3$quantile$T,H3$quantile$q,type='l',lwd=2,ylim=c(mini,maxi),log='x',
       xlab='T',ylab='Q(T)',col=curve_color)
  graphics::lines(H3$quantile$T,H3$quantile$IC.low,lty=2,col=curve_color)
  graphics::lines(H3$quantile$T,H3$quantile$IC.high,lty=2,col=curve_color)
  graphics::segments(H3$empirical$T,H3$y[H3$ixEmpirical,1],
                     H3$empirical$T,H3$y[H3$ixEmpirical,2],
                     col=H3$SystErrorIndex[H3$ixEmpirical]+2)
  graphics::points(H3$empirical$T,H3$empirical$y,
                   pch=19,col=H3$SystErrorIndex[H3$ixEmpirical]+2)
  graphics::legend('topleft',paste('Systematic error',0:max(H3$SystErrorIndex)),
         pch=19,col = 2+(0:max(H3$SystErrorIndex)))
}

#****************************
# Fonctions privées ----
#****************************

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
      z=matrix(y[indx,],ncol=2) # make sure this stays a 2-column matrix
      dt=dataType[indx]
      if(k==0) {g=1} else {g=gamma[k]}
      distPar=c(theta[1]/g,theta[2]/g)
      if(dist=='GEV'){distPar=c(distPar,theta[3])}
      # Perfectly-known data
      w=matrix(z[dt==0,],ncol=2)
      if(NROW(w)>0){
        lkh=lkh+GetLogLikelihood(distPar,w[,1],dist)
      }
      # Interval data
      w=matrix(z[dt==1,],ncol=2)
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

readOnePrior <- function(threeLines){
  dist=as.character(threeLines[2])
  if(dist=='Gaussian') dist='Normal'
  if(dist %in% c('FlatPrior+','FlatPrior-')) dist='FlatPrior'
  param=strsplit(as.character(threeLines[3]),',')[[1]]
  if(length(param)==0){
    pars=NULL
  } else {
    pars=as.numeric(param)
  }
  return(list(dist=dist,par=pars))
}
