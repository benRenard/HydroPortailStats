#****************************
# HydRoStat v1.0
# Copyright 2017 Irstea, IDDN.FR.001.460013.000.S.C.2017.000.20700
# Author: Benjamin Renard

# This program is free software: you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the Free 
# Software Foundation, either version 3 of the License, or (at your option) 
# any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
# more details.

# A copy of the GNU General Public License is provided within this 
# distribution.
# See also <https://www.gnu.org/licenses/>.
#****************************

#~******************************************************************************
#~* OBJET: Calculs statistiques de la Banque HYDRO3
#~******************************************************************************
#~* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
#~******************************************************************************
#~* CREE/MODIFIE: 09/07/2015
#~******************************************************************************
#~* PRINCIPALES FONCTIONS
#~*    1. Hydro3_Estimation: calcule tous les éléments nécessaires à HYDRO3
#~*    2. Hydro3_Plot: trace les résultats de la fonction ci-dessus
#~******************************************************************************
#~* REF.: XXX
#~******************************************************************************
#~* A FAIRE:
#~******************************************************************************
#~* COMMENTAIRES: XXX
#~******************************************************************************

#****************************
# Déclaration de variables globales ----
#****************************

# Méthodes d'estimation par défaut
Emeth_def="LMOM" # "MOM","LMOM","BAY", "ML"
Umeth_def="PBOOT" # "BOOT","PBOOT","BAY","ML", NONE"

# objet "H3" par défaut (contient tous les résultats nécessaires à HYDRO3)
empirical_def=data.frame(y=NA,freq=NA,T=NA,u=NA)
pcdf_def=data.frame(x=NA,pdf=NA,cdf=NA)
quantile_def=data.frame(T=NA,p=NA,u=NA,q=NA,IC.low=NA,IC.high=NA)
par_def=data.frame(index=NA,name=NA,estimate=NA,IC.low=NA,IC.high=NA,mean=NA,median=NA,sdev=NA)
test_def=list(pval=NA,stat=NA,xtra=NA)
H3_def_fail=list(empirical=empirical_def,pcdf=pcdf_def,quantile=quantile_def,
                 par=par_def,KS=test_def,MK=test_def,Pettitt=test_def,
                 u=Uncertainty_fail,ok=FALSE,err=666,message="fatal",dist=NA)
H3_def_success=list(empirical=empirical_def,pcdf=pcdf_def,quantile=quantile_def,
                 par=par_def,KS=test_def,MK=test_def,Pettitt=test_def,
                 u=Uncertainty_success,ok=TRUE,err=0,message="ok",dist=NA)

# objets "options" par défaut
options_def=list(FreqFormula="Hazen",pgrid=seq(0.001,0.999,0.001),
                 Tgrid=c(seq(1.1,1.9,0.1),seq(2,9,1),seq(10,90,10),seq(100,1000,100)),
                 IClevel=0.9,p2T=1,invertT=FALSE,splitZeros=FALSE,lang="fr",
                 nsim=nsim_def)
mcmcoptions_def=list(mult=0.1,eps=0.1,batch.length=100,batch.n=100,
                     moverate.min=0.1,moverate.max=0.5,mult.down=0.9, mult.up=1.1,
                     burn=0.5,slim=5)

#****************************
# Fonctions principales ----
#****************************

#' Hydro3 estimation
#'
#' Main estimation function used in the HydroPortail.
#' In short, this function estimates a distribution and the associated uncertainty, 
#' and returns all needed information to display and plot the results 
#' (parameter estimates, quantile curves, etc.)
#'
#' @param y numeric vector, data.
#' @param dist character, distribution name. See dataset distInfo for a description of available distributions.
#'     In particular, type names(distInfo) for the list of available distributions, and 
#'     distInfo[['GEV']] for more information on a particular distribution (here, GEV).
#' @param Emeth character, estimation method. Default is 'LMOM' (L-Moments), available: 
#'     'MOM' (Moments), 'ML' (Maximum Likelihood), 'BAY' (Bayesian).
#' @param Umeth character, uncertainty quantification method. Default is 'PBOOT' (Parametric bootstrap), available: 
#'     'BOOT' (Bootstrap, not recommended), 'NONE', 'ML' (only usable when Emeth='ML' as well), 
#'     and 'BAY' (the only usable method when Emeth='BAY').
#' @param options list, options, see details below.
#' @param mcmcoptions list, MCMC options, see details below.
#' @param prior list, prior distributions, only used when Emeth='BAY'. See ?GetEstimate_BAY for details.
#' @param do.KS,do.MK,do.Pettitt logical, perform KS/MK/Pettitt tests?
#' @return A list with the following components:
#'     \item{dist}{character, estimated distribution.}
#'     \item{ok}{logical, did estimation succeed?}
#'     \item{err}{integer, error code (0 if ok).}
#'     \item{message}{error message.}
#'     \item{empirical}{data frame, sorted data and empirical estimates 
#'         (nonexceedance frequency, return period and reduced variate) }
#'     \item{pcdf}{data frame, estimated pdf and cdf}
#'     \item{quantile}{data frame, estimated quantiles and uncertainty intervals}
#'     \item{par}{data frame, estimated parameters and uncertainty intervals}
#'     \item{KS}{list, result of the Kolmogorov-Smirnov test, see ?KS}
#'     \item{MK}{list, result of the Mann-Kendall test, see ?MK}
#'     \item{Pettitt}{list, result of the Pettitt test, see ?Pettitt}
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
#'         p=1-1/(p2T*T). In general p2T=1 but for a peak-over-threshold approach leading to say 3 events
#'         per year on average, p2T=3.}
#'     \item{invertT, logical, when invertT=TRUE, LARGE return periods correspond to SMALL data values.
#'         This is typically used for low-flow statistics.}
#'     \item{splitZeros, logical, when splitZeros=TRUE zero and negative values are removed from the data y before 
#'         estimating the distribution,and are used to estimate the probability of zeros p0. This is 
#'         typically used for low-flow statistics to estimate the probability of zero streamflow.}
#'     \item{lang, chanracter, language ('fr' or 'en').}
#'     \item{nsim, integer, number of replicated parameters representing uncertainty.}
#'     }
#'     The argument 'mcmcoptions' is only used when Emeth='BAY' and is a list controlling MCMC properties: 
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
#' y=stats::rnorm(50)
#' H3=Hydro3_Estimation(y,'Normal')
#' H3=Hydro3_Estimation(y,'GEV',Emeth='ML',Umeth='ML')
#' @importFrom stats sd median quantile
#' @export
Hydro3_Estimation<-function(y,dist,Emeth=Emeth_def,Umeth=Umeth_def,
                            options=options_def,mcmcoptions=mcmcoptions_def,
                            prior=GetDefaultPrior(GetParNumber(dist)),
                            do.KS=TRUE,do.MK=TRUE,do.Pettitt=TRUE){
  #^******************************************************************************
  #^* OBJET: Principale fonction d'estimation d'HYDRO3
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 09/07/2015
  #^******************************************************************************
  #^* IN
  #^*    1. [real vector] y, données 
  #^*    2. [character] dist, distribution supposée 
  #^*    3. [character] Emeth, methode d'estimation 
  #^*    4. [character] Umeth, methode de quantification de l'incertitude 
  #^*    5. [list] options, objet "options"
  #^*    5. [list] mcmcmoptions, options pour les simulations MCMC
  #^*    5. [list] prior, a prioris pour chaque paramètre
  #^*    6. [logical] do.KS, perform KS test?
  #^*    7. [logical] do.MK, perform MK test? (a bit long)
  #^*    8. [logical] do.Pettitt, perform Pettitt test? (a bit long)
  #^* OUT
  #^*    1. objet "H3" contenant tous les résultats utiles  
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************
  par.ncol=8
  empirical.ncol=4
  pcdf.ncol=3
  quantile.ncol=6
  
  out=H3_def_success
  
  ###########################
  # Handle values <= 0
  ###########################
  ny=length(y)
  if(options$splitZeros){
    z=y[y>0]
    n0=sum(y<=0)
    p0=n0/ny
  } else {
    z=y;
    n0=0
    p0=0
  }

  ###########################
  # Estimation des paramètres
  ###########################
  if(Emeth=="BAY"){
    par0=GetEstimate_ROUGH(z,dist)
    if(par0$ok==FALSE) {out=H3_def_fail;out$message="estimation:echec";return(out)}
    mcmc=GetEstimate_BAY(z,dist,prior,par0$par,
                         mult=mcmcoptions$mult,eps=mcmcoptions$eps,
                         batch.length=mcmcoptions$batch.length,batch.n=mcmcoptions$batch.n,
                         moverate.min=mcmcoptions$moverate.min,moverate.max=mcmcoptions$moverate.max,
                         mult.down=mcmcoptions$mult.down,mult.up=mcmcoptions$mult.up)
    w=GetMode(mcmc)
  } else {
    w=GetEstimate_OneSample(z,dist,Emeth)
  }
  if(w$ok==FALSE) {out=H3_def_fail;out$message="estimation:echec";return(out)}
  out$dist=dist
  # fill in $par
  npar=GetParNumber(dist)
  out$par=data.frame(matrix(NA,npar+1,par.ncol))
  names(out$par)<-names(H3_def_success$par)
  out$par$index=1:(npar+1)
  out$par$name=c("p0",GetParName(dist,options$lang))
  out$par$estimate=c(p0,w$par)
  # fill in $empirical
  ny=length(y)
  out$empirical=data.frame(matrix(NA,ny,empirical.ncol))
  names(out$empirical)<-names(H3_def_success$empirical)
  out$empirical$y=sort(y)
  out$empirical$freq=sapply(1:ny,GetEmpFreq,ny,options$FreqFormula)
  out$empirical$T=sapply(out$empirical$freq,p2T,options$p2T,options$invertT)
  out$empirical$u=sapply(out$empirical$freq,GetReducedVariate,dist)
  # fill in $pcdf
  mask=(options$pgrid-p0<=sqrt(.Machine$double.eps)) # a complicated alternative to options$pgrid<=p0 to avoid rounding errors...
  xgrid=c(rep(0,sum(mask)),sapply((options$pgrid[!mask]-p0)/(1-p0),GetQuantile,dist,w$par))
  if(options$splitZeros) {xgrid[xgrid<=0]=0}
  nx=length(xgrid)
  out$pcdf=data.frame(matrix(NA,nx,pcdf.ncol))
  names(out$pcdf)<-names(H3_def_success$pcdf)
  out$pcdf$x=xgrid
  out$pcdf$pdf=c(rep(p0,sum(mask)),(1-p0)*sapply(xgrid[!mask],GetPdf,dist,w$par))
  out$pcdf$cdf=c(rep(p0,sum(mask)),p0+(1-p0)*sapply(xgrid[!mask],GetCdf,dist,w$par))  
  # fill in $quant
  pgrid=sapply(options$Tgrid,T2p,options$p2T,options$invertT)
  mask=(pgrid<=p0)
  np=length(pgrid)
  out$quantile=data.frame(matrix(NA,np,quantile.ncol))
  names(out$quantile)<-names(H3_def_success$quantile)
  out$quantile$T=options$Tgrid
  out$quantile$p=pgrid
  out$quantile$u=sapply(pgrid,GetReducedVariate,dist)
  out$quantile$q=0
  out$quantile$q[!mask]=sapply((pgrid[!mask]-p0)/(1-p0),GetQuantile,dist,w$par)
  if(options$splitZeros) {out$quantile$q[out$quantile$q<=0]=0}
  # Tests
  if(do.KS){out$KS=KS(y=z,dist=dist,par=w$par)}
  if(do.MK){out$MK=MK(y=y)}
  if(do.Pettitt){out$Pettitt=Pettitt(y=y)}
  
  ##################################
  # quantification de l'incertitude
  ##################################
  # combinaisons non autorisées
  if( Umeth=="NONE"|(Umeth=="BAY" & Emeth!="BAY" )|(Umeth=="ML" & Emeth!="ML" )) {
    out$u=Uncertainty_fail
    out$u$message=paste("Fatal: combination Emeth/Umeth non autorisee:",Emeth,"/",Umeth,sep="")
    return(out)
  }
  u=switch(Umeth,
           ML=GetUncertainty_ML(y,dist,w$par,nsim=options$nsim),
           BOOT=GetUncertainty_Bootstrap(y,dist,w$par,Emeth,type="NP",nsim=options$nsim),
           PBOOT=GetUncertainty_Bootstrap(y,dist,w$par,Emeth,type="P",nsim=options$nsim),
           BAY=GetUncertainty(mcmc,burn=mcmcoptions$burn,slim=mcmcoptions$slim),
           Uncertainty_fail)
  if(u$ok==TRUE){
    # actually not ok if simulated pars encompass infinity
    if( any(is.infinite(u$sim)) ) {u$ok=FALSE}
  }
  if(u$ok==FALSE) {out$err=u$err
               out$message=paste("incertitude:echec:",u$message,sep="")
               out$u=Uncertainty_fail
               return(out)}
  # fill in $par
  out$u=u
  out$par$mean=c(p0,apply(u$sim,2,mean,na.rm=TRUE))
  out$par$median=c(p0,apply(u$sim,2,stats::median,na.rm=TRUE))
  out$par$sdev=c(0,apply(u$sim,2,stats::sd,na.rm=TRUE))
  out$par$IC.low=c(p0,apply(u$sim,2,stats::quantile,probs=0.5*(1-options$IClevel),na.rm=TRUE))
  out$par$IC.high=c(p0,apply(u$sim,2,stats::quantile,probs=1-0.5*(1-options$IClevel),na.rm=TRUE))
  # fill in $quantile
  Q=matrix(NA,options$nsim,length(pgrid))
  for (i in 1:options$nsim) {
    Q[i,mask]=0
    Q[i,!mask]=sapply((pgrid[!mask]-p0)/(1-p0),GetQuantile,dist,u$sim[i,])
    if(options$splitZeros) {Q[i,Q[i,]<=0]=0}
  }
  out$quantile$IC.low=apply(Q,2,stats::quantile,probs=0.5*(1-options$IClevel),na.rm=TRUE)
  out$quantile$IC.high=apply(Q,2,stats::quantile,probs=1-0.5*(1-options$IClevel),na.rm=TRUE)
  
  return(out)
}

#' Hydro3 plot
#'
#' Plot summarizing the results of Hydro3_Estimation()
#'
#' @param H3 list, resulting from a call to Hydro3_Estimation()
#' @param useU logical, use reduced variate u rather than return period T in plots?
#' @param lwd,cex.lab,cex.axis,pch numeric, graphical parameters, see ?graphics::par
#' @param col character, graphical parameter (points color)
#' @return nothing (just creates a plot)
#' @examples
#' y=stats::rnorm(50)
#' H3=Hydro3_Estimation(y,'Normal')
#' Hydro3_Plot(H3)
#' @importFrom graphics plot lines points rug par text
#' @export
Hydro3_Plot<-function(H3,useU=FALSE,lwd=2,cex.lab=2,cex.axis=1.3,pch=19,col="red"){
  #^******************************************************************************
  #^* OBJET: Tracé des figures HYDRO3
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [list] H3, objet contenant tous les résultats de l'estimation 
  #^*    2. [logical] useU, utiliser la variable réduite u plutôt que T?
  #^* OUT
  #^*    1. [plot] 
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************
  oldpar=graphics::par(no.readonly=TRUE) # save original par settings
  on.exit(graphics::par(oldpar)) # reverts to original settings on exit
  graphics::par(mfrow=c(2,2),mar=c(4,5,2,1)+0.2)
  # par plot
  graphics::plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",main=H3$dist)
  npar=nrow(H3$par)
  for (i in 1:npar){
    x=0;y=1-i/npar+(0.5/npar)
    graphics::text(x,y,paste(H3$par$name[i]," = ",format(H3$par$estimate[i],digits=2,nsmall=2),
                   " [",format(H3$par$IC.low[i],digits=2,nsmall=2),";",
                   format(H3$par$IC.high[i],digits=2,nsmall=2),"]",sep=""),pos=4,cex=cex.axis)
  }
  # pdf plot
  graphics::plot(H3$pcdf$x,H3$pcdf$pdf,type="l",xlab="y",ylab="f(y)",lwd=lwd,cex.lab=cex.lab,cex.axis=cex.axis)
  graphics::rug(H3$empirical$y,col=col)
  # cdf plot
  graphics::plot(H3$pcdf$x,H3$pcdf$cdf,type="l",xlab="y",ylab="F(y)",lwd=lwd,cex.lab=cex.lab,cex.axis=cex.axis,ylim=c(0,1))
  graphics::points(H3$empirical$y,H3$empirical$freq,pch=pch,col=col)
  # quantile plot
  if(H3$u$ok) {
    yl=c(min(H3$empirical$y,H3$quantile$IC.low),max(H3$empirical$y,H3$quantile$IC.high))
  } else{
    yl=c(min(H3$empirical$y,H3$quantile$q),max(H3$empirical$y,H3$quantile$q))
  }
  if(useU) {
    x=H3$quantile$u
    dolog=""
    xla="u"
    yla="q(u)"
  } else {
    x=H3$quantile$T
    dolog="x"
    xla="T"
    yla="q(T)"
  }
  graphics::plot(x,H3$quantile$q,type="l",log=dolog,xlab=xla,ylab=yla,
       lwd=lwd,cex.lab=cex.lab,cex.axis=cex.axis,ylim=yl)
  if(H3$u$ok) graphics::lines(x,H3$quantile$IC.low,lty=2)
  if(H3$u$ok) graphics::lines(x,H3$quantile$IC.high,lty=2)
  if(useU) {x=H3$empirical$u} else {x=H3$empirical$T}
  graphics::points(x,H3$empirical$y,pch=pch,col=col)
}

#' Get quantile from return period
#'
#' Compute the T-quantile from the results of Hydro3_Estimation() 
#'
#' @param RP numeric, return period
#' @param H3 list, resulting from a call to Hydro3_Estimation()
#' @param options list, see ?Hydro3_Estimation
#' @return A list with the following components:
#'     \item{q}{numeric, quantile}
#'     \item{IC}{numeric vector, uncertainty interval}
#' @examples
#' y=stats::rnorm(50)
#' H3=Hydro3_Estimation(y,'Normal')
#' GetQfromT(100,H3)
#' @export
GetQfromT<-function(RP,H3,options=options_def){
  #^******************************************************************************
  #^* OBJET: Calcul d'un quantile pour une période de retour donnée
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 28/08/2017
  #^******************************************************************************
  #^* IN
  #^*    1. [real] RP, période de retour
  #^*    2. [list] H3, objet H3 résultant d'un appel de Hydro3_Estimation
  #^*    3. [list] options, objet "options"
  #^* OUT
  #^*    1. liste contenant le quantile et son intervalle d'incertitude
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************
  fail=list(q=NA,IC=c(NA,NA))
  # check whether RP is valid
  if(RP<=0) {return(fail)}
  p=T2p(RP,factor=options$p2T,invert=options$invertT)
  if(is.na(p)){return(fail)}
  # compute q
  p0=H3$par$estimate[1]
  q=GetQ_engine(par=H3$par$estimate[-1],p=p,p0=p0,dist=H3$dist,options=options)
  # compute IC
  uq=apply(H3$u$sim,1,GetQ_engine,p=p,p0=p0,dist=H3$dist,options=options)
  IC=as.numeric(stats::quantile(uq,probs=c(0.5*(1-options$IClevel),1-0.5*(1-options$IClevel)),na.rm=TRUE))
  return(list(q=q,IC=IC))
}

#' Get return period from value
#'
#' Compute the return period associated with a value from the results of Hydro3_Estimation() 
#'
#' @param q numeric, value
#' @param H3 list, resulting from a call to Hydro3_Estimation()
#' @param options list, see ?Hydro3_Estimation
#' @return A list with the following components:
#'     \item{RP}{numeric, return period}
#'     \item{IC}{numeric vector, uncertainty interval}
#' @examples
#' y=stats::rnorm(50)
#' H3=Hydro3_Estimation(y,'Normal')
#' GetTfromQ(3,H3)
#' @export
GetTfromQ<-function(q,H3,options=options_def){
  #^******************************************************************************
  #^* OBJET: Calcul d'une période de retour pour une valeur donnée
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 28/08/2017
  #^******************************************************************************
  #^* IN
  #^*    1. [real] q, valeur
  #^*    2. [list] H3, objet H3 résultant d'un appel de Hydro3_Estimation
  #^*    3. [list] options, objet "options"
  #^* OUT
  #^*    1. liste contenant la période de retour et son intervalle d'incertitude
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************
  # compute q
  p0=H3$par$estimate[1]
  RP=GetT_engine(par=H3$par$estimate[-1],q=q,p0=p0,dist=H3$dist,options=options)
  # compute IC
  uT=apply(H3$u$sim,1,GetT_engine,q=q,p0=p0,dist=H3$dist,options=options)
  IC=as.numeric(stats::quantile(uT,probs=c(0.5*(1-options$IClevel),1-0.5*(1-options$IClevel)),na.rm=TRUE))
  return(list(RP=RP,IC=IC))
}

#****************************
# Fonctions privées ----
#****************************

p2T<-function(p,factor=1,invert=FALSE){
  if(p>=1) {if(invert) {return(-Inf)} else {return(Inf)}}
  if(p<=0) {if(invert) {return(Inf)} else {return(-Inf)}}
  if(invert) {
    RP=1/(factor*p)
  } else {
    RP=1/(factor*(1-p))
  }
  return(RP)
}

T2p<-function(RP,factor=1,invert=FALSE){
  if(RP*factor<1) {return(NA)}
  if(invert) {
    p=1/(factor*RP)
  } else {
    p=1-1/(factor*RP)
  }
  return(p)
}

GetQ_engine<-function(par,p,p0,dist,options){
  if(options$splitZeros) {
    if(p<p0){q=0} else{
      q=GetQuantile(p=(p-p0)/1-p0,dist=dist,par=par)
    }
    q=max(q,0)
  } else {
    q=GetQuantile(p=p,dist=dist,par=par)
  }
  return(q)
}

GetT_engine<-function(par,q,p0,dist,options){
  p=GetCdf(y=q,dist=dist,par=par)
  RP=p2T(max(p,p0),options$p2T,options$invertT)
  return(RP)
}
