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
#~* OBJET: Estimation Bayesienne
#~******************************************************************************
#~* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
#~******************************************************************************
#~* CREE/MODIFIE: 09/07/2015
#~******************************************************************************
#~* PRINCIPALES FONCTIONS
#~*    1. XXX
#~*    2. XXX
#~******************************************************************************
#~* REF.: XXX
#~******************************************************************************
#~* A FAIRE: XXX
#~******************************************************************************
#~* COMMENTAIRES: XXX
#~******************************************************************************

#****************************
# Fonctions principales ----
#****************************

#' Bayesian estimation of a distribution
#'
#' Returns MCMC samples from the posterior distribution.
#'
#' @param y numeric vector, data
#' @param dist character, distribution name
#' @param prior list of lists, prior distributions. For each parameter to be estimated, the prior
#'    is a list of the form pr=list(dist=..., par=...). See example below.
#' @param par0 numeric vector, initial parameter guess. You may use GetEstimate_ROUGH().
#' @param mult numeric, initial jump standard deviations are set to mult * abs(par0)
#' @param eps numeric, where par0 is zero, initial jump standard deviations are set to eps (to avoid jumps of size zero)
#' @param batch.length integer, MCMC parameter: length of each non-adaptive batch 
#' @param batch.n integer, MCMC parameter: number of batches (= adaptation period). Total number of simulations is nsim=batch.n*batch.length
#' @param moverate.min numeric in (0;1), MCMC parameter: lower bound for the desired move rate interval
#' @param moverate.max numeric in (0;1), MCMC parameter: upper bound for the desired move rate interval
#' @param mult.down numeric in (0;1), MCMC parameter: multiplication factor used to decrease jump size when move rate is too low.
#' @param mult.up numeric (>1, avoid 1/mult.down), MCMC parameter: multiplication factor used to increase jump size when move rate is too high.
#' @return A list with the following components:
#'     \item{x}{numeric matrix nsim*length(x0), MCMC simulations}
#'     \item{fx}{numeric vector, corresponding values f(x)}
#' @examples
#' y=c(9.2,9.5,11.4,9.5,9.4,9.6,10.5,11.1,10.5,10.4)
#' prior1=list(dist='FlatPrior',par=NULL)
#' prior2=list(dist='LogNormal',par=c(1,1))
#' prior3=list(dist='Normal',par=c(0,0.25))
#' prior=list(prior1,prior2,prior3)
#' par0=GetEstimate_ROUGH(y,'GEV')$par
#' mcmc=GetEstimate_BAY(y,'GEV',prior,par0,batch.length=50,batch.n=50)
#' graphicalpar=par(mfrow=c(2,3))
#' plot(mcmc$x[,1],type='l'); plot(mcmc$x[,2],type='l'); plot(mcmc$x[,3],type='l')
#' hist(mcmc$x[,1]); hist(mcmc$x[,2]); hist(mcmc$x[,3])
#' par(graphicalpar)
#' @export
GetEstimate_BAY<-function(y,dist,prior,par0,
                          mult=0.1,eps=0.1,
                          batch.length=100,batch.n=100,
                          moverate.min=0.1,moverate.max=0.5,
                          mult.down=0.9, mult.up=1.1){
  #^******************************************************************************
  #^* OBJET: Retourne le resultat d'une estimation BAY+MCMC
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 09/07/2015
  #^******************************************************************************
  #^* IN
  #^*    1. [real] y, vecteur des données 
  #^*    2. [character] dist, nom de la distribution 
  #^*    3. [object] prior, liste des a prioris pour chaque paramètre 
  #^*    4. [real] par0, point de départ pour l'algo MCMC 
  #^* OUT
  #^*    1. [list] Une liste comprenant: 
  #^*         $x [real matrix]: simulations MCMC
  #^*         $fx [real vector]: valeur de f correspondant à chaque simu MCMC
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************  
  
  sdjump=mult*par0;sdjump[sdjump==0]=eps
  out=Metropolis_OAAT_adaptive(GetLogPost,par0,sdjump,
                               y,dist,prior,
                               batch.length=batch.length,batch.n=batch.n,
                               moverate.min=moverate.min,moverate.max=moverate.max,
                               mult.down=mult.down, mult.up=mult.up)
  return(out)
}

#****************************
# Fonctions privées ----
#****************************

GetLogPrior<-function(par,prior){
  #^******************************************************************************
  #^* OBJET: Retourne la densité a priori 
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 09/07/2015
  #^******************************************************************************
  #^* IN
  #^*    1. [real] par, vecteur de paramètres 
  #^*    2. [object] prior, liste des a prioris pour chaque paramètre 
  #^* OUT
  #^*    1. [real] la log-densité a priori
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: si les paramètres sont impossibles, retourne NA
  #^*               si l'a priori est nul, retourne -Inf
  #^****************************************************************************** 
  n=length(par)
  LP=0
  for(i in 1:n){
    w=GetPdf(par[i],prior[[i]]$dist,prior[[i]]$par,log=TRUE)
    if(is.na(w)) return(NA)
    if(w==-Inf) return(-Inf)
    LP=LP+w
  }  
  return(LP)
}

GetLogPost<-function(par,y,dist,prior){
  #^******************************************************************************
  #^* OBJET: Retourne la distribution a posteriori pour la distribution 'dist' de  
  #^*        paramètres 'par'
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 09/07/2015
  #^******************************************************************************
  #^* IN
  #^*    1. [real] par, vecteur de paramètres 
  #^*    2. [real] y, vecteur des données 
  #^*    3. [character] dist, nom de la distribution 
  #^*    4. [object] prior, liste des a prioris pour chaque paramètre 
  #^* OUT
  #^*    1. [real] la log-densité a posteriori
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: si les paramètres sont impossibles, retourne NA
  #^*               si la vraisemblance ou l'a priori est nulle, retourne -Inf
  #^****************************************************************************** 
  LL=GetLogLikelihood(par,y,dist)
  if(is.na(LL)) {return(NA)}
  if(LL==-Inf) {return(-Inf)}
  LP=GetLogPrior(par,prior)
  if(is.na(LP)) {return(NA)}
  if(LP==-Inf) {return(-Inf)}
  return(LL+LP)  
}

GetDefaultPrior<-function(n){
  p=list()
  for(i in 1:n){
    p[[i]]=list(dist="FlatPrior",par=NULL)    
  }
  return(p)
}
