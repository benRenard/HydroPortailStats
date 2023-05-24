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

#' Adaptive One-At-A-Time Metropolis sampler
#'
#' Performs nsim iterations of the Adaptive version of the OAAT Metropolis sampler
#' (see ?Metropolis_OAAT).
#' Adaptation is performed by monitoring move rates every batch.length iterations, and 
#' increasing / decreasing the jump standard deviation if the move rate is not within specified bounds. 
#'
#' @param f function, log-pdf of the target distribution
#' @param x0 numeric vector, starting point
#' @param sdjump numeric vector, initial standard deviation of the Gaussian jump for each component
#' @param ... other arguments passed to f
#' @param batch.length integer, length of each non-adaptive batch 
#' @param batch.n integer, number of batches (= adaptation period). Total number of simulations is nsim=batch.n*batch.length
#' @param moverate.min numeric in (0;1), lower bound for the desired move rate interval
#' @param moverate.max numeric in (0;1), upper bound for the desired move rate interval
#' @param mult.down numeric in (0;1), multiplication factor used to decrease jump size when move rate is too low.
#' @param mult.up numeric (>1, avoid 1/mult.down) multiplication factor used to increase jump size when move rate is too high.
#' @return A list with the following components:
#'     \item{x}{numeric matrix nsim*length(x0), MCMC simulations}
#'     \item{fx}{numeric vector, corresponding values f(x)}
#' @examples
#' # Bivariate target distribution: beta(0.8,0.4) X exp(1)
#' f=function(x){stats::dbeta(x[1],0.8,0.4,log=TRUE)+stats::dexp(x[2],log=TRUE)}
#' x0=c(0.5,2)
#' sdjump=c(0.5,1)
#' mcmc=Metropolis_OAAT_adaptive(f,x0,sdjump)
#' graphicalpar=par(mfrow=c(1,3))
#' plot(mcmc$x);hist(mcmc$x[,1]); hist(mcmc$x[,2])
#' par(graphicalpar)
#' @export
Metropolis_OAAT_adaptive<-function(f,x0,sdjump,...,
                                   batch.length=100,batch.n=100,
                                   moverate.min=0.1,moverate.max=0.5,
                                   mult.down=0.9, mult.up=1.1){
  #^******************************************************************************
  #^* OBJET: Algorithme MCMC Metropolis_OAAT (One At A Time) dans sa version adaptative
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 09/07/2015
  #^******************************************************************************
  #^* IN
  #^*    1. [function] f, log-densité de la distribution cible
  #^*    2. [real vector] x0, point de départ
  #^*    3. [real vector] sdjump, écart-type des sauts pour chaque composante
  #^*    5. [whatever] ... extra arguments passed to f
  #^* OUT
  #^*    1. [list] Une liste comprenant: 
  #^*         $x [real matrix]: simulations MCMC
  #^*         $fx [real vector]: valeur de f correspondant à chaque simu MCMC
  #^******************************************************************************
  #^* REF.: Renard, B., V. Garreta, and M. Lang (2006), An application of Bayesian 
  #^*       analysis and MCMC methods to the estimation of a regional trend in 
  #^*       annual maxima, Water Resources Research, 42(12).
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************
  n=length(x0)
  x=matrix(NA,batch.n*batch.length,n)
  fx=matrix(NA,batch.n*batch.length,1)
  sdj=sdjump
  xini=x0
  for(i in 1:batch.n){
    w<-Metropolis_OAAT(f=f,x0=xini,nsim=batch.length,sdjump=sdj,...)
    indx=( (i-1)*batch.length+1 ) : (i*batch.length) 
    x[indx,]=w$x
    fx[indx]=w$fx
    xini=w$x[batch.length,]
    for(j in 1:n){
      if(w$moverate[j]<moverate.min){
        sdj[j]=sdj[j]*mult.down
      } else if(w$moverate[j]>moverate.max){
        sdj[j]=sdj[j]*mult.up
      }
    }
  }
  return(list(x=x,fx=fx))
}

#' One-At-A-Time Metropolis sampler
#'
#' Performs nsim iterations of the OAAT Metropolis sampler
#' (simulated vector is updated one component at a time).
#'  a.k.a block Metropolis sampler with blocks of length one.
#'  Sometimes also called 'Metropolis-within-Gibbs'.
#'
#' @param f function, log-pdf of the target distribution
#' @param x0 numeric vector, starting point
#' @param nsim integer, number of simulations
#' @param sdjump numeric vector, standard deviation of the Gaussian jump for each component
#' @param ... other arguments passed to f
#' @return A list with the following components:
#'     \item{x}{numeric matrix nsim*length(x0), MCMC simulations}
#'     \item{fx}{numeric vector, corresponding values f(x)}
#'     \item{moverate}{numeric vector, move rate associated with each component}
#' @examples
#' # Bivariate target distribution: beta(0.8,0.4) X exp(1)
#' f=function(x){stats::dbeta(x[1],0.8,0.4,log=TRUE)+stats::dexp(x[2],log=TRUE)}
#' x0=c(0.5,2)
#' sdjump=c(0.5,1)
#' mcmc=Metropolis_OAAT(f,x0,1000,sdjump)
#' graphicalpar=par(mfrow=c(1,3))
#' plot(mcmc$x);hist(mcmc$x[,1]); hist(mcmc$x[,2])
#' par(graphicalpar)
#' @export
Metropolis_OAAT<-function(f,x0,nsim,sdjump,...){
  #^******************************************************************************
  #^* OBJET: nsim itérations pour l'algorithme MCMC Metropolis_OAAT (One At A Time)
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 09/07/2015
  #^******************************************************************************
  #^* IN
  #^*    1. [function] f, log-densité de la distribution cible
  #^*    2. [real vector] x0, point de départ
  #^*    3. [integer] nsim, nombre de simulations
  #^*    4. [real vector] sdjump, écart-type des sauts pour chaque composante
  #^*    5. [whatever] ... extra arguments passed to f
  #^* OUT
  #^*    1. [list] Une liste comprenant: 
  #^*         $x [real matrix]: simulations MCMC
  #^*         $fx [real vector]: valeur de f correspondant à chaque simu MCMC
  #^*         $moverate [real vector]: taux d'acceptation pour chaque composante
  #^******************************************************************************
  #^* REF.: Renard, B., V. Garreta, and M. Lang (2006), An application of Bayesian 
  #^*       analysis and MCMC methods to the estimation of a regional trend in 
  #^*       annual maxima, Water Resources Research, 42(12).
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************
  n=length(x0)
  x=matrix(NA,nsim,n)
  fx=matrix(NA,nsim,1)
  moverate=matrix(0,1,n)
  z=x0
  fz=f(x0,...)
  # abandon if starting point is unfeasible
  if(is.na(fz) | fz==-Inf){return(list(x=x,fx=fx,moverate=moverate/nsim))}
  for(i in 1:nsim){    
    w<-Metropolis_OAAT_jump(f,z,fz,sdjump,...)
    x[i,]=w$x
    fx[i]=w$fx
    moverate=moverate+w$move
    z=w$x
    fz=w$fx
  }
  return(list(x=x,fx=fx,moverate=moverate/nsim))
}

#' One-At-A-Time Metropolis sampler
#'
#' Performs a single iteration of the OAAT Metropolis sampler
#' (simulated vector is updated one component at a time).
#'  a.k.a block Metropolis sampler with blocks of length one.
#'  Sometimes also called 'Metropolis-within-Gibbs'.
#'
#' @param f function, log-pdf of the target distribution
#' @param x0 numeric vector, starting point
#' @param fx0 numeric, f(x0)
#' @param sdjump numeric vector, standard deviation of the Gaussian jump for each component
#' @param ... other arguments passed to f
#' @return A list with the following components:
#'     \item{x}{numeric vector, updated point after the iteration}
#'     \item{fx}{numeric, updated value f(x)}
#'     \item{move}{logical vector, TRUE for components of the vector x that changed}
#' @examples
#' # Bivariate target distribution: beta(2,10) X exp(1)
#' f=function(x){stats::dbeta(x[1],2,10,log=TRUE)+stats::dexp(x[2],log=TRUE)}
#' x0=c(0.5,0.5)
#' fx0=f(x0)
#' sdjump=c(0.1,0.1)
#' Metropolis_OAAT_jump(f,x0,fx0,sdjump)
#' @export
Metropolis_OAAT_jump<-function(f,x0,fx0,sdjump,...){
  #^******************************************************************************
  #^* OBJET: Une itération pour l'algorithme MCMC Metropolis_OAAT (One At A Time)
  #^*        
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 09/07/2015
  #^******************************************************************************
  #^* IN
  #^*    1. [function] f, log-densité de la distribution cible
  #^*    2. [real vector] x0, point de départ
  #^*    3. [real] fx0, f(x0) 
  #^*    4. [real vector] sdjump, écart-type des sauts pour chaque composante
  #^*    5. [whatever] ... extra arguments passed to f
  #^* OUT
  #^*    1. [list] Une liste comprenant: 
  #^*         $x [real vector]: nouveau point à l'issue de l'itération
  #^*         $fx [real vector]: nouvelle valeur de f à l'issue de l'itération
  #^*         $move [logical vector]: TRUE pour les composantes ayant bougé
  #^******************************************************************************
  #^* REF.: Renard, B., V. Garreta, and M. Lang (2006), An application of Bayesian 
  #^*       analysis and MCMC methods to the estimation of a regional trend in 
  #^*       annual maxima, Water Resources Research, 42(12).
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************
  
  #initialize
  fx=fx0
  x=x0
  n=length(x0)
  move=matrix(FALSE,1,n)
  for (i in 1:n){ # loop on each component of x
    # generate candidate
    candid=x
    candid[i]=candid[i]+sdjump[i]*stats::rnorm(1)
    # evaluate f(candid)
    fcandid=f(candid,...)
    # if NA or -Inf, reject candid
    if(is.na(fcandid) | fcandid==-Inf){next}
    # Otherwise apply Metropolis rule
    ratio=exp(fcandid-fx);
    u=runif(1) # throw the dice
    if(u<=ratio){ # accept
      x=candid
      fx=fcandid
      move[i]=TRUE
    }    
  }
  return(list(x=x,fx=fx,move=move))
}

GetMode<-function(mcmc){
  #^******************************************************************************
  #^* OBJET: get modal estimate from mcmc simulation
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 09/07/2015
  #^******************************************************************************
  #^* IN
  #^*    1. [object] mcmc
  #^* OUT
  #^*    1. [object] object "estimate" (cf. Estimation_SharedTools)
  #^******************************************************************************
  #^* REF.:
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************
  w=Estimate_success
  i=which.max(mcmc$fx)
  w$par=mcmc$x[i,]
  w$obj=mcmc$fx[i]
  return(w)
}

GetUncertainty<-function(mcmc,burn=0.5,slim=5){
#^******************************************************************************
#^* OBJET: get uncertainty estimate from mcmc simulation
#^******************************************************************************
#^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
#^******************************************************************************
#^* CREE/MODIFIE: 09/07/2015
#^******************************************************************************
#^* IN
#^*    1. [object] mcmc
#^* OUT
#^*    1. [object] object "uncertainty" (cf. Estimation_SharedTools)
#^******************************************************************************
#^* REF.:
#^******************************************************************************
#^* A FAIRE: 
#^******************************************************************************
#^* COMMENTAIRES: 
#^******************************************************************************
  out<-Uncertainty_success
  nsim=length(mcmc$fx)
  ix=seq(trunc(burn*nsim),nsim,slim)
  out$cov=stats::cov(as.matrix(mcmc$x[ix,]),use="na.or.complete")
  out$sim=as.matrix(mcmc$x[ix,])
  return(out)
}