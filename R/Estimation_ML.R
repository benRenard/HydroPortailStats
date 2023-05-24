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
#~* OBJET: Estimation par maximum de vraisemblance
#~******************************************************************************
#~* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
#~******************************************************************************
#~* CREE/MODIFIE: XXX
#~******************************************************************************
#~* PRINCIPALES FONCTIONS
#~*    1. GetEstimate_ML, estimateur du maximum de vraisemblance
#~*    2. GetUncertainty_ML, quantification de l'incertitude basée sur le hessien
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

#' Maximum-likelihood estimate of a distribution
#'
#' Returns an estimate of a distribution using the method of maximum likelihood.
#'
#' @param y numeric vector, data
#' @param dist character, distribution name
#' @param par0 numeric vector, initial parameter guess. You may use GetEstimate_ROUGH().
#' @param method character, method used to maximize likelihood, see ?optim
#' @param lower numeric vector, lower bounds, see ?optim
#' @param upper numeric vector, upper bounds, see ?optim
#' @return A list with the following components:
#'     \item{par}{numeric vector, estimated parameter vector.}
#'     \item{obj}{numeric, objective fonction (maximum log-likelihood)}
#'     \item{ok}{logical, did computation succeed?}
#'     \item{err}{integer, error code (0 if ok)}
#'     \item{message}{error message}
#' @examples
#' y=c(9.2,9.5,11.4,9.5,9.4,9.6,10.5,11.1,10.5,10.4)
#' GetEstimate_ML(y,'Normal')
#' GetEstimate_ML(y,'LogNormal')
#' GetEstimate_ML(y,'Gumbel')
#' GetEstimate_ML(y,'Gumbel',par0=GetEstimate_ROUGH(y,'Gumbel')$par)
#' GetEstimate_ML(y,'GEV',par0=GetEstimate_ROUGH(y,'GEV')$par)
#' GetEstimate_ML(y,'Poisson')
#' @export
GetEstimate_ML<-function(y,dist,par0=NULL,method=optim_method_def,
                         lower = -Inf, upper = Inf){
  #^******************************************************************************
  #^* OBJET: Retourne l'estimateur du maximum de vraisemblance
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [real] y, vecteur des données 
  #^*    2. [character] dist, nom de la distribution 
  #^*    3. [real] par0, point de départ pour l'optimiseur 
  #^*    4. [character] method, GetEstimate_ML_optim
  #^*    5-6. [real] lower-upper, vecteur de contraintes, utilisé si method = "L-BFGS-B"
  #^* OUT
  #^*    1. [list] Une liste comprenant: 
  #^*         $par: paramètres du max. de vraisemblance
  #^*         $obj: vraisemblance maximisée
  #^*         $ok: TRUE si ok, FALSE si pb lors de l'optimisation
  #^*         $err: code d'erreur (0 = pas d'erreur)
  #^*         $message: message
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************  
  
  out=switch(dist,
             # Normal distribution
             Normal={ReturnEstimate(rough_Normal(y),y,dist)},
             # LogNormal distribution
             LogNormal={ReturnEstimate(rough_Normal(log(y)),y,dist)},   
             # 1-parameter exponential distribution
             Exponential1={ReturnEstimate(mean(y),y,dist)},   
             # 2-parameter exponential distribution
             Exponential2={ReturnEstimate(rough_Exponential(y),y,dist)},   
             # Poisson
             Poisson={if(mean(y)>0){ReturnEstimate(mean(y),y,dist)} else {ReturnEstimate(NA,y,dist)}},   
             # Autres distributions - pas de formules explicites
{
  if(is.null(par0)) {
    fail=Estimate_fail
    fail$message="[par0=NULL] not allowed for this distribution"
    return(fail)}
  w<-GetEstimate_ML_optim(y,dist,par0,method,lower,upper,do.hessian=FALSE)
  if(w$convergence==0){
    ReturnEstimate(w$par,y,dist)
  } else {
    fail=Estimate_fail
    fail$err=w$convergence
    fail$message=w$message
    return(fail)}
}
  )
return(out)
}

#' Maximum-likelihood estimation of uncertainty
#'
#' Returns an estimate of the uncertainty around the maximum-likelihood estimate, 
#' in the form of a covariance matrix and some simulations from the corresponding 
#' Gaussian distribution.
#'
#' @param y numeric vector, data
#' @param dist character, distribution name
#' @param par numeric vector, estimated parameter (using GetEstimate_ML()).
#' @param nsim integer, number of simulated parameter replicates.
#' @return A list with the following components:
#'     \item{cov}{numeric matrix npar*npar, covariance matrix.}
#'     \item{sim}{numeric matrix nsim*npar, simulated parameter replicates.}
#'     \item{ok}{logical, did computation succeed?}
#'     \item{err}{integer, error code (0 if ok)}
#'     \item{message}{error message}
#' @examples
#' y=c(9.2,9.5,11.4,9.5,9.4,9.6,10.5,11.1,10.5,10.4)
#' estim=GetEstimate_ML(y,'Gumbel',par0=GetEstimate_ROUGH(y,'Gumbel')$par)
#' GetUncertainty_ML(y,'Gumbel',par=estim$par)
#' @importFrom numDeriv hessian
#' @importFrom mvtnorm rmvnorm
#' @export
GetUncertainty_ML<-function(y,dist,par,nsim=nsim_def){
  #^******************************************************************************
  #^* OBJET: Retourne la matrice de covariance des paramètres + des simulations
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [real] y, vecteur des données 
  #^*    2. [character] dist, nom de la distribution 
  #^*    3. [real] par, vecteur de paramètres estimés 
  #^*    4. [integer] nsim, nombre de simulations
  #^* OUT
  #^*    1. [list] Une liste comprenant: 
  #^*         $cov: matrice de covariance (issue du hessien)
  #^*         $sim: paramètres simulés
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************  
  # calcul de la covariance
  hess=tryCatch(numDeriv::hessian(GetLogLikelihood,x=par,y=y,dist=dist),error=function(e) NA)
  if(any(is.na(hess))) {
    fail=Uncertainty_fail
    fail$message="problem computing hessian"
    return(fail)}
  cov=tryCatch(solve(-1*hess),error=function(e) NA)
  if(any(is.na(cov))) {
    fail=Uncertainty_fail
    fail$message="problem computing covariance from hessian"
    return(fail)}
  # simulations
  sim=tryCatch(mvtnorm::rmvnorm(nsim,mean=par,sigma=cov),error=function(e) NA)
  if(any(is.na(sim))) {
    fail=Uncertainty_fail
    fail$message="problem simulating from multivariate normal"
    return(fail)}
  
  out<-Uncertainty_success
  out$cov=cov;out$sim=sim
  return(out)
}

#****************************
# Fonctions privées ----
#****************************

GetLogLikelihood<-function(par,y,dist){
  #^******************************************************************************
  #^* OBJET: Retourne la vraisemblance des données x pour la distribution 'dist' de  
  #^*        paramètres 'par'
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [real] par, vecteur de paramètres 
  #^*    2. [real] y, vecteur des données 
  #^*    3. [character] dist, nom de la distribution 
  #^* OUT
  #^*    1. [real] la log-vraisemblance
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: si les paramètres sont impossibles, retourne NA
  #^*               si la vraisemblance est nulle, retourne -Inf
  #^******************************************************************************  
  p=sapply(y,GetPdf,dist,par,TRUE)
  if(any(is.na(p))) {return(NA)}
  if(any(p==-Inf)) {return(-Inf)}
  return(sum(p))
}

GetEstimate_ML_optim<-function(y,dist,par0,method=optim_method_def,
                               lower = -Inf, upper = Inf,do.hessian=FALSE){
  #^******************************************************************************
  #^* OBJET: Retourne l'estimateur du maximum de vraisemblance par optimisation
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [real] y, vecteur des données 
  #^*    2. [character] dist, nom de la distribution 
  #^*    3. [real] par0, point de départ pour l'optimiseur 
  #^*    4. [character] method, voir ci-dessous et ?optim
  #^*    5-6. [real] lower-upper, vecteur de contraintes, utilisé si method = "L-BFGS-B"
  #^*    7. [logical] do.hessian: calcul du hessien? 
  #^*                 (utile pour calculer la matrice de covariance des paramètres
  #^* OUT
  #^*    1. [list] Une liste comprenant (voir ?optim): 
  #^*         $par: paramètres du max. de vraisemblance
  #^*         $value: vraisemblance maximisée
  #^*         $counts: nombre d'évaluations de la vraisemblance et de son gradient
  #^*         $convergence: 0=ok, sinon: voir ?optim
  #^*         $message: message renvoyé par la fonction optim()
  #^*         $hessian: la matrice hessienne (si demandée)
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: method = "Nelder-Mead": compromis vitesse - fiabilité.
  #^*               method = "BFGS": méthode de Newton, rapide mais peu fiable si 
  #^*                         la vraisemblance n'est pas lisse.
  #^*               method = "L-BFGS-B": comme ci-dessus, mais possibilité d'ajouter
  #^*                         des contraintes.
  #^*               method = "SANN": "Simulated annealing", fonctionne si la
  #^*                         vraisemblance n'est pas lisse, mais lente.
  #^******************************************************************************  
  if(length(par0)==1) {method="BFGS"} # "Nelder-Mead" non fiable en dimension 1 - voir ?optim  
  w<-tryCatch(stats::optim(par0, GetLogLikelihood, gr = NULL, y,dist,
           method = method,
           lower = lower, upper = upper,
           control = list(fnscale=-1), hessian = do.hessian),
           error=function(e) list(par=NA,convergence=666,value=NA,counts=0))
           
  return(w)
}

