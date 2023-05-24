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
#~* OBJET: Estimation grossière
#~******************************************************************************
#~* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
#~******************************************************************************
#~* CREE/MODIFIE: XXX
#~******************************************************************************
#~* PRINCIPALES FONCTIONS
#~*    1. GetEstimate_ROUGH, estimateur grossier
#~******************************************************************************
#~* REF.: XXX
#~******************************************************************************
#~* A FAIRE: XXX
#~******************************************************************************
#~* COMMENTAIRES: Le but de cet estimateur est simplement de fournir un 
#~*               point de départ pour les estimateurs plus complexes 
#~*               (e.g. ML ou BAY). L'estimateur grossier doit simplement être 
#~*               facile à calculer et être "passe-partout" (éviter les NA etc.)
#~******************************************************************************

#****************************
# Fonctions principales ----
#****************************

#' Rough estimate of a distribution
#'
#' Returns a rough first-guess estimate of a distribution.
#' This estimate may be poor but it solely aims at being used as a starting point 
#' for more advanced estimation approaches (e.g. max-likelihood or Bayesian).
#' It is therefore chosen as an easy-to-compute explicit formula, robust and error-proof.
#'
#' @param y numeric vector, data
#' @param dist character, distribution name
#' @return A list with the following components:
#'     \item{par}{numeric vector, estimated parameter vector.}
#'     \item{obj}{numeric, objective fonction (NA for this estimate)}
#'     \item{ok}{logical, did computation succeed?}
#'     \item{err}{integer, error code (0 if ok)}
#'     \item{message}{error message}
#' @examples
#' y=c(9.2,9.5,11.4,9.5,9.4,9.6,10.5,11.1,10.5,10.4)
#' GetEstimate_ROUGH(y,'Normal')
#' GetEstimate_ROUGH(y,'LogNormal')
#' GetEstimate_ROUGH(y,'Gumbel')
#' GetEstimate_ROUGH(y,'GEV')
#' GetEstimate_ROUGH(y,'Poisson')
#' @export
GetEstimate_ROUGH<-function(y,dist){
  #^******************************************************************************
  #^* OBJET: Retourne l'estimateur grossier pour la distribution "dist"
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [real] y, vecteur des données 
  #^*    2. [character] dist, nom de la distribution 
  #^* OUT
  #^*    1. [list] Une liste comprenant: 
  #^*         $par: paramètres estimés
  #^*         $obj: fonction objectif (NA pour cet estimateur)
  #^*         $ok: TRUE si ok, FALSE si pb lors du calcul
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
             Normal={ReturnEstimate(rough_Normal(y))},
             # LogNormal distribution
             LogNormal={ReturnEstimate(rough_Normal(log(y)))},   
             # 1-parameter exponential distribution
             Exponential1={ReturnEstimate(mean(y))},   
             # 2-parameter exponential distribution
             Exponential2={ReturnEstimate(rough_Exponential(y))},   
             # 2-parameter GPD
             GPD2={ReturnEstimate(c(mean(y),0))},   
             # 3-parameter GPD
             GPD3={ReturnEstimate(c(rough_Exponential(y),0))},   
             # Gumbel
             Gumbel={ReturnEstimate(rough_Gumbel(y))},   
             # GEV
             GEV={ReturnEstimate(c(rough_Gumbel(y),0))},   
             # PearsonIII
             PearsonIII={ReturnEstimate(rough_PearsonIII(y))},   
             # LogPearsonIII
             LogPearsonIII={ReturnEstimate(rough_PearsonIII(log(y)))},   
             # Poisson
             Poisson={if(mean(y)>0){ReturnEstimate(mean(y))} else {ReturnEstimate(NA)}},
             # Gumbel for minima
             Gumbel_min={ReturnEstimate(rough_Gumbel_min(y))},   
             # GEV for minima
             GEV_min={ReturnEstimate(c(rough_Gumbel_min(y),0))},   
             # GEV for minima with lower bound at zero
             GEV_min_pos={ReturnEstimate(rough_Gumbel_min(y))}   
  )
  return(out)
}

#****************************
# Fonctions privées ----
#****************************

rough_Normal<-function(y){
  n=length(y)
  m=mean(y)
  s=stats::sd(y)*sqrt((n-1)/n)
  par=c(m,s)
  return(par)
}

rough_Exponential<-function(y){
  s=min(y)
  m=mean(y)-s
  par=c(s,m)
  return(par)
}

rough_Gumbel<-function(y){
  n=length(y)
  m=mean(y)
  s=stats::sd(y)*sqrt((n-1)/n)
  scale=(sqrt(6)/pi)*s
  loc=m-gamma_gumbel*scale
  par=c(loc,scale)
  return(par)
}

rough_Gumbel_min<-function(y){
  par=rough_Gumbel(-1*y)
  par[1]=-1*par[1]
  return(par)
}

rough_PearsonIII<-function(y){
  n=length(y)
  m=mean(y)
  mini=min(y)
  s=stats::sd(y)*sqrt((n-1)/n)
  m3=(1/n)*sum((y-m)^3);skew=m3/(s^3)
  if(skew==0) {par3=4} else {par3=4/(skew^2)}
  par2=sign(skew)*s/sqrt(par3)
  par1=m-par2*par3
  par=c(par1,par2,par3)
  return(par)
}
