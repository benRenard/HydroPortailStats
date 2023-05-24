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
#~* OBJET: Estimation par la méthode des L-moments
#~******************************************************************************
#~* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
#~******************************************************************************
#~* CREE/MODIFIE: XXX
#~******************************************************************************
#~* PRINCIPALES FONCTIONS
#~*    1. GetEstimate_LMOM, estimateur des moments
#~******************************************************************************
#~* REF.: Hosking and Wallis, Regional Frequency Analysis: an approach based on
#~*       L-moments, Cambridge University Press, 1997
#~*       see also http://tig.dvo.ru/_ld/0/25_WatRes4_10Gubar.pdf
#~******************************************************************************
#~* A FAIRE: 
#~******************************************************************************
#~* COMMENTAIRES: 1/ Pas à strictement parler l'estimateur des Lmoments pour les
#~*    distributions LogNormal et LogPearsonIII,puisqu'on calcule l'estimateur 
#~*    des L-moments de log(y)
#~******************************************************************************

#****************************
# Fonctions principales ----
#****************************

#' L-Moment estimate of a distribution
#'
#' Returns an estimate of a distribution using the method of L-moments.
#' Note that for some distributions, this is not strictly speaking the L-moment estimate:
#' For LogNormal and LogPearsonIII, the L-moment estimate of log(data) is used.
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
#' GetEstimate_LMOM(y,'Normal')
#' GetEstimate_LMOM(y,'LogNormal')
#' GetEstimate_LMOM(y,'Gumbel')
#' GetEstimate_LMOM(y,'GEV')
#' GetEstimate_LMOM(y,'Poisson')
#' @export
GetEstimate_LMOM<-function(y,dist){
  #^******************************************************************************
  #^* OBJET: Retourne l'estimateur des L-moments
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
  #^*         $obj: NA
  #^*         $ok: TRUE si ok, FALSE si pb
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
             Normal={ReturnEstimate(LMOM_Normal(y))},
             # LogNormal distribution
             LogNormal={ReturnEstimate(LMOM_Normal(log(y)))},   
             # 1-parameter exponential distribution
             Exponential1={ReturnEstimate(LMOM_Exponential1(y))},   
             # 2-parameter exponential distribution
             Exponential2={ReturnEstimate(LMOM_Exponential2(y))},   
             # 2-parameter GPD
             GPD2={ReturnEstimate(LMOM_GPD2(y))},   
             # 3-parameter GPD
             GPD3={ReturnEstimate(LMOM_GPD3(y))},   
             # Gumbel
             Gumbel={ReturnEstimate(LMOM_Gumbel(y))},   
             # GEV
             GEV={ReturnEstimate(LMOM_GEV(y))},   
             # PearsonIII
             PearsonIII={ReturnEstimate(LMOM_PearsonIII(y))},   
             # LogPearsonIII
             LogPearsonIII={ReturnEstimate(LMOM_PearsonIII(log(y)))},
             # Poisson
             Poisson={ReturnEstimate(LMOM_Poisson(y))}, 
             # Gumbel for minima
             Gumbel_min={ReturnEstimate(LMOM_Gumbel_min(y))},   
             # GEV for minima
             GEV_min={ReturnEstimate(LMOM_GEV_min(y))},   
             # GEV for minima with lower bound at zero
             GEV_min_pos={ReturnEstimate(c(NA,NA))}   
  )
  return(out)
}

#****************************
# Fonctions privées ----
#****************************

pstar<-function(r,k){
  if(r==k){w=choose(r,k)*choose(r+k,k)
  } else {
    w=((-1)^(r-k))*choose(r,k)*choose(r+k,k)
  }
  return(w)
}

Get3LMoments<-function(y){
  w=rep(NA,3)
  if(any(is.na(y))){return(c(NA,NA,NA))}
  sy=sort(y);n=length(y)
  b0=sum(sy)/n
  w[1]=b0
  b1=sum( (((1:n)-1)/(n-1))*sy )/n
  w[2]=b0*pstar(1,0)+b1*pstar(1,1)
  if(is.na(w[2])){return(w)}
  if(w[2]==0){return(w)}
  b2=sum( ( (((1:n)-1)*((1:n)-2))/((n-1)*(n-2)) )*sy )/n
  w[3]=(b0*pstar(2,0)+b1*pstar(2,1)+b2*pstar(2,2))/w[2]
  return(w)
}

LMOM_Normal<-function(y){
  w=Get3LMoments(y)
  if(any(is.na(w))){return(c(NA,NA))}
  par1=w[1]
  par2=w[2]*sqrt(pi)
  if(par2<=0){return(c(NA,NA))}
  par=c(par1,par2)
  return(par)
}

LMOM_Gumbel<-function(y){
  w=Get3LMoments(y)
  if(any(is.na(w))){return(c(NA,NA))}
  par2=w[2]/log(2)
  if(par2<=0){return(c(NA,NA))}
  par1=w[1]-gamma_gumbel*par2
  par=c(par1,par2)
  return(par)
}

LMOM_Gumbel_min<-function(y){
  par=LMOM_Gumbel(-1*y)
  par[1]=-1*par[1]
  return(par)
}

LMOM_GEV<-function(y){
  w=Get3LMoments(y)
  if(any(is.na(w))){return(c(NA,NA,NA))}
  c=2/(3+w[3])-log(2)/log(3)
  par3=7.8590*c+2.9554*c^2
  par2=(w[2]*par3)/((1-2^(-1*par3))*exp(lgamma(1+par3)))
  if(par2<=0){return(c(NA,NA,NA))}
  par1=w[1]-(par2/par3)*(1-exp(lgamma(1+par3)))
  par=c(par1,par2,par3)
  return(par)
}

LMOM_GEV_min<-function(y){
  par=LMOM_GEV(-1*y)
  par[1]=-1*par[1]
  return(par)
}

LMOM_PearsonIII<-function(y){
  w=Get3LMoments(y)
  if(any(is.na(w))){return(c(NA,NA,NA))}
  if(abs(w[3])<1/3){
    z=3*pi*w[3]^2
    a=(1+0.2906*z)/(z+0.1882*z^2+0.0442*z^3)
  } else {
    z=1-abs(w[3])
    a=(0.36067*z-0.59567*z^2+0.25361*z^3)/(1-2.78861*z+2.56096*z^2-0.77045*z^3)
  }
  if(a<=0){return(c(NA,NA,NA))}
  gama=2*sign(w[3])/sqrt(a)
  sigma=w[2]*sqrt(pi*a)*exp(lgamma(a)-lgamma(0.5+a))
  mu=w[1]
  par1=mu-2*(sigma/gama)
  par2=0.5*sigma*gama
  par3=4/gama^2
  if( abs(par3)<=1e-05 | abs(par3)>=1e+05 ){return(c(NA,NA,NA))} # avoid numerical issues
  par=c(par1,par2,par3)
  return(par)  
}

LMOM_Exponential1<-function(y){
  w=Get3LMoments(y)
  if(any(is.na(w))){return(c(NA))}
  par=w[1]
  if(par<=0){return(c(NA))}
  return(par)
}

LMOM_Exponential2<-function(y){
  w=Get3LMoments(y)
  if(any(is.na(w))){return(c(NA,NA))}
  par2=2*w[2]
  if(par2<=0){return(c(NA,NA))}
  par1=w[1]-par2
  par=c(par1,par2)
  return(par)
}

LMOM_GPD2<-function(y){
  w=Get3LMoments(y)
  if(any(is.na(w))){return(c(NA,NA))}
  if(w[2]<=0){return(c(NA,NA))}
  par2=w[1]/w[2]-2
  if(par2<=0){return(c(NA,NA))}
  par1=(1+par2)*w[1]
  par=c(par1,par2)
  return(par)
}

LMOM_GPD3<-function(y){
  w=Get3LMoments(y)
  if(any(is.na(w))){return(c(NA,NA,NA))}
  if(w[3]==-1){return(c(NA,NA,NA))}
  par3=(1-3*w[3])/(1+w[3])       
  if(abs(par3)>=1e+06){return(c(NA,NA,NA))} # numerical issues
  par2=(1+par3)*(2+par3)*w[2]
  if(par2<=0){return(c(NA,NA,NA))}
  par1=w[1]-(2+par3)*w[2]
  par=c(par1,par2,par3)
  return(par)
}

LMOM_Poisson<-function(y){
  par=mean(y)
  if(par<=0){return(NA)}
  return(par)
}

