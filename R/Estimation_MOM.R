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
#~* OBJET: Estimation par la méthode des moments
#~******************************************************************************
#~* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
#~******************************************************************************
#~* CREE/MODIFIE: XXX
#~******************************************************************************
#~* PRINCIPALES FONCTIONS
#~*    1. GetEstimate_MOM, estimateur des moments
#~******************************************************************************
#~* REF.: XXX
#~******************************************************************************
#~* A FAIRE: XXX
#~******************************************************************************
#~* COMMENTAIRES: 1/ Pas à strictement parler l'estimateur des moments pour la GPD
#~*                  à seuil inconnu, puisque le seuil est estimé comme min(y)
#~* 2/ De même pas à proprement parler l'estimateur des moments pour la LogPearsonIII, 
#~*    puisqu'on calcule l'estimateur des moments de log(y)
#~******************************************************************************

#****************************
# Fonctions principales ----
#****************************

#' Moment estimate of a distribution
#'
#' Returns an estimate of a distribution using the method of moments.
#' Note that for some distributions, this is not strictly speaking the moment estimate.
#' For LogPearsonIII for instance, the moment estimate of log(data) is used.
#' Also for GPD3, the threshold is estimated as min(data).
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
#' GetEstimate_MOM(y,'Normal')
#' GetEstimate_MOM(y,'LogNormal')
#' GetEstimate_MOM(y,'Gumbel')
#' GetEstimate_MOM(y,'GEV')
#' GetEstimate_MOM(y,'Poisson')
#' @export
GetEstimate_MOM<-function(y,dist){
  #^******************************************************************************
  #^* OBJET: Retourne l'estimateur des moments
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [real] y, vecteur des données 
  #^*    2. [character] dist, nom de la distribution 
  #^* OUT
  #^*    1. [list] Une liste comprenant (voir ?optim): 
  #^*         $par: paramètres estimés
  #^*         $obj: NA
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
             Normal={ReturnEstimate(rough_Normal(y))},
             # LogNormal distribution
             LogNormal={ReturnEstimate(MOM_LogNormal(y))},   
             # 1-parameter exponential distribution
             Exponential1={ReturnEstimate(mean(y))},   
             # 2-parameter exponential distribution
             Exponential2={ReturnEstimate(MOM_Exponential2(y))},   
             # 2-parameter GPD
             GPD2={ReturnEstimate(MOM_GPD2(y))},   
             # 3-parameter GPD
             GPD3={ReturnEstimate(MOM_GPD3(y))},   
             # Gumbel
             Gumbel={ReturnEstimate(rough_Gumbel(y))},   
             # GEV
             GEV={ReturnEstimate(MOM_GEV(y))},   
             # PearsonIII
             PearsonIII={ReturnEstimate(MOM_PearsonIII(y))},   
             # LogPearsonIII
             LogPearsonIII={ReturnEstimate(MOM_PearsonIII(log(y)))},
             # Poisson
             Poisson={if(mean(y)>0){ReturnEstimate(mean(y))} else {ReturnEstimate(NA)}}, 
             # Gumbel for minima
             Gumbel_min={ReturnEstimate(rough_Gumbel_min(y))},   
             # GEV for minima
             GEV_min={ReturnEstimate(MOM_GEV_min(y))},   
             # GEV for minima with lower bound at zero
             GEV_min_pos={ReturnEstimate(c(NA,NA))} 
  )
return(out)
}

#****************************
# Fonctions privées ----
#****************************

MOM_Exponential2<-function(y){
  n=length(y)
  m=mean(y)
  s=stats::sd(y)*sqrt((n-1)/n)
  par2=s
  if(par2<=0){return(c(NA,NA))}
  par1=m-s
  par=c(par1,par2)
  return(par)
}

MOM_LogNormal<-function(y){
  n=length(y)
  m=mean(y)
  s=stats::sd(y)*sqrt((n-1)/n)
  x=log(1+s^2/m^2)
  if(x<=0){return(c(NA,NA))}
  par1=log(m)-0.5*x
  par2=sqrt(x)
  par=c(par1,par2)
  return(par)
}

MOM_GPD2<-function(y){
  n=length(y)
  m=mean(y)
  s=stats::sd(y)*sqrt((n-1)/n)
  if(s<=0){return(c(NA,NA))}
  x=m^2/s^2
  par1=0.5*m*(x+1)
  if(par1<=0){return(c(NA,NA))}
  par2=0.5*(x-1)
  par=c(par1,par2)
  return(par)
}

MOM_GPD3<-function(y){
  n=length(y)
  m=mean(y)
  s=stats::sd(y)*sqrt((n-1)/n)
  if(s<=0){return(c(NA,NA,NA))}
  y0=min(y)
  x=(m-y0)^2/s^2
  par1=0.5*(m-y0)*(x+1)
  if(par1<=0){return(c(NA,NA,NA))}
  par2=0.5*(x-1)
  par=c(y0,par1,par2)
  return(par)
}

MOM_GEV<-function(y,lower=MOM_lower,upper=MOM_upper,epsilon=MOM_epsilon){
  n=length(y)
  m=mean(y)
  s=stats::sd(y)*sqrt((n-1)/n)
  m3=(1/n)*sum((y-m)^3);skew=m3/(s^3)
  w=tryCatch(stats::uniroot(fgev,c(lower,upper),skew,epsilon)$root,error=function(e) NA)
  if(is.na(w)) {return(c(NA,NA,NA))}
  par3=w
  if(abs(par3)<epsilon){#Gumbel by continuity
    return(c(rough_Gumbel(y),0))
  } else{
    g1=gamma(par3+1)
    g2=gamma(2*par3+1)
    par2=abs(par3)*s*(g2-g1^2)^(-0.5)
    if(par2<=0){return(c(NA,NA,NA))}
    par1=m-(par2/par3)*(1-g1)    
    par=c(par1,par2,par3)
    return(par)
  }
}

MOM_GEV_min<-function(y){
  par=MOM_GEV(-1*y)
  par[1]=-1*par[1]
  return(par)
}

MOM_PearsonIII<-function(y){
  n=length(y)
  m=mean(y)
  s=stats::sd(y)*sqrt((n-1)/n)
  m3=(1/n)*sum((y-m)^3);skew=m3/(s^3)
  if(skew==0){return(c(NA,NA,NA))}
  par3=4/(skew^2)
  par2=sign(skew)*s/sqrt(par3)
  if(par2<=0){return(c(NA,NA,NA))}
  par1=m-par2*par3
  par=c(par1,par2,par3)
  return(par)  
}

fgev<-function(z,skew,epsilon){
  # fonction à annuler pour obtenir l'estimateur du paramètre de forma de la GEV
  if(abs(z)<=epsilon) {return(skew+skew_gumbel)} # loi de Gumbel par continuité
  g1=gamma(z+1)
  g2=gamma(2*z+1)
  g3=gamma(3*z+1)
  w=skew+(z/abs(z))*(g3-3*g1*g2+2*g1^3)/((g2-g1^2)^1.5)
  return(w)
}
