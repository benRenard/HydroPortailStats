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
#~* OBJET: Outils utiles pour tous les estimateurs
#~******************************************************************************
#~* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
#~******************************************************************************
#~* CREE/MODIFIE: XXX
#~******************************************************************************
#~* PRINCIPALES FONCTIONS
#~*    1. XXX
#~******************************************************************************
#~* REF.: XXX
#~******************************************************************************
#~* A FAIRE: XXX
#~******************************************************************************
#~* COMMENTAIRES: XXX
#~******************************************************************************

#****************************
# Déclaration de variables globales
#****************************

# objets "estimate" et "uncertainty" par défaut
Estimate_fail=list(par=NA,obj=NA,ok=FALSE,err=666,message="fatal")
Estimate_success=list(par=NA,obj=NA,ok=TRUE,err=0,message="ok")
Uncertainty_fail=list(cov=NA,sim=NA,ok=FALSE,err=666,message="fatal")
Uncertainty_success=list(cov=NA,sim=NA,ok=TRUE,err=0,message="ok")

# Constantes utilisée pour la loi de Gumbel
gamma_gumbel=0.57721566490153
skew_gumbel=-1.13954709940465

# paramètres par défaut de l'optimiseur
optim_method_def="Nelder-Mead"

# paramètres par défaut pour les simulations Monte Carlo
nsim_def=1000

# paramètres par défaut de la méthode des moments
MOM_lower=-1/3+0.0001
MOM_upper=10
MOM_epsilon=0.0001

#****************************
# Funk
#****************************

ReturnEstimate<-function(par,y=NULL,dist=NULL){
  # return estimate
  out<-Estimate_success
  if(any(is.na(par))){return(Estimate_fail)}
  out$par=par
  if(!is.null(y)) {out$obj=GetLogLikelihood(par,y,dist)}
  return(out)
}

GetUncertainty_Bootstrap<-function(y,dist,par,estimator,type="P",
                                   nsim=nsim_def,method=optim_method_def,
                                   lower = -Inf, upper = Inf){
  n=length(y);npar=length(par)
  sim=matrix(NA,nsim,npar)
  if(type=="P") { # Parametric Bootstrap
    z=Generate(dist,par,n*nsim)
    x=matrix(z,nsim,n)
  } else { # "standard" bootstrap
    x=matrix(NA,nsim,n)
    for(i in 1:nsim){x[i,]=sample(y,length(y),replace=TRUE)}
  }
  v=apply(x,1,GetEstimate_OneSample_vector,dist,estimator,method,lower,upper)
  if(npar>1) {sim=t(v)} else {sim=as.matrix(v)}
  out<-Uncertainty_success
  out$cov=stats::cov(sim,use="na.or.complete")
  out$sim=sim
  if(any(is.na(sim))){
    out$err=sum(is.na(sim[,1]))
    message="warning: certains parametres sont NA"
  }
  return(out)
}

GetEstimate_OneSample<-function(y,dist,estimator,
                                method=optim_method_def,
                                lower = -Inf, upper = Inf){
  w=switch(estimator,
           MOM=GetEstimate_MOM(y,dist),
           ROUGH=GetEstimate_ROUGH(y,dist),
           ML={
             w=GetEstimate_ROUGH(y,dist)
             if(w$ok) {
               par0=w$par
               w=GetEstimate_ML(y,dist,par0,method,lower,upper)
             }               
           },
           LMOM=GetEstimate_LMOM(y,dist),
           BAY=Estimate_fail,
           Estimate_fail)
  return(w)
}

GetEstimate_OneSample_vector<-function(y,dist,estimator,
                                       method=optim_method_def,
                                       lower = -Inf, upper = Inf){
  w=GetEstimate_OneSample(y,dist,estimator,method,lower,upper)
  if(w$ok) {return(w$par)} else {return(rep(NA,GetParNumber(dist)))}
}

PlotPdf<-function(y,dist,par,x=y,lang='fr',lwd=3,col="black",cex.lab=1.5,cex.axis=1.25){
  #^******************************************************************************
  #^* OBJET: trace la densité de probabilité de la distribution 'dist' de  
  #^*        paramètres 'par', ainsi que les données 'y'  
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [real] y, valeur où la densité est évaluée 
  #^*    2. [character] dist, nom de la distribution 
  #^*    3. [real] par, vecteur de paramètres 
  #^*    4. [real] x, valeurs où la densité est calculée 
  #^*    5. [character] lang, langue utilisée pour les légendes
  #^* OUT
  #^*    1. [plot]
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES:
  #^******************************************************************************
  x0=sort(x)
  fx=sapply(x0,GetPdf,dist,par)
  ylab=switch(lang,fr="densite",en="density","")
  graphics::plot(x0,fx,type="l",lwd=lwd,col=col,xlab="y",ylab=ylab,cex.lab=cex.lab,cex.axis=cex.axis)
  graphics::rug(y)
}

PlotCdf<-function(y,dist,par,x=y,lang='fr',lwd=3,col="black",cex.lab=1.5,cex.axis=1.25,pch=19){
  #^******************************************************************************
  #^* OBJET: trace la fonction de répartition de la distribution 'dist' de  
  #^*        paramètres 'par', ainsi que les données 'y'  
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [real] y, valeur où la fdr est évaluée 
  #^*    2. [character] dist, nom de la distribution 
  #^*    3. [real] par, vecteur de paramètres 
  #^*    4. [real] x, valeurs où la fdr 
  #^*    5. [character] lang, langue utilisée pour les légendes
  #^* OUT
  #^*    1. [plot]
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES:
  #^******************************************************************************
  x0=sort(x)
  n=length(y)
  pp=((1:n)-0.5)/n
  fx=sapply(x0,GetCdf,dist,par)
  ylab=switch(lang,fr="probabilite au non-depassement",en="cdf","")
  graphics::plot(x0,fx,type="l",lwd=lwd,col=col,xlab="y",ylab=ylab,cex.lab=cex.lab,
       cex.axis=cex.axis,ylim=c(0,1))
  graphics::points(sort(y),pp,pch=pch)
}
