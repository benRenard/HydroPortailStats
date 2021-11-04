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
#~* OBJET: Tests statistiques pour la Banque HYDRO3
#~******************************************************************************
#~* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
#~******************************************************************************
#~* CREE/MODIFIE: 09/08/2016
#~******************************************************************************
#~* PRINCIPALES FONCTIONS
#~*    1. KS: test d'adéquation de Kolmogorov-Smirnov
#~*    2. MK: test de tendance de Mann-Kendall
#~*    2. Pettitt: test de rupture de Pettitt
#~******************************************************************************
#~* REF.: XXX
#~******************************************************************************
#~* A FAIRE:
#~******************************************************************************
#~* COMMENTAIRES: Attention aux nombreuses hypothèses derrière ces tests...
#~******************************************************************************

GetCdf_vect<-function(y,dist,par){
  #^******************************************************************************
  #^* OBJET: a vectorial version of the GetCdf function, needed by the ks.test function
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 15/03/2017
  #^******************************************************************************
  #^* IN
  #^*    1. [real vector] y, données 
  #^*    2. [character] dist, distribution 
  #^*    3. [real vector] par, paramètres
  #^* OUT
  #^*    1. vecteur contenant la cdf évaluée en chaque point de y  
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************
  return(sapply(y,GetCdf,dist,par))
}

#' Kolmogorov-Smirnov Test
#'
#' Applies a one-sample Kolmogorov-Smirnov test (see ?stats::ks.test)
#'
#' @param y numeric vector, data
#' @param dist character, distribution name
#' @param par numeric vector, parameter vector
#' @return A list with the following components:
#'     \item{pval}{numeric, p-value of the test}
#'     \item{stat}{numeric, test statistics}
#'     \item{xtra}{numeric, xtra information: empty for this test}
#' @examples
#' y=stats::rnorm(20)
#' KS(y,'Normal',c(0,1))
#' KS(y,'Normal',c(1,1))
#' KS(y,'Gumbel',c(0,1))
#' @importFrom stats ks.test
#' @export
KS<-function(y,dist,par){
  #^******************************************************************************
  #^* OBJET: test d'adéquation de Kolmogorov-Smirnov
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 09/08/2016
  #^******************************************************************************
  #^* IN
  #^*    1. [real vector] y, données 
  #^*    2. [character] dist, distribution supposée 
  #^*    3. [real vector] par, paramètres estimés
  #^* OUT
  #^*    1. liste contenant les éléments suivants: 
  #^*       $pval, p-valeur du test  
  #^*       $stat, statistique du test  
  #^*       $xtra, rien pour ce test  
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: Application abusive dans Hydro3, car il n'y a pas de correction
  #^*               pour prendre en compte le fait que les paramètres sont estimés.
  #^*               De plus hypothèse d'indépendance surement douteuse pour 
  #^*               certaines séries.
  #^******************************************************************************
  
  # initialisation 
  OUT=list(pval=NA,stat=NA,xtra=NA)
  # use R's ks.test function
  w=stats::ks.test(y, "GetCdf_vect",dist,par)
  # return
  OUT$pval=w$p.value
  OUT$stat=as.numeric(w$statistic)
  return(OUT)
}

#' Mann-Kemdall Test
#'
#' Applies the Mann-Kendall trend test
#'
#' @param y numeric vector, data
#' @return A list with the following components:
#'     \item{pval}{numeric, p-value of the test}
#'     \item{stat}{numeric, test statistics}
#'     \item{xtra}{numeric, xtra information: empty for this test}
#' @examples
#' y=stats::rnorm(50)
#' MK(y)
#' y=y+0.1*(1:length(y))
#' MK(y)
#' @export
MK<-function(y){
  #^******************************************************************************
  #^* OBJET: test de tendance de Mann-Kendall
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 09/08/2016
  #^******************************************************************************
  #^* IN
  #^*    1. [real vector] y, données 
  #^* OUT
  #^*    1. liste contenant les éléments suivants: 
  #^*       $pval, p-valeur du test  
  #^*       $stat, statistique du test  
  #^*       $xtra, rien pour ce test  
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: Hypothèse d'indépendance surement douteuse pour 
  #^*               certaines séries.
  #^******************************************************************************  
  
  # initialisation 
  OUT=list(pval=NA,stat=NA,xtra=NA)
  
  # Calcul de la statistique de test basique 
  n=length(y);count.p=0;count.m=0;
  for(j in 2:n){
    for(i in 1:(j-1)){
      if(y[j]>y[i]){count.p=count.p+1}
      else if (y[j]<y[i]){count.m=count.m+1}
    }
  }
  mk=count.p-count.m
  
  # Calcul de la variance de mk
  var0=((n*(n-1)*(2*n+5)))/18
  # Correction pour les ex-aequos
  w=matrix(NA,n,1);tie=matrix(NA,n,1);v=matrix(NA,n,1)
  for(i in 1:n){w[i]=sum(y==y[i])} # compte combien de fois chaque valeur est dupliquée
  for(i in 1:n){
    tie[i]=sum(w==i)/i # vecteur contenant le nombre d'ex-aequo d'étendue i
    v[i]=tie[i]*i*(i-1)*(2*i+5) # contribution de l'étendue i à la correction
  }
  tie.correction=sum(v)/18
  MKvar=var0-tie.correction
  if(is.na(MKvar)){return(OUT)}
  if(MKvar<=0){return(OUT)}

  # statistique finale et p-valeur
  if(mk>0){stat=(mk-1)/sqrt(MKvar)}
  else if(mk<0){stat=(mk+1)/sqrt(MKvar)}
  else{stat=mk/sqrt(MKvar)}
  OUT$stat=stat
  # p-val (2-sided test)
  OUT$pval=2*stats::pnorm(-1*abs(stat),mean=0,sd=1)
  
  return(OUT)
}

#' Pettitt Test
#'
#' Applies the Pettitt step-change test
#'
#' @param y numeric vector, data
#' @return A list with the following components:
#'     \item{pval}{numeric, p-value of the test}
#'     \item{stat}{numeric, test statistics}
#'     \item{xtra}{numeric, xtra information: position of the step change}
#' @examples
#' y=stats::rnorm(50)
#' Pettitt(y)
#' y[26:50]=y[26:50]+2
#' Pettitt(y)
#' @export
Pettitt<-function(y){
  #^******************************************************************************
  #^* OBJET: test de rupture de Pettitt
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: 09/08/2016
  #^******************************************************************************
  #^* IN
  #^*    1. [real vector] y, données 
  #^* OUT
  #^*    1. liste contenant les éléments suivants: 
  #^*       $pval, p-valeur du test  
  #^*       $stat, statistique du test  
  #^*       $xtra, rien pour ce test  
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: Hypothèse d'indépendance surement douteuse pour 
  #^*               certaines séries + pas de corrections pour les ex-aequo.
  #^******************************************************************************  
  
  # initialisation 
  OUT=list(pval=NA,stat=NA,xtra=NA)
  
  # Calcul de la statistique de test 
  n=length(y);U=0*(1:(n-1))
  for(k in 1:(n-1)){
    count.p=0;count.m=0;
    for(i in 1:k){
      for(j in (k+1):n){
        if(y[j]>y[i]){count.p=count.p+1} else if (y[j]<y[i]){count.m=count.m+1}
      }
    }
    U[k]=count.p-count.m
  }
  # return
  OUT$stat=max(abs(U))
  OUT$xtra=which.max(abs(U))
  OUT$pval=min(c(1,2*exp((-6*OUT$stat^2)/(n^2+n^3))))
  return(OUT)
}