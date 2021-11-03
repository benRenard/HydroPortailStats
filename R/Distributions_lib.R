# --------------------------
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
# --------------------------

#~******************************************************************************
#~* OBJET: Librairie de distributions & outils de base associés
#~******************************************************************************
#~* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
#~******************************************************************************
#~* CREE/MODIFIE: modifié le 26/07/2016: ajout des distributions pour minima
#~******************************************************************************
#~* PRINCIPALES FONCTIONS
#~*    1. GetParNumber, nombre de paramètres
#~*    2. GetParName, nom des paramètres
#~*    3. GetParFeas, faisabilité des paramètres
#~*    4. GetPdf, calcul de la densité
#~*    5. GetCdf, calcul de la fonction de répartition
#~*    6. GetQuantile, calcul du quantile
#~*    7. Generate, simulation
#~*    8. GetReducedVariate, calcul de la variable réduite
#~******************************************************************************
#~* REF.: XXX
#~******************************************************************************
#~* A FAIRE: XXX
#~******************************************************************************
#~* COMMENTAIRES: XXX
#~******************************************************************************

GetParNumber<-function(dist){
  #^******************************************************************************
  #^* OBJET: Retourne le nombre de paramètres de la distribution dist 
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [character] dist, nom de la distribution 
  #^* OUT
  #^*    1. [integer] nombre de paramètres de dist, NA si dist n'est pas reconnue 
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************
  Npar=switch(dist,
              FlatPrior=0,
              Uniform=2,
              Normal=2,
              LogNormal=2,
              Gumbel=2,
              Exponential1=1, 
              Exponential2=2, 
              GEV=3,
              GPD2=2, 
              GPD3=3, 
              Poisson=1,
              PearsonIII=3,
              LogPearsonIII=3,
              Gumbel_min=2,
              GEV_min=3,
              GEV_min_pos=2,
              NA)
  return(Npar)}

GetParName<-function(dist,lang='fr'){
  #^******************************************************************************
  #^* OBJET: Retourne le nom des paramètres de la distribution dist 
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [character] dist, nom de la distribution 
  #^*    2. [character] lang, langue ('fr' ou 'en')
  #^* OUT
  #^*    1. [character] nom des paramètres, NA si dist n'est pas reconnue 
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************
  name=switch(dist,
              FlatPrior="",
              Uniform=switch(lang,fr=c('borne_inf','borne_sup'),en=c('lower-bound','higher_bound'),NA), 
              Normal=switch(lang,fr=c('moyenne','ecart_type'),en=c('mean','standard_deviation'),NA),
              LogNormal=switch(lang,fr=c('moyenne_log','ecart_type_log'),en=c('mean_log','standard_deviation_log'),NA),
              Gumbel=switch(lang,fr=c('position','echelle'),en=c('location','scale'),NA),
              Exponential1=switch(lang,fr=c('echelle'),en=c('scale'),NA), 
              Exponential2=switch(lang,fr=c('seuil','echelle'),en=c('threshold','scale'),NA), 
              GEV=switch(lang,fr=c('position','echelle','forme'),en=c('location','scale','shape'),NA),
              GPD2=switch(lang,fr=c('echelle','forme'),en=c('scale','shape'),NA),
              GPD3=switch(lang,fr=c('seuil','echelle','forme'),en=c('threshold','scale','shape'),NA), 
              Poisson=switch(lang,fr=c('taux'),en=c('rate'),NA),
              PearsonIII=switch(lang,fr=c('position','echelle','forme'),en=c('location','scale','shape'),NA),
              LogPearsonIII=switch(lang,fr=c('position_log','echelle_log','forme_log'),en=c('location_log','scale_log','shape_log'),NA),
              Gumbel_min=switch(lang,fr=c('position','echelle'),en=c('location','scale'),NA),
              GEV_min=switch(lang,fr=c('position','echelle','forme'),en=c('location','scale','shape'),NA),
              GEV_min_pos=switch(lang,fr=c('position','echelle'),en=c('location','scale'),NA),
              NA)
  return(name)}

GetParFeas<-function(dist,par){
  #^******************************************************************************
  #^* OBJET: Retourne la faisabilité du paramètre 'par' de la distribution 'dist'  
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [character] dist, nom de la distribution 
  #^*    2. [real] par, vecteur de paramètres 
  #^* OUT
  #^*    1. [logical] TRUE si 'par' est faisable, FALSE sinon
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************
  
  if(!is.null(par)) {if(any(is.na(par))){return(NA)}}
  
  feas=switch(dist,
              FlatPrior=T,
              Uniform={if(par[2]<=par[1]){F} else {T}},
              Normal={if(par[2]<=0){F} else {T}},
              LogNormal={if(par[2]<=0){F} else {T}},
              Gumbel={if(par[2]<=0){F} else {T}},
              Exponential1={if(par[1]<=0){F} else {T}}, 
              Exponential2={if(par[2]<=0){F} else {T}}, 
              GEV={if(par[2]<=0){F} else {T}},
              GPD2={if(par[1]<=0){F} else {T}}, 
              GPD3={if(par[2]<=0){F} else {T}}, 
              Poisson={if(par[1]<=0){F} else {T}},
              PearsonIII={if( (par[2]==0) | (par[3]<=0) ){F} else {T}},
              LogPearsonIII={if( (par[2]==0) | (par[3]<=0) ){F} else {T}},
              Gumbel_min={if(par[2]<=0){F} else {T}},
              GEV_min={if(par[2]<=0){F} else {T}},
              GEV_min_pos={if( (par[1]<=0) | (par[2]<=0)){F} else {T}},
              NA)
  return(feas)}

GetPdf<-function(y,dist,par,log=F){
  #^******************************************************************************
  #^* OBJET: Retourne la densité de probabilité de la distribution 'dist' de  
  #^*        paramètres 'par', évaluée en 'y'  
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [real] y, valeur où la densité est évaluée 
  #^*    2. [character] dist, nom de la distribution 
  #^*    3. [real] par, vecteur de paramètres 
  #^*    4. [logical] log, si TRUE retourne le log de la densité 
  #^* OUT
  #^*    1. [real] la densité évaluée en y
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: si les paramètres sont impossibles, retourne NA
  #^*               si la densité est nulle et log=TRUE, retourne -Inf
  #^******************************************************************************
  
  if(is.na(y)){return(NA)}
  if(!GetParFeas(dist,par)){return(NA)}
  
  pdf=switch(dist,
             FlatPrior=1,
             Uniform=dunif(y,min=par[1],max=par[2],log=log),
             Normal=dnorm(y,mean=par[1],sd=par[2],log=log),
             LogNormal=dlnorm(y,meanlog=par[1],sdlog=par[2],log=log),
             Gumbel=dgumbel(y,loc=par[1],scale=par[2],log=log),
             Exponential1=dexp(y,rate=1/par[1],log=log), 
             Exponential2=dexp(y-par[1],rate=1/par[2],log=log), 
             GEV=dgev(y,loc=par[1],scale=par[2],shape=-1*par[3],log=log),
             GPD2={if(y==0){ if(log==F) {1/par[1]} else {-log(par[1])}
             } else {dgpd(y,loc=0,scale=par[1],shape=-1*par[2],log=log)}}, 
             GPD3={if(y==par[1]){if(log==F) {1/par[2]} else {-log(par[2])}
             } else {dgpd(y,loc=par[1],scale=par[2],shape=-1*par[3],log=log)}}, 
             Poisson=dpois(y,lambda=par[1],log=log),
             PearsonIII={ if(par[2]>0){dgamma(y-par[1],shape=par[3],scale=par[2],log=log)
             } else {dgamma(par[1]-y,shape=par[3],scale=-1*par[2],log=log)}},
             LogPearsonIII={ if(par[2]>0){fy=dgamma(log(y)-par[1],shape=par[3],scale=par[2],log=T)-log(y)
             } else {fy=dgamma(par[1]-log(y),shape=par[3],scale=-1*par[2],log=T)-log(y)}
             if(log==T) {fy
             } else {exp(fy)}},
             Gumbel_min=dgumbel(-1*y,loc=-1*par[1],scale=par[2],log=log),
             GEV_min=dgev(-1*y,loc=-1*par[1],scale=par[2],shape=-1*par[3],log=log),
             GEV_min_pos=dgev(-1*y,loc=-1*par[1],scale=par[2],shape=-1*par[2]/par[1],log=log),
             NA)
  return(pdf)}

GetCdf<-function(y,dist,par){
  #^******************************************************************************
  #^* OBJET: Retourne la fonction de répartition de la distribution 'dist' de  
  #^*        paramètres 'par', évaluée en 'y'  
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [real] y, valeur où la densité est évaluée 
  #^*    2. [character] dist, nom de la distribution 
  #^*    3. [real] par, vecteur de paramètres 
  #^* OUT
  #^*    1. [real] la fonction de répartition évaluée en y
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: si les paramètres sont impossibles, retourne NA
  #^******************************************************************************
  
  if(is.na(y)){return(NA)}
  if(!GetParFeas(dist,par)){return(NA)}
  
  cdf=switch(dist,
             Uniform=punif(y,min=par[1],max=par[2]),
             Normal=pnorm(y,mean=par[1],sd=par[2]),
             LogNormal=plnorm(y,meanlog=par[1],sdlog=par[2]),
             Gumbel=pgumbel(y,loc=par[1],scale=par[2]),
             Exponential1=pexp(y,rate=1/par[1]), 
             Exponential2=pexp(y-par[1],rate=1/par[2]), 
             GEV=pgev(y,loc=par[1],scale=par[2],shape=-1*par[3]),
             GPD2=pgpd(y,loc=0,scale=par[1],shape=-1*par[2]), 
             GPD3=pgpd(y,loc=par[1],scale=par[2],shape=-1*par[3]), 
             Poisson=ppois(y,lambda=par[1]),
             PearsonIII={ if(par[2]>0){pgamma(y-par[1],shape=par[3],scale=par[2])
             } else {1-pgamma(par[1]-y,shape=par[3],scale=-1*par[2])}},
             LogPearsonIII={ if(par[2]>0){pgamma(log(y)-par[1],shape=par[3],scale=par[2])
             } else {1-pgamma(par[1]-log(y),shape=par[3],scale=-1*par[2])}},
             Gumbel_min=1-pgumbel(-1*y,loc=-1*par[1],scale=par[2]),
             GEV_min=1-pgev(-1*y,loc=-1*par[1],scale=par[2],shape=-1*par[3]),
             GEV_min_pos=1-pgev(-1*y,loc=-1*par[1],scale=par[2],shape=-1*par[2]/par[1]),
             NA)
  return(cdf)}

GetQuantile<-function(p,dist,par){
  #^******************************************************************************
  #^* OBJET: Retourne le p-quantile de la distribution 'dist' de paramètres 'par'  
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [real] p, probabilité entre 0 et 1 (strictement) 
  #^*    2. [character] dist, nom de la distribution 
  #^*    3. [real] par, vecteur de paramètres 
  #^* OUT
  #^*    1. [real] le p-quantile
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: si les paramètres sont impossibles, retourne NA
  #^******************************************************************************
  
  if(is.na(p)) {return(NA)}
  if(any(is.na(par))) {return(NA)}
  if(!GetParFeas(dist,par) | p<0 | p>1){return(NA)}
  if(p==0) {return(-Inf)} else if(p==1) {return(Inf)}
  q=switch(dist,
           Uniform=qunif(p,min=par[1],max=par[2]),
           Normal=qnorm(p,mean=par[1],sd=par[2]),
           LogNormal=qlnorm(p,meanlog=par[1],sdlog=par[2]),
           Gumbel=qgumbel(p,loc=par[1],scale=par[2]),
           Exponential1=qexp(p,rate=1/par[1]), 
           Exponential2=par[1]+qexp(p,rate=1/par[2]), 
           GEV=qgev(p,loc=par[1],scale=par[2],shape=-1*par[3]),
           GPD2=qgpd(p,loc=0,scale=par[1],shape=-1*par[2]), 
           GPD3=qgpd(p,loc=par[1],scale=par[2],shape=-1*par[3]), 
           Poisson=qpois(p,lambda=par[1]),
           PearsonIII={ if(par[2]>0){par[1]+qgamma(p,shape=par[3],scale=par[2])
           } else {par[1]-qgamma(1-p,shape=par[3],scale=-1*par[2])}},
           LogPearsonIII={ if(par[2]>0){exp(par[1]+qgamma(p,shape=par[3],scale=par[2]))
           } else {exp(par[1]-qgamma(1-p,shape=par[3],scale=-1*par[2]))}},
           Gumbel_min=-1*qgumbel(1-p,loc=-1*par[1],scale=par[2]),
           GEV_min=-1*qgev(1-p,loc=-1*par[1],scale=par[2],shape=-1*par[3]),
           GEV_min_pos=-1*qgev(1-p,loc=-1*par[1],scale=par[2],shape=-1*par[2]/par[1]),
           NA)
  return(q)}

Generate<-function(dist,par,n=1){
  #^******************************************************************************
  #^* OBJET: simule une réalisation de la distribution "dist" de paramètres "par"  
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [character] dist, nom de la distribution 
  #^*    2. [real] par, vecteur de paramètres 
  #^*    3. [integer] n, nombre de réalisations 
  #^* OUT
  #^*    1. [real] la (ou les) valeurs simulées
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: si les paramètres sont impossibles, retourne NA
  #^******************************************************************************
  
  if(!GetParFeas(dist,par)){return(NA)}
  r=switch(dist,
           Uniform=runif(n,min=par[1],max=par[2]),
           Normal=rnorm(n,mean=par[1],sd=par[2]),
           LogNormal=rlnorm(n,meanlog=par[1],sdlog=par[2]),
           Gumbel=rgumbel(n,loc=par[1],scale=par[2]),
           Exponential1=rexp(n,rate=1/par[1]), 
           Exponential2=par[1]+rexp(n,rate=1/par[2]), 
           GEV=rgev(n,loc=par[1],scale=par[2],shape=-1*par[3]),
           GPD2=rgpd(n,loc=0,scale=par[1],shape=-1*par[2]), 
           GPD3=rgpd(n,loc=par[1],scale=par[2],shape=-1*par[3]), 
           Poisson=rpois(n,lambda=par[1]),
           PearsonIII={ if(par[2]>0){par[1]+rgamma(n,shape=par[3],scale=par[2])
           } else {par[1]-rgamma(n,shape=par[3],scale=-1*par[2])}},
           LogPearsonIII={ if(par[2]>0){exp(par[1]+rgamma(n,shape=par[3],scale=par[2]))
           } else {exp(par[1]-rgamma(n,shape=par[3],scale=-1*par[2]))}},
           Gumbel_min=-1*rgumbel(n,loc=-1*par[1],scale=par[2]),
           GEV_min=-1*rgev(n,loc=-1*par[1],scale=par[2],shape=-1*par[3]),
           GEV_min_pos=-1*rgev(n,loc=-1*par[1],scale=par[2],shape=-1*par[2]/par[1]),
           NA)
  return(r)}

GetReducedVariate<-function(p,dist){
  #^******************************************************************************
  #^* OBJET: Retourne la variable réduite "canonique" pour la distribution 'dist'  
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: Modifie le 17/11/2020 pour la loi de Poisson: utilisation 
  #^*               d'une N(0,1) a la place d'une Poi(1) pour eviter des artefacts
  #^*               graphiques lies a la nature discrete de la distribution.
  #^*               Modifie le 05/01/2021 pour les lois XXX_min ('-1*' manquant).
  #^******************************************************************************
  #^* IN
  #^*    1. [real] p, probabilité entre 0 et 1 (strictement) 
  #^*    2. [character] dist, nom de la distribution 
  #^* OUT
  #^*    1. [real] la variable réduite
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: si les paramètres sont impossibles, retourne NA
  #^******************************************************************************
  
  if(is.na(p)) {return(NA)}
  if(p==0) {return(-Inf)} else if(p==1) {return(Inf)}
  q=switch(dist,
           Uniform=qunif(p,min=0,max=1),
           Normal=qnorm(p,mean=0,sd=1),
           LogNormal=qnorm(p,mean=0,sd=1),
           Gumbel=qgumbel(p,loc=0,scale=1),
           Exponential1=qexp(p,rate=1), 
           Exponential2=qexp(p,rate=1), 
           GEV=qgumbel(p,loc=0,scale=1),
           GPD2=qexp(p,rate=1), 
           GPD3=qexp(p,rate=1), 
           Poisson=qnorm(p,mean=0,sd=1),
           PearsonIII=qgamma(p,shape=1,scale=1),
           LogPearsonIII=qgamma(p,shape=1,scale=1),
           Gumbel_min=-1*qgumbel(1-p,loc=0,scale=1),
           GEV_min=-1*qgumbel(1-p,loc=0,scale=1),
           GEV_min_pos=-1*qgumbel(1-p,loc=0,scale=1),
           NA)
  return(q)}