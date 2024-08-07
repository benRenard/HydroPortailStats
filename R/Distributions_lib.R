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

#' Number of parameters.
#'
#' Returns the number of parameters of a distribution.
#'
#' @param dist character, distribution name
#' @return An integer.
#' @examples
#' GetParNumber('Normal')
#' GetParNumber('GEV')
#' @export
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
              Triangle=3,
              NA)
  return(Npar)}

#' Parameter names.
#'
#' Returns the names of the parameters of a distribution, 
#' in French (default) or English.
#'
#' @param dist character, distribution name
#' @param lang character, language ('en' or 'fr')
#' @return A character vector.
#' @examples
#' GetParName('Normal')
#' GetParName('GEV')
#' GetParName('GEV',lang='en')
#' @export
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
              Uniform=switch(lang,fr=c('borne_inf','borne_sup'),en=c('lower_bound','higher_bound'),NA), 
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
              Triangle=switch(lang,fr=c('pic','borne_inf','borne_sup'),en=c('peak','lower_bound','higher_bound'),NA),
              NA)
  return(name)}

#' Parameter feasibility
#'
#' Evaluates whether a parameter vector is feasible 
#' (for instance, are scale parameters >0 ?)
#'
#' @param dist character, distribution name
#' @param par numeric vector, parameter vector
#' @return A logical.
#' @examples
#' # Feasible
#' GetParFeas('Normal',c(0,1))
#' # Not feasible because second parameter (standard deviation) is negative
#' GetParFeas('Normal',c(0,-1))
#' @export
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
              FlatPrior=TRUE,
              Uniform={if(par[2]<=par[1]){FALSE} else {TRUE}},
              Normal={if(par[2]<=0){FALSE} else {TRUE}},
              LogNormal={if(par[2]<=0){FALSE} else {TRUE}},
              Gumbel={if(par[2]<=0){FALSE} else {TRUE}},
              Exponential1={if(par[1]<=0){FALSE} else {TRUE}}, 
              Exponential2={if(par[2]<=0){FALSE} else {TRUE}}, 
              GEV={if(par[2]<=0){FALSE} else {TRUE}},
              GPD2={if(par[1]<=0){FALSE} else {TRUE}}, 
              GPD3={if(par[2]<=0){FALSE} else {TRUE}}, 
              Poisson={if(par[1]<=0){FALSE} else {TRUE}},
              PearsonIII={if( (par[2]==0) | (par[3]<=0) ){FALSE} else {TRUE}},
              LogPearsonIII={if( (par[2]==0) | (par[3]<=0) ){FALSE} else {TRUE}},
              Gumbel_min={if(par[2]<=0){FALSE} else {TRUE}},
              GEV_min={if(par[2]<=0){FALSE} else {TRUE}},
              GEV_min_pos={if( (par[1]<=0) | (par[2]<=0)){FALSE} else {TRUE}},
              Triangle={if( (par[3]<=par[2]) | (par[3]<=par[1]) | (par[1]<=par[2])){FALSE} else {TRUE}},
              NA)
  return(feas)}

#' Probability Density Function (pdf)
#'
#' Evaluates the pdf of a distribution
#'
#' @param y numeric, value at which the pdf is evaluated
#' @param dist character, distribution name
#' @param par numeric vector, parameter vector
#' @param log logical, returns log-pdf if TRUE
#' @return The pdf or the log-pdf as a numeric.
#' @examples
#' GetPdf(0,'Normal',c(0,1))
#' GetPdf(200,'GEV',c(100,25,-0.2))
#' GetPdf(200,'GEV',c(100,25,0.2))
#' GetPdf(3,'Poisson',0.75)
#' @importFrom stats dunif dnorm dlnorm dexp dgamma dpois
#' @importFrom evd dgev dgpd dgumbel
#' @export
GetPdf<-function(y,dist,par,log=FALSE){
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
             Uniform=stats::dunif(y,min=par[1],max=par[2],log=log),
             Normal=stats::dnorm(y,mean=par[1],sd=par[2],log=log),
             LogNormal=stats::dlnorm(y,meanlog=par[1],sdlog=par[2],log=log),
             Gumbel=evd::dgumbel(y,loc=par[1],scale=par[2],log=log),
             Exponential1=stats::dexp(y,rate=1/par[1],log=log), 
             Exponential2=stats::dexp(y-par[1],rate=1/par[2],log=log), 
             GEV=evd::dgev(y,loc=par[1],scale=par[2],shape=-1*par[3],log=log),
             GPD2={if(y==0){ if(log==FALSE) {1/par[1]} else {-log(par[1])}
             } else {evd::dgpd(y,loc=0,scale=par[1],shape=-1*par[2],log=log)}}, 
             GPD3={if(y==par[1]){if(log==FALSE) {1/par[2]} else {-log(par[2])}
             } else {evd::dgpd(y,loc=par[1],scale=par[2],shape=-1*par[3],log=log)}}, 
             Poisson=stats::dpois(y,lambda=par[1],log=log),
             PearsonIII={ if(par[2]>0){stats::dgamma(y-par[1],shape=par[3],scale=par[2],log=log)
             } else {stats::dgamma(par[1]-y,shape=par[3],scale=-1*par[2],log=log)}},
             LogPearsonIII={ if(par[2]>0){fy=stats::dgamma(log(y)-par[1],shape=par[3],scale=par[2],log=TRUE)-log(y)
             } else {fy=stats::dgamma(par[1]-log(y),shape=par[3],scale=-1*par[2],log=TRUE)-log(y)}
             if(log==TRUE) {fy
             } else {exp(fy)}},
             Gumbel_min=evd::dgumbel(-1*y,loc=-1*par[1],scale=par[2],log=log),
             GEV_min=evd::dgev(-1*y,loc=-1*par[1],scale=par[2],shape=-1*par[3],log=log),
             GEV_min_pos=evd::dgev(-1*y,loc=-1*par[1],scale=par[2],shape=-1*par[2]/par[1],log=log),
             Triangle=dtriangle(y,peak=par[1],min=par[2],max=par[3],log=log),
             NA)
  return(pdf)}

#' Cumulative Distribution Function (cdf)
#'
#' Evaluates the cdf of a distribution
#'
#' @param y numeric, value at which the cdf is evaluated
#' @param dist character, distribution name
#' @param par numeric vector, parameter vector
#' @return The cdf as a numeric.
#' @examples
#' GetCdf(0,'Normal',c(0,1))
#' GetCdf(200,'GEV',c(100,25,-0.2))
#' GetCdf(200,'GEV',c(100,25,0.2))
#' GetCdf(3,'Poisson',0.75)
#' @importFrom stats punif pnorm plnorm pexp pgamma ppois
#' @importFrom evd pgev pgpd pgumbel
#' @export
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
             Uniform=stats::punif(y,min=par[1],max=par[2]),
             Normal=stats::pnorm(y,mean=par[1],sd=par[2]),
             LogNormal=stats::plnorm(y,meanlog=par[1],sdlog=par[2]),
             Gumbel=evd::pgumbel(y,loc=par[1],scale=par[2]),
             Exponential1=stats::pexp(y,rate=1/par[1]), 
             Exponential2=stats::pexp(y-par[1],rate=1/par[2]), 
             GEV=evd::pgev(y,loc=par[1],scale=par[2],shape=-1*par[3]),
             GPD2=evd::pgpd(y,loc=0,scale=par[1],shape=-1*par[2]), 
             GPD3=evd::pgpd(y,loc=par[1],scale=par[2],shape=-1*par[3]), 
             Poisson=stats::ppois(y,lambda=par[1]),
             PearsonIII={ if(par[2]>0){stats::pgamma(y-par[1],shape=par[3],scale=par[2])
             } else {1-stats::pgamma(par[1]-y,shape=par[3],scale=-1*par[2])}},
             LogPearsonIII={ if(par[2]>0){stats::pgamma(log(y)-par[1],shape=par[3],scale=par[2])
             } else {1-stats::pgamma(par[1]-log(y),shape=par[3],scale=-1*par[2])}},
             Gumbel_min=1-evd::pgumbel(-1*y,loc=-1*par[1],scale=par[2]),
             GEV_min=1-evd::pgev(-1*y,loc=-1*par[1],scale=par[2],shape=-1*par[3]),
             GEV_min_pos=1-evd::pgev(-1*y,loc=-1*par[1],scale=par[2],shape=-1*par[2]/par[1]),
             Triangle=ptriangle(y,peak=par[1],min=par[2],max=par[3]),
             NA)
  return(cdf)}

#' Quantile Function
#'
#' Evaluates the quantiles of a distribution
#'
#' @param p numeric in (0;1), nonexceedance probability
#' @param dist character, distribution name
#' @param par numeric vector, parameter vector
#' @return The p-quantile as a numeric.
#' @examples
#' GetQuantile(0.99,'Normal',c(0,1))
#' GetQuantile(0.99,'GEV',c(100,25,-0.2))
#' GetQuantile(0.99,'GEV',c(100,25,0.2))
#' GetQuantile(0.99,'Poisson',0.75)
#' @importFrom stats qunif qnorm qlnorm qexp qgamma qpois
#' @importFrom evd qgev qgpd qgumbel
#' @export
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
           Uniform=stats::qunif(p,min=par[1],max=par[2]),
           Normal=stats::qnorm(p,mean=par[1],sd=par[2]),
           LogNormal=stats::qlnorm(p,meanlog=par[1],sdlog=par[2]),
           Gumbel=evd::qgumbel(p,loc=par[1],scale=par[2]),
           Exponential1=stats::qexp(p,rate=1/par[1]), 
           Exponential2=par[1]+stats::qexp(p,rate=1/par[2]), 
           GEV=evd::qgev(p,loc=par[1],scale=par[2],shape=-1*par[3]),
           GPD2=evd::qgpd(p,loc=0,scale=par[1],shape=-1*par[2]), 
           GPD3=evd::qgpd(p,loc=par[1],scale=par[2],shape=-1*par[3]), 
           Poisson=stats::qpois(p,lambda=par[1]),
           PearsonIII={ if(par[2]>0){par[1]+stats::qgamma(p,shape=par[3],scale=par[2])
           } else {par[1]-stats::qgamma(1-p,shape=par[3],scale=-1*par[2])}},
           LogPearsonIII={ if(par[2]>0){exp(par[1]+stats::qgamma(p,shape=par[3],scale=par[2]))
           } else {exp(par[1]-stats::qgamma(1-p,shape=par[3],scale=-1*par[2]))}},
           Gumbel_min=-1*evd::qgumbel(1-p,loc=-1*par[1],scale=par[2]),
           GEV_min=-1*evd::qgev(1-p,loc=-1*par[1],scale=par[2],shape=-1*par[3]),
           GEV_min_pos=-1*evd::qgev(1-p,loc=-1*par[1],scale=par[2],shape=-1*par[2]/par[1]),
           Triangle=qtriangle(p,peak=par[1],min=par[2],max=par[3]),
           NA)
  return(q)}

#' Random numbers generator
#'
#' Generate random realizations from a distribution
#'
#' @param dist character, distribution name
#' @param par numeric vector, parameter vector
#' @param n integer, number of values to generate
#' @return The generated values as a numeric vector.
#' @examples
#' Generate('Normal',c(0,1),10)
#' Generate('GEV',c(100,25,-0.2),10)
#' Generate('GEV',c(100,25,0.2),10)
#' Generate('Poisson',0.75,10)
#' @importFrom stats runif rnorm rlnorm rexp rgamma rpois
#' @importFrom evd rgev rgpd rgumbel
#' @export
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
           Uniform=stats::runif(n,min=par[1],max=par[2]),
           Normal=stats::rnorm(n,mean=par[1],sd=par[2]),
           LogNormal=stats::rlnorm(n,meanlog=par[1],sdlog=par[2]),
           Gumbel=evd::rgumbel(n,loc=par[1],scale=par[2]),
           Exponential1=stats::rexp(n,rate=1/par[1]), 
           Exponential2=par[1]+stats::rexp(n,rate=1/par[2]), 
           GEV=evd::rgev(n,loc=par[1],scale=par[2],shape=-1*par[3]),
           GPD2=evd::rgpd(n,loc=0,scale=par[1],shape=-1*par[2]), 
           GPD3=evd::rgpd(n,loc=par[1],scale=par[2],shape=-1*par[3]), 
           Poisson=stats::rpois(n,lambda=par[1]),
           PearsonIII={ if(par[2]>0){par[1]+stats::rgamma(n,shape=par[3],scale=par[2])
           } else {par[1]-stats::rgamma(n,shape=par[3],scale=-1*par[2])}},
           LogPearsonIII={ if(par[2]>0){exp(par[1]+stats::rgamma(n,shape=par[3],scale=par[2]))
           } else {exp(par[1]-stats::rgamma(n,shape=par[3],scale=-1*par[2]))}},
           Gumbel_min=-1*evd::rgumbel(n,loc=-1*par[1],scale=par[2]),
           GEV_min=-1*evd::rgev(n,loc=-1*par[1],scale=par[2],shape=-1*par[3]),
           GEV_min_pos=-1*evd::rgev(n,loc=-1*par[1],scale=par[2],shape=-1*par[2]/par[1]),
           Triangle=rtriangle(n,peak=par[1],min=par[2],max=par[3]),
           NA)
  return(r)}

#' Constrained random numbers generator
#'
#' Generate random realizations from a distribution, 
#' constraining these realizations to stay within bounds.
#'
#' @param dist character, distribution name
#' @param par numeric vector, parameter vector
#' @param n integer, number of values to generate
#' @param lowerBound Numeric, lower bound
#' @param higherBound Numeric, higher bound, should be strictly larger than the lower bound
#' @return The generated values as a numeric vector.
#' @examples
#' set.seed(123456)
#' y0=GenerateWithinBounds(dist='GEV',par=c(0,1,-0.2),n=1000)
#' y1=GenerateWithinBounds(dist='GEV',par=c(0,1,-0.2),n=1000,lowerBound=0,higherBound=5)
#' plot(y0);points(y1,col='red')
#' @export
GenerateWithinBounds<-function(dist,par,n=1,lowerBound=-Inf,higherBound=Inf){
  if(!GetParFeas(dist,par)){return(NA)}
  if(higherBound<=lowerBound){
    message('Fatal: bounds are not in increasing order')
    return(NA)
  }
  u0=runif(n) # uniform numbers in (0;1)
  pLow=GetCdf(lowerBound,dist,par) # nonexceedance prob associated with lower bound
  pHigh=GetCdf(higherBound,dist,par) # nonexceedance prob associated with higher bound
  u=pLow+u0*(pHigh-pLow) # rescale u0 between plow and phigh
  out=sapply(u,GetQuantile,dist=dist,par=par) # transform u into deviates from dist
  return(out)
}

#' Reduced variate
#'
#' Returns the 'reduced variate' that is used in some quantile plots 
#' (see e.g. quantile curve on Gumbel paper)
#'
#' @param p numeric in (0;1), nonexceedance probability
#' @param dist character, distribution name
#' @return The reduced variate with nonexceedance probability p.
#' @examples
#' GetReducedVariate(0.99,'Normal')
#' GetReducedVariate(0.99,'Gumbel')
#' GetReducedVariate(0.99,'GEV')
#' GetReducedVariate(0.99,'Poisson')
#' @importFrom stats qunif qnorm qexp qgamma
#' @importFrom evd qgumbel
#' @export
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
           Uniform=stats::qunif(p,min=0,max=1),
           Normal=stats::qnorm(p,mean=0,sd=1),
           LogNormal=stats::qnorm(p,mean=0,sd=1),
           Gumbel=evd::qgumbel(p,loc=0,scale=1),
           Exponential1=stats::qexp(p,rate=1), 
           Exponential2=stats::qexp(p,rate=1), 
           GEV=evd::qgumbel(p,loc=0,scale=1),
           GPD2=stats::qexp(p,rate=1), 
           GPD3=stats::qexp(p,rate=1), 
           Poisson=stats::qnorm(p,mean=0,sd=1),
           PearsonIII=stats::qgamma(p,shape=1,scale=1),
           LogPearsonIII=stats::qgamma(p,shape=1,scale=1),
           Gumbel_min=-1*evd::qgumbel(1-p,loc=0,scale=1),
           GEV_min=-1*evd::qgumbel(1-p,loc=0,scale=1),
           GEV_min_pos=-1*evd::qgumbel(1-p,loc=0,scale=1),
           Triangle=stats::qunif(p,min=0,max=1),
           NA)
  return(q)}


#****************************
# Private functions ----
#****************************

dtriangle <- function(x,peak=0,min=-1,max=1,log=FALSE){
  out=rep(0,length(x))
  out[ x<=peak & x>min ] = (2*(x[x<=peak & x>min]-min))/((max-min)*(peak-min))
  out[ x>peak & x<max ]  = (2*(max-x[x>peak & x<max]))/((max-min)*(max-peak))
  if(log==TRUE) {out=log(out)}
  return(out)
}

ptriangle <- function(x,peak=0,min=-1,max=1){
  out=rep(NA,length(x))
  out[x<=min]=0
  out[x>=max]=1
  out[ x<=peak & x>min ] = ((x[x<=peak & x>min]-min)**2)/((max-min)*(peak-min))
  out[ x>peak & x<max ]  = 1-((max-x[x>peak & x<max])**2)/((max-min)*(max-peak))
  return(out)
}

qtriangle <- function(p,peak=0,min=-1,max=1){
  r=(peak-min)/(max-min)
  out=rep(NA,length(p))
  out[p==0]=min
  out[p==1]=max
  out[ p>0 & p<=r ] = min+sqrt(p[p>0 & p<=r]*(max-min)*(peak-min))
  out[ p<1 & p>r ]  = max-sqrt((1-p[p<1 & p>r])*(max-min)*(max-peak))
  return(out)
}

rtriangle <- function(n,peak=0,min=-1,max=1){
  u=runif(n)
  out=qtriangle(u,peak,min,max)
  return(out)
}
