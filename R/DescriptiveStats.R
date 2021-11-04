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
#~* OBJET: Outils pour le calcul de stat. descriptives (empiriques)
#~******************************************************************************
#~* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
#~******************************************************************************
#~* CREE/MODIFIE: XXX
#~******************************************************************************
#~* PRINCIPALES FONCTIONS
#~*    1. GetEmpFreq, calcul de la fréquence empirique au non-dépassement
#~*    2. 
#~*    3. 
#~*    4. 
#~*    5. 
#~*    6. 
#~*    7. 
#~*    8. XXX
#~******************************************************************************
#~* REF.: XXX
#~******************************************************************************
#~* A FAIRE: XXX
#~******************************************************************************
#~* COMMENTAIRES: XXX
#~******************************************************************************

#' Empirical nonexceedance frequency
#'
#' Computes the empirical nonexceedance frequency of the ith sorted value amongst n
#'
#' @param i integer or integer vector, observation rank(s)
#' @param n integer, number of observations
#' @param formula character, formula, available: 'Hazen', 'Standard', 'MinusOne', 'Weibull',
#'     'Benard', 'Cunnane', 'Beard', 'Blom', 'Gringorten', 'Landwehr', 'Tukey'.
#' @return The nonexceedance frequency.
#' @examples
#' GetEmpFreq(i=1:10,n=10)
#' GetEmpFreq(i=1:10,n=10,formula='Standard')
#' GetEmpFreq(i=1:10,n=10,formula='MinusOne')
#' GetEmpFreq(i=1:10,n=10,formula='Cunnane')
#' @export
GetEmpFreq<-function(i,n,formula="Hazen"){
  #^******************************************************************************
  #^* OBJET: Retourne la fréquence empirique au non-dépassement 
  #^******************************************************************************
  #^* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREE/MODIFIE: XXX
  #^******************************************************************************
  #^* IN
  #^*    1. [integer] i, rang 
  #^*    2. [integer] n, taille de l'échantillon 
  #^*    3. [character] formula, formule de calcul 
  #^* OUT
  #^*    1. [real] fréquence empirique au non-dépassement  
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* A FAIRE: 
  #^******************************************************************************
  #^* COMMENTAIRES: 
  #^******************************************************************************
  f=switch(formula,
           Hazen=(i-0.5)/n,
           Standard=i/n,
           MinusOne=(i-1)/n,
           Weibull=i/(n+1),
           Benard=(i-0.3)/(n+0.4),
           Cunnane=(i-0.4)/(n+0.2),
           Beard=(i-0.31)/(n+0.38),
           Blom=(i-0.375)/(n+0.25),
           Gringorten=(i-0.44)/(n+0.12),
           Landwehr=(i-0.35)/n,
           Tukey=(3*i-1)/(3*n+1),
           NA)
  return(f)
}
