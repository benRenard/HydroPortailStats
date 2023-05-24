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
#~* OBJET: Estimation comme dans HYDRO2
#~******************************************************************************
#~* PROGRAMMEUR: Benjamin Renard, Irstea Lyon
#~******************************************************************************
#~* CREE/MODIFIE: XXX
#~******************************************************************************
#~* PRINCIPALES FONCTIONS
#~*    1. GetEstimate_HYDRO2, estimateur des moments
#~******************************************************************************
#~* REF.: XXX
#~******************************************************************************
#~* A FAIRE: XXX
#~******************************************************************************
#~* COMMENTAIRES: 1/ Les seules distributions disponibles sont Normal, LogNormal, 
#~*                  et Gumbel
#~******************************************************************************

#****************************
# Fonctions principales ----
#****************************

#' Hydro2 estimate of a distribution
#'
#' Returns an estimate of a distribution as it was computed in the old HYDRO2 software.
#' Only available for distributions 'Normal', 'LogNormal', and 'Gumbel'.
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
#' GetEstimate_HYDRO2(y,'Normal')
#' GetEstimate_HYDRO2(y,'LogNormal')
#' GetEstimate_HYDRO2(y,'Gumbel')
#' GetEstimate_HYDRO2(y,'GEV')
#' GetEstimate_HYDRO2(y,'Poisson')
#' @export
GetEstimate_HYDRO2<-function(y,dist){
  #^******************************************************************************
  #^* OBJET: Retourne l'estimateur HYDRO2
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
             LogNormal={ReturnEstimate(rough_Normal(log(y)))},   
              # Gumbel
             Gumbel={ReturnEstimate(rough_Gumbel(y))},
            {w=Estimate_fail;w$message="fatal:distribution non disponible avec HYDRO2";w}
  )
  return(out)
}
