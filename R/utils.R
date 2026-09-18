#~******************************************************************************
#~* OBJET: Outils divers
#~******************************************************************************
#~* PROGRAMMEUR: Benjamin Renard, INRAE Aix-En-Provence
#~******************************************************************************
#~* CREE/MODIFIE: 18/09/2026
#~******************************************************************************
#~* PRINCIPALES FONCTIONS
#~*    1. Filtres
#~******************************************************************************
#~* REF.: XXX
#~******************************************************************************
#~* A FAIRE:
#~******************************************************************************
#~* COMMENTAIRES: 
#~******************************************************************************

#' Moving-statistics filter
#'
#' Apply a statistics (typically the mean) on a moving window.
#'
#' @param x numeric vector, values. 
#' @param window numeric, size of the moving window
#' @param stat function, function to apply on the moving window
#' @param ... other arguments passed to stat
#' @examples
#' x=rnorm(100)
#' y1=applyMovingStat(x,9)
#' y2=applyMovingStat(x,9,stat=min)
#' y3=applyMovingStat(x,9,stat=quantile,probs=0.75)
#' plot(x);lines(y1,col='red');lines(y2,col='blue');lines(y3,col='green')
#' @export
applyMovingStat <- function(x,window,stat=mean,...){
  out=x*NA
  n=length(x)
  ix=1:n
  for (i in 1:n){
    mask= (ix >= (i-0.5*window)) & (ix <= (i+0.5*window))
    out[i]=stat(x[mask],...)
  }
  return(out)
}
