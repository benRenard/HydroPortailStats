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

#' Baseflow separation
#'
#' Apply the baseflow separation (BFS) algorithm described in
#' Tallaksen and Van Lanenś book (2004):
#' Hydrological Drought: Processes and Estimation Methods for Streamflow and Groundwater. Elsevier. 
#'
#' @param x numeric vector, values. 
#' @param d integer, bloc size
#' @param w numeric, smoothing parameter
#' @examples
#' x=rnorm(1001)
#' bf=applyBFS(x)
#' plot(x,type='l')
#' lines(bf,col='red')
#' @export
#' @importFrom stats approx
applyBFS <- function(x,d=5,w=0.9){
  n=length(x)
  ix=1:n
  nr=floor(n/d)
  M=matrix(x[1:(nr*d)],ncol=d,byrow=TRUE)
  mins=apply(M,1,min)
  imins=rep(NA,length(mins))
  imins[!is.na(mins)]=apply(M[!is.na(mins),],1,which.min)
  isPivot=rep(FALSE,n)
  for(j in 2:(NROW(M)-1)){
    if(!is.na(mins[j]+mins[j-1]+mins[j+1])){
      if(w*mins[j]<min(mins[j-1],mins[j+1])){
        isPivot[(j-1)*d+imins[j]]=TRUE
      }
    }
  }
  out=approx(x=ix[isPivot],y=x[isPivot],xout=1:n,na.rm=FALSE)
  return(out$y)
}

#' Sequent Peak Algorithm (SPA)
#'
#' Apply the Sequent Peak Algorithm.
#'
#' @param x numeric vector, values. 
#' @param threshold numeric vector, threshold
#' @examples
#' x=rnorm(100)
#' spa=applySPA(x,-1)
#' plot(x,type='l')
#' lines(rep(-1,100),col='gray')
#' lines(spa,col='red')
#' @export
#' @importFrom stats quantile
applySPA <- function(x,threshold){
  n=length(x)
  if(length(threshold)==1){threshold=rep(threshold,n)}
  out=rep(0,n)
  for(j in 2:length(out)){
    out[j]=out[j-1]+threshold[j]-x[j]
    if(out[j]<0){out[j]=0}
  }
  return(out)
}

