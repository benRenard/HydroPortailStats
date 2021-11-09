library(HydroPortailStats)

#------------------------------------------------------------------------------------------
# PROPRIETES DES DONNEES SIMULEES
#------------------------------------------------------------------------------------------
n=30 # taille de l'échantillon
dist="GEV" #distribution: "Normal", "LogNormal", "Exponential1", "Exponential2", "Gumbel", "GPD2", "GPD3", "GEV", "PearsonIII", "LogPearsonIII", "Poisson"
param=c(100,50,-0.2) # paramètres de la distribution
y<-Generate(dist=dist,par=param,n=n) # simulation des données

#------------------------------------------------------------------------------------------
# EXEMPLE 1: Appel minimaliste
#------------------------------------------------------------------------------------------
#ooooooooooooooooooo
# Estimation en utilisant les propriétés par défaut (définies dans Hydro3Stats.R)
h3=Hydro3_Estimation(y=y,dist=dist) # objet "hydro3" contenant tous les résultats
X11(title="exemple 1: appel minimaliste");Hydro3_Plot(h3) # graphique résumant l'inférence
# utisation de l'estimation interactive Q=f(T) ou T=f(Q)
GetQfromT(T=100,H3=h3)
GetTfromQ(q=250,H3=h3)

#------------------------------------------------------------------------------------------
# EXEMPLE 2: Autres méthodes d'estimation & de quantification des incertitudes
#------------------------------------------------------------------------------------------
#ooooooooooooooooooo
# Estimation par la méthode des moments, avec du Bootstrap pour les incertitudes
h3=Hydro3_Estimation(y=y,dist=dist,Emeth="MOM",Umeth="BOOT") # objet "hydro3" contenant tous les résultats
X11(title="exemple 2: MOM + BOOT");Hydro3_Plot(h3) # graphique résumant l'inférence
#ooooooooooooooooooo
# Idem mais Bootstrap paramétrique à la place de Bootstrap (recommandé pour les extrêmes)
h3=Hydro3_Estimation(y=y,dist=dist,Emeth="MOM",Umeth="PBOOT") # objet "hydro3" contenant tous les résultats
X11(title="exemple 2: MOM + BOOT Parametrique");Hydro3_Plot(h3) # graphique résumant l'inférence
#ooooooooooooooooooo
# Maximum de vraisemblance
h3=Hydro3_Estimation(y=y,dist=dist,Emeth="ML",Umeth="ML") # objet "hydro3" contenant tous les résultats
X11(title="exemple 2: ML");Hydro3_Plot(h3) # graphique résumant l'inférence
#ooooooooooooooooooo
# On peut zapper la quantification des incertitudes (mais c'est pas bien) 
h3=Hydro3_Estimation(y=y,dist=dist,Emeth="ML",Umeth="NONE") # objet "hydro3" contenant tous les r?sultats
X11(title="exemple 2: ML sans incertitudes");Hydro3_Plot(h3) # graphique r?sumant l'inf?rence
#ooooooooooooooooooo
# Certaines combinaisons sont interdites; par exemple, la quantification des incertitudes "ML" ne marche qu'avec la méthode d'estimation "ML" 
h3=Hydro3_Estimation(y=y,dist=dist,Emeth="LMOM",Umeth="ML") # objet "hydro3" contenant tous les résultats
h3$u$message # message d'erreur, la quantification des incertitudes est ignorée
X11(title="exemple 2: Combinaison interdite Emeth=LMOM / Umeth=ML");Hydro3_Plot(h3) # graphique résumant l'inférence

#------------------------------------------------------------------------------------------
# EXEMPLE 3: Un peu de Bayesien... (temps de calcul un peu plus long)
#------------------------------------------------------------------------------------------
#ooooooooooooooooooo
# Estimation bayesienne avec les a priori par défaut (plats)
h3=Hydro3_Estimation(y=y,dist=dist,Emeth="BAY",Umeth="BAY") # objet "hydro3" contenant tous les résultats
X11(title="exemple 3: BAY");Hydro3_Plot(h3) # graphique résumant l'inférence
#ooooooooooooooooooo
# Spécification d'a prioris informatif
prior=list() # initialisation
prior[[1]]=list(dist="Uniform",par=c(0,500)) # a priori pour le 1er paramètre (paramètre de position): U[0;500]
prior[[2]]=list(dist="FlatPrior",par=NULL) # a priori pour le 2nd paramètre (paramètre d'échelle): "FlatPrior" (aucune info)
prior[[3]]=list(dist="Normal",par=c(-0.13, 0.05)) # a priori pour le 3e paramètre (paramètre de forme): N(-0.13,0.05)
h3=Hydro3_Estimation(y=y,dist=dist,Emeth="BAY",Umeth="BAY",prior=prior) # objet "hydro3" contenant tous les résultats
X11(title="exemple 3: BAY avec a prioris informatifs");Hydro3_Plot(h3) # graphique résumant l'inférence

#------------------------------------------------------------------------------------------
# EXEMPLE 4: Quelques options utiles
#------------------------------------------------------------------------------------------
#ooooooooooooooooooo
# On s'interesse aux extrêmes bas => les grandes périodes de retour sont pour les petites valeurs!
options=options_def # options par défaut, définies dans Hydro3Stats.R
options$invertT=TRUE # les grandes périodes de retour sont pour les petites valeurs
h3=Hydro3_Estimation(y=y,dist=dist,Emeth="LMOM",Umeth="PBOOT",options=options) # objet "hydro3" contenant tous les résultats
X11(title="exemple 4: options invertT=TRUE");Hydro3_Plot(h3) # graphique résumant l'inférence
#ooooooooooooooooooo
# Données issues d'un échantillonnage SUP-SEUIL avec 3 valeurs par an => il faut corriger la période de retour!
y.supseuil<-Generate(dist="GPD3",par=param,n=n) # simulation des données dans une GPD
options=options_def # options par défaut, définies dans Hydro3Stats.R
h3=Hydro3_Estimation(y=y.supseuil,dist="GPD3",Emeth="LMOM",Umeth="PBOOT",options=options) # sans correction
X11(title="exemple 4: sans correction de la periode de retour par le nombre moyen d'evenements annuels");Hydro3_Plot(h3) # graphique résumant l'inférence
options$p2T=3 # en moyenne 3 événements par an
h3=Hydro3_Estimation(y=y.supseuil,dist="GPD3",Emeth="LMOM",Umeth="PBOOT",options=options) # avec correction
X11(title="exemple 4: avec correction de la periode de retour par le nombre moyen d'evenements annuels");Hydro3_Plot(h3) # graphique résumant l'inférence
#ooooooooooooooooooo
# Présence de zéros dans les données: on peut les traiter à part en estimant la probabilité d'une valeur nulle - utile pour les basses eaux?
y<-c(y,0,0) # on rajoute deux zéros aux données
options=options_def # options par défaut, définies dans Hydro3Stats.R
options$invertT=TRUE # on s'intéresse aux extrêmes bas => les grandes périodes de retour sont pour les petites valeurs
h3=Hydro3_Estimation(y=y,dist="GEV",Emeth="LMOM",Umeth="PBOOT",options=options) # on traite les zéros comme les autres données
X11(title="exemple 4: pas de traitement particulier des zeros");Hydro3_Plot(h3) # graphique résumant l'inférence
options$splitZeros=TRUE # activation de l'option pour traiter les zéros à part
h3=Hydro3_Estimation(y=y,dist="GEV",Emeth="LMOM",Umeth="PBOOT",options=options) # les zéros sont traités à part
X11(title="exemple 4: traitement a part des zeros");Hydro3_Plot(h3) # graphique résumant l'inférence
#ooooooooooooooooooo
# Une pensée pour nos collègues non-francophones
options$lang="en" # language = anglais
h3=Hydro3_Estimation(y=y,dist="GEV",Emeth="LMOM",Umeth="PBOOT",options=options) # objet "hydro3" contenant tous les résultats
X11(title="exemple 4: nom des parametres en anglais");Hydro3_Plot(h3) # graphique résumant l'inférence

#------------------------------------------------------------------------------------------
# EXEMPLE 5: distributions pour les minima
#------------------------------------------------------------------------------------------
#ooooooooooooooooooo
# simulation des données
dist="GEV_min"
y<-Generate(dist=dist,par=c(1,0.5,0.5),n=n)
# On s'interesse aux extrêmes bas => les grandes périodes de retour sont pour les petites valeurs!
options=options_def # options par défaut, définies dans Hydro3Stats.R
options$invertT=TRUE # les grandes périodes de retour sont pour les petites valeurs
# Estimation en utilisant les propriétés par défaut (définies dans Hydro3Stats.R)
h3=Hydro3_Estimation(y=y,dist=dist,options=options) # objet "hydro3" contenant tous les résultats
X11(title="exemple 5: GEV_min");Hydro3_Plot(h3) # graphique résumant l'inférence
# Estimation par maximum de vraisemblance (la seule disponible) pour GEV_min_pos
h3=Hydro3_Estimation(y=y,dist="GEV_min_pos",Emeth="ML",Umeth="ML",options=options) # objet "hydro3" contenant tous les résultats
X11(title="exemple 5: GEV_min_pos");Hydro3_Plot(h3) # graphique résumant l'inférence
