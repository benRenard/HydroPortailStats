# Packages utilises pour realiser les figures
library(ggplot2);library(gridExtra);library(HydroPortailStats)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# REGLAGES ----
# Liste des distributions a illustrer
distList=c('Normal','LogNormal','Exponential1','Exponential2',
           'GPD2','GPD3','Gumbel','GEV','PearsonIII','LogPearsonIII',
           'Gumbel_min','GEV_min','GEV_min_pos','Poisson')
# Parametres graphiques
cols=c('#2E3D7C','#FFA73C','#BA292E') # jeu de 3 couleurs 
ptSize=2;ptAlpha=0.7 # points (taille et opacite)
lineSize=1;barSize=1 # courbes ou barres (pour Poisson)
aSize=6 # annotations sur la figure
fWidth=12;fHeight=4 # taille de la figure finale
# Parametres divers
ngrid=1000 # numbre de valeurs x ou la densite f(x) est calculee
pmax=0.995 # proba. pour calculer le min/max de la grille de x
nsim=100 # nombre de valeurs simulees
# Reproductibilite
set.seed(3)
# Fonction pour la mise en forme cummune aux 2 figures
formatFigure <- function(g){
  g=g+theme_bw()
  g=g+theme(axis.title=element_text(size=14),axis.text=element_text(size=12),
            plot.margin=margin(t=10,l=10,r=10,b=5),
            legend.text=element_text(size=12,face='bold'),legend.title=element_blank(),
            legend.background=element_rect(colour='black'),legend.text.align=0)
  return(g)
}
for(dist in distList){
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # CALCULS ----
  # Choix du (ou des) vecteur(s) de parametres.
  # Un unique vecteur pour les lois n'ayant pas de parametre de forme,
  # 3 vecteurs pour celles en ayant un
  par=switch(dist,
             Normal=list(c(50,10)),
             LogNormal=list(c(4,0.5)),
             Gumbel=list(c(50,10)),
             Exponential1=list(c(10)), 
             Exponential2=list(c(50,10)), 
             GPD2=list(c(10,-0.3),c(10,0),c(10,0.3)),
             GPD3=list(c(50,10,-0.3),c(50,10,0),c(50,10,0.3)), 
             GEV=list(c(50,10,-0.3),c(50,10,0),c(50,10,0.3)),
             PearsonIII=list(c(50,10,2),c(50,10,1.5),c(50,10,1)),
             LogPearsonIII=list(c(4,0.1,1.5),c(4,0.1,1.25),c(4,0.1,1)),
             Gumbel_min=list(c(50,10)),
             GEV_min=list(c(50,10,-0.3),c(50,10,0),c(50,10,0.3)),
             GEV_min_pos=list(c(50,40)),
             Poisson=list(c(0.5),c(2),c(3)),
             NA)
  # Symboles pour chaque parametre
  symb=switch(dist,
              Normal=c('mu','sigma'),
              LogNormal=c('mu','sigma'),
              Gumbel=c('mu','sigma'),
              Exponential1=c('sigma'), 
              Exponential2=c('mu','sigma'), 
              GPD2=c('sigma','xi'),
              GPD3=c('mu','sigma','xi'), 
              GEV=c('mu','sigma','xi'),
              PearsonIII=c('mu','sigma','xi'),
              LogPearsonIII=c('mu','sigma','xi'),
              Gumbel_min=c('mu','sigma'),
              GEV_min=c('mu','sigma','xi'),
              GEV_min_pos=c('mu','sigma'),
              Poisson=c('lambda'),
              NA)
  # Texte a utiliser pour la legende des figures (nom de la distribution et parametres)
  dname=c()
  for(i in 1:length(par)){
    txt=c()
    for(j in 1:length(par[[i]])){
      txt=c(txt,paste0(symb[j],'==',par[[i]][j]))
    }
    dname=c(dname,paste0(dist,'(',paste0(txt,collapse=','),')'))
  }
  # Calcul des valeurs x ou la densite est evaluee
  xhigh=GetQuantile(pmax,dist,par[[1]])
  xlow=GetQuantile(1-pmax,dist,par[[1]])
  xwidth=xhigh-xlow
  xlow=xlow-0.1*xwidth;xhigh=xhigh+0.1*xwidth
  x=seq(xlow,xhigh,length.out=ngrid)
  # Cas particulier de la loi de Poisson (discrete)
  if(dist == 'Poisson'){x=0:8}
  # Calul de la densite et des donnees simulees
  densite=data.frame();sim=data.frame()
  for(i in 1:length(par)){
    param=par[[i]]
    y=GetPdf(x,dist,param)
    densite=rbind(densite,data.frame(x=x,y=y,label=dname[i]))
    y=Generate(dist,param,nsim)
    sim=rbind(sim,data.frame(x=1:nsim,y=y,label=dname[i]))
  }
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # FIGURE 1: DENSITE ----
  g=ggplot(densite,aes(x=x,y=y))
  # texte de la legende
  legendLabs=parse(text=as.character(unique(densite$label)))
  if(dist=='Poisson'){ # barres verticales
    g=g+geom_col(aes(fill=label),position='dodge',size=barSize)
    g=g+scale_fill_manual(values=cols,labels=legendLabs)
    g=g+labs(x='x',y='probabilité f(x)')
  } else { # courbe
    g=g+geom_line(aes(col=label),size=lineSize)
    g=g+scale_color_manual(values=cols,labels=legendLabs)
    g=g+labs(x='x',y='densité f(x)')
  }
  # Annotations: mu
  if(dist %in% c('Normal','LogNormal','Exponential1','Exponential2',
                 'Gumbel','Gumbel_min','GPD2','GPD3','PearsonIII',
                 'LogPearsonIII','GEV','GEV_min')){
    x0=par[[1]][1];txt='mu'
    if(dist %in% c('LogNormal','LogPearsonIII')){x0=exp(x0);txt='e^{mu}'}
    if(dist%in% c('Exponential1','GPD2')){x0=0;txt=''}
    g=g+geom_vline(xintercept=x0,col='black',linetype='dashed')
    g=g+annotate('text',x=x0,y=0,label=parse(text=txt),size=aSize,hjust=-0.5,vjust=0)
    if(dist=='LogNormal'){
      g=g+annotate('text',x=x0,y=0,label='(mediane)',size=aSize,hjust=-0.4,vjust=0)
    }
  }
  # Annotations: mu+sigma
  if(dist %in% c('Normal','Exponential1','Exponential2','Gumbel','Gumbel_min')){
    if(dist=='Exponential1'){
      x1=par[[1]][1]
    } else {
      x1=par[[1]][1]+par[[1]][2]
    }
    yArrow=GetPdf(x1,dist,param)
    g=g+geom_vline(xintercept=x1,col='black',linetype='dotted')
    g=g+annotate('segment',x=x0,y=yArrow,xend=x1,yend=yArrow,arrow=arrow(ends='both',length=unit(2,'mm')))
    g=g+annotate('text',x=0.5*(x0+x1),y=yArrow,label=parse(text='sigma'),size=aSize,hjust=0.5,vjust=1.5)
  }
  # Annotations: mu-sigma
  if(dist %in% c('Normal','Gumbel','Gumbel_min')){
    x2=par[[1]][1]-par[[1]][2]
    g=g+geom_vline(xintercept=x2,col='black',linetype='dotted')
  }
  # position de la legende, mise en forme
  if(dist %in% c('GEV_min','Gumbel_min')){ # en haut a gauche
    lpos=c(0.01,0.99);ljust=c(0,1)
  } else { # en haut a droite
    lpos=0.99*c(1,1);ljust=c(1,1) 
  }
  g1=formatFigure(g)+theme(legend.justification=ljust,legend.position=lpos)
  if(dist=='Poisson'){g1=g1+scale_x_continuous(breaks=x)}
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # FIGURE 2: DONNEES SIMULEES ----
  g=ggplot(sim,aes(x=x,y=y))
  g=g+geom_point(aes(col=label),size=ptSize,alpha=ptAlpha)
  g=g+scale_color_manual(values=cols)
  g=g+labs(x='index',y='données simulées')
  # mise en forme
  g2=formatFigure(g)+theme(legend.position='none')
  if(dist=='Poisson'){g2=g2+scale_y_continuous(breaks=0:max(sim$y))}
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Sauvegarde du fichier ----
  pdf(file=paste0(dist,'.pdf'),width=fWidth,height=fHeight,useDingbats=FALSE)
  grid.arrange(g1,g2,nrow=1)
  dev.off()
}
