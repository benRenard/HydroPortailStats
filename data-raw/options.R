options_def=list(FreqFormula="Hazen",pgrid=seq(0.001,0.999,0.001),
                 Tgrid=c(seq(1.1,1.9,0.1),seq(2,9,1),seq(10,90,10),seq(100,1000,100)),
                 IClevel=0.9,p2T=1,invertT=F,splitZeros=F,lang="fr",
                 nsim=1000)
save(options_def,file='../data/options_def.RData')

mcmcoptions_def=list(mult=0.1,eps=0.1,batch.length=100,batch.n=100,
                     moverate.min=0.1,moverate.max=0.5,mult.down=0.9, mult.up=1.1,
                     burn=0.5,slim=5)
save(mcmcoptions_def,file='../data/mcmcoptions_def.RData')
