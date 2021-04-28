library(ggplot2)
library(coda)
library("colorblindr")
library(rjson)


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


logs <- list.files(path="./combined/", pattern=paste("*all*.log", sep=""), full.names = TRUE)

dat = data.frame()
dat.ratio = data.frame()

human.covs = c("229E", "OC43", "NL63")

for (i in seq(1, length(logs))){
  if (!grepl("sars2",logs[[i]])){
    t = read.table(logs[[i]], header=T, sep="\t")
    tmp = strsplit(logs[[i]], split="[/\\.]")[[1]]
    seq.length = 29000#l.dat[which(l.dat$virus==tmp[[length(tmp)-1]]), "length"]
    virus = toupper(gsub("_all","",tmp[[length(tmp)-1]]))
    virus = gsub("LIKE", "like", virus)
    
    mean.vals = (t$relativeRecombinationRate1+t$relativeRecombinationRate2+t$relativeRecombinationRate3)/3

    
    dat = rbind(dat, data.frame(recombination.rate=t$recombinationRate*exp(mean.vals), 
                                Ne=t$popSize, 
                                clock=t$mut1, 
                                virus=virus, 
                                length=seq.length))
    
    
    dat.ratio = rbind(dat.ratio, data.frame(recombination.rate=t$relativeRecombinationRate1, 
                                            mean.rate = (t$relativeRecombinationRate2 + t$relativeRecombinationRate3)/2,
                                            virus=virus, 
                                            protein="5\' of spike"))
    dat.ratio = rbind(dat.ratio, data.frame(recombination.rate=t$relativeRecombinationRate2, 
                                            mean.rate = (t$relativeRecombinationRate1 + t$relativeRecombinationRate3)/2,
                                            virus=virus, 
                                            protein="spike"))
    dat.ratio = rbind(dat.ratio, data.frame(recombination.rate=t$relativeRecombinationRate3, 
                                            mean.rate = (t$relativeRecombinationRate2 + t$relativeRecombinationRate1)/2,
                                            virus=virus, 
                                            protein="3\' of spike"))
  }
}

dat$virus <- factor(dat$virus, levels = c("SARS-like","MERS", "229E", "OC43", "NL63"))
dat.ratio$virus <- factor(dat.ratio$virus, levels = c("SARS-like","MERS", "229E", "OC43", "NL63"))

use.covs = c("MERS", "229E", "OC43", "NL63")

for (i in seq(1,length(use.covs))){
  ind = which(dat$virus==use.covs[[i]])
  print(use.covs[[i]])
  interval = HPDinterval(as.mcmc(dat[ind,"clock"]))
  print(paste("median", median(dat[ind,"clock"]), "CI:", interval[1,"lower"], interval[1,"upper"]))
}

p.rec = ggplot(dat[which(dat$virus!="SARS-like"),]) +
  geom_violin(aes(x=virus, y=recombination.rate)) +
  scale_y_log10() +
  theme_minimal() +
  ylab("recombination rate pre site and year") +
  xlab("")

plot(p.rec)
ggsave(plot=p.rec, paste("../../Recombination-Text/Figures/RecombinationRates.pdf", sep=""), height=3, width=9)


dat.withflu = dat[which(dat$virus!="SARS-like"),c("recombination.rate", "virus", "length")]
dat.withflu$type = "Corona virus"
t = read.table("../../Reassortment-Material/Applications/H3N2/combined/h3n2_sub0.log", header=T, sep="\t")
dat.withflu = rbind(dat.withflu, data.frame(recombination.rate=t$reassortmentRate, virus="Influenza A/H3N2", length=1, type="Influenza virus"))
t = read.table("../../Reassortment-Material/Applications/InfB/combined/infB_sub0.log", header=T, sep="\t")
dat.withflu = rbind(dat.withflu, data.frame(recombination.rate=t$reassortmentRate, virus="Influenza B", length=1, type="Influenza virus"))


dat.plot = data.frame()
u = unique(dat.withflu$virus)
dat.withflu$yearrate = dat.withflu$recombination.rate*dat.withflu$length
for (i in seq(1,length(u))){
  interval = HPDinterval(as.mcmc(dat.withflu[which(dat.withflu$virus==u[[i]]), "yearrate"]))
  dat.plot = rbind(dat.plot, data.frame(virus=u[[i]], median=median(dat.withflu[which(dat.withflu$virus==u[[i]]), "yearrate"]),
                                        l=interval[1,"lower"], u=interval[1,'upper']))
}


p.rec.lin = ggplot(dat.plot) +
  geom_point(aes(x=virus, y=median), size=2)+
  geom_linerange(aes(x=virus, ymin=l, ymax=u)) +
  # geom_boxplot(aes(x=virus, y=recombination.rate*length), size=1) +
  scale_y_log10() +
  theme_minimal() +
  # scale_color_OkabeIto() +
  # scale_fill_OkabeIto() +
  ylab("recombination/reassortment events\nper lineage and year") +
  xlab("")
plot(p.rec.lin)

ggsave(plot=p.rec.lin, paste("../../Recombination-Text/Figures/networks/RecombinationYear.pdf", sep=""), height=3, width=9)

p.Ne = ggplot(dat[which(dat$virus!="SARS-like"),]) +
  geom_violin(aes(x=virus, y=Ne)) +
  scale_y_log10() +
  theme_minimal() +
  ylab("effective population size") +
  xlab("")
plot(p.Ne)

ggsave(plot=p.Ne, paste("../../Recombination-Text/Figures/Ne.pdf", sep=""), height=3, width=5)

p.clock = ggplot(dat[which(dat$virus!="SARS-like"),]) +
  geom_violin(aes(x=virus, y=clock)) +
  # scale_y_log10(limits=c(0.0001,0.001)) +
  theme_minimal() +
  ylab("evolutionary rate per site and year") +
  xlab("")
plot(p.clock)
ggsave(plot=p.clock, paste("../../Recombination-Text/Figures/Clock.pdf", sep=""), height=3, width=5)


fnames = data.frame(virus="OC43",protein="5\' of spike",section="replicase1ab")
fnames = rbind(fnames, data.frame(virus="OC43",protein="5\' of spike",section="he"))
fnames = rbind(fnames, data.frame(virus="OC43",protein="5\' of spike",section="nonstructural2a"))
fnames = rbind(fnames, data.frame(virus="229E",protein="5\' of spike",section="replicase1ab"))
fnames = rbind(fnames, data.frame(virus="NL63",protein="5\' of spike",section="replicase1ab"))
fnames = rbind(fnames, data.frame(virus="OC43",protein="spike",section="spike"))
fnames = rbind(fnames, data.frame(virus="229E",protein="spike",section="spike"))
fnames = rbind(fnames, data.frame(virus="NL63",protein="spike",section="spike"))
fnames = rbind(fnames, data.frame(virus="OC43",protein="3\' of spike",section="nonstructural2"))
fnames = rbind(fnames, data.frame(virus="OC43",protein="3\' of spike",section="envelope"))
fnames = rbind(fnames, data.frame(virus="OC43",protein="3\' of spike",section="membrane"))
fnames = rbind(fnames, data.frame(virus="OC43",protein="3\' of spike",section="nucleocapsid"))
fnames = rbind(fnames, data.frame(virus="OC43",protein="3\' of spike",section="n2protein"))
fnames = rbind(fnames, data.frame(virus="229E",protein="3\' of spike",section="envelope"))
fnames = rbind(fnames, data.frame(virus="229E",protein="3\' of spike",section="membrane"))
fnames = rbind(fnames, data.frame(virus="229E",protein="3\' of spike",section="nucleocapsid"))
fnames = rbind(fnames, data.frame(virus="229E",protein="3\' of spike",section="protein4a"))
fnames = rbind(fnames, data.frame(virus="229E",protein="3\' of spike",section="protein4b"))
fnames = rbind(fnames, data.frame(virus="NL63",protein="3\' of spike",section="protein3"))
fnames = rbind(fnames, data.frame(virus="NL63",protein="3\' of spike",section="envelope"))
fnames = rbind(fnames, data.frame(virus="NL63",protein="3\' of spike",section="membrane"))
fnames = rbind(fnames, data.frame(virus="NL63",protein="3\' of spike",section="nucleocapsid"))

adapt.vals = data.frame()
u = unique(fnames$protein)
require(seqinr)
for (i in seq(1,length(human.covs))){
  prot.vals = data.frame();
  for (j in seq(1,length(u))){
    ind = which(fnames$virus==human.covs[[i]] & fnames$protein==u[[j]])  
    totLength = 0
    vals = rep(0,100)
    for (k in seq(1,length(ind))){
      fname = paste("../../seasonal-cov/antigenic_evolution/bhatt_results/", fnames[ind[[k]], "virus"], 
      "_", fnames[ind[[k]], "section"],"_bhatt_analysis_bootstrapped.json", sep="")
      fname2 = paste("../../seasonal-cov/", fnames[ind[[k]], "virus"], 
                    "/data/",fnames[ind[[k]], "virus"], "_", fnames[ind[[k]], "section"],".fasta", sep="")
      
      adapt = fromJSON(file=fname)
      fas = read.fasta(file=fname2, as.string=T,seqonly=T)
      length = nchar(fas[[1]][[1]])
      totLength = totLength+length
      vals = vals + adapt$bootstrap_rate_of_adaptation * length
    }
    prot.vals = rbind(prot.vals, data.frame(protein=u[[j]], vals=vals/totLength))
  }
  
  for (j in seq(1,length(u))){
    a=seq(1,3)
    a=a[a!=j]

    vals1 = prot.vals[which(prot.vals$protein==u[[j]]), "vals"]
    vals2 = prot.vals[which(prot.vals$protein==u[[a[[1]]]]), "vals"]
    vals3 = prot.vals[which(prot.vals$protein==u[[a[[2]]]]), "vals"]

    vals4 = vals1/(0.5*(vals2+vals3))
    interval = HPDinterval(as.mcmc(vals4), prob=0.95)
    interval.8 = HPDinterval(as.mcmc(vals4), prob=0.8)
    interval.5 = HPDinterval(as.mcmc(vals4), prob=0.5)

    adapt.vals = rbind(adapt.vals, data.frame(u.95 = interval[1,"upper"], l.95 = interval[1,"lower"],
                                                u.8 = interval.8[1,"upper"], l.8 = interval.8[1,"lower"],
                                                u.5 = interval.5[1,"upper"], l.5 = interval.5[1,"lower"],
                                                median = median(vals4),
                                                virus = human.covs[[i]],
                                                protein=u[[j]]))
  }
}


# read in the rates of adaptation
require("rjson")

logs <- list.files(path="./combined/", pattern=paste("*subunit.*.log", sep=""), full.names = TRUE)

dat = data.frame()
adatation = data.frame()
adapt.ratio = data.frame()
for (i in seq(1, length(logs))){
  t = read.table(logs[[i]], header=T, sep="\t")
  tmp = strsplit(logs[[i]], split="[/\\.]")[[1]]
  tmp2 = strsplit(tmp[[5]], split="_")[[1]]

  fname = paste("../../seasonal-cov/antigenic_evolution/bhatt_results/", tmp2[[1]], "_s1_bhatt_analysis_bootstrapped.json", sep="")
  
  dat = rbind(dat, data.frame(rate=exp(t$relativeRecombinationRate1-t$relativeRecombinationRate2),
                              protein="recombination",
                              virus=gsub("LIKE","like",toupper(tmp2[[1]]))))

  
  if (file.exists(fname)){
    adapt = fromJSON(file=fname)
    adapt2 = fromJSON(file=gsub("s1","s2",fname))
    
    adapt2$bootstrap_rate_of_adaptation[which(adapt2$bootstrap_rate_of_adaptation<=0)] = min(adapt2$bootstrap_rate_of_adaptation[which(adapt2$bootstrap_rate_of_adaptation>0)])
    ratio = adapt$bootstrap_rate_of_adaptation/mean(adapt2$bootstrap_rate_of_adaptation)
    
    adatation = rbind(adatation, data.frame(s1=adapt$bootstrap_rate_of_adaptation,
                                          s2=adapt2$bootstrap_rate_of_adaptation,
                              virus=gsub("LIKE","like",toupper(tmp2[[1]])), protein="adaptation"))
    
    rec.interval = HPDinterval(as.mcmc(exp(t$relativeRecombinationRate1-t$relativeRecombinationRate2)), prob=0.95)
    rec.interval.8 = HPDinterval(as.mcmc(exp(t$relativeRecombinationRate1-t$relativeRecombinationRate2)), prob=0.8)
    rec.interval.5 = HPDinterval(as.mcmc(exp(t$relativeRecombinationRate1-t$relativeRecombinationRate2)), prob=0.5)
    
    
    interval = HPDinterval(as.mcmc(adapt$bootstrap_rate_of_adaptation/adapt2$bootstrap_rate_of_adaptation), prob=0.95)
    interval.8 = HPDinterval(as.mcmc(adapt$bootstrap_rate_of_adaptation/adapt2$bootstrap_rate_of_adaptation), prob=0.8)
    interval.5 = HPDinterval(as.mcmc(adapt$bootstrap_rate_of_adaptation/adapt2$bootstrap_rate_of_adaptation), prob=0.5)
    
    adapt.ratio = rbind(adapt.ratio, data.frame(median=median(adapt$bootstrap_rate_of_adaptation/adapt2$bootstrap_rate_of_adaptation),
                                                u.95 = interval[1,"upper"], l.95 = interval[1,"lower"],
                                                u.8 = interval.8[1,"upper"], l.8 = interval.8[1,"lower"],
                                                u.5 = interval.5[1,"upper"], l.5 = interval.5[1,"lower"],
                                                median = median(adapt$bootstrap_rate_of_adaptation/adapt2$bootstrap_rate_of_adaptation),
                                                rec.u.95 = rec.interval[1,"upper"], rec.l.95 = rec.interval[1,"lower"],
                                                rec.u.8 = rec.interval.8[1,"upper"], rec.l.8 = rec.interval.8[1,"lower"],
                                                rec.u.5 = rec.interval.5[1,"upper"], rec.l.5 = rec.interval.5[1,"lower"],
                                                rec.median = median(exp(t$relativeRecombinationRate1-t$relativeRecombinationRate2)),
                                                virus = toupper(tmp2[[1]]), protein="adaptation"))
    
    
    
  }
}


order_prots = c("5\' of spike", "spike", "3\' of spike", "s1 over s2")
dat.ratio$genome="full genome"
dat$genome="spike only"
p.sarslike = ggplot(dat.ratio[which(dat.ratio$virus=="SARS-like"),]) +
  geom_violin(aes(x=protein, y=exp(recombination.rate-mean.rate), group=protein)) +
  geom_violin(data=dat[which(dat$virus=="SARS-like"),], aes(x="s1 over s2", y=rate,group=2)) +
  facet_grid(.~genome, scales="free_x", space="free")+
  theme_minimal() +
  coord_cartesian(ylim=c(0.05,20))+
  # scale_x_discrete(limits=order_prots)+
  ylab("recombination rate ratio") +
  scale_y_log10()+
  xlab("")

plot(p.sarslike)
ggsave(plot=p.sarslike, paste("../../Recombination-Text/Figures/SarsLikeRelRates.pdf", sep=""), width=7, height=3)


order_prots = c("5\' of spike", "spike", "3\' of spike", "s1 over s2")
dat.ratio$genome="full genome"
dat$genome="spike only"
p.mers = ggplot(dat.ratio[which(dat.ratio$virus=="MERS"),]) +
  geom_violin(aes(x=protein, y=exp(recombination.rate-mean.rate), group=protein)) +
  geom_violin(data=dat[which(dat$virus=="MERS"),], aes(x="s1 over s2", y=rate,group=2)) +
  facet_grid(.~genome, scales="free_x", space="free")+
  theme_minimal() +
  coord_cartesian(ylim=c(0.05,20))+
  # scale_x_discrete(limits=order_prots)+
  ylab("recombination rate ratio") +
  scale_y_log10()+
  xlab("")

plot(p.mers)
ggsave(plot=p.mers, paste("../../Recombination-Text/Figures/MERSRelRates.pdf", sep=""), width=7, height=3)




dat$virus <- factor(dat$virus, levels = c("229E", "OC43", "NL63"))

adapt.ratio$virus <- factor(adapt.ratio$virus, levels = c("229E", "OC43", "NL63"))

p.rec = ggplot(dat) +
  geom_violin(data=adatation, aes(x=virus, y=s1, color="s1")) +
  geom_violin(data=adatation, aes(x=virus, y=s2, color="s2")) +
  # scale_y_log10() +
  theme_minimal() +
  # facet_wrap(.~protein, ncol=3)+
  ylab("rate on subunit 1\nover rate on subunit 2") +
  xlab("") +
  # scale_y_log10()+
  # coord_cartesian(ylim=c(0.1,100))+
  scale_color_OkabeIto() +
  scale_fill_OkabeIto(guide="none")
plot(p.rec)



for (i in seq(1,length(human.covs))){
  for (j in seq(1,length(u))){
    vals = exp(dat.ratio[which(dat.ratio$virus==human.covs[[i]] & dat.ratio$protein==u[[j]]), "recombination.rate"]-
                 dat.ratio[which(dat.ratio$virus==human.covs[[i]] & dat.ratio$protein==u[[j]]), "mean.rate"])
    interval = HPDinterval(as.mcmc(vals), prob=0.95)
    interval.8 = HPDinterval(as.mcmc(vals), prob=0.8)
    interval.5 = HPDinterval(as.mcmc(vals), prob=0.5)
    
    
    
    adapt.vals[which(adapt.vals$virus==human.covs[[i]] & adapt.vals$protein==u[[j]]), "rec.u.95"] = interval[1,"upper"]
    adapt.vals[which(adapt.vals$virus==human.covs[[i]] & adapt.vals$protein==u[[j]]), "rec.l.95"] = interval[1,"lower"]

    adapt.vals[which(adapt.vals$virus==human.covs[[i]] & adapt.vals$protein==u[[j]]), "rec.u.8"] = interval.8[1,"upper"]
    adapt.vals[which(adapt.vals$virus==human.covs[[i]] & adapt.vals$protein==u[[j]]), "rec.l.8"] = interval.8[1,"lower"]

    adapt.vals[which(adapt.vals$virus==human.covs[[i]] & adapt.vals$protein==u[[j]]), "rec.u.5"] = interval.5[1,"upper"]
    adapt.vals[which(adapt.vals$virus==human.covs[[i]] & adapt.vals$protein==u[[j]]), "rec.l.5"] = interval.5[1,"lower"]
    adapt.vals[which(adapt.vals$virus==human.covs[[i]] & adapt.vals$protein==u[[j]]), "rec.median"] = median(vals)
    
    }
}


# p.rec = ggplot(dat.ratio[which(dat.ratio$virus!="SARS-like" & dat.ratio$virus!="MERS"),]) +
#   geom_linerange(data=adapt.vals, aes(x=protein, ymin=l.5, ymax=u.5, color="adaptation rate ratio",size="50% HPD")) +
#   geom_linerange(data=adapt.vals, aes(x=protein, ymin=l.8, ymax=u.8, color="adaptation rate ratio",size="80% HPD")) +
#   geom_linerange(data=adapt.vals, aes(x=protein, ymin=l.95, ymax=u.95, color="adaptation rate ratio",size="95% HPD")) +
#   geom_violin(aes(x=protein, y=exp(recombination.rate-mean.rate), color="recombination rate ratio"), size=0.6, fill=NA) +
#   scale_size_manual(name="",values=c(5,2,1))+
#   scale_fill_OkabeIto(guide="none")+
#   scale_color_OkabeIto(name="")+
#   theme_minimal() +
#   facet_wrap(.~virus, ncol=1)+
#   coord_cartesian(ylim=c(0.05,20))+
#   ylab("recombination rate") +
#   scale_y_log10()+
#   ylab("rate ratio")+
#   xlab("")
# 
# plot(p.rec)
# ggsave(plot=p.rec, paste("../../Recombination-Text/Figures/RecombinationRateProtein.pdf", sep=""), width=6, height=4)



adapt.ratio$protein ="s1 over s2"
adapt.ratio$genome = "spike only"

dat$protein ="s1 over s2"
dat$genome = "spike only"
 
adapt.vals$genome = "full genome"
dat.ratio$genome = "full genome"

require("ggstance")
p.rec = ggplot() +
  geom_abline(linetype="dashed")+
  geom_point(data=adapt.vals, aes(x=rec.median, y=median, color=paste(protein, "vs rest"))) +
  # geom_linerange(data=adapt.vals, aes(x=rec.median, ymin=l.5, ymax=u.5, color=paste(protein, "vs rest"),size="50% HPD"), alpha=0.8) +
  # geom_linerange(data=adapt.vals, aes(x=rec.median, ymin=l.8, ymax=u.8, color=paste(protein, "vs rest"),size="80% HPD"), alpha=0.8) +
  # geom_linerange(data=adapt.vals, aes(x=rec.median, ymin=l.95, ymax=u.95, color=paste(protein, "vs rest"),size="95% HPD"), alpha=0.8) +
  geom_linerange(data=adapt.vals, aes(x=rec.median, ymin=l.95, ymax=u.95, color=paste(protein, "vs rest"))) +
  # geom_linerangeh(data=adapt.vals, aes(y=median, xmin=rec.l.5, xmax=rec.u.5, color=paste(protein, "vs rest"),size="50% HPD"), alpha=0.8) +
  # geom_linerangeh(data=adapt.vals, aes(y=median, xmin=rec.l.8, xmax=rec.u.8, color=paste(protein, "vs rest"),size="80% HPD"), alpha=0.8) +
  # geom_linerangeh(data=adapt.vals, aes(y=median, xmin=rec.l.95, xmax=rec.u.95, color=paste(protein, "vs rest"),size="95% HPD"), alpha=0.8) +
  geom_linerangeh(data=adapt.vals, aes(y=median, xmin=rec.l.95, xmax=rec.u.95, color=paste(protein, "vs rest"))) +
  
  
  geom_point(data=adapt.ratio, aes(x=rec.median, y=median, color="subunit 1 over subunit 2")) +
  # geom_linerange(data=adapt.ratio, aes(x=rec.median, ymin=l.5, ymax=u.5, color="subunit 1 over subunit 2",size="50% HPD")) +
  # geom_linerange(data=adapt.ratio, aes(x=rec.median, ymin=l.8, ymax=u.8, color="subunit 1 over subunit 2",size="80% HPD")) +
  # geom_linerange(data=adapt.ratio, aes(x=rec.median, ymin=l.95, ymax=u.95, color="subunit 1 over subunit 2",size="95% HPD")) +
  geom_linerange(data=adapt.ratio, aes(x=rec.median, ymin=l.95, ymax=u.95, color="subunit 1 over subunit 2")) +
  # geom_linerangeh(data=adapt.ratio, aes(y=median, xmin=rec.l.5, xmax=rec.u.5, color="subunit 1 over subunit 2",size="50% HPD")) +
  # geom_linerangeh(data=adapt.ratio, aes(y=median, xmin=rec.l.8, xmax=rec.u.8, color="subunit 1 over subunit 2",size="80% HPD")) +
  # geom_linerangeh(data=adapt.ratio, aes(y=median, xmin=rec.l.95, xmax=rec.u.95, color="subunit 1 over subunit 2",size="95% HPD")) +
  geom_linerangeh(data=adapt.ratio, aes(y=median, xmin=rec.l.95, xmax=rec.u.95, color="subunit 1 over subunit 2")) +
  xlab("relative recombination rates")+
  ylab("relative adaptation rates")+
  
  
  # scale_size_manual(name="",values=c(1.5,1,0.5)) +
  scale_y_log10()+
  scale_x_log10() +
  theme_minimal()+
  facet_grid(.~virus) +
  scale_color_OkabeIto(name="section of the genome")

  
  
plot(p.rec)
  
ggsave(plot=p.rec, paste("../../Recombination-Text/Figures/RecombinationRateRatio.pdf", sep=""), width=9, height=3)



p.rec = ggplot(dat.ratio[which(dat.ratio$virus!="SARS-like" & dat.ratio$virus!="MERS"),]) +
  
  geom_linerange(data=adapt.vals, aes(x=protein, ymin=l.5, ymax=u.5, color="adaptation rate ratio",size="50% HPD")) +
  geom_linerange(data=adapt.vals, aes(x=protein, ymin=l.8, ymax=u.8, color="adaptation rate ratio",size="80% HPD")) +
  geom_linerange(data=adapt.vals, aes(x=protein, ymin=l.95, ymax=u.95, color="adaptation rate ratio",size="95% HPD")) +
  
  # geom_linerange(data=adapt.vals, aes(x=protein, ymin=rec.l.5, ymax=rec.u.5, color="recombination rate ratio",size="50% HPD")) +
  # geom_linerange(data=adapt.vals, aes(x=protein, ymin=rec.l.8, ymax=rec.u.8, color="recombination rate ratio",size="80% HPD")) +
  geom_linerange(data=adapt.vals, aes(x=protein, ymin=rec.l.95, ymax=rec.u.95, color="recombination rate ratio",size="95% HPD")) +
  
  
  geom_violin(aes(x=protein, y=exp(recombination.rate-mean.rate), color="recombination rate ratio", fill=NA), size=0.6) +
  geom_linerange(data=adapt.ratio, aes(x=protein, ymin=l.5, ymax=u.5, color="adaptation rate ratio",size="50% HPD")) +
  geom_linerange(data=adapt.ratio, aes(x=protein, ymin=l.8, ymax=u.8, color="adaptation rate ratio",size="80% HPD")) +
  geom_linerange(data=adapt.ratio, aes(x=protein, ymin=l.95, ymax=u.95, color="adaptation rate ratio",size="95% HPD")) +
  
  geom_violin(data=dat[which(dat$virus=="229E"|dat$virus=="OC43" | dat$virus=="NL63"),],
              aes(x=protein, y=rate, color="recombination rate ratio", fill=NA), size=0.6) +

  theme_minimal() +
  scale_size_manual(name="",values=c(5,2,1))+
  facet_grid(virus~genome, scales="free_x",space="free")+
  ylab("rate ratio") +
  xlab("") +
  coord_cartesian(ylim=c(0.05,1000))+
  scale_y_log10()+
  scale_color_discrete_qualitative(name="") +
  scale_fill_brewer(guide="none")
plot(p.rec)

ggsave(plot=p.rec, paste("../../Recombination-Text/Figures/RecombinationRateProtein.pdf", sep=""), width=9, height=6)



dat[which(dat$virus!="SARS-like"),]