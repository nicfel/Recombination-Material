######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(coda)


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# read in the true rates
true.rates <- read.table("rates.csv", header=TRUE, sep=",")
true.evol <- read.table("rates_sequences.csv", header=TRUE, sep=",")

first = T
first.evol = T
first.freq = T


# get the names of all output files of the first replicate
log <- list.files(path="./out/", pattern=paste("inf", ".*rep1.*..log", sep=""), full.names = TRUE)

for (i in seq(1,length(log))){
  fname1 = log[[i]]
  print(fname1)
  
  # read in the log file
  t.1 <- read.table(fname1, header=TRUE, sep="\t")
  t.2 <- read.table(gsub("rep1","rep2",fname1), header=TRUE, sep="\t")
  t.3 <- read.table(gsub("rep1","rep3",fname1), header=TRUE, sep="\t")
  
  # take a 10% burnin
  t.1 <- t.1[-seq(1,ceiling(length(t.1$posterior)/5)), ]
  t.2 <- t.3[-seq(1,ceiling(length(t.2$posterior)/5)), ]
  t.3 <- t.3[-seq(1,ceiling(length(t.3$posterior)/5)), ]

  
  # combine the replictes
  # t = t.1
  t = rbind(t.1,t.2,t.3)

  # calculate ess values
  if (length(t$Sample)>10){
    ess <- effectiveSize(t)
  
    if (min(ess[2:5])>30){
      
      # find the correct run number
      tmp = strsplit(log[[i]], '_')[[1]][[3]]
      lam = as.numeric(strsplit(strsplit(log[[i]], '_')[[1]][[4]], split="\\.")[[1]][[2]])
      
      used.rates = true.rates[which(true.rates$run==as.numeric(tmp) & true.rates$lambda==lam),];
    
      # combine with the other replicates
      new.rea = data.frame(true=used.rates$recombination*lam*exp(used.rates$relRate1), 
                          estimated=mean(t$recombinationRate*exp(t$relativeRecombinationRate1)), 
                          upper=quantile(t$recombinationRate*exp(t$relativeRecombinationRate1),0.975),
                          lower=quantile(t$recombinationRate*exp(t$relativeRecombinationRate1),0.025),
                          mean = as.character(lam)
                          )
      
      new.rea = rbind(new.rea,data.frame(true=used.rates$recombination*lam*exp(used.rates$relRate2), 
                           estimated=mean(t$recombinationRate*exp(t$relativeRecombinationRate2)), 
                           upper=quantile(t$recombinationRate*exp(t$relativeRecombinationRate2),0.975),
                           lower=quantile(t$recombinationRate*exp(t$relativeRecombinationRate2),0.025),
                           mean = as.character(lam)
      ))
      
      new.rea = rbind(new.rea,data.frame(true=used.rates$recombination*lam*exp(used.rates$relRate3), 
                                         estimated=mean(t$recombinationRate*exp(t$relativeRecombinationRate3)), 
                                         upper=quantile(t$recombinationRate*exp(t$relativeRecombinationRate3),0.975),
                                         lower=quantile(t$recombinationRate*exp(t$relativeRecombinationRate3),0.025),
                                         mean = as.character(lam)
      ))
      
      
      new.Ne = data.frame(true=used.rates$Ne, 
                          estimated=mean(t$popSize),
                          upper=quantile(t$popSize,0.975), 
                          lower=quantile(t$popSize,0.025),
                          mean = as.character(lam)
                          )
      
      used.evol = true.evol[which(true.evol$run==as.numeric(tmp) & true.rates$lambda==lam),];
      j=1
      kappa = paste("kappa",j,sep="")
      mut = paste("mut",j,sep="")
      gamma = paste("gam",j,sep="")
      
      new.kappa = data.frame(true=used.evol$k, 
                           estimated=mean(t[,kappa]), 
                           upper=quantile(t[,kappa],0.975),
                           lower=quantile(t[,kappa],0.025),
                           mean = as.character(lam))
      
      new.m = data.frame(true=used.evol$m, 
                         estimated=mean(t[,mut]), 
                         upper=quantile(t[,mut],0.975),
                         lower=quantile(t[,mut],0.025),
                         mean = as.character(lam))
      
      new.g = data.frame(true=used.evol$g, 
                         estimated=mean(t[,gamma]), 
                         upper=quantile(t[,gamma],0.975),
                         lower=quantile(t[,gamma],0.025),
                         mean = as.character(lam))
      
      for (k in seq(1,4)){
        freq = paste("freq",j,k,sep="")
        freq.name = paste("f",k,sep="")
        new.f = data.frame(true=used.evol[,freq.name], 
                           estimated=mean(t[,freq]), 
                           upper=quantile(t[,freq],0.975),
                           lower=quantile(t[,freq],0.025),
                           mean = as.character(lam))
        
        if (first.freq){
          f = new.f
          first.freq = F
        }else{
          f = rbind(f, new.f)
        }
        
      }
      
      if (first.evol){
        kap = new.kappa
        m = new.m
        g = new.g
        first.evol = F
      }else{
        kap = rbind(kap, new.kappa)
        m = rbind(m, new.m)
        g = rbind(g, new.g)
      }
  
      
      
      
      if (first){
        rea = new.rea
        Ne = new.Ne
        first = F
      }else{
        rea = rbind(rea, new.rea)
        Ne = rbind(Ne, new.Ne)
      }
    }else{
      print("ess low")
    }
  }
}




p.Ne <- ggplot(Ne)+
  geom_abline(intercept = 0, color="red")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper,color=mean), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated,color=mean), size=2) + 
  theme_minimal() +
  scale_y_log10() +
  scale_x_log10() +
  xlab("true effective population size")+
  ylab("estimated effective population size")


p.rea <- ggplot(rea)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper, color=mean), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated, color=mean), size=2) + 
  theme_minimal() +
  scale_y_log10() +
  scale_x_log10() +
  xlab("true recombination rate * mean of Poisson distribution")+
  ylab("estimated recombination rate")

p.kappa <- ggplot(kap)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) + 
  theme_minimal() +
  scale_y_log10() +
  scale_x_log10() +
  xlab("true transition/transversion rate")+
  ylab("estimated transition/transversion rate")


p.mut <- ggplot(m)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) + 
  theme_minimal() +
  scale_y_log10() +
  scale_x_log10() +
  xlab("true evolutionary rate")+
  ylab("estimated evolutionary rate")

p.gam <- ggplot(g)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) + 
  theme_minimal() +
  scale_y_log10() +
  scale_x_log10() +
  xlab("true gamma")+
  ylab("estimated gamma")

p.freq <- ggplot(f) +
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) + 
  theme_minimal() +
  xlab("true nucleotide frequency")+
  ylab("estimated nucleotide frequency")


plot(p.Ne)
plot(p.rea)
plot(p.kappa)
plot(p.mut)
plot(p.gam)
plot(p.freq)

# ggsave(plot=p.Ne,paste("../../../Recombination-Text/Figures/simulations/hky_Ne.pdf", sep=""),width=4, height=4)
# ggsave(plot=p.rea,paste("../../../Recombination-Text/Figures/simulations/hky_rea.pdf", sep=""),width=4, height=4)
# ggsave(plot=p.kappa,paste("../../../Recombination-Text/Figures/simulations/hky_kappa.pdf", sep=""),width=4, height=4)
# ggsave(plot=p.mut,paste("../../../Recombination-Text/Figures/simulations/hky_mut.pdf", sep=""),width=4, height=4)
# ggsave(plot=p.gam,paste("../../../Recombination-Text/Figures/simulations/hky_gam.pdf", sep=""),width=4, height=4)
# ggsave(plot=p.freq,paste("../../../Recombination-Text/Figures/simulations/hky_freq.pdf", sep=""),width=4, height=4)
