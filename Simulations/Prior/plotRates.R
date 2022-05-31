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
log <- list.files(path="./out/", pattern=paste("inf", ".*.log", sep=""), full.names = TRUE)

for (i in seq(1,length(log))){
  fname1 = log[[i]]
  
  # read in the log file
  t.1 <- read.table(fname1, header=TRUE, sep="\t")

  # take a 10% burnin
  t <- t.1[-seq(1,ceiling(length(t.1$posterior)/5)), ]

  # combine the replictes
  
  # calculate ess values
  if (length(t$Sample)>10){
    ess <- effectiveSize(t)
    
    if (min(ess[2:5])>50){
      print(fname1)
      
      
      # find the correct run number
      tmp = strsplit(log[[i]], '_')[[1]][[3]]
      
      
      scale_val = 0.995
      if (grepl( "rep1", log[[i]])){
        prior="low narrow"
        offset=scale_val*scale_val
      }else if (grepl( "rep2", log[[i]])){
        prior="centered"
        offset=1
      }else if (grepl( "rep3", log[[i]])){
        prior="high narrow"
        offset=1/scale_val/scale_val
      }else if (grepl( "rep4", log[[i]])){
        prior="low wide"
        offset=scale_val
      }else if (grepl( "rep5", log[[i]])){
        prior="high wide"
        offset=1/scale_val
        
      }else{
        dsadadsawsa
      }
      
      used.rates = true.rates[which(true.rates$run==as.numeric(tmp)),];
      # combine with the other replicates
      new.rea = data.frame(true=used.rates$recombination*exp(used.rates$relRate1), 
                           estimated=mean(t$recombinationRate), 
                           upper=quantile(t$recombinationRate,0.975),
                           lower=quantile(t$recombinationRate,0.025),
                           prior=prior,offset=offset,est="estimates",run=tmp
      )
      
      # das

      new.Ne = data.frame(true=used.rates$Ne, 
                          estimated=mean(t$popSize),
                          upper=quantile(t$popSize,0.975), 
                          lower=quantile(t$popSize,0.025),
                          prior=prior,offset=offset,run=tmp
      )
      
      used.evol = true.evol[which(true.evol$run==as.numeric(tmp)),];
      j=1
      kappa = paste("kappa",j,sep="")
      mut = paste("mut",j,sep="")
      gamma = paste("gam",j,sep="")
      

      new.m = data.frame(true=used.evol$m, 
                         estimated=mean(t[,mut]), 
                         upper=quantile(t[,mut],0.975),
                         lower=quantile(t[,mut],0.025),
                         prior=prior,offset=offset,run=tmp)
      


      if (first.evol){
        m = new.m
        first.evol = F
      }else{
        m = rbind(m, new.m)
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
      # print("ess low")
    }
  }else{
    # print("not enough its")
  }
}

require(pracma)

val_x = logseq(min(rea$true),max(rea$true),5)


# combine with the other replicates
prior = data.frame(true="low narrow",
                   x=val_x[[1]],
                   vals = rlnorm(n=10000,  -12.7365, 0.5),
                     prior="low narrow"
)
prior = rbind(prior, data.frame(true="low wide",
                   x=val_x[[2]],
                   vals = rlnorm(n=10000, -12.7365, 2),
                   prior="low wide"
))

prior = rbind(prior, data.frame(true="centered",
                   vals = rlnorm(n=10000, -11.1271, 0.5),
                   x=val_x[[3]],
                   prior="centered"
))

prior = rbind(prior, data.frame(true="high wide",
                                vals = rlnorm(n=10000,-9.5177, 2),
                                x=val_x[[4]],
                                prior="high wide"
))

prior = rbind(prior, data.frame(true="high narrow",
                                vals = rlnorm(n=10000,-9.5177, 0.5),
                                x=val_x[[5]],
                                prior="high narrow"
))


diff_priors = c("low narrow"="#A6CEE3", "low wide"="#1F78B4", "centered"="#B2DF8A", "high wide"="#33A02C", "high narrow"="#FB9A99")

p.Ne <- ggplot(Ne)+
  geom_abline(intercept = 0, color="red")+
  geom_errorbar(aes(x=true*offset, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true*offset, y=estimated, color=prior), size=2) + 
  geom_line(aes(x=true*offset, y=estimated, group=true))+
  theme_minimal() +
  scale_y_log10() +
  scale_x_log10() +
  xlab("true effective population size")+
  ylab("estimated effective population size")

p.rea = list()
for (i in seq(1,length(diff_priors))){
  p.rea[[i]] <- ggplot(rea[which(rea$prior==labels(diff_priors)[[i]]),])+
    geom_abline(intercept = 0, color="red", linetype="dashed")+
    geom_errorbar(aes(x=true, ymin=lower, ymax=upper),color=diff_priors[[i]], width=0.01) +
    geom_point(aes(x=true, y=estimated),color=diff_priors[[i]], size=2) + 
    geom_line(aes(x=true, y=estimated),color=diff_priors[[i]])+
    theme_minimal() +
    scale_y_log10() +
    scale_x_log10() +
    # scale_color_manual(name="prior distributions\nused for inference", values=diff_priors)+
    # scale_fill_manual(name="prior distributions",values=diff_priors)+
    coord_cartesian(ylim=c(0.0000005,0.001))+
    ggtitle(labels(diff_priors)[[i]])+
    xlab("true recombination rate")+
    ylab("estimated recombination rate")
  plot(p.rea[[i]])
}

p.rea[[length(diff_priors)+1]] <- ggplot(rea)+
  geom_violin(data=prior,aes(x=prior, y=vals, fill=prior),color=NA) +
  theme_minimal() +
  scale_y_log10() +
  scale_color_manual(name="prior distributions\nused for inference", values=diff_priors)+
  scale_fill_manual(name="prior distributions",values=diff_priors)+
  coord_cartesian(ylim=c(0.0000005,0.001))+
  xlab(" ")+
  ylab("recombination rate") +
  ggtitle("log-normal prior distributions")+
  theme(legend.position = "none") + theme(axis.text.x=element_text(angle=30, hjust=1))

plot(p.rea[[length(diff_priors)+1]] )

library("ggpubr")

p_combined1 = ggarrange(p.rea[[1]],p.rea[[2]],p.rea[[3]], ncol=3, labels=c("A","B","C"))
p_combined2 = ggarrange(p.rea[[4]],p.rea[[5]],p.rea[[6]], ncol=3, labels=c("D","E","F"))
p_combined  = ggarrange(p_combined1,p_combined2,ncol=1)
plot(p_combined)

p.mut <- ggplot(m)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true*offset, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true*offset, y=estimated, color=prior), size=2) + 
  geom_line(aes(x=true*offset, y=estimated, group=true))+
  theme_minimal() +
  scale_y_log10() +
  scale_x_log10() +
  xlab("true evolutionary rate")+
  ylab("estimated evolutionary rate")



plot(p.Ne)

plot(p.mut)

ggsave(plot=p.Ne,paste("../../../Recombination-Text/Figures/simulations/prior_Ne.pdf", sep=""),width=4, height=4)
ggsave(plot=p_combined,paste("../../../Recombination-Text/Figures/prior_rea.pdf", sep=""),width=9, height=6)
# ggsave(plot=p.kappa,paste("../../../Recombination-Text/Figures/simulations/hky_kappa.pdf", sep=""),width=4, height=4)
ggsave(plot=p.mut,paste("../../../Recombination-Text/Figures/simulations/prior_mut.pdf", sep=""),width=4, height=4)
# ggsave(plot=p.gam,paste("../../../Recombination-Text/Figures/simulations/hky_gam.pdf", sep=""),width=4, height=4)
# ggsave(plot=p.freq,paste("../../../Recombination-Text/Figures/simulations/hky_freq.pdf", sep=""),width=4, height=4)
