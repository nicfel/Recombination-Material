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

log <- list.files(path=paste("./out/",sep=""), pattern="*.log", full.names = TRUE)

dat = data.frame()
for (i in seq(1,length(log))){
  t = read.table(log[[i]], header=T,sep="\t")
  t = t[-seq(1,length(t$Sample)/10),]
  
  name = strsplit(log[[i]], split="_")[[1]][[2]]
  name = strsplit(name, split="\\.")[[1]]
  
  ess = effectiveSize(as.mcmc(t))
  
  print(ess[length(ess)-11])
  
  dat = rbind(dat, data.frame(run=name[[2]], type=name[[1]], ess=ess[2], values="posterior"))
  dat = rbind(dat, data.frame(run=name[[2]], type=name[[1]], ess=ess[3], values="likelihood"))
  dat = rbind(dat, data.frame(run=name[[2]], type=name[[1]], ess=ess[5], values="height"))
  dat = rbind(dat, data.frame(run=name[[2]], type=name[[1]], ess=ess[length(ess)-11], values="effective population size"))
}

p = ggplot(dat, aes(x=type)) +
  geom_point(aes(y=ess)) +
  theme_minimal()+
  facet_wrap(.~values) +
  xlab("")+
  ylab("effective sample size") +
  coord_cartesian(ylim=c(0,2500))

plot(p)

ggsave(plot=p, paste("../../../Recombination-Text/Figures/ESS.pdf", sep=""), height=3, width=6)
