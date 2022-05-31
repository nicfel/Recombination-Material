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
log <- list.files(path="./out/", pattern=paste("inf", ".*rep1.log", sep=""), full.names = TRUE)

ess.dat = data.frame()

for (i in seq(1,length(log))){
  fname1 = log[[i]]
  print(fname1)
  
  # read in the log file
  t.1 <- read.table(fname1, header=TRUE, sep="\t")

  # take a 10% burnin
  t <- t.1[-seq(1,ceiling(length(t.1$posterior)/2)), ]
  
  time = c()

  fileName <- gsub(".log", ".out", fname1)
  conn <- file(fileName,open="r")
  linn <-readLines(conn)
  for (j in 1:length(linn)){
    tmp=strsplit(linn[j], split="\\s+")[[1]]
    if (length(tmp)>1){
      if (grepl("Msamples", tmp[length(tmp)])){
        tmp3 = strsplit(gsub("/Msamples","",tmp[length(tmp)]), split="[hms]")[[1]]
        if (length(tmp3)==3){
          time = c(time, as.numeric(tmp3[[1]]) + as.numeric(tmp3[[2]])/60 + as.numeric(tmp3[[1]])/3600)
        }else if (length(tmp3)==2){
          time = c(time, as.numeric(tmp3[[2]])/60 + as.numeric(tmp3[[1]])/3600)
        }else{
          time = c(time, as.numeric(tmp3[[1]])/3600)
        }
      }
    }
  }
  close(conn)

  # calculate ess values
  ess <- effectiveSize(t)
    
  tmp = strsplit(log[[i]], '_')[[1]][[3]]
  used.rates = true.rates[which(true.rates$run==as.numeric(tmp)),];
  
  iterations = (max(t$Sample)-min(t$Sample))/1000000
  
  ess.dat = rbind(ess.dat, 
                  data.frame(samples.name=paste(used.rates$nr_samples, "samples"), 
                             samples=used.rates$nr_samples, 
                             ess.post = ess[2], 
                             ess.likelihood = ess[3], 
                             iterations=iterations,
                             mean.rea = mean(t$network.recombinationNodeCount),
                             name = fname1,
                             clock=mean(t$mut1),
                             recomb=mean(t$recombinationRate),
                             mstime=mean(time)))
}





# 
# p.ess <- ggplot(ess.dat)+
#   # geom_abline(intercept = 0, color="red")+
#   # geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
#   geom_point(aes(x=mean.rea, y=ess.post/iterations, color=samples.name), size=2) + 
#   geom_smooth(aes(x=mean.rea, y=ess.post/iterations, group=samples.name, color=samples.name), method='lm', formula= y~x)+
#   
#   theme_minimal() +
#   scale_y_log10() +
#   scale_x_log10() +
#   xlab("number of recombination events")+
#   ylab("ESS")
# plot(p.ess)
# 
# p.ess <- ggplot(ess.dat)+
#   # geom_abline(intercept = 0, color="red")+
#   # geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
#   geom_boxplot(aes(x=samples, y=ess.post/iterations, color=mean.rea, group=samples)) + 
#   geom_point(aes(x=samples, y=ess.post/iterations, color=mean.rea), size=2) + 
#   # geom_smooth(aes(x=samples, y=ess.post))+
# 
#   theme_minimal() +
#   scale_y_log10() +
#   scale_x_log10() +
#   xlab("number of samples")+
#   ylab("ESS")
# plot(p.ess)

ess.dat.tmp = ess.dat[which(ess.dat$ess.post>40 & ess.dat$iterations>4.5),]

max=300
p.ess.rea <- ggplot(ess.dat.tmp)+
  # geom_abline(intercept = 0, color="red")+
  # geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_smooth(aes(x=mean.rea, y=ess.post/iterations/mstime), method='lm', formula= y~x)+
  geom_point(aes(x=mean.rea, y=ess.post/iterations/mstime))+
  theme_minimal() +
  scale_y_log10(limits=c(30,1000)) +
  xlim(c(1,max)) +
  xlab("number of recombination events")+
  ylab("ESS per hour")
plot(p.ess.rea)

p.ess.samp <- ggplot(ess.dat.tmp)+
  # geom_abline(intercept = 0, color="red")+
  # geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_smooth(aes(x=samples, y=ess.post/iterations/mstime), method='lm', formula= y~x)+
  geom_point(aes(x=samples, y=ess.post/iterations/mstime))+
  theme_minimal() +
  scale_y_log10(limits=c(30,1000)) +
  xlim(c(1,max)) +
  xlab("number of samples")+
  ylab("ESS per hour")
plot(p.ess.samp)

p.ess.clock <- ggplot(ess.dat.tmp)+
  # geom_abline(intercept = 0, color="red")+
  # geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_smooth(aes(x=recomb, y=ess.post/iterations/mstime), method='lm', formula= y~x)+
  geom_point(aes(x=recomb, y=ess.post/iterations/mstime))+
  theme_minimal() +
  scale_y_log10(limits=c(5,100)) +
  scale_y_log10(limits=c(30,1000)) +
  xlab("clock rate")+
  ylab("ESS per hour")
plot(p.ess.clock)

p.ess.tot <- ggplot(ess.dat.tmp)+
  # geom_abline(intercept = 0, color="red")+
  # geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_smooth(aes(x=mean.rea+samples, y=ess.post/iterations/mstime), method='lm', formula= y~x)+
  geom_point(aes(x=mean.rea+samples, y=ess.post/iterations/mstime))+
  theme_minimal() +
  scale_y_log10(limits=c(30,1000)) +
  xlim(c(1,max)) +
  xlab("number of recombination + sampling events")+
  ylab("ESS per hour")
plot(p.ess.tot)

require("gridExtra")
require("ggpubr")
p_ess = ggarrange(p.ess.rea,p.ess.samp,p.ess.clock,
                   ncol = 3)
# p_ess2 = ggarrange(p.ess.tot,p.ess.clock,
#                   ncol = 2)
# p_ess = ggarrange(p_ess1,p_ess2,
#                   ncol = 1)
plot(p_ess)
summary(glm(log(ess.post) ~ samples+mean.rea+1, family = poisson, data=ess.dat.tmp))



ggsave(plot=p_ess,paste("../../../Recombination-Text/Figures/ess_sims.pdf", sep=""),width=8, height=3)
