library(ggplot2)
library(coda)
library(seewave)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# get the names of all output files of the first replicate
log <- list.files(path="./out/", pattern=paste("test", ".*log", sep=""), full.names = TRUE)

sims <- read.table("./out/simulate_serial.log", header=TRUE, sep="\t")
dat = data.frame()
for (i in seq(1,length(log))){
  print(log[[i]])
  t <- read.table(log[[i]], header=TRUE, sep="\t")
  t = t[-seq(1,length(t$Sample)/10),]
  lab =  labels(sims)[[2]]
  for (k in seq(2,7)){
    m = mean(sims[, lab[[k]]])
    v = var(sims[, lab[[k]]])
    y = vector(,0)
    y2 = vector(,0)
    for (j in seq(1000,length(t$Sample),1000)){
      # y = append(y, ks.test(sims[, lab[[k]]], t[seq(1,j),lab[[k]]])[1]$statistic)
      y = append(y, mean(t[seq(1,j),lab[[k]]])/m)
      y2 = append(y2, var(t[seq(1,j),lab[[k]]])/v)
    }
    label = gsub("./out//test", "", log[[i]])
    label = gsub(".log", "", label)
    tmp = strsplit(label, split="\\.")
    dat = rbind(dat,data.frame(x=seq(1000,length(t$Sample), 1000), mean=y, var=y2, label=tmp[[1]][[1]], group=tmp[[1]][[2]], summary=gsub("network.","",lab[[k]])))
  }
}

dat$summary = as.character(dat$summary)
dat[!grepl("obs",dat$summary),"summary"] = paste("above root\n",dat[!grepl("obs",dat$summary),"summary"], sep="")


dat$summary = gsub("obs","below root\n", dat$summary)
dat$summary = gsub("recombinationNodeCount","recombination\nNode Count", dat$summary)
dat$summary = gsub("ReassortmentNodeCount","recombination\nNode Count", dat$summary)
dat$summary = gsub("totalLength","Length", dat$summary)
dat$summary = gsub("TotalLength","Length", dat$summary)

p = ggplot(data=dat) +
  geom_line(aes(x=x,y=mean, group=interaction(group,summary))) +
  # scale_y_log10()+
  scale_x_log10()+
  geom_hline(yintercept=1, linetype=2) +
  # theme_minimal()+
  facet_grid(summary~label)+xlab("") +ylab("mean sampled/mean simulated")
  
  
plot(p)
ggsave(plot=p, paste("../../Recombination-Text/Figures/Validation.pdf", sep=""), width=15, height=10)



