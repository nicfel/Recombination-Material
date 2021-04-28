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


# get the names of all output files of the first replicate
# trees <- list.files(path="./out/", pattern=paste(".*rep1.*.trees", sep=""), full.names = TRUE)
# for (i in seq(1, length(trees))){
#   in_command <- " -b 80 -log"
#   for (j in seq(1,3)){
#     in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), trees[i]), sep="")
#   }
# 
#   out_command = gsub("_rep1", "", trees[i])
#   out_command = gsub("out", "combined", out_command)
# 
#   combined_command = gsub(".trees",".trees", out_command)
#   combined_command = paste(" -o ", gsub("_rep1", "",combined_command), sep="")
# 
#   # combine the trees
#   system(paste("/Applications/BEAST\\ 2.6.3/bin/logcombiner", in_command, combined_command, sep=" "))
# }

# get the names of all output files of the first replicate
# trees <- list.files(path="./out/", pattern=paste("*rep1.*.log", sep=""), full.names = TRUE)
# for (i in seq(1, length(trees))){
#   in_command <- " -b 20 -log"
#   for (j in seq(1,3)){
#     in_command = paste(in_command, " ", gsub("rep1", paste("rep", j,sep=""), trees[i]), sep="")
#   }
# 
#   out_command = gsub("_rep1", "", trees[i])
#   out_command = gsub("out", "combined", out_command)
# 
#   combined_command = gsub(".log",".log", out_command)
#   combined_command = paste(" -o ", gsub("_rep1", "",combined_command), sep="")
# 
#   # combine the trees
#   system(paste("/Applications/BEAST\\ 2.6.3/bin/logcombiner", in_command, combined_command, sep=" "))
# }


# get the names of all output files of the first replicate
# trees <- list.files(path="./combined/", pattern=paste("*.trees", sep=""), full.names = TRUE)
# for (i in seq(1, length(trees))){
#   system(paste("java -jar ../../Software/RateAfterRecombination.jar -burnin 0 -maxDist 100", trees[i], gsub("*.trees", ".recrate", trees[i]), sep=" "))
# }


# get the names of all output files of the first replicate
log <- list.files(path="./combined/", pattern="*recrate", full.names = TRUE)
# log.files <- list.files(path="./combined/", pattern="*log", full.names = TRUE)

axis_labels = c()
axis_pos = c()

rec.rates = data.frame()
require(coda)
for (i in seq(1,length(log))){
  # read in the log file
  t <- read.table(log[[i]], header=F, sep="\t")
  # t.log <- read.table(log.files[[i]], header=T, sep="\t")
  runnr = gsub(".recrate","",strsplit(log[[i]], split="_")[[1]][[3]])
  lambda = strsplit(log[[i]], split="\\.")[[1]][[3]]
  
  vals = (t$V4/t$V5)/(t$V2/t$V3)
  vals = vals[!is.na(vals)]
  i = HPDinterval(as.mcmc(vals))
  
  runval = as.numeric(strsplit(runnr, split="\\.")[[1]][[1]])
  if (grepl("1", lambda)){
    xval = 12 + runval
  }else{
    xval = 24 + runval
  }
  
  axis_labels = c(axis_labels, paste("replicate", runval))
  axis_pos = c(axis_pos, xval)
  
  

  rec.rates = rbind(rec.rates, data.frame(run=runnr,
                                          rate=vals,
                                          l=i[1,"lower"],
                                          u=i[1,"upper"],
                                          method=paste("lambda=",lambda,sep=""),
                                          xval=xval))

}


# log <- list.files(path="./../../Applications/combined/", pattern=paste(".*\\.recrate", sep=""), full.names = TRUE)
log <- list.files(path="./../../Applications/combined/", pattern=paste("*.subunits.*recrate", sep=""), full.names = TRUE)
xval1 = 0
xval2 = 7
for (i in seq(1,length(log))){
  if (!grepl("sars2", log[[i]])){
    t = read.table(log[i], header=F, sep="\t")
    # t.log = read.table(gsub("recrate","log",log[i]), header=F, sep="\t")
  
  
    tmp = strsplit(log[[i]], split="[/\\.]")[[1]]
    virus = strsplit(tmp[[length(tmp)-1]], split="_")
    sim = "inferred"
    if (grepl("sim",log[i])){
      sim = "posterior predictive\nsimulations"
      xval2 = xval2+1
      xval = xval2
    }else{
      xval1 = xval1+1
      xval = xval1
    }
    axis_labels = c(axis_labels, toupper(virus[[1]][[1]]))
    axis_pos = c(axis_pos, xval)
    
    if (grepl("all",virus[[1]][[2]])){
      dataset="full genome"
    }else{
      dataset="spike only"
    }
    
    runnr = paste(virus[[1]][[1]], dataset, sim, sep="\n")
    
    vals = (t$V4/t$V5)/(t$V2/t$V3)
    vals = vals[!is.na(vals)]
    i = HPDinterval(as.mcmc(vals), prob=0.8)
    
    
    rec.rates = rbind(rec.rates, data.frame(run=runnr,
                                            rate=vals,
                                            l=i[1,"lower"],
                                            u=i[1,"upper"],
                                            method=sim,
                                            xval=xval))
  }
}

rec.rates$method <- factor(rec.rates$method, levels = c("inferred","posterior predictive\nsimulations", "lambda=1", "lambda=2"))

p1 = ggplot(rec.rates[which(rec.rates$method=="inferred"),]) +
  geom_boxplot(aes(x=xval,group=run, y=rate), color="black") +
  scale_y_log10() +
  theme_minimal()+
  ylab("rate of recombination\nafter recombination event\nover after coalescent event")+
  xlab("") +
  scale_x_continuous(breaks = axis_pos, labels=axis_labels) +
  ggtitle("Rates from inferred\nrecombination networks\nfor different coronaviruses") +
  coord_cartesian(ylim=c(0.3,30)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(p1)

p2 = ggplot(rec.rates[which(rec.rates$method=="posterior predictive\nsimulations"),]) +
  geom_boxplot(aes(x=xval,group=run, y=rate), color="black") +
  scale_y_log10() +
  theme_minimal()+
  ylab("")+
  xlab("") +
  scale_x_continuous(breaks = axis_pos, labels=axis_labels) +
  ggtitle("Rates from posterior predictive\nrecombination networks\nfor different coronaviruses") +
  coord_cartesian(ylim=c(0.3,30))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(p2)

p3 = ggplot(rec.rates[which(rec.rates$method=="lambda=1"),]) +
  geom_boxplot(aes(x=xval,group=run, y=rate), color="black") +
  scale_y_log10() +
  theme_minimal()+
  ylab("")+
  xlab("") +
  scale_x_continuous(breaks = axis_pos, labels=axis_labels) +
  ggtitle("Rates from inferred\nrecombination networks\nfor simulated data using a mean = 1") +
  coord_cartesian(ylim=c(0.3,30))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(p3)

p4 = ggplot(rec.rates[which(rec.rates$method=="lambda=2"),]) +
  geom_boxplot(aes(x=xval,group=run, y=rate), color="black") +
  scale_y_log10() +
  theme_minimal()+
  ylab("")+
  xlab("") +
  scale_x_continuous(breaks = axis_pos, labels=axis_labels) +
  ggtitle("Rates from inferred\nrecombination networks\nfor simulated data using a mean = 2") +
  coord_cartesian(ylim=c(0.3,30))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(p4)

require(ggpubr)
p = ggarrange(p1+ylab(""),p2,p3,p4,labels=c("A", "B", "C", "D"), ncol=4)  
plot(p)
ggsave(paste("../../../Recombination-Text/Figures/RecombinationRateDiffSimulations.pdf", sep=""), plot=p, width=15, height=4)






# p = ggplot(rec.rates) +
#   # geom_point(aes(x=run, y=rate)) +
#   # geom_linerange(aes(x=run, ymin=l,ymax=u, color=method)) +
#   geom_boxplot(aes(x=xval,group=run, y=rate, color=method, fill=method)) +
#   scale_y_log10() +
#   theme_minimal()+
#   facet_wrap(.~method, strip.position = "bottom", scales = "free_x",  nrow=1)+
#   # theme(legend.title = element_blank())+
#   theme(panel.spacing = unit(0, "lines"),
#         # strip.background = element_blank(),
#         # axis.line = element_line(colour = "grey"),
#         # panel.grid.major.y =element_line(colour = "grey"),
#         strip.placement = "outside",
#         axis.text.x = element_text(angle = 90, hjust = 1),
#         # panel.background = element_rect(fill = 'white', colour = 'white'),
#         legend.position = "none"
#   )+
#   ylab("rate of recombination\nafter recombination event\nover after coalescent event")+
#   xlab("") +
#   scale_x_continuous(breaks = axis_pos, labels=axis_labels)

# p = ggplot(rec.rates) +
#   geom_count(aes(x=total, y=after.recombination/after.coalescence, color=run)) +
#   # geom_violin(aes(x=run,group=run, y=after.recombination/(after.coalescence))) +
#   scale_y_log10() +
#   theme_minimal()+
#   theme(legend.title = element_blank())+
#   ylab("rate of recombination\nafter recombination event\nover after coalescent event")+
#   xlab("")
# plot(p)


