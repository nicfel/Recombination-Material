library(ggplot2)
library(coda)
library("colorblindr")
library(rjson)


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


logs <- list.files(path="./combined/", pattern=paste("*all*.length", sep=""), full.names = TRUE)

dat = data.frame()

for (i in seq(1, length(logs))){
  if (!grepl("sars2",logs[[i]])){
    t = read.table(logs[[i]], header=T, sep="\t")
    tmp = strsplit(logs[[i]], split="[/\\.]")[[1]]
    seq.length = 29000#l.dat[which(l.dat$virus==tmp[[length(tmp)-1]]), "length"]
    virus = toupper(gsub("_all","",tmp[[length(tmp)-1]]))
    virus = gsub("LIKE", "like", virus)
    

      dat = rbind(dat, data.frame(events=t$events,
                                  segmentLength=t$mean, 
                                virus=virus, 
                                stat="average segment length"))

    
    
  }
}

dat$virus <- factor(dat$virus, levels = c("SARS-like","MERS", "229E", "OC43", "NL63"))
# 
p1 = ggplot(dat) +
  geom_violin(aes(x=virus, y=segmentLength)) +
  scale_y_log10() +
  theme_minimal() +
  # ggtitle("segment length without recombination event") +
  ylab("number of basepairs")+
  xlab("")

p2 = ggplot(dat) +
  geom_violin(aes(x=virus, y=events)) +
  scale_y_log10() +
  theme_minimal() +
  # ggtitle("Number of recombination events that affected a sampled genome") +
  ylab("number of recombination events")+
  xlab("")

require(ggpubr)
p.rec = ggarrange(p2,p1, labels=c("A", "B"), ncol=1)
plot(p.rec)
ggsave(plot=p.rec, paste("../../Recombination-Text/Figures/RecombinationLength.pdf", sep=""), height=6, width=9)





logs <- list.files(path="./combined/", pattern=paste("*all\\..*.trunk", sep=""), full.names = TRUE)

dat.rate = data.frame()

human.covs = c("229E", "OC43", "NL63", "MERS")

for (i in seq(1, length(logs))){
  if (!grepl("sars",logs[[i]])){
    t = read.table(logs[[i]], header=F, sep="\t")
    tmp = strsplit(logs[[i]], split="[/\\.]")[[1]]
    virus = gsub("_all","",tmp[[length(tmp)-2]])
    virus = gsub("_subunits","",virus)
    years = paste(tmp[[length(tmp)-1]], "years")
    
    # dat.rate = rbind(dat.rate, data.frame(rate=t$V1/t$V3-t$V2/t$V4,
    #                                       virus=toupper(virus),
    #                                       stat="diff",
    #                                       years=years,
    #                                       method="Inferred"))
    
    dat.rate = rbind(dat.rate, data.frame(rate=t$V1/t$V3,
                                virus=toupper(virus),
                                stat="fit",
                                years=years,
                                method="Inferred"))
    dat.rate = rbind(dat.rate, data.frame(rate=t$V2/t$V4,
                                virus=toupper(virus),
                                stat="unfit",
                                years=years,
                                method="Inferred"))
    
    
  }
}

logs <- list.files(path="./combined/", pattern=paste("*all\\..*.trunksim", sep=""), full.names = TRUE)
for (i in seq(1, length(logs))){
  if (!grepl("sars",logs[[i]])){
    t = read.table(logs[[i]], header=F, sep="\t")
    tmp = strsplit(logs[[i]], split="[/\\.]")[[1]]
    virus = gsub("_all","",tmp[[length(tmp)-2]])
    virus = gsub("_subunits","",virus)
    years = paste(tmp[[length(tmp)-1]], "years")

    # dat.rate = rbind(dat.rate, data.frame(rate=t$V1/t$V3-t$V2/t$V4,
    #                                       virus=toupper(virus),
    #                                       stat="diff",
    #                                       years=years,
    #                                       method="Posterior predictive simulations"))
    
    
    dat.rate = rbind(dat.rate, data.frame(rate=t$V1/t$V3,
                                          virus=toupper(virus),
                                          stat="fit",
                                          years=years,
                                          method="Posterior predictive simulations"))
    dat.rate = rbind(dat.rate, data.frame(rate=t$V2/t$V4,
                                          virus=toupper(virus),
                                          stat="unfit",
                                          years=years,
                                          method="Posterior predictive simulations"))
    
  }
}

dat.rate$years <- factor(dat.rate$years, levels = c("1 years","2 years", "5 years", "10 years"))
# human.covs = c("mers")
p.fit = ggplot(dat.rate[is.element(dat.rate$virus,human.covs), ]) +
  geom_violin(aes(x=virus, y=rate, color=stat, fill=stat)) +
  scale_y_log10() +
  theme_minimal() +
  facet_grid(method~years)+
  ylab("recombination rate") +
  xlab("")

plot(p.fit)
ggsave(plot=p.fit, paste("../../Recombination-Text/Figures/RecombinationFitness.pdf", sep=""), height=6, width=9)


# dat.rate$years <- factor(dat.rate$years, levels = c("1 years","2 years", "5 years", "10 years"))
# # human.covs = c("mers")
# p.fit = ggplot(dat.rate[is.element(dat.rate$virus,human.covs), ]) +
#   geom_violin(aes(x=virus, y=rate)) +
#   # scale_y_log10() +
#   theme_minimal() +
#   facet_grid(method.~years)+
#   # facet_grid(method~years)+
#   ylab("recombination rate") +
#   xlab("")
# 
# plot(p.fit)
# ggsave(plot=p.fit, paste("../../Recombination-Text/Figures/RecombinationFitness.pdf", sep=""), height=6, width=9)
