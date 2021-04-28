library(ggplot2)
library(coda)


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

trunk.logs <- list.files(path="./combined/", pattern=paste("*pruned.trunk.txt", sep=""), full.names = TRUE)
dat = data.frame()
diff = data.frame()
for (i in seq(1, length(trunk.logs)-2)){
  t = read.table(trunk.logs[[i]], sep="\t")
  tmp = gsub(".trunk.txt", "", strsplit(trunk.logs[[i]], split="/")[[1]][[4]])
  tmp2 = strsplit(tmp, split="\\.")[[1]]
  if (!grepl("simulated", trunk.logs[[i]])){
    dat = rbind(dat, data.frame(events=t$V1, time=t$V3, virus=tmp2[[1]], simulation="inferred", trunk="is trunk"))
    dat = rbind(dat, data.frame(events=t$V2, time=t$V4, virus=tmp2[[1]], simulation="inferred", trunk="not trunk"))
    diff = rbind(diff, data.frame(diff=t$V1/t$V3 -t$V2/t$V4, virus=tmp2[[1]], simulation="inferred"))
  }else{
    dat = rbind(dat, data.frame(events=t$V1, time=t$V3, virus=tmp2[[1]], simulation="simulated", trunk="is trunk"))
    dat = rbind(dat, data.frame(events=t$V2, time=t$V4, virus=tmp2[[1]], simulation="simulated", trunk="not trunk"))
    diff = rbind(diff, data.frame(diff=t$V1/t$V3 -t$V2/t$V4, virus=tmp2[[1]], simulation="simulated"))
  }
}

p = ggplot(dat, aes(x=interaction(virus), y=events/time))+
  geom_violin(aes(color=trunk)) +
  facet_grid(.~simulation) +
  scale_y_log10() +
  coord_cartesian(ylim=c(0.01,2)) +
  theme_minimal()

plot(p)
diff = diff[-which(diff$virus=="mers"),]
p = ggplot(diff, aes(x=virus, y=diff))+
  geom_violin(aes(color=simulation)) +
  # coord_cartesian(ylim=c(-0.25,0.25)) +
  theme_minimal() +
  ylab("difference between recombination rates\non fit and unfit edges")

plot(p)
ggsave(plot=p, paste("../../Recombination-Text/Figures/trunk_rate.pdf", sep=""), width=6, height=4)




trunk.logs <- list.files(path="./combined_spike/", pattern=paste("*.trunk.txt", sep=""), full.names = TRUE)
dat = data.frame()
diff = data.frame()
for (i in seq(1, length(trunk.logs))){
  t = read.table(trunk.logs[[i]], sep="\t")
  tmp = gsub(".trunk.txt", "", strsplit(trunk.logs[[i]], split="/")[[1]][[4]])
  tmp2 = strsplit(tmp, split="\\.")[[1]]
  if (!grepl("simulated", trunk.logs[[i]])){
    dat = rbind(dat, data.frame(events=t$V1, time=t$V3, virus=tmp2[[1]], simulation="inferred", trunk="is trunk"))
    dat = rbind(dat, data.frame(events=t$V2, time=t$V4, virus=tmp2[[1]], simulation="inferred", trunk="not trunk"))
    diff = rbind(diff, data.frame(diff=t$V1/t$V3 -t$V2/t$V4, virus=tmp2[[1]], simulation="inferred"))
  }else{
    dat = rbind(dat, data.frame(events=t$V1, time=t$V3, virus=tmp2[[1]], simulation="simulated", trunk="is trunk"))
    dat = rbind(dat, data.frame(events=t$V2, time=t$V4, virus=tmp2[[1]], simulation="simulated", trunk="not trunk"))
    diff = rbind(diff, data.frame(diff=t$V1/t$V3 -t$V2/t$V4, virus=tmp2[[1]], simulation="simulated"))
  }
}

p = ggplot(dat, aes(x=interaction(virus), y=events/time))+
  geom_violin(aes(color=trunk)) +
  facet_grid(.~simulation) +
  scale_y_log10() +
  coord_cartesian(ylim=c(0.01,2)) +
  theme_minimal()

plot(p)
diff = diff[-which(diff$virus=="mers"),]
p = ggplot(diff, aes(x=virus, y=diff))+
  geom_violin(aes(color=simulation)) +
  # coord_cartesian(ylim=c(-0.25,0.25)) +
  theme_minimal() +
  ylab("difference between recombination rates\non fit and unfit edges")

plot(p)
ggsave(plot=p, paste("../../Recombination-Text/Figures/trunk_rate_spike.pdf", sep=""), width=6, height=4)


