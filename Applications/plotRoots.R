library(ggplot2)
library(coda)


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

root.logs <- list.files(path="./combined/", pattern=paste("*all\\.root.tsv", sep=""), full.names = TRUE)
dat = data.frame()

for (i in seq(1, length(root.logs)-1)){
  t = read.table(root.logs[[i]], sep="\t", header=T)
  tmp = gsub(".trunk.txt", "", strsplit(root.logs[[i]], split="/")[[1]][[4]])
  tmp2 = strsplit(tmp, split="\\.")[[1]]
  dat = rbind(dat, data.frame(position=t$position, mean=t$mean, median=t$median, lower=t$lower, upper=t$upper, virus=gsub("_all","",tmp2[[1]])))
}

dat$virus <- factor(dat$virus, levels = c("229e", "oc43", "nl63"))

p = ggplot(dat[which(dat$virus=="229e" | dat$virus=="oc43" | dat$virus=="nl63"), ], aes(x=position, y=2019-mean, fill=virus, color=virus))+
  geom_line()+
  geom_ribbon(aes(ymin=2019-lower, ymax=2019-upper), alpha=0.2) +
  facet_wrap(.~virus, ncol=4, scales="free") +
  theme_minimal()
plot(p)


