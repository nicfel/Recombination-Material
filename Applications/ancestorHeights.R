library(ggplot2)
library(coda)
library("colorblindr")



# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

trees <- list.files(path="./combined/", pattern=paste("*.pruned", sep=""), full.names = TRUE)

# system(paste("java -Xmx16g -jar ../Software/CommonAncestorHeights.jar -burnin 0 -sequence Wuhan-Hu-1\\/2019\\|2019", trees[9], gsub("*.pruned", ".wuhan.tsv", trees[9]), sep=" "))
# system(paste("java -Xmx16g -jar ../Software/CommonAncestorHeights.jar -burnin 0 -sequence bat\\/Yunnan\\/RmYN02\\/2019\\|2019", trees[9], gsub("*.pruned", ".rmyn02.tsv", trees[9]), sep=" "))
# system(paste("java -Xmx16g -jar ../Software/CommonAncestorHeights.jar -burnin 0 -sequence pangolin\\/Guangdong\\/P2S\\/2019\\|2019", trees[9], gsub("*.pruned", ".p2s.tsv", trees[9]), sep=" "))
# system(paste("java -Xmx16g -jar ../Software/CommonAncestorHeights.jar -burnin 0 -sequence bat\\/Yunnan\\/RaTG13\\/2013\\|2013", trees[9], gsub("*.pruned", ".ratg13.tsv", trees[9]), sep=" "))
# system(paste("java -Xmx16g -jar ../Software/CommonAncestorHeights.jar -burnin 0 -sequence Sin842\\|2003", trees[9], gsub("*.pruned", ".sin842.tsv", trees[9]), sep=" "))

compareto = c("Sin842|2003", 
              "bat/Yunnan/RmYN02/2019|2019", 
              "pangolin/Guangdong/P2S/2019|2019", 
              "bat/Yunnan/RaTG13/2013|2013",
              "Wuhan-Hu-1/2019|2019")

name = c("SARS-CoV-1",
         "RmYN02",
         "P2S",
         "RaTG13",
         "SARS-CoV-2")

cols = c("SARS-CoV-1"="#1B9E77",
         "RmYN02"="#D95F02",
         "P2S"="#7570B3",
         "RaTG13"="#E7298A",
         "SARS-CoV-2"="#66A61E")


t = read.table("./combined/sars-like_all.wuhan.tsv", header=T, sep="\t")

# read genbak file for sars-like
require("genbankr")
require("gggenes")

gb = readGenBank(paste(this.dir, "/reference/sars-like.gb", sep="/"))
gen_data = data.frame(min=as.numeric(gb@cds@ranges@start), width=as.numeric(gb@cds@ranges@width), gene=gb@cds$gene)

t = t[is.element(t$name, compareto), ]
t$plotname = ""
for(i in seq(1,length(name))){
  t[which(t$name== compareto[[i]]), "plotname"] = name[[i]]
}

t = t[which(t$lower>0),]

mrsi=2019
breaks = c(0,19,49,119,219,419,1019)



gen_data$min = gen_data$min-100

p = ggplot(t[which(t$position!=-1),]) +
  geom_ribbon(aes(x=position, ymin=lower, ymax=upper, fill=plotname, group=plotname), alpha=0.2)+
  geom_line(aes(x=position, y=median, color=plotname, group=plotname))+
  theme_minimal() +
  scale_fill_manual(name="",values=cols)+
  scale_color_manual(name="",values=cols)+
  theme(legend.position = "top", legend.direction = "vertical", legend.title = element_blank()) +
  ylab("most recent common ancestor time\nwith SARS-CoV-2")+
  xlab("")+
  coord_cartesian(ylim = c(20,1200), xlim=c(0,29900))+
  scale_y_log10(breaks=breaks, labels=mrsi-breaks) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  geom_gene_arrow(data=gen_data, aes(xmin = min, xmax = min+width, y = 20), arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"))+
  geom_gene_label(data=gen_data, aes(xmin = min, xmax = min+width, y = 20, label = gene), align = "left")
plot(p)

ggsave(paste("../../Recombination-Text/Figures/networks/commonAncestorTime.sars2.pdf", sep=""), plot=p, width=6 ,heigh=3)

# breaks = c(0,19,69,119,169,219,269,319)

name = c(
         "RmYN02",
         "P2S",
         "RaTG13",
         "SARS-CoV-2")

p = ggplot(t[which(t$position==-1),]) +
  geom_point(aes(x=plotname, y=median, color=plotname))+
  geom_linerange(aes(x=plotname, ymin=lower, ymax=upper, color=plotname))+
  theme_minimal() +
  scale_color_manual(name="",values=cols)+
  theme(legend.position = "top", legend.direction = "vertical", legend.title = element_blank()) +
  ylab("most recent\nMRCA time\nwith SARS-CoV-2")+
  xlab("")+
  scale_y_log10(breaks=breaks, labels=mrsi-breaks) +
  scale_x_discrete(limits=c("RmYN02","RaTG13","P2S","SARS-CoV-1"))+
  coord_cartesian(ylim = c(10,350))+
  theme(legend.position = "none")
plot(p)

ggsave(paste("../../Recombination-Text/Figures/networks/mrmrca.sars2.pdf", sep=""), plot=p, width=4 ,heigh=2)

t = read.table("./combined/sars-like_all.rmyn02.tsv", header=T, sep="\t")

t = t[is.element(t$name, compareto), ]
t$plotname = ""
for(i in seq(1,length(name))){
  t[which(t$name== compareto[[i]]), "plotname"] = name[[i]]
}

t = t[which(t$lower>0),]

mrsi=2019
breaks = c(0,19,49,119,1019)

p = ggplot(t) +
  geom_ribbon(aes(x=position, ymin=lower, ymax=upper, fill=plotname, group=plotname), alpha=0.2)+
  geom_line(aes(x=position, y=median, color=plotname, group=plotname))+
  theme_minimal() +
  scale_fill_manual(name="",values=cols)+
  scale_color_manual(name="",values=cols)+
  theme(legend.position = "top", legend.direction = "vertical", legend.title = element_blank()) +
  ylab("common ancestor time\nwith RmYN02")+
  xlab("")+
  coord_cartesian(ylim = c(20,1200))+
  scale_y_log10(breaks=breaks, labels=mrsi-breaks) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  geom_gene_arrow(data=gen_data, aes(xmin = min, xmax = min+width, y = 20), arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"))+
  geom_gene_label(data=gen_data, aes(xmin = min, xmax = min+width, y = 20, label = gene), align = "left")
plot(p)

ggsave(paste("../../Recombination-Text/Figures/commonAncestorTime.rmyn02.pdf", sep=""), plot=p, width=8 ,heigh=4)


t = read.table("./combined/sars-like_all.ratg13.tsv", header=T, sep="\t")

t = t[is.element(t$name, compareto), ]
t$plotname = ""
for(i in seq(1,length(name))){
  t[which(t$name== compareto[[i]]), "plotname"] = name[[i]]
}

t = t[which(t$lower>0),]

mrsi=2019
breaks = c(0,19,49,119,1019)

p = ggplot(t) +
  geom_ribbon(aes(x=position, ymin=lower, ymax=upper, fill=plotname, group=plotname), alpha=0.2)+
  geom_line(aes(x=position, y=median, color=plotname, group=plotname))+
  theme_minimal() +
  scale_fill_manual(name="",values=cols)+
  scale_color_manual(name="",values=cols)+
  theme(legend.position = "top", legend.direction = "vertical", legend.title = element_blank()) +
  ylab("common ancestor time\nwith RaTG13")+
  xlab("")+
  coord_cartesian(ylim = c(20,1200))+
  scale_y_log10(breaks=breaks, labels=mrsi-breaks) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  geom_gene_arrow(data=gen_data, aes(xmin = min, xmax = min+width, y = 20), arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"))+
  geom_gene_label(data=gen_data, aes(xmin = min, xmax = min+width, y = 20, label = gene), align = "left")
plot(p)

ggsave(paste("../../Recombination-Text/Figures/commonAncestorTime.ratg13.pdf", sep=""), plot=p, width=8 ,heigh=4)


t = read.table("./combined/sars-like_all.p2s.tsv", header=T, sep="\t")

t = t[is.element(t$name, compareto), ]
t$plotname = ""
for(i in seq(1,length(name))){
  t[which(t$name== compareto[[i]]), "plotname"] = name[[i]]
}

t = t[which(t$lower>0),]
p = ggplot(t) +
  geom_ribbon(aes(x=position, ymin=lower, ymax=upper, fill=plotname, group=plotname), alpha=0.2)+
  geom_line(aes(x=position, y=median, color=plotname, group=plotname))+
  theme_minimal() +
  scale_fill_manual(name="",values=cols)+
  scale_color_manual(name="",values=cols)+
  theme(legend.position = "top", legend.direction = "vertical", legend.title = element_blank()) +
  ylab("common ancestor time\nwith P2S")+
  xlab("")+
  coord_cartesian(ylim = c(20,1200))+
  scale_y_log10(breaks=breaks, labels=mrsi-breaks) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  geom_gene_arrow(data=gen_data, aes(xmin = min, xmax = min+width, y = 20), arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"))+
  geom_gene_label(data=gen_data, aes(xmin = min, xmax = min+width, y = 20, label = gene), align = "left")
plot(p)

ggsave(paste("../../Recombination-Text/Figures/commonAncestorTime.p2s.pdf", sep=""), plot=p, width=8 ,heigh=4)




t = read.table("./combined/sars-like_all.sin842.tsv", header=T, sep="\t")

t = t[is.element(t$name, compareto), ]
t$plotname = ""
for(i in seq(1,length(name))){
  t[which(t$name== compareto[[i]]), "plotname"] = name[[i]]
}

t = t[which(t$lower>0),]
p = ggplot(t) +
  geom_ribbon(aes(x=position, ymin=lower, ymax=upper, fill=plotname, group=plotname), alpha=0.2)+
  geom_line(aes(x=position, y=median, color=plotname, group=plotname))+
  theme_minimal() +
  scale_fill_manual(name="",values=cols)+
  scale_color_manual(name="",values=cols)+
  theme(legend.position = "top", legend.direction = "vertical", legend.title = element_blank()) +
  ylab("common ancestor time\nwith SARS-CoV-1")+
  xlab("")+
  coord_cartesian(ylim = c(20,1200))+
  scale_y_log10(breaks=breaks, labels=mrsi-breaks) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  geom_gene_arrow(data=gen_data, aes(xmin = min, xmax = min+width, y = 20), arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"))+
  geom_gene_label(data=gen_data, aes(xmin = min, xmax = min+width, y = 20, label = gene), align = "left")
plot(p)

ggsave(paste("../../Recombination-Text/Figures/commonAncestorTime.sars1.pdf", sep=""), plot=p, width=8 ,heigh=4)



system(paste("java -Xmx16g -jar ../Software/CommonAncestorHeights.jar -burnin 0 -sequence bat\\_SL\\_CoVZXC21\\|2015", trees[9], gsub("*.pruned", ".zxc21.tsv", trees[9]), sep=" "))
system(paste("java -Xmx16g -jar ../Software/CommonAncestorHeights.jar -burnin 0 -sequence bat\\_SL\\_CoVZC45\\|2017", trees[9], gsub("*.pruned", ".zc45.tsv", trees[9]), sep=" "))


t = read.table("./combined/sars-like_all.zxc21.tsv", header=T, sep="\t")
name = c("SARS-CoV-1",
         "RmYN02",
         "P2S",
         "RaTG13",
         "SARS-CoV-2")

t = t[is.element(t$name, compareto), ]
t$plotname = ""
for(i in seq(1,length(name))){
  t[which(t$name== compareto[[i]]), "plotname"] = name[[i]]
}

t = t[which(t$lower>0),]
p = ggplot(t) +
  geom_ribbon(aes(x=position, ymin=lower, ymax=upper, fill=plotname, group=plotname), alpha=0.2)+
  geom_line(aes(x=position, y=median, color=plotname, group=plotname))+
  theme_minimal() +
  scale_fill_manual(name="",values=cols)+
  scale_color_manual(name="",values=cols)+
  theme(legend.position = "top", legend.direction = "vertical", legend.title = element_blank()) +
  ylab("common ancestor time\nwith ZXC21")+
  xlab("")+
  coord_cartesian(ylim = c(20,1200))+
  scale_y_log10(breaks=breaks, labels=mrsi-breaks) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  geom_gene_arrow(data=gen_data, aes(xmin = min, xmax = min+width, y = 20), arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"))+
  geom_gene_label(data=gen_data, aes(xmin = min, xmax = min+width, y = 20, label = gene), align = "left")
plot(p)

ggsave(paste("../../Recombination-Text/Figures/commonAncestorTime.zxc21.png", sep=""), plot=p, width=8 ,heigh=4)



t = read.table("./combined/sars-like_all.zc45.tsv", header=T, sep="\t")

t = t[is.element(t$name, compareto), ]
t$plotname = ""
for(i in seq(1,length(name))){
  t[which(t$name== compareto[[i]]), "plotname"] = name[[i]]
}

t = t[which(t$lower>0),]
p = ggplot(t) +
  geom_ribbon(aes(x=position, ymin=lower, ymax=upper, fill=plotname, group=plotname), alpha=0.2)+
  geom_line(aes(x=position, y=median, color=plotname, group=plotname))+
  theme_minimal() +
  scale_fill_manual(name="",values=cols)+
  scale_color_manual(name="",values=cols)+
  theme(legend.position = "top", legend.direction = "vertical", legend.title = element_blank()) +
  ylab("common ancestor time\nwith ZC45")+
  xlab("")+
  coord_cartesian(ylim = c(20,1200))+
  scale_y_log10(breaks=breaks, labels=mrsi-breaks) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  geom_gene_arrow(data=gen_data, aes(xmin = min, xmax = min+width, y = 20), arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"))+
  geom_gene_label(data=gen_data, aes(xmin = min, xmax = min+width, y = 20, label = gene), align = "left")
plot(p)

ggsave(paste("../../Recombination-Text/Figures/commonAncestorTime.zc45.png", sep=""), plot=p, width=8 ,heigh=4)

