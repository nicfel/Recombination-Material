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
log <- list.files(path="./combined/", pattern=".*all.trees", full.names = TRUE)

# read genbak file for sars-like
require("genbankr")
require("gridExtra")
require("ggpubr")
require("gggenes")




plot_list = list()

ymax = c(20, 5, 30,50,500, 100)

name = c("229E", "MERS", "NL63","OC43","SARS-like")

for (i in seq(1,4)){
  system(paste("java -jar ../Software/RecombinationDistance.jar -burnin 0 ", log[[i]], gsub(".trees",".dist.txt", log[[i]])))
  
  gb = readGenBank(paste(this.dir, "/reference/",tolower(name[[i]]),".gb", sep=""))
  gen_data = data.frame(min=as.numeric(gb@cds@ranges@start), width=as.numeric(gb@cds@ranges@width), gene=gb@cds$product)
  gen_data$min = gen_data$min-100
  
  t.est = read.table(gsub(".trees",".dist.txt", log[[i]]), header=F, sep="\t")

  p1 <- ggplot() + 
    geom_density_2d_filled(data=t.est, aes(x=V1, y=V3), adjust=c(0.05, 1),contour_var="ndensity")+
    geom_point(data=t.est, aes(x=V1, y=V3), alpha = 0.025, size=0.1, color="red") +
    ylab("time in the past")+
    xlab("position on the genome") +
    theme_minimal() +
    ylim(c(0,ymax[[i]]))+
    # scale_size_continuous(name="Information content", limits=c(1,max((t.true$V2-t.true$V3)*t.true$V4)), range=c(0,3)) +
    theme(legend.position="none")+
    scale_fill_brewer()+
    ggtitle(name[[i]]) +   
    geom_gene_arrow(data=gen_data, aes(xmin = min, xmax = min+width, y = 0),
                    arrow_body_height=unit(5, "mm"), 
                    arrowhead_height = unit(5, "mm"), 
                    arrowhead_width = unit(5, "mm"))+
    geom_gene_label(data=gen_data, aes(xmin = min, xmax = min+width, y = 0, label = gene), align = "left")
  plot(p1)
  plot_list[[i]]=p1
  
}
p_all1 = ggarrange(plot_list[[1]],plot_list[[2]],
                  ncol = 2)
p_all2 = ggarrange(plot_list[[3]],plot_list[[4]],
                  ncol = 2)

p_all= ggarrange(p_all1,p_all2,ncol=1)

plot(p_all)

ggsave(p_all, filename="../../Recombination-Text/Figures/cov_node_support.png", height=8,width=16)


system(paste("java -jar ../Software/RecombinationDistance.jar -burnin 0 -sequence Wuhan-Hu-1/2019\\|2019", log[[5]], gsub(".trees",".sars2dist.txt", log[[5]])))
t.est = read.table(gsub(".trees",".sars2dist.txt", log[[5]]), header=F, sep="\t")
gb = readGenBank(paste(this.dir, "/reference/",tolower(name[[5]]),".gb", sep=""))
gen_data = data.frame(min=as.numeric(gb@cds@ranges@start), width=as.numeric(gb@cds@ranges@width), gene=gb@cds$gene)
gen_data$min = gen_data$min-100


p1 <- ggplot(t.est) + 
  
  geom_density_2d_filled(data=t.est, aes(x=V1, y=V3), adjust=c(0.025, 1),contour_var="ndensity")+
  geom_point(data=t.est, aes(x=V1, y=V3), alpha = 0.25, size=0.1, color="red") +
  ylab("time in the past")+
  xlab("position on the genome") +
  theme_minimal() +
  # coord_cartesian(ylim=c(0,100))+
  theme(legend.position="none") +
  geom_gene_arrow(data=gen_data, aes(xmin = min, xmax = min+width, y = 0),arrow_body_height=unit(5, "mm"), arrowhead_height = unit(5, "mm"), arrowhead_width = unit(5, "mm"))+
  geom_gene_label(data=gen_data, aes(xmin = min, xmax = min+width, y = 0, label = gene), align = "left") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks=c(0,25,50,75,100),labels=c("2020", "1995","1970","1945","1920"), limits=c(-2,100)) +
  scale_fill_brewer()+
  ylab("")

plot(p1)
print(sum(t.est$V3<100)/length(unique(t.est$V5)))
print(1/length(unique(t.est$V5)))

ggsave(p1, filename="../../Recombination-Text/Figures/sars2_node_support.pdf", height=4,width=8)



plot_list = list()

ymax = c(20, 5, 30,50,500, 100)

name = c("229E", "MERS", "NL63","OC43","SARS-like")

for (i in seq(1,4)){
  gb = readGenBank(paste(this.dir, "/reference/",tolower(name[[i]]),".gb", sep=""))
  gen_data = data.frame(min=as.numeric(gb@cds@ranges@start), width=as.numeric(gb@cds@ranges@width), gene=gb@cds$product)
  gen_data$min = gen_data$min-100
  t.est = read.table(gsub(".trees",".dist.txt", log[[i]]), header=F, sep="\t")
  p1 <- ggplot() + 
    geom_density(data=t.est, aes(x=V1), adjust=c(0.05, 1))+
    ylab("Probability density")+
    xlab("position on the genome") +
    theme_minimal() +
    # ylim(c(0,ymax[[i]]))+
    # scale_size_continuous(name="Information content", limits=c(1,max((t.true$V2-t.true$V3)*t.true$V4)), range=c(0,3)) +
    theme(legend.position="none")+
    scale_fill_brewer()+
    ggtitle(name[[i]]) +   
    geom_gene_arrow(data=gen_data, aes(xmin = min, xmax = min+width, y = 0),
                    arrow_body_height=unit(5, "mm"), 
                    arrowhead_height = unit(5, "mm"), 
                    arrowhead_width = unit(5, "mm"))+
    geom_gene_label(data=gen_data, aes(xmin = min, xmax = min+width, y = 0, label = gene), align = "left")
  plot(p1)
  plot_list[[i]]=p1
  
}
p_all1 = ggarrange(plot_list[[1]],plot_list[[2]],
                   ncol = 2)
p_all2 = ggarrange(plot_list[[3]],plot_list[[4]],
                   ncol = 2)

p_all= ggarrange(p_all1,p_all2,ncol=1)
plot(p_all)

ggsave(p_all, filename="../../Recombination-Text/Figures/cov_node_density.pdf", height=8,width=16)

system(paste("java -jar ../Software/RecombinationDistance.jar -burnin 0", log[[5]], gsub(".trees",".sars2dist.all.txt", log[[5]])))


t.est = read.table(gsub(".trees",".sars2dist.all.txt", log[[5]]), header=F, sep="\t")
gb = readGenBank(paste(this.dir, "/reference/",tolower(name[[5]]),".gb", sep=""))
gen_data = data.frame(min=as.numeric(gb@cds@ranges@start), width=as.numeric(gb@cds@ranges@width), gene=gb@cds$gene)
gen_data$min = gen_data$min-100


p1 <- ggplot(t.est) + 
  geom_density(data=t.est, aes(x=V1), adjust=c(0.05, 1))+
  ylab("Probability density")+
  xlab("position on the genome") +
  theme_minimal() +
  # coord_cartesian(ylim=c(0,100))+
  theme(legend.position="none") +
  geom_gene_arrow(data=gen_data, aes(xmin = min, xmax = min+width, y = 0),arrow_body_height=unit(5, "mm"), arrowhead_height = unit(5, "mm"), arrowhead_width = unit(5, "mm"))+
  geom_gene_label(data=gen_data, aes(xmin = min, xmax = min+width, y = 0, label = gene), align = "left") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # scale_x_continuous(breaks = NULL) +
  # scale_y_continuous(breaks=c(0,25,50,75,100),labels=c("2020", "1995","1970","1945","1920"), limits=c(-2,100)) +
  scale_fill_brewer()+
  ylab("")

plot(p1)


ggsave(p1, filename="../../Recombination-Text/Figures/sars2_node_density.pdf", height=4,width=8)
