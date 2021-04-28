library(ggplot2)
library(coda)
library(phangorn)
library(ggtree)
library(phytools)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

Nnet <- read.nexus("./combined/sars-like_all.loci.tree")
par(mar=c(1,1,1,1))



tree.plot <- Nnet[seq(1,length(Nnet),10)]

# tree.plot = list()
color=list()
cols = c("#E69F00", "#009E73","#56B4E9")
  
for (i in seq(1,length(tree.plot))){
  for (j in seq(1,length(tree.plot[[1]]$tip.label))){
    tree.plot[[i]]$tip.label[[j]] = strsplit(tree.plot[[i]]$tip.label[[j]], split="\\|")[[1]][[1]]
  }
  tmp = as.numeric(strsplit(labels(tree.plot)[[i]], split="_")[[1]][[2]])
  if (tmp <21494){
    color[i] = "#E69F00"
  }else if (tmp < 25261){
    color[i] = "#009E73"
  }else{
    color[i] = "#56B4E9"
  }
}
pdf("../../Recombination-Text/Figures/sars_like_densi.pdf", width=15, height=8, onefile=FALSE)

end_point = 3000

densiTree(tree.plot, col=color, alpha=0.02, width=4, scale.bar=T, underscore=T)
axis(1,at=c((end_point-2020)/end_point,(end_point-1520)/end_point,(end_point-1020)/end_point,(end_point-520)/end_point,end_point/end_point), 
     labels=c(0,500,1000,1500,2020),col.axis="red", las=2)
legend("bottom", inset=.05, title="",border="white",box.col="white",
       c("5' of spike","spike","3' of spike"), fill=cols, horiz=F)
# axisPhylo(side = 1)
dev.off()


# p = ggdensitree(tree.plot, alpha=0.01) +scale_color_OkabeIto() + theme_tree2()
# plot(p)
# ggsave(plot=p, "../../Recombination-Text/Figures/sars_like_densi.png", width=10, height=5)
# 
