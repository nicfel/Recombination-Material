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
log <- list.files(path="./out/", pattern=paste("inf", ".*rep1.trees", sep=""), full.names = TRUE)


plot_list = list()
# for (i in seq(1,length(log))){
jj=1
for (i in seq(1,9)){
    
  in_command <- " -b 20 -log"
  for (j in seq(1,3)){
    in_command = paste(in_command, " ", gsub("rep1", paste("rep", j,sep=""), log[i]), sep="")
  }
  
  out_command = gsub("_rep1", "", log[i])
  out_command = gsub("out", "combined", out_command)

  combined_command = gsub(".trees",".trees", out_command)
  combined_command = paste(" -o ", combined_command, sep="")
  
  true_file = gsub("combined", "simout", out_command)
  true_file = gsub("inf_hky", "sim", true_file)
  true_file = gsub(".trees", ".tree", true_file)
  
  a = readLines(true_file)
  write(paste("tree","=",a), file=gsub("simout", "combined", true_file))

  # system(paste("/Applications/BEAST\\ 2.6.3/bin/logcombiner", in_command, combined_command, "", sep=" "), intern=TRUE, show.output.on.console = F, ignore.stderr=T)
  # system(paste("java -jar ../../Software/RecombinationDistance.jar ", out_command, gsub(".trees",".est.txt", out_command)))
  # system(paste("java -jar ../../Software/RecombinationDistance.jar -burnin 0 ", gsub("simout", "combined", true_file), gsub(".trees",".true.txt", out_command)))
  
  t.est = read.table( gsub(".trees",".est.txt", out_command), header=F, sep="\t")
  t.true = read.table( gsub(".trees",".true.txt", out_command), header=F, sep="\t")
  its = length(unique(t.est$V5))
  
  p1 <- ggplot(t.est, aes(x=V1, y=V3)) + 
    # geom_point(alpha = 0.1, size=0.1) +
    geom_density_2d_filled(, adjust=c(.5, 0.5),contour_var="ndensity")+
    # geom_density_2d(adjust=c(.5, 0.5),colour = "black") +
    # geom_point(alpha = 0.1, size=0.1, color="red") +
    
    geom_point(data=t.true, aes(x=V1, y=V3, size=(V2-V3)*V4), color="red", shape=4)+
    ylab("time in the past")+
    xlab("position on the genome") +
    theme_minimal() +
    # ylim(c(0,max(t.true$V3)*1.2))+
    ylim(c(0,30))+
    # ggtitle(log[[i]])+
    scale_size_continuous(name="Information content", limits=c(1,max((t.true$V2-t.true$V3)*t.true$V4)), range=c(0,3)) +
    theme(legend.position="none") +
    scale_fill_brewer()
  
  # dsa
  plot_list[[jj]]=p1
  jj=jj+1
}

require("gridExtra")
p_all1 = ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],
                  ncol = 3)
p_all2 = ggarrange(plot_list[[4]],plot_list[[5]],plot_list[[6]],
                  ncol = 3)

p_all3 = ggarrange(plot_list[[7]],plot_list[[8]],plot_list[[9]],
                  ncol = 3)

p_all= ggarrange(p_all1,p_all2,p_all3,ncol=1)

plot(p_all)

ggsave(p_all, filename="../../../Recombination-Text/Figures/recomb_support.pdf", height=8,width=16)
