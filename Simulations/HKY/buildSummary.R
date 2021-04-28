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

# remove folders
system("rm -r combined")
system("rm -r mcc")
system("rm -r target")

system("mkdir combined")
system("mkdir mcc")
system("mkdir target")

trees <- list.files(path=paste("./treeout/",sep=""), pattern="*rep1.trees$", full.names = TRUE)

for (i in seq(1,length(trees))){

  in_command <- " -b 20 -log"
  for (j in seq(1,3)){
    in_command = paste(in_command, " ", gsub("rep1", paste("rep", j,sep=""), trees[i]), sep="")
  }

  out_command = gsub("rep1_", "", trees[i])
  out_command = gsub("treeout", "combined", out_command)

  combined_command = gsub(".trees",".trees", out_command)
  outfile_name = gsub("_rep1", "",combined_command)
  outfile = paste(" -o ", outfile_name, sep="")
  
  true_network = gsub(".trees", ".trueNetwork.tree", outfile_name)
  true_network = gsub("inf_hky_", "sim_", true_network)

  # combine the trees
  system(paste("/Applications/BEAST\\ 2.5.2/bin/logcombiner", in_command, outfile, sep=" "))
  # build the mcc tree
  system(paste("java -jar ./../../Software/ReassortmentNetworkSummarizer.jar -burnin 0", outfile_name,  gsub("/combined/","/mcc/", outfile_name), sep=" "))
  # add tree STATE_0 = to simulated file
  tree <- readLines(gsub("/combined/", "/simout/", true_network))
  if (length(tree)==1){
    fileConn<-file(gsub("/combined/", "/simout/", true_network))
    writeLines("#nexus", fileConn)
    close(fileConn)
    write("begin trees;", gsub("/combined/", "/simout/", true_network), append=TRUE)
    write(paste("tree STATE_0 =", tree), gsub("/combined/", "/simout/", true_network), append=TRUE)
    write("END;", gsub("/combined/", "/simout/", true_network), append=TRUE)
  }
  
  outfile_for_target = gsub("/combined/","/target/", outfile_name)
  
  # build the mcc tree with the true network as the target
  system(paste("java -jar ./../../Software/ReassortmentNetworkSummarizer.jar -burnin 0 -target ",  gsub("/combined/", "/simout/", true_network), 
               outfile_name,  outfile_for_target, sep=" "))

    # build the mcc tree with the true network as the target
  system(paste("java -jar ./../../Software/ReassortmentDistance.jar -burnin 0 ", 
               gsub("/combined/", "/simout/", true_network),  gsub(".trees", ".txt", outfile_for_target), sep=" "))
  
  # build the mcc tree with the true network as the target for only the first 2 segments
  system(paste("java -jar ./../../Software/ReassortmentNetworkSummarizer.jar -removeSegments 2,3 -burnin 0 -target ",  gsub("/combined/", "/simout/", true_network), 
               outfile_name,  gsub(".trees",".2seg.trees", outfile_for_target), sep=" "))
  
  # build the mcc tree with the true network as the target
  system(paste("java -jar ./../../Software/ReassortmentDistance.jar -burnin 0 -removeSegments 2,3 ", 
               gsub("/combined/", "/simout/", true_network),  gsub(".trees", ".2seg.txt", outfile_for_target), sep=" "))
  
  # build the mcc tree with the true network as the target for only the first 3 segments
  system(paste("java -jar ./../../Software/ReassortmentNetworkSummarizer.jar -removeSegments 3 -burnin 0 -target ",  gsub("/combined/", "/simout/", true_network), 
               outfile_name,  gsub(".trees",".3seg.trees", outfile_for_target), sep=" "))
  
  # build the mcc tree with the true network as the target
  system(paste("java -jar ./../../Software/ReassortmentDistance.jar -burnin 0 -removeSegments 3 ", 
               gsub("/combined/", "/simout/", true_network),  gsub(".trees", ".3seg.txt", outfile_for_target), sep=" "))
  
}

