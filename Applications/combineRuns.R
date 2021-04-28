library(ggplot2)
library(coda)


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# # get the names of all output files of the first replicate
# log <- list.files(path="./out/", pattern=paste("*rep0.log", sep=""), full.names = TRUE)
# for (i in seq(1, length(log))){
#   in_command <- " -b 20 -resample 50000 -log"
#   for (j in seq(0,2)){
#     in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), log[i]), sep="")
#   }
# 
#   out_command = gsub("rep0_", "", log[i])
#   out_command = gsub("out", "combined", out_command)
# 
#   combined_command = gsub(".trees",".trees", out_command)
#   combined_command = paste(" -o ", gsub("_rep0", "",combined_command), sep="")
#   # combine the trees
#   system(paste("/Applications/BEAST\\ 2.6.3/bin/logcombiner", gsub(".network.trees",".log", in_command), gsub(".network.trees",".log", combined_command), sep=" "))
# }
# 
# 
# # get the names of all output files of the first replicate
# trees <- list.files(path="./out/", pattern=paste("*rep0.trees", sep=""), full.names = TRUE)
# for (i in seq(1, length(trees))){
#   in_command <- " -b 20 -resample 100000 -log"
#   for (j in seq(0,2)){
#     in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), trees[i]), sep="")
#   }
# 
#   out_command = gsub("rep0_", "", trees[i])
#   out_command = gsub("out", "combined", out_command)
# 
#   combined_command = gsub(".trees",".trees", out_command)
#   combined_command = paste(" -o ", gsub("_rep0", "",combined_command), sep="")
#   # combine the trees
#   system(paste("/Applications/BEAST\\ 2.6.3/bin/logcombiner", gsub(".network.trees",".log", in_command), gsub(".network.trees",".log", combined_command), sep=" "))
# 
# }
# 
# # get the names of all output files of the first replicate
# trees <- list.files(path="./combined/", pattern=paste("*.trees", sep=""), full.names = TRUE)
# for (i in seq(1, length(trees))){
# 
#   system(paste("java -jar ../Software/ReSim.jar -log", gsub(".trees", ".log", trees[i]), "-xmlFile", gsub("combined","xmls",gsub("*.trees", "_rep0.xml", trees[i])),  trees[i], gsub("*.trees", ".fullsim", trees[i]), sep=" "))
# }
# 
# 
# trees <- list.files(path="./combined/", pattern=paste("*.trees", sep=""), full.names = TRUE)
# for (i in seq(1, length(trees))){
#   if (!grepl("simulated", trees[[i]])){
#     system(paste("java -jar ../Software/NetworkPruner.jar -xmlFile", gsub("combined","xmls",gsub("*.trees", "_rep0.xml", trees[i])), trees[i], gsub("*.trees", ".pruned", trees[i]), sep=" "))
#   }
# }
# trees <- list.files(path="./combined/", pattern=paste("*.fullsim", sep=""), full.names = TRUE)
# for (i in seq(1, length(trees))){
#   if (!grepl("simulated", trees[[i]])){
#     system(paste("java -jar ../Software/NetworkPruner.jar -xmlFile", gsub("combined","xmls",gsub("*.fullsim", "_rep0.xml", trees[i])), trees[i], gsub("*.fullsim", ".simulated", trees[i]), sep=" "))
#   }
# }

# get the names of all output files of the first replicate
trees <- list.files(path="./combined/", pattern=paste("*.pruned", sep=""), full.names = TRUE)
for (i in seq(1, length(trees))){
  if (!grepl("simulated", trees[[i]])){
      # system(paste("java -jar ../Software/RecombinationSummarizer.jar -burnin 0 -positions NONE", trees[i], gsub("*.pruned", ".tree", trees[i]), sep=" "))
    if (grepl("sars-like_all", trees[[i]])){
      # system(paste("java -jar ../Software/RecombinationSummarizer.jar -burnin 0 -positions NONE", trees[i], gsub("*.pruned", ".tree", trees[i]), sep=" "))
      # system(paste("java -jar ../Software/RecombinationSummarizer.jar -burnin 0 -positions NONE -useEveryTree 10", trees[i], gsub("*.pruned", ".loci.tree", trees[i]), sep=" "))
    }
    # }else{
    #   system(paste("java -jar ../Software/RecombinationSummarizer.jar -burnin 0 ", trees[i], gsub("*.pruned", ".tree", trees[i]), sep=" "))
    #
    # }
  }
}

# get the names of all output files of the first replicate
trees <- list.files(path="./combined/", pattern=paste("*trees", sep=""), full.names = TRUE)
for (i in seq(1, length(trees))){
  if (!grepl("simulated", trees[[i]])){
    # system(paste("java -jar ../Software/RecombinationSummarizer.jar -burnin 0 -positions NONE", trees[i], gsub("*.pruned", ".tree", trees[i]), sep=" "))
    if (grepl("sars-like_all", trees[[i]])){
      # system(paste("java -jar ../Software/RecombinationSummarizer.jar -burnin 0 -positions NONE", trees[i], gsub("*.pruned", ".tree", trees[i]), sep=" "))
      system(paste("java -jar ../Software/RecombinationSummarizer.jar -burnin 0 -positions NONE -useEveryTree 10", trees[i], gsub("*.trees", ".loci.tree", trees[i]), sep=" "))
    }
    # }else{
    #   system(paste("java -jar ../Software/RecombinationSummarizer.jar -burnin 0 ", trees[i], gsub("*.pruned", ".tree", trees[i]), sep=" "))
    #
    # }
  }
}


# trees <- list.files(path="./combined/", pattern=paste("*.pruned", sep=""), full.names = TRUE)
# for (i in seq(1, length(trees))){
#   system(paste("java -jar ../Software/SameTree.jar -burnin 0", trees[i], gsub("*.pruned", ".length", trees[i]), sep=" "))
# }

# years = c(1,2,5,10)
# for (j in seq(1,length(years))){
#   trees <- list.files(path="./combined/", pattern=paste("*.trees", sep=""), full.names = TRUE)
#   for (i in seq(1, length(trees))){
#     if (!grepl("sars", trees[[i]]) && !grepl("mers", trees[[i]]) && grepl("subunits", trees[[i]])){
#       system(paste("java -jar ../Software/TrunkRecombination.jar -burnin 0 -trunkDefinition minTipDistance -minTipDistance", years[j], trees[i], gsub("*.trees", paste(".",years[j],".trunk", sep=""), trees[i]), sep=" "))
#       # system(paste("java -jar ../Software/TrunkRecombination.jar -burnin 0", trees[i], gsub("*.trees", ".trunk", trees[i]), sep=" "))
#     }
#   }
#
#   trees <- list.files(path="./combined/", pattern=paste("*.fullsim", sep=""), full.names = TRUE)
#   for (i in seq(1, length(trees))){
#     if (!grepl("sars", trees[[i]]) && !grepl("mers", trees[[i]]) && !grepl("subunits", trees[[i]])){
#       system(paste("java -jar ../Software/TrunkRecombination.jar -burnin 0 -trunkDefinition minTipDistance -minTipDistance", years[j], trees[i], gsub("*.fullsim",  paste(".",years[j],".trunksim", sep=""), trees[i]), sep=" "))
#       # system(paste("java -jar ../Software/TrunkRecombination.jar -burnin 0 /", trees[i], gsub("*.fullsim", ".trunksim", trees[i]), sep=" "))
#     }
#   }
# }


# # # get the names of all output files of the first replicate
# trees <- list.files(path="./combined/", pattern=paste("*.pruned", sep=""), full.names = TRUE)
# for (i in seq(1, length(trees))){
#   if (!grepl("sars", trees[[i]])){
#     system(paste("java -jar ../Software/RateAfterRecombination.jar -burnin 0 -maxDist 10", trees[i], gsub("*.pruned", ".recrate", trees[i]), sep=" "))
#     system(paste("java -jar ../Software/RateAfterRecombination.jar -burnin 0 -maxDist 10", gsub("*.pruned", ".simulated", trees[i]), gsub("*.pruned", ".simrecrate", trees[i]), sep=" "))
#   }
# }
# get the names of all output files of the first replicate
trees <- list.files(path="./combined/", pattern=paste("*pruned", sep=""), full.names = TRUE)
for (i in seq(1, length(trees))){
  if (!grepl("sarss", trees[[i]])){
    system(paste("java -jar ../Software/RateAfterRecombination.jar -burnin 0 -maxDist 100", trees[i], gsub("*.pruned", ".recrate", trees[i]), sep=" "))
    system(paste("java -jar ../Software/RateAfterRecombination.jar -burnin 0 -maxDist 100", gsub("*.pruned", ".simulated", trees[i]), gsub("*.pruned", ".simrecrate", trees[i]), sep=" "))
  }
}

das
# # get the names of all output files of the first replicate
# trees <- list.files(path="./combined/", pattern=paste("*pruned.trees", sep=""), full.names = TRUE)
# for (i in seq(1, length(trees))){
#   system(paste("java -jar ../Software/RecombinationDistance.jar -burnin 0", trees[i], gsub("*.trees", ".distance.txt", trees[i]), sep=" "))
# }

# get the names of all output files of the first replicate
trees <- list.files(path="./combined/", pattern=paste("*.trees", sep=""), full.names = TRUE)
for (i in seq(1, length(trees))){
  if (!grepl("simulated", trees[[i]])){
    system(paste("java -jar ../Software/CommonAncestorHeights.jar -burnin 0 -rootOnly true", trees[i], gsub("*.trees", ".root.tsv", trees[i]), sep=" "))
  }
}

