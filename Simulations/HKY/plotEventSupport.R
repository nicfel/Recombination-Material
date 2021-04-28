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

# read in the true rates
support.vals <- read.table("event_support.csv", header=TRUE, sep=",")


p.dist <- ggplot(support.vals, aes(dist, support)) +
  # geom_density_2d() +
  geom_point(size=2, alpha=0.05) +
  geom_smooth(span = 0.5)+
  theme_light() +
  xlab("reassortment distance in years")+
  ylab("posterior support")+
  facet_grid(.~nrSegments) +
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,10))
plot(p.dist)
ggsave(plot=p.dist,"../../../Reassortment-Text/Figures/simulation/hky_reassortmentsupport.pdf",width=6, height=2.5)
  


