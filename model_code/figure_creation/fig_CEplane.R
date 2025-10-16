###############################################
# TITLE: CE plane
### this script is for plotting the cost-effectiveness plane
### inputs: outputs/psa_outputs.csv
### outputs: plot.ce.plane.png

### Clear environment
rm(list=ls())

## Load required packages
library(ggplot2)
library(data.table)

psa.outputs.dt <- read.csv("outputs/psa_outputs.csv")
psa.outputs.dt <- as.data.table(psa.outputs.dt)



options(scipen=10000)

plot.ce.plane <- function(results, wtp){
  ### FUNCTION: PLOTTING THE COST-EFFECTIVENESS PLANE (for one comparator, using incremetnal costs and outcomes)
  ### INPUTS: a results data frame that has the columns "inc.qol" and "inc.cost"
  ### Need ggplot called
  ### OUTPUTS: cost-effectiveness plane plot
  require(ggplot2)
  xlabel = "Incremental QALYs"
  ylabel = "Incremental Cost"
  
  stopifnot(all(c("inc.qol","inc.cost") %in% names(results)))
  
  plot = ggplot(results) + 
    geom_point(shape = 21, size = 2, colour = "black", fill = NA, alpha = 0.5,
               aes(x=inc.qol, y=inc.cost)) + 
    labs(x = xlabel, text = element_text(size=10)) + 
    labs (y = ylabel, text = element_text(size=10)) + theme_classic() +
    theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.key.width=unit(1.8,"line"), text = element_text(size=12),
          plot.margin=unit(c(1.2,0.5,0,1.2),"cm"))+
    geom_vline(xintercept=0,linetype=3)+
    geom_hline(yintercept=0,linetype=3)+
    geom_abline(slope = wtp,
                intercept = 0,
                color="blue")+ theme_linedraw()
  
  return(plot)
  
}
plot.ce.plane(psa.outputs.dt,20000)
ggsave("outputs/plot.ce.plane.png", width = 100, height = 100, units='mm',dpi=1000)


