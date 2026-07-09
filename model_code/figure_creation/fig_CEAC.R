###############################################
# TITLE: CEAC plotting
### This script is for plotting the cost-effectiveness acceptability curve (CEAC)

### Clear environment
rm(list=ls())

## Load required packages
library(ggplot2)
library(data.table)

## load the PSA outputs
psa.outputs.dt <- read.csv("outputs/psa_outputs.csv")
psa.outputs.dt <- as.data.table(psa.outputs.dt)



options(scipen=10000)

# Generate CEAC table
WTP.values <- seq(from = 0, to = 50000, by = 1000) ## use the seq() function to get a vector of specified numeric values

CEAC <- data.frame(WTP = WTP.values, 
                   `no intervention`= rep(as.numeric(NA),length(WTP.values)),
                   intervention= rep(as.numeric(NA),length(WTP.values)))


pCE.3b <-function(WTP, psa.outputs.dt) {
  ## a function that estimates the probability of the new intervention
  # being cost-effective 
  # OUTPUTS: A numeric value specifying the probability of cost-effectiveness 
  
  need <- c("standard_qol","intervention_qol","standard_costs","intervention_costs")
  if (!all(need %in% names(psa.outputs.dt))) stop("CEAC: missing columns: ", paste(setdiff(need, names(psa.outputs.dt)), collapse=", "))
  
  
  nmb <- ((psa.outputs.dt[,c("standard_qol",
                             "intervention_qol")])*WTP) - 
    (psa.outputs.dt[,c("standard_costs",
                       "intervention_costs")] )
  
  max.nmb <- apply(nmb, 1, max) # selecting max value indication by row within nmb
  
  ## creating an indication of TRUE/FALSE as to whether each treatment column == that max value:
  CE <- nmb[1:nrow(psa.outputs.dt),] == max.nmb[1:nrow(psa.outputs.dt)] 
  probCE<- apply(CE, 2, mean) ## averaging over TRUE (=1) and FALE (=0) for each column
  
  return(probCE)
  
}


pb = txtProgressBar(min = 0, max = length(WTP.values), initial = 0, style = 3)

for (i in 1:length(WTP.values)) {
  setTxtProgressBar(pb,i)
  CEAC[i,"WTP"] <- WTP.values[i]
  CEAC[i,2:3]<- pCE.3b(WTP.values[i], psa.outputs.dt)
  
}


# Look at the results
head(CEAC)  
tail(CEAC)
CEAC <- as.data.table(CEAC)
CEAC.long <- melt(CEAC, id.vars = c("WTP"))
colnames(CEAC.long) <- c("WTP", "group", "pCE")
head(CEAC.long)

write.csv(CEAC.long, file="outputs/CEAC_long.csv")

CEAC.long$Group<-CEAC.long$group


plot.ceac.all <- function(results){
  ## FUNCTION: Cost-effectiveness acceptability curve (CEAC) for multiple comparators
  ## INPUTS: a results data frame that has the columns "WTP", "group" & "pCE"
  ## OUTPUT: CEAC curve where the color of the CEAC differs for each comparator
  
  xlabel = "Willingness To Pay Threshold"
  ylabel = "Probability Cost-Effective"
  
  plot = ggplot(results) + 
    geom_line(aes(x=WTP, y=pCE, color=Group), size=1) + 
    labs(x = xlabel, y = ylabel) + 
    scale_color_discrete(labels = c("no.intervention" = "Standard Care", 
                                    "intervention" = "B360S"))+
    theme_classic() +
    theme(
      legend.title = element_blank(),
      axis.title = element_text(face="bold", size = 25), # Increase axis title font size
      axis.text.x = element_text(size = 25), # Increase x-axis text font size
      axis.text.y = element_text(size = 25), # Increase y-axis text font size
      axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
      axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      legend.key.width = unit(1.8, "line"),
      text = element_text(size = 14) # Adjust other text elements if necessary
    ) + 
    scale_x_continuous(expand = c(0, 0.1)) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), expand = c(0, 0)) + 
    theme_linedraw()
  
  return(plot)
  
}

plot.ceac.all(CEAC.long)
ggsave("outputs/plot_ceac_all.png", width = 120, height = 100, units='mm',dpi=1000)

#### for narrative text
# CEAC thresholds from PSA
ceac_int <- CEAC.long[group == "intervention"]

# WTP where pCE first exceeds 90% and 95%
wtp_90 <- ceac_int[pCE >= 0.90, min(WTP)]
wtp_95 <- ceac_int[pCE >= 0.95, min(WTP)]

cat(sprintf("WTP at which pCE >= 90%%: £%s\n", comma(wtp_90)))
cat(sprintf("WTP at which pCE >= 95%%: £%s\n", comma(wtp_95)))

# pCE at NICE thresholds
pCE_20k <- ceac_int[WTP == 20000, pCE]
pCE_30k <- ceac_int[WTP == 30000, pCE]
cat(sprintf("pCE at £20,000/QALY: %.1f%%\n", pCE_20k * 100))
cat(sprintf("pCE at £30,000/QALY: %.1f%%\n", pCE_30k * 100))