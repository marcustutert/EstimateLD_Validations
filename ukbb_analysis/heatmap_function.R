#Function to generate the heatmaps to show proportion of weights changing over time
#To do this I will need to have labelled samples (by super population and by sub-population, and an output of the weighst over time)
#Going to do this in ggplot2 for now...but hopefully transition to plotly
#install.packages("tidyverse")
#install.packages("tidyr")
#install.packages("reshape2")
library(tidyverse)
library(tidyr)
library(reshape2)
library(data.table)
#Get some of the data to play around with
gibbs_array = readRDS("~/Desktop/gibbs_array_heatmap_test")
#First thing we can do is subset this down to a matrix across the Nhaps and nSamples (full)
sample_ids  = read.table("~/Desktop/heatmap_id_input", header = T) #Pops and super pops/Samples
#Because weights do not change along the genome
data = gibbs_array[,1,]
#Rbind in the sample names/pops
data = cbind(sample_ids$pop,data)
data = cbind(sample_ids$super_pop,data)
data = cbind(sample_ids$sample,data)
#data[, lapply(.SD, sum, na.rm=TRUE), .SDcols=colsToSum]
data = as.data.table(data)
mine.long = melt(data,id.vars=c("V1","V2","V3"))
setDT(mine.long)[, ("value") := lapply(.SD, as.numeric), .SDcols = "value"]
setDT(mine.long)[, ("V2") := lapply(.SD, as.character), .SDcols = "V2"]
x = mine.long[, lapply(.SD, sum, na.rm=TRUE), by=c("V2","variable"), .SDcols=c("value") ]
afr = x[which(x$V2 == "AFR"),]$value/600
eur = x[which(x$V2 == "EUR"),]$value/500
fig <- plot_ly(data, x = ~seq(1:5760))
fig <- fig %>% add_trace(y = ~afr/(eur+afr), name = 'AFR',mode = 'lines')
fig <- fig %>% add_trace(y = ~eur/(eur+afr), name = 'EUR', mode = 'lines')
fig
#Sum up by column
mine.heatmap <- ggplot(data = mine.long, mapping = aes(x = V1,
                                                       y = variable,
                                                       fill = value)) +
  geom_tile() +
  xlab(label = "Update Number")

mine.heatmap

fig <- plot_ly(data, x = ~seq(1:5760))
fig <- fig %>% add_trace(y = ~ x[which(x$V2 == "AFR"),]$value/600, name = 'AFR',mode = 'lines')
fig <- fig %>% add_trace(y = ~x[which(x$V2 == "EUR"),]$value/500, name = 'EUR', mode = 'lines+markers')


