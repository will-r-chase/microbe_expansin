library(tidyverse)
library(ggthemes)
library(ape)
library(ggtree)
library(readxl)
library(gridExtra)
library(viridis)
library(colormap)

setwd("~/Dropbox/AncientDNA/Expansin-endoglucanase-phylogeny/Microbial_Expansins")

## Read tree and find nodes to collapse, only monophyletic clades can be collapsed
tree <- read.tree("fileS1_msa_fixed_rooted.tree")

## Read data, clean up columns, attach to tree
data <- read_xlsx("ecology_metadata.xlsx")

data <- data[, 1:6]
data
colnames(data) <- c("organism", "group", "group2", "plant_path", "plant_associated", "ecology")
data$plant_associated[which(data$plant_associated=="No" & data$plant_path=="Yes")] <- "Yes"

## Attach to tree
p <- ggtree(tree)
p <- p %<+% data

colors<-c("#a6cee3", "#744f9f", "#cbcb3b", "black", "#1f36a0ff", "#52ae5aff", "#193f29ff", "#bf6b05ff", "#6a060fff", "black")
# "#6a060fff", "gray21" "#744f9f"

## Group tips by ecology
eco_split <- split(p$data$label, p$data$ecology)
eco_tree <- groupOTU(tree, eco_split, group_name = "ecology")

## Plot tree with legend
ecology_tree2 <-  ggtree(eco_tree, size = .8, aes(color=ecology)) + 
  xlim(0, 8) +
  scale_color_manual(values = colors) +
  theme(legend.position = "right") +
  geom_strip(296, 315, label = "Xanthomonads", barsize = 2, color = "black", align = T, fontsize = 6, offset = 1) +
  geom_strip(317, 347, label = "Firmicutes", barsize = 2, color = "black", align = T, fontsize = 6, offset = 1) +
  geom_strip(349, 360, label = "Enterobacteria", barsize = 2, color = "black", align = T, fontsize = 6, offset = 1) +
  geom_strip(363, 374, label = "Firmicutes", barsize = 2, color = "black", align = T, fontsize = 6, offset = 1) +
  geom_strip(380, 394, label = "Myxobacteria", barsize = 2, color = "black", align = T, fontsize = 6, offset = 1) +
  geom_strip(426, 500, label = "Actinobacteria", barsize = 2, color = "black", align = T, fontsize = 6, offset = 1) +
  geom_strip(510, 520, label = "B-Proteobacteria", barsize = 2, color = "black", align = T, fontsize = 6, offset = 1) +
  geom_strip(527, 602, label = "Actinobacteria", barsize = 2, color = "black", align = T, fontsize = 6, offset = 1) +
  geom_strip(2, 161, label = "Ascomycetes", barsize = 2, color = "black", align = T, fontsize = 6, offset = 1) +
  geom_strip(163, 174, label = "Basidiomycetes", barsize = 2, color = "black", align = T, fontsize = 6, offset = 1) +
  geom_strip(178, 207, label = "Amoebozoa", barsize = 2, color = "black", align = T, fontsize = 6, offset = 1) +
  geom_strip(225, 231, label = "Archaeplastida", barsize = 2, color = "black", align = T, fontsize = 6, offset = 1) +
  geom_strip(259, 291, label = "Stramenopiles", barsize = 2, color = "black", align = T, fontsize = 6, offset = 1)
ecology_tree2 + theme(legend.text = element_text(size = 10, face = "bold")) 


ggsave("fig2.trial.pdf", height=20, width=20)

## Make ecology counts table with colors defined
color_df <- data.frame(ecology = unique(data$ecology), color = "NA")
color_df <- color_df[order(color_df$ecology), ]
color_df$color <- colors

#######################################
######## important part ###############
#######################################
ecology_counts <- data %>%   # Change as.factor to domain to get Fungi vs Bacteria plots
  mutate(group2 = as.factor(group2), ecology = as.factor(ecology)) %>% #have to make a factor to use .drop = FALSE in group_by
  group_by(group2, ecology, .drop = FALSE) %>% #important!!! as of dplyr 0.8, group_by with .drop = FALSE will preserve 0-length groups
  tally() %>%
  mutate(percent = round((n/sum(n))*100, digits = 1)) %>%  ## Here have it divide by sum within the group, not by total number in the dataset
  inner_join(., color_df, by="ecology")

ecology_split <- split(ecology_counts, ecology_counts$group2)

## Make ecology breakdown barplots
ecology_plots <- imap(ecology_split, ~
  ggplot(.x, aes(x = ecology, y = percent)) + 
    geom_col(fill = .x$color) + 
    ylim(0, 115) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
    theme(legend.position = "none") +
    scale_fill_manual(values = colors) +
    geom_text(aes(label = paste0(percent, "%"), hjust = 0)) +
    labs(title = .y)
  )

ecology_plots 

## save as 5x5

## Make plant association barplots 
plant_associated_counts <- data %>% 
  group_by(group, plant_associated) %>%
  tally() %>%
  group_by(group) %>% 
  mutate(percent = round((n/sum(n))*100, digits = 1)) 

plant_associated_split <- split(plant_associated_counts, plant_associated_counts$group)

plant_associated_plots <- lapply(plant_associated_split, function(x){
  ggplot(x, aes(x = group, y = percent, fill = plant_associated)) +
    scale_fill_manual(values = c("white", "springgreen4", "black")) +
    geom_col(color = "black", size = 0.5) + 
    theme_inset() +
    labs(y = paste0(x$percent[2], "%"))
}
)
plant_associated_plots

## Supplemental Figure 2

ecology_tree2 <-
  ggtree(eco_tree, size = 1.25, aes(color=ecology)) +  # size = 1 controls line thickness
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))>=70 & as.numeric(sub(".*/", "", label))>=95 & !isTip), color="black", size=.01) +
  xlim(0, 6) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme(legend.position = "left") +  # Add legend at the top of the tree
  geom_treescale(width = 0.5, linesize = 0.5, fontsize = 5, y = 0, x = 0) +
 ecology_tree2 


## Define colors for ecology coloring
# Coloring guide https://bhaskarvk.github.io/colormap/
scales::show_col(colormap(colormap=colormaps$temperature, nshades=50))
scales::show_col(colormap(colormap=colormaps$phase, nshades=50))
scales::show_col(colormap(colormap=colormaps$cubehelix, nshades=50))
scales::show_col(colormap(colormap=colormaps$turbidity, nshades=50))
#scales::show_col(colormap(colormap=colormaps$salinity, nshades=100))
#scales::show_col(colormap(colormap=colormaps$oxygen, nshades=100))
#scales::show_col(colormap(colormap=colormaps$chlorophyll, nshades=100))

