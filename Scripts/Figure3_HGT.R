library(purrr)
library(ape)
library(tidytree)
library(ggplot2)
library(ggtree)
library(plyr)
library(dplyr)
library(tidyverse)
library(magrittr)
library(readxl)

setwd("~/Dropbox/AncientDNA/Expansin-endoglucanase-phylogeny/Microbial_Expansins")

## Read tree and find nodes to collapse, only monophyletic clades can be collapsed
tree <- read.tree("fileS1_msa_fixed_rooted.tree")

## Read in tree and look at node numbers
nodes_tree <- ggtree(tree) + 
  geom_tiplab(size = 2, hjust = -0.5) + 
  geom_text2(aes(label = node), size = 3, hjust = 1) +
  xlim(0,3)
nodes_tree
#ggsave("nodes.pdf", height=50, width=10, `limitsize = FALSE`)

nodes_to_collapse <- c( 614, # Ascomycota
                        774, # Basidiomycetes
                        788, # Amoebozoa
                        823, # Ascomycetes; 2nd HGT
                        837, # Viridiplantae and Haptophyta
                        848, # Stramenopiles
                        923, # Xanthomonads
                        938, # Xanthomonads 
                        943, # Paenibacillus spp Firmicutes
                        951, # Paenibacillus spp Firmicutes
                        983, # Bacillus spp Firmicutes
                        994, # Myxococcus
                        1001, # Myxococcus
                        1020, # Actinobacteria
                        1134) # Actinobacteria
                        


labels = c("Ascomycota", 
           "Basidiomycetes", 
           "Amoebozoa", 
           "Ascomycetes", 
           "Viridiplantae", 
           "Stramenopiles", 
           "Xanthomonads", 
           "Xanthomonads", 
           "Firmicutes", 
           "Firmicutes", 
           "Firmicutes", 
           "Myxobacteria", 
           "Myxobacteria", 
           "Actinobacteria", 
           "Actinobacteria")

## Read data, clean up columns, attach to tree
data <- read_delim("fileS4_microbe-data.txt", delim = "\t")
data <- data[, 1:3]
data
colnames(data) <- c("organism", "group", "group2")
p <- ggtree(tree)
p <- p %<+% data

## Group tips by taxonomy
taxonomy_tree <- split(p$data$label, p$data$group2)
taxonomy_tree <- groupOTU(tree, taxonomy_tree, group_name = "taxonomy")
taxonomy_tree$tip.label <- stringr::str_replace(taxonomy_tree$tip.label, "-", " ")

## Define colors for taxonomy coloring
colors1<-c("#a5c8e1", "#1c78b3", "#b1d689", "#b15928", "#f49596", "#2e9b47", "#4c3571", "#2a6b55", "#cbcb3b", "#fcbc6d", "#e11f26", "#c9afd4", "#f37c20", "#c9afd4", "#744f9f", "#72242b", "#744f9f", "black")
colors2<-c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black")

full_taxonomy_tree <- ggtree(taxonomy_tree, aes(color=taxonomy)) + 
  geom_tiplab(color = "black", size = 1.5) +
  geom_treescale(width = 0.5, linesize = 1, fontsize = 5, y = -30, x = 0) +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))>=70 & as.numeric(sub(".*/", "", label))>=95 & !isTip), color = "black", size = .1) +
  scale_color_manual(values = colors1) +
  theme(legend.position = "right")
full_taxonomy_tree

#ggsave("Supp_file2.pdf", height=25, width=15)

## Get number of children in each clade to print with the collapsed nodes
## This will produce a warning that can be ignored (Warning message:Unknown or uninitialised column: 'node'. )
tree_dat <- taxonomy_tree %>% 
  as.treedata() %>% 
  tidytree::as_tibble()

children<-list()

for(i in 1:length(nodes_to_collapse)){
  children[[i]]<-offspring(tree_dat, .node = nodes_to_collapse[i])
}

names(children)<-labels

children <- map(children, ~.[-(which(grepl("/", .$label))), ])

number_collapsed <- map_int(children, ~nrow(.))

## Make labels for collapsed clades
tree_labels <- paste(labels, paste0("(", paste(number_collapsed, "taxa)", sep = " ")), sep = " ")

## Plot tree
collapse_tree <-
  ggtree(taxonomy_tree, aes(color=taxonomy), size = 1.5) + 
  geom_treescale(width = 0.5, linesize = 1, fontsize = 5, y = -5, x = 0) +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))>=70 & as.numeric(sub(".*/", "", label))>=95 & !isTip), color = "black", size = 4) +
  geom_cladelabel(node = 903, label = "Bacteria", color = "black", barsize = 2, align = T, fontsize = 10, offset = 1) +
  geom_cladelabel(node = 610, label = "Eukaryotes", color = "black", barsize = 2, align = T, fontsize = 10, offset = 1) +
  geom_strip(608, 604, label = "Mixed Bacteria/Eukaryotes", color = "black", barsize = 2, align = T, fontsize = 10, offset = 1) +
  xlim(0, 8) +
  scale_color_manual(values = colors1) +
  scale_size_manual(values=c(1, 3)) + 
  geom_tiplab(color = "black", size = 2)
collapse_tree

## Collapse nodes
for(i in 1:length(nodes_to_collapse)){
  collapse_tree <- collapse(collapse_tree, node = nodes_to_collapse[i], clade_name=tree_labels[i])
}

## Plot collapsed tree w/ points
p2<-collapse_tree + 
  geom_point2(aes(subset = (node %in% nodes_to_collapse), fill=taxonomy), size = 8, shape = 23) + 
  geom_nodelab(aes(subset = (node %in% nodes_to_collapse)), color = "black", hjust = -0.05, size = 6) +
  scale_fill_manual(values = colors2)
p2
ggsave("trial_hgt.pdf", height = 30, width = 30)
### Edited for aesthetics in Illustrator ###

