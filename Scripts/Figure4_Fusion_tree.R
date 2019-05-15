library(tidyverse)
library(ape)
library(tidytree)
library(ggtree)
library(readxl)
library(RColorBrewer)

setwd("~/Dropbox/AncientDNA/Expansin-endoglucanase-phylogeny/Microbial_Expansins")

## Read tree and find nodes to collapse, only monophyletic clades can be collapsed
tree <- read.tree("fileS1_msa_fixed_rooted.tree")

## Assign a list of nodes to collapse
nodes_to_collapse <- c(616,  # "Random Fungi"
                       738, # "Penicillium"
                       774, # "Puccinia"
                       788, # "Amoebozoa"
                       836, # "Single-celled Eukaryotes"
                       846, # "Oomycetes"
                       943, # "Paenibacillus"
                       951, # "Bacillus"
                       970, # "Enterobacteria"
                       988, # "Ruminococcus"
                       983, # "Paenibacillus"
                       990, # "Other Bacteria"
                       1001, # "Myxos"
                       759, # Aspergillus
                       768, # Aspergillus
                       1020, # "Streptomyces"
                       1135) # "Streptomyces"
 
labels = c("Random Fungi", "Penicillium", "Puccinia", "Amoebozoa", "Single-celled Eukaryotes", "Oomycetes", 
           "Paenibacillus", "Bacillus", "Enterobacteria", "Ruminococcus", "Paenibacillus", "Other Bacteria", 
           "Myxos", "Aspergillus", "Aspergillus", "Streptomyces", "Streptomyces")

## Read data, clean up columns, attach to tree
data_fusion <- read_xlsx("fusions_for_tree.xlsx")
p <- ggtree(tree)
p <- p %<+% data_fusion

data_tax <- read_xlsx("fileS4_microbe-data.xlsx")


data_tax <- data_tax[, 1:5]
colnames(data_tax) <- c("organism", "group", "plant_path", "plant_associate", "ecology")
p <- p %<+% data_tax

## Group tips by taxonomy
taxonomy_tree <- split(p$data$label, p$data$group)
taxonomy_tree <- groupOTU(tree, taxonomy_tree, group_name = "taxonomy")

## Group tips by fusion
fusion_tree <- split(p$data$label, p$data$`Domain Architecture`)
fusion_tree <- groupOTU(taxonomy_tree, fusion_tree, group_name = "fusion")

## Define colors for taxonomy coloring
colors1<-c("#a6cee3", "#1f78b4", "#b15928", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#cccc3d")
shapes <- c(32, 15, 15, 15, 15, 15, 15, 15, 15)

treedat <- tidytree::as_tibble(fusion_tree)

fusion_nodes <-
  treedat %>%
  filter(fusion != 0) %>%
  pull(node)

normal_nodes <-
  treedat %>%
  filter(fusion == 0) %>%
  pull(node)


### Collapse nodes
## This will produce a warning that can be ignored (Warning message:Unknown or uninitialised column: 'node'. )
tree_dat <- fusion_tree %>% 
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
  ggtree(fusion_tree, aes(color=taxonomy)) + 
  geom_tiplab(aes(subset = node %in% fusion_nodes & isTip, label = paste(label, fusion, sep = ", ")), color = "black", size = 2.5, fontface = "bold", hjust = -0.05) +
  geom_tiplab(aes(subset = node %in% normal_nodes & isTip, label = label), color = "black", size = 1.5, hjust = -0.05) +
  geom_point2(aes(subset = isTip, shape = fusion), color = "black", size = 2) +
  scale_color_manual(values = colors1) +
  geom_treescale(width = 0.5, linesize = 1, fontsize = 5, y = -5, x = 0) +
  scale_shape_manual(values = shapes) +
  xlim(0, 8) +
  theme(legend.position = "right")

## Collapse nodes
for(i in 1:length(nodes_to_collapse)){
  collapse_tree <- collapse(collapse_tree, node = nodes_to_collapse[i], clade_name=tree_labels[i])
}

## Plot collapsed tree w/ points
p2 <- collapse_tree + 
  geom_point2(aes(subset = (node %in% nodes_to_collapse)), size = 3, shape = 23) + 
  geom_nodelab(aes(subset = (node %in% nodes_to_collapse)), hjust = -0.05, size = 3) 
p2

## Edited for colors and aesthetics in Illustrator ##

