# Create qgraph to demonstrate shared clones between subsets 

install.packages("qgraph")
install.packages("DescTools")
install.packages("ggpattern")

library(qgraph)
library(dplyr)
library(readr)
library(ggplot2)
library(alakazam)
library(RColorBrewer)  # Load the RColorBrewer package
library(igraph)
library(shazam)
library(DescTools)
library(ggpattern)

# read in file 
HCtsv_DUPCOUNT <- read_tsv("Data/Reprocessed_HC_Combined_filtered_dist_ham_clone-pass_germ-pass_DUPCOUNT.tsv")

# identify and remove chimeric sequences 
is_chimeric <- slideWindowDb(
  HCtsv_DUPCOUNT,
  sequenceColumn = "sequence_alignment", 
  germlineColumn = "germline_alignment_d_mask", 
  mutThres = 6, 
  windowSize = 10, 
  nproc = 8
)

is_chimeric <- data.frame(is_chimeric) # change chimeric indexes to dataframe 
HCtsv_DUPCOUNT_non_chimeric <- HCtsv_DUPCOUNT[!is_chimeric, ] # filter to non-chimeric sequences 

# count clones 
clones <- countClones(HCtsv_DUPCOUNT_non_chimeric, group=c("subset")) # group subsets separately 
clones_clone_id <- subset(clones, select = -c(seq_freq)) # remove seq_freq column

# subset clones by cell type 
## ASC clones
ASC_clones <- subset(clones_clone_id, subset == "ASC", select = -c(subset)) 
colnames(ASC_clones)[-1] <- paste0("ASC_", colnames(ASC_clones)[-1])
## CD21lo clones 
CD21lo_clones <- subset(clones_clone_id, subset == "CD21lo", select = -c(subset))
colnames(CD21lo_clones)[-1] <- paste0("CD21lo_", colnames(CD21lo_clones)[-1])
## FOB clones 
FOB_clones <- subset(clones_clone_id, subset == "FOB", select = -c(subset))
colnames(FOB_clones)[-1] <- paste0("FOB_", colnames(FOB_clones)[-1])
## GC clones 
GC_clones <- subset(clones_clone_id, subset == "GC", select = -c(subset))
colnames(GC_clones)[-1] <- paste0("GC_", colnames(GC_clones)[-1])

# merge subset clone tables to generate matrix
subset_merge <- merge(ASC_clones, CD21lo_clones, by = "clone_id", all = TRUE)
subset_merge <- merge(subset_merge, FOB_clones, by = "clone_id", all = TRUE)
subset_merge <- merge(subset_merge, GC_clones, by = "clone_id", all = TRUE)

colnames(subset_merge)[-1] <- c("ASC", "CD21lo", "FOB", "GC") # add column names 

subset_merge <- replace(subset_merge, is.na(subset_merge), 0) # replace any NA values by 0

# remove row names 
rownames(subset_merge) <- subset_merge$clone_id 
subset_merge <- subset_merge[, -1]

# selecting 20 seq as the threshold for presenting a positive connection between clones 
subset_merge <- subset_merge %>%
  mutate(across(everything(), ~ ifelse(. > 20, 1, 0)))

# convert to adjacency matrix 
adj_matrix <- as.matrix(subset_merge)

# Identify the non-zero elements in the adjacency matrix
indices <- which(adj_matrix != 0, arr.ind = TRUE)

# Create the edgelist from the indices
edgelist <- data.frame(from = rownames(adj_matrix)[indices[, 1]],
                       to = colnames(adj_matrix)[indices[, 2]])
edgelist_shared <- edgelist %>% # filter to shared clones 
  group_by(from) %>% 
  filter(n_distinct(to) >1 ) %>% 
  ungroup() %>%
  data.frame()

# set node sizes
to_column_nodes <- unique(c(edgelist_shared[,1], edgelist_shared[,2])) # create index of unique values  
node_sizes <- rep(1, length(to_column_nodes)) # set default node size to 1 
node_sizes[which(to_column_nodes %in% c("ASC", "CD21lo", "GC", "FOB"))] <- 15 # set subset nodes to size 15
for(i in which(node_sizes==1)) {
  node_sizes[i]=Gmean(clones_clone_id$seq_count[which(clones_clone_id$clone_id==to_column_nodes[i])])
  node_sizes[i]=node_sizes[i]/20
} # scale remaining node sizes to clone sizes  

# set node label sizes
label_sizes <- rep(0.00000001, length(to_column_nodes))
label_sizes[which(to_column_nodes %in% c("ASC", "CD21lo", "GC", "FOB"))] <- 1

# set node colors 
node_colors <- rep("white", length(to_column_nodes))
node_colors[which(to_column_nodes %in% c("ASC"))] <- "#69BED8"
node_colors[which(to_column_nodes %in% c("FOB"))] <- "#E7899E"
node_colors[which(to_column_nodes %in% c("CD21lo"))] <- "#CF76D6"
node_colors[which(to_column_nodes %in% c("GC"))] <- "#F8AC6D"

## pull clone_id's that are shared by three subsets 
shared_by_three <- edgelist_shared %>% 
  group_by(from) %>% 
  filter(n()==3) %>% 
  pull(from) %>%
  unique()

## pull clone_id's that are shared by four subsets 
shared_by_four <- edgelist_shared %>% 
  group_by(from) %>% 
  filter(n()==4) %>% 
  pull(from) %>%
  unique()

## set different node colors for shared clones  
node_colors[which(to_column_nodes %in% shared_by_three)] <- "#CBCBCB"
node_colors[which(to_column_nodes %in% shared_by_four)] <- "#757575"

# set edge colors 
edge_colors <- rep("black", length(edgelist_shared$to))
edge_colors[which(edgelist_shared$to %in% c("ASC"))] <- "#08768A"
edge_colors[which(edgelist_shared$to %in% c("FOB"))] <- "#E7899E"
edge_colors[which(edgelist_shared$to %in% c("CD21lo"))] <- "#7B4C96"
edge_colors[which(edgelist_shared$to %in% c("GC"))] <- "#E8952A"

#graph
qgraph(edgelist_shared, vsize = node_sizes, directed = FALSE, open = TRUE, color = node_colors, labels = TRUE, label.cex = label_sizes, repulsion = 20, esize = 1, edge.color = edge_colors, edge.width = 1.2)





                                           