# The purpose of this script is to calculate the evolutionary rates of samples based on their branch lengths in a phylogenetic tree.
## import the libraries
library(phangorn)
library(ape)

## midpoint-root tree preparation 
tree = read.tree(file='/Users/lizewei/Desktop/lab.data/7.sample_ecology/final_file/submission/data/analysis_data/20230921_31marker_filt10_drop_wrong.fa.treefile')
tree_mid_rooted = midpoint(tree)

## calculate the branch length to the root
root_distances <- node.depth.edgelength(tree_mid_rooted)

all_node_numbers <- 1:(length(tree_mid_rooted$tip.label) + tree_mid_rooted$Nnode)

node_names <- ifelse(
  all_node_numbers <= length(tree_mid_rooted$tip.label),
  tree_mid_rooted$tip.label[all_node_numbers],
  paste0("Node_", all_node_numbers)
)

tip_length <- data.frame(
  Node = node_names,
  DistanceToRoot = root_distances[all_node_numbers]
)

## length of the leaves
tip_length_leaf <- tip_length[grepl(".pro.fa", tip_length$Node), ]
branch_length_map <- setNames(tip_length_leaf$DistanceToRoot, tip_length_leaf$Node)

## load the 3132 leaves rMAG otu
otu_as_matrix_df = as.data.frame(t(read.delim('/Users/lizewei/Desktop/lab.data/7.sample_ecology/final_file/submission/data/analysis_data/20250628_3132abundance_otu.csv',sep = ',',header = T,row.names = 1)))

## community evolutionary rate calculation
otu_as_matrix_df$evo_rate <- apply(otu_as_matrix_df, 1, function(row) {
  positive_count <- sum(row > 0)   
  if (positive_count == 0) return(0)   
  sum(row * branch_length_map[names(row)]) / sum(row) / positive_count  
})

evo_rate_df = data.frame(
  sample = rownames(otu_as_matrix_df),
  evolutionary_rate = otu_as_matrix_df$evo_rate
)
