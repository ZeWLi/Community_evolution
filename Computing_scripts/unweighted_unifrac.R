# The purpose of this script is to calculate unweighted unifrac distance and perform NMDS and ANOSIM analysis.
## load the libraries
library(picante)
library(phyloseq)
library(vegan)
library(phytools)

## load the data
### rooted phologenetic trees
tree = read.tree(file='/Users/lizewei/Desktop/lab.data/7.sample_ecology/final_file/submission/data/analysis_data/20230705_31marker_filt10_for_ecology.treefile')
rooted_tree = midpoint.root(tree)

### otu table
otu_matrix = as.matrix(read.delim('/Users/lizewei/Desktop/lab.data/7.sample_ecology/final_file/submission/data/analysis_data/otu_matrix.txt',sep = '\t'))

### group
group = read.delim('/Users/lizewei/Desktop/lab.data/7.sample_ecology/final_file/submission/data/analysis_data/20231212_sample_group_for_nst.txt',sep = '\t')

## unifrac calculation
### unweighted
out.unifrac = unifrac(otu_matrix,rooted_tree)
out.unifrac.matrix = as.matrix(out.unifrac)

### unifrac NMDS
dist_unifrac = as.dist(out.unifrac.matrix)
nmds_result <- metaMDS(dist_unifrac, k = 2)

nmds_group = merge(nmds_points,group,by='sample',sort=F)

nmds_result$stress # 0.1258345

### ANOSIM test
anosim_result = anosim(dist_unifrac,nmds_group$group)
plot(anosim_result)

anosim_df <- data.frame(
  Distance = c(anosim_result$dis.rank),
  Group = anosim_result$class.vec
)

