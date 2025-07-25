# The purpose of the script is to calculate GRAVY value of community-level amino acid composition
library(ggplot2)

aa_comp = read.delim('/Users/lizewei/Desktop/lab.data/7.sample_ecology/final_file/submission/data/analysis_data/20240406_pro_usage_proportion.txt', header = TRUE, sep = '\t')

# Kyte-Doolittle scale
hydropathy_values <- c(
  A = 1.8, R = -4.5, N = -3.5, D = -3.5, C = 2.5, 
  Q = -3.5, E = -3.5, G = -0.4, H = -3.2, I = 4.5, 
  L = 3.8, K = -3.9, M = 1.9, F = 2.8, P = -1.6, 
  S = -0.8, T = -0.7, W = -0.9, Y = -1.3, V = 4.2
)

# define the calculation method
gravy_values <- apply(aa_comp, 1, function(row) {
  total_aa <- sum(as.numeric(row[names(hydropathy_values)]), na.rm = TRUE)
  gravy_value <- sum(sapply(names(hydropathy_values), function(aa) {
    if (aa %in% colnames(aa_comp)) {
      return((as.numeric(row[[aa]]) / total_aa) * hydropathy_values[aa])
    } else {
      return(0)
    }
  }), na.rm = TRUE)
  
  return(gravy_value)  
})

## calculate GRAVY value
output_df <- data.frame(sample = aa_comp$sample, GRAVY = gravy_values)

output_group = merge(output_df,group,by = 'sample')
