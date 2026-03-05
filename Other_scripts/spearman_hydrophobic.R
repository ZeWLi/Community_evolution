# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript script.R <temperature_file> <directory> <output_file>")
}

temp_file <- args[1]
input_dir <- args[2]
output_file <- args[3]

# Function to read FASTA file
read_fasta <- function(file_path) {
  lines <- readLines(file_path)
  
  # Find header lines (starting with >)
  header_indices <- grep("^>", lines)
  
  if (length(header_indices) == 0) {
    return(data.frame())
  }
  
  records <- data.frame()
  
  for (i in 1:length(header_indices)) {
    # Get header
    header <- sub("^>", "", lines[header_indices[i]])
    
    # Get sequence
    if (i < length(header_indices)) {
      seq_lines <- lines[(header_indices[i] + 1):(header_indices[i + 1] - 1)]
    } else {
      seq_lines <- lines[(header_indices[i] + 1):length(lines)]
    }
    
    # Remove empty lines and concatenate sequence
    seq_lines <- seq_lines[seq_lines != ""]
    sequence <- paste(seq_lines, collapse = "")
    
    records <- rbind(records, 
                     data.frame(id = header, 
                                seq = sequence, 
                                stringsAsFactors = FALSE))
  }
  
  return(records)
}

## Parse the sample temperature
temp_data <- read.table(temp_file, header = TRUE, sep = "\t", na.strings = 'n.a.')
sample_temp <- setNames(as.numeric(temp_data[, ncol(temp_data)]), temp_data[, 1])

## Define hydropathy values (Kyte-Doolittle scale)
hydropathy_values <- c(
  'A' = 1.8, 'R' = -4.5, 'N' = -3.5, 'D' = -3.5, 'C' = 2.5, 
  'Q' = -3.5, 'E' = -3.5, 'G' = -0.4, 'H' = -3.2, 'I' = 4.5, 
  'L' = 3.8, 'K' = -3.9, 'M' = 1.9, 'F' = 2.8, 'P' = -1.6, 
  'S' = -0.8, 'T' = -0.7, 'W' = -0.9, 'Y' = -1.3, 'V' = 4.2
)

## Parse the directory
pwd <- getwd()

# Check if the input directory exists
if (!file.exists(input_dir) || !file.info(input_dir)$isdir) {
  stop("Input directory does not exist or is not a directory: ", input_dir)
}

# Get the list of items in the directory
all_items <- list.files(input_dir, full.names = FALSE)
dir_list <- c()

# Filter for actual directories
for (item in all_items) {
  full_path <- file.path(input_dir, item)
  if (file.exists(full_path) && file.info(full_path)$isdir) {
    dir_list <- c(dir_list, item)
  }
}

if (length(dir_list) == 0) {
  stop("No directories found in: ", input_dir)
}

setwd(input_dir)

results <- data.frame()

for (directory in dir_list) {
  # Check if the directory exists
  if (!file.exists(directory) || !file.info(directory)$isdir) {
    next
  }
  
  tryCatch({
    setwd(directory)
    
    ps_site_path <- "ps_site"
    if (!file.exists(ps_site_path) || !file.info(ps_site_path)$isdir) {
      setwd("..")
      next
    }
    
    setwd(ps_site_path)
    
    # Find FASTA files (common extensions)
    all_files <- list.files()
    fasta_files <- all_files[grepl("\\.(fasta|fa|fas|faa)$", all_files, ignore.case = TRUE)]
    
    if (length(fasta_files) == 0) {
      setwd("../..")
      next
    }
    
    for (file in fasta_files) {
      # Read FASTA file using custom function
      records <- read_fasta(file)
      
      if (nrow(records) == 0) next
      
      # Extract records matching sample names
      extracted_records <- data.frame()
      for (record_idx in 1:nrow(records)) {
        for (sample in names(sample_temp)) {
          if (grepl(sample, records$id[record_idx], fixed = TRUE)) {
            extracted_records <- rbind(extracted_records, 
                                       data.frame(sample = sample, 
                                                  seq = records$seq[record_idx],
                                                  stringsAsFactors = FALSE))
            break
          }
        }
      }
      
      # Remove duplicated records
      extracted_records <- unique(extracted_records)
      
      if (nrow(extracted_records) > 0) {
        seq_length <- nchar(extracted_records$seq[1])
        
        # Calculate GRAVY value for each sample
        sample_gravy_data <- data.frame()
        
        for (seq_idx in 1:nrow(extracted_records)) {
          sample_name <- extracted_records$sample[seq_idx]
          sequence <- extracted_records$seq[seq_idx]
          temp_val <- sample_temp[sample_name]
          
          if (!is.na(temp_val)) {
            # Calculate GRAVY value for the entire sequence
            hydropathy_sum <- 0
            valid_aa_count <- 0
            
            for (pos in 1:seq_length) {
              aa <- substr(sequence, pos, pos)
              if (aa != '-' && aa != 'X' && aa %in% names(hydropathy_values)) {
                hydropathy_sum <- hydropathy_sum + hydropathy_values[aa]
                valid_aa_count <- valid_aa_count + 1
              }
            }
            
            if (valid_aa_count > 0) {
              gravy <- hydropathy_sum / valid_aa_count
              sample_gravy_data <- rbind(sample_gravy_data,
                                       data.frame(temperature = temp_val,
                                                gravy = gravy,
                                                sample = sample_name,
                                                valid_positions = valid_aa_count))
            }
          }
        }
        
        # If multiple samples have the same temperature, take the average GRAVY value
        if (nrow(sample_gravy_data) > 0) {
          # Aggregate GRAVY values by temperature
          aggregated_data <- aggregate(gravy ~ temperature, 
                                     data = sample_gravy_data, 
                                     FUN = mean)
          
          # Count the number of samples for each temperature
          sample_counts <- aggregate(gravy ~ temperature, 
                                   data = sample_gravy_data, 
                                   FUN = length)
          names(sample_counts)[2] <- "sample_count"
          
          # Merge data
          final_data <- merge(aggregated_data, sample_counts, by = "temperature")
          
          if (nrow(final_data) > 2 && var(final_data$gravy) > 0) {
            tryCatch({
              # Analyze the relationship between GRAVY and temperature using Spearman correlation
              cor_result <- cor.test(final_data$temperature, final_data$gravy, method = 'spearman')
              
              rho_val <- cor_result$estimate  # Spearman correlation coefficient
              p_val <- cor_result$p.value     # p-value
              
              results <- rbind(results, 
                             data.frame(cluster = directory,
                                        gene_cluster = file,
                                        unique_temperatures = nrow(final_data),
                                        mean_gravy = mean(final_data$gravy),
                                        gravy_range = max(final_data$gravy) - min(final_data$gravy),
                                        spearman_rho = rho_val,
                                        p_value = p_val,
                                        total_samples = sum(final_data$sample_count)))
            }, error = function(e) {
              results <<- rbind(results, 
                              data.frame(cluster = directory,
                                         gene_cluster = file,
                                         unique_temperatures = nrow(final_data),
                                         mean_gravy = mean(final_data$gravy),
                                         gravy_range = max(final_data$gravy) - min(final_data$gravy),
                                         spearman_rho = NA,
                                         p_value = NA,
                                         total_samples = sum(final_data$sample_count)))
            })
          } else {
            # Case where data is insufficient
            results <- rbind(results, 
                           data.frame(cluster = directory,
                                      gene_cluster = file,
                                      unique_temperatures = nrow(final_data),
                                      mean_gravy = ifelse(nrow(final_data) > 0, mean(final_data$gravy), NA),
                                      gravy_range = ifelse(nrow(final_data) > 0, max(final_data$gravy) - min(final_data$gravy), NA),
                                      spearman_rho = NA,
                                      p_value = NA,
                                      total_samples = ifelse(nrow(final_data) > 0, sum(final_data$sample_count), 0)))
          }
        }
      }
    }
    
    setwd("../..")
    
  }, error = function(e) {
    # Ensure returning to the correct directory
    tryCatch(setwd(input_dir), error = function(e2) {})
  })
}

# Ensure returning to the original directory
setwd(pwd)

# Output the results
write.table(results, file = output_file, 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Analysis completed. Results saved to:", output_file, "\n")
