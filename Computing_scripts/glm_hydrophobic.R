## The purpose of this script is to test the response of ps site for hydrophobic amino acids conversion to temperature changes in protein sequences.
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

print(sample_temp)

## Define the hydrophobic and hydrophilic amino acids
hydrophobic <- c('A', 'C', 'I', 'L', 'M', 'F', 'V')
hydrophilic <- c('R', 'N', 'D', 'Q', 'E', 'G', 'H', 'K', 'P', 'S', 'T', 'W', 'Y')

# Create binary hydropathy values
hydropathy_binary <- c()
for (aa in hydrophobic) {
  hydropathy_binary[aa] <- 1
}
for (aa in hydrophilic) {
  hydropathy_binary[aa] <- 0
}

## Parse the directory
pwd <- getwd()
dir_list <- list.files(input_dir, full.names = FALSE)
dir_list <- dir_list[file.info(file.path(input_dir, dir_list))$isdir]

setwd(input_dir)

results <- data.frame()

for (directory in dir_list) {
  if (file.exists(directory) && file.info(directory)$isdir) {
    setwd(directory)
    
    ps_site_path <- "ps_site"
    if (file.exists(ps_site_path) && file.info(ps_site_path)$isdir) {
      setwd(ps_site_path)
      
      # Find FASTA files (common extensions)
      all_files <- list.files()
      fasta_files <- all_files[grepl("\\.(fasta|fa|fas|faa)$", all_files, ignore.case = TRUE)]
      
      for (file in fasta_files) {
        # Read FASTA file using custom function
        records <- read_fasta(file)
        
        if (nrow(records) == 0) next
        
        # Extract records matching sample names
        extracted_records <- data.frame()
        for (i in 1:nrow(records)) {
          for (sample in names(sample_temp)) {
            if (grepl(sample, records$id[i], fixed = TRUE)) {
              extracted_records <- rbind(extracted_records, 
                                         data.frame(sample = sample, 
                                                    seq = records$seq[i],
                                                    stringsAsFactors = FALSE))
              break
            }
          }
        }
        
        # Remove duplicated records
        extracted_records <- unique(extracted_records)
        
        if (nrow(extracted_records) > 0) {
          seq_length <- nchar(extracted_records$seq[1])
          
          for (pos in 1:seq_length) {
            test_data <- data.frame()
            
            for (j in 1:nrow(extracted_records)) {
              aa <- substr(extracted_records$seq[j], pos, pos)
              sample_name <- extracted_records$sample[j]
              
              if (aa != '-' && aa != 'X' && aa %in% names(hydropathy_binary)) {
                temp_val <- sample_temp[sample_name]
                if (!is.na(temp_val)) {
                  test_data <- rbind(test_data, 
                                     data.frame(temperature = temp_val,
                                                hydrophobic = hydropathy_binary[aa]))
                }
              }
            }
            
            # Remove rows with NA values
           # test_data <- test_data[complete.cases(test_data), ]
            
            if (nrow(test_data) > 2 && length(unique(test_data$hydrophobic)) > 1) {
              tryCatch({
                # Fit GLM with binomial family
                model <- glm(hydrophobic ~ temperature, 
                             data = test_data, 
                             family = binomial())
                print(test_data)
                
                coef_val <- coef(model)[2]  # coefficient for temperature
                p_val <- summary(model)$coefficients[2, 4]  # p-value for temperature
                
                results <- rbind(results, 
                                 data.frame(cluster = directory,
                                            gene_cluster = file,
                                            position = pos,
                                            coefficient = coef_val,
                                            p_value = p_val,
                                            data_point = nrow(test_data)))
              }, error = function(e) {
                results <<- rbind(results, 
                                  data.frame(cluster = directory,
                                             gene_cluster = file,
                                             position = pos,
                                             coefficient = NA,
                                             p_value = NA,
                                             data_point = nrow(test_data)))
              })
            } else {
              results <- rbind(results, 
                               data.frame(cluster = directory,
                                          gene_cluster = file,
                                          position = pos,
                                          coefficient = NA,
                                          p_value = NA,
                                          data_point = nrow(test_data)))
            }
          }
        }
      }
      setwd("..")
    }
    setwd("..")
  }
}

# Output the results
setwd(pwd)

write.table(results, file = output_file, 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Analysis completed. Results saved to:", output_file, "\n")