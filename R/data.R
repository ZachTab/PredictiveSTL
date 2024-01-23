library(Reach)

download_data <- function() {
  Reach::downloadOSFdata(repository='j67bv',
                         filelist=list('data'=c('performance.zip')), 
                         folder='data', 
                         unzip=TRUE, 
                         removezips=TRUE)
}


preprocess_file <- function(id) {
  
  # Read file for participant id
  filename <- sprintf('data/%s_performance.csv', id) 
  df <- read.csv(filename, stringsAsFactors = FALSE)
  
  # Extract relevant data 
  idx <- which(df$label == 'stl-target-rotation')
  response <- df$reachdeviation_deg[idx + 1] - df$reachdeviation_deg[idx - 1]
  rotation <- df$rotation[idx]
  
  # Store into data frame
  output <- data.frame(rotation, response, 'trial' = idx)
  output$participant <- id
  
  # Normalize data
  output$response[which(output$rotation > 0)] <- -1 *  output$response[which(output$rotation > 0)]
  output$rotation[which(output$rotation < 0)] <- -1 *  output$rotation[which(output$rotation < 0)]
  return(output)
  
}

# Setup a vector of id strings
csv_files <- list.files("data", pattern = "_performance.csv", full.names = TRUE)

# Create an empty list to store identifiers
identifiers_list <- list()

# Loop through each CSV file and extract the identifier
for (file in csv_files) {
  # Extract the identifier (xxxxxx) from the file name
  identifier <- gsub("_performance.csv", "", basename(file))
  
  # Add the identifier to the list
  identifiers_list[[identifier]] <- identifier
}


# Concatenate their output data into one data frame
preprocess_all_files <- function(identifiers_list) {
  # Create an empty list to store processed data frames
  processed_data_list <- list()
  
  # Loop through each identifier and apply the preprocess_file function
  for (identifier in identifiers_list) {
    # Apply the preprocess_file function and store the result in the processed_data_list
    processed_data <- preprocess_file(identifier)
    processed_data_list[[identifier]] <- processed_data
  }
  
  # Combine all the individual data frames into one
  combined_data <- do.call(rbind, processed_data_list)
  
  # Now, combined_data contains all the results in one singular data frame
  return(combined_data)
}

# Call the function with the identifiers list obtained earlier
combined_data <- preprocess_all_files(identifiers_list)


# Plot all of the combined data
plot(combined_data$rotation, combined_data$response)

