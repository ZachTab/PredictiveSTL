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
  
  # Use aggregate to calculate the mean for each rotation value
  result <- aggregate(response ~ rotation, data = output, mean)
  
  # Store both "output" and "result" data frames in a list
  processed_data <- list(output = output, result = result)
  
  return(processed_data)
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

# Function to automate and preprocess all participant data
preprocess_all_files <- function(identifiers_list) {
  # Loop through each identifier and apply the preprocess_file function
  for (identifier in identifiers_list) {
    preprocess_file(identifier)
  }
}

# Function to plot all participant raw/averaged data
plot_combined_data <- function(identifiers_list) {
  rows <- ceiling(length(identifiers_list) / 2)
  par(mfrow=c(rows, 2), mar=c(2, 2, 2, 2))
  
  for (i in seq_along(identifiers_list)) {
    identifier <- identifiers_list[[i]]
    processed_data <- preprocess_file(identifier)
    
  
    # Plot "output" data points
    plot(processed_data$output$rotation, processed_data$output$response,
           type = "p", pch = 16, col = "black",
           xlab = "Rotation in Degrees", ylab = "Average Response",
           main = identifier)
      
      
    # Plot "result" data points with lines connecting them
    lines(processed_data$result$rotation, processed_data$result$response,
              type = "l", col = "red")
    
    
  }
}

# To call preprocess_all_files and plot_combined_data:
# preprocess_all_files(identifiers_list)
# plot_combined_data(identifiers_list)