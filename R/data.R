library(Reach)

download_data <- function() {
  Reach::downloadOSFdata(repository='j67bv',
                         filelist=list('data'=c('performanceurpp.zip')), 
                         folder='data', 
                         unzip=TRUE, 
                         removezips=TRUE)
}

# Preprocess participant data
preprocess_file <- function(id) {
  # Read file for participant id
  filename <- sprintf('data/%s_performance.csv', id) 
  df <- read.csv(filename, stringsAsFactors = FALSE)
  
  # Extract relevant data 
  idx <- which(df$label == 'stl-target-rotation')
  response <- df$reachdeviation_deg[idx + 1] - df$reachdeviation_deg[idx - 1]
  rotation <- df$rotation[idx]
  
  # Exclude outliers and filter corresponding elements
  filtered_idx <- idx[response >= -60 & response <= 60]
  filtered_response <- response[response >= -60 & response <= 60]
  filtered_rotation <- rotation[response >= -60 & response <= 60]
  
  # Store into data frame
  output <- data.frame(rotation = filtered_rotation, 
                       response = filtered_response, 
                       trial = filtered_idx)
  output$participant <- id
  
  # Normalize data
  output$response[which(output$rotation > 0)] <- -1 * output$response[which(output$rotation > 0)]
  output$rotation[which(output$rotation < 0)] <- -1 * output$rotation[which(output$rotation < 0)]
  
  # Use aggregate to calculate the mean for each rotation value
  result_mean <- aggregate(response ~ rotation, data = output, mean)
  result_median <- aggregate(response ~ rotation, data = output, median)
  
  # Store both "output" and "result" data frames in a list
  processed_data <- list(output = output, result_mean = result_mean, result_median = result_median)
  
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

preprocess_all_files <- function(identifiers_list) {
  processed_data_list <- list()
  
  for (identifier in identifiers_list) {
    processed_data <- preprocess_file(identifier)
    processed_data_list[[identifier]] <- processed_data
  }
  
  return(processed_data_list)
}

plot_combined_data <- function(identifiers_list, pdf_filename) {
  pdf(pdf_filename, width = 8.5, height = 11)  # Open PDF device
  
  layout(mat=matrix(c(1:6),byrow=TRUE,nrow=3))
  
  for (identifier in identifiers_list) {
    processed_data <- preprocess_file(identifier)
    
    # Set up the plot
    plot(processed_data$output$rotation, processed_data$output$response,
         type = "p", pch = 16, col = "black",
         xlab = "Rotation in Degrees", ylab = "Average Response",
         main = identifier)
    
    # Add a dotted grey line at y = 0
    abline(h = 0, lty = 2, col = "grey")
    
    # Plot mean data points with red line
    lines(processed_data$result_mean$rotation, processed_data$result_mean$response,
          type = "l", col = "red")
    
    # Plot median data points with blue line
    lines(processed_data$result_median$rotation, processed_data$result_median$response,
          type = "l", col = "blue")
    
    # Close PDF device if this is the last participant
    if (identifier == identifiers_list[length(identifiers_list)]) {
      dev.off()  # Close PDF device
    }
  }
}


download_data()
preprocess_all_files(identifiers_list)
# Specify the PDF file name and path
pdf_filename <- "~/Desktop/participant_plots.pdf"  
plot_combined_data(identifiers_list, pdf_filename)

