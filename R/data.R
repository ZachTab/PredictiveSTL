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


# STL Prediction Function
STLpredict <- function(slope, width, rotations) {
  
  return((rotations * slope) * (dnorm(rotations, mean=0, sd= width) / dnorm(0, mean=0, sd= width)))
}


# STL Error Function
STLerrors <- function(rotations, deviations, slope, width) {
  
  # Generate model predictions using STLpredict function
  model_predictions <- STLpredict(slope, width, rotations)
  
  # Calculate squared errors between model predictions and actual deviations
  squared_errors <- (deviations - model_predictions)^2
  
  # Compute the mean squared error (MSE)
  MSE <- mean(squared_errors)
  
  return(MSE)
}

# STL Grid Search Function
STLgridsearch <- function(rotations, deviations) {
  # Define parameter grids
  s_values <- seq(0.1, 1, length.out = 10)  # Example range for parameter s
  w_values <- seq(0, 60, length.out = 60)  # Example range for parameter w
  
  # Generate all combinations of parameters
  parameter_combinations <- expand.grid(s = s_values, w = w_values)
  
  # Initialize a data frame to store MSEs for each parameter set
  MSE_df <- data.frame(parameters = NA, MSE = NA)
  
  # Perform grid search and store MSEs in a data frame
  for (i in 1:nrow(parameter_combinations)) {
    s <- parameter_combinations$s[i]
    w <- parameter_combinations$w[i]
    MSE <- STLerrors(rotations, deviations, s, w)
    MSE_df <- rbind(MSE_df, data.frame(parameters = toString(c(s = s, w = w)), MSE = MSE))
  }
  
  # Sort the MSEs and select top 10 parameter sets
  top_10_parameter_sets <- MSE_df[order(MSE_df$MSE), ][1:10, ]
  
  return(top_10_parameter_sets)
}



