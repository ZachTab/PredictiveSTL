library(Reach)
library(optimx)

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

# Function to preprocess data and write to CSV
RawData_CSV <- function() {
  # Initialize a data frame to store all processed data
  all_data <- data.frame()
  
  # Setup a vector of id strings
  csv_files <- list.files("data", pattern = "_performance.csv", full.names = TRUE)
  
  # Loop through each CSV file and preprocess data
  for (file in csv_files) {
    # Extract the identifier (xxxxxx) from the file name
    identifier <- gsub("_performance.csv", "", basename(file))
    
    # Read file for participant id
    df <- read.csv(file, stringsAsFactors = FALSE)
    
    # Extract relevant data
    idx <- which(df$label == 'stl-target-rotation')
    response <- df$reachdeviation_deg[idx + 1] - df$reachdeviation_deg[idx - 1]
    rotation <- df$rotation[idx]
    
    # Exclude outliers and filter corresponding elements
    filtered_idx <- idx[response >= -60 & response <= 60]
    filtered_response <- response[response >= -60 & response <= 60]
    filtered_rotation <- rotation[response >= -60 & response <= 60]
    
    # Store into data frame
    output <- data.frame(Participant_ID = rep(identifier, length(filtered_idx)),
                         Rotation = filtered_rotation,
                         Response = filtered_response)
    
    # Normalize data
    output$Response[output$Rotation > 0] <- -1 * output$Response[output$Rotation > 0]
    output$Rotation[output$Rotation < 0] <- -1 * output$Rotation[output$Rotation < 0]
    
    # Add data for the current participant to the all_data data frame
    all_data <- rbind(all_data, output)
  }
  
  # Write processed data to CSV file
  write.csv(all_data, "raw_data.csv", row.names = FALSE)
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


# Specify the PDF file name and path
pdf_filename <- "~/Desktop/participant_plots.pdf"  
plot_combined_data(identifiers_list, pdf_filename)



# STL Prediction Function
STLpredict <- function(par, rotations) {
  
  return((rotations * par['s']) * (dnorm(rotations, mean=0, sd= par['w']) / dnorm(0, mean=0, sd= par['w'])))
}


# STL Error Function
STLerrors <- function(par, rotations, deviations) {
  
  # Generate model predictions using STLpredict function
  model_predictions <- STLpredict(par, rotations)
  
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
  
  MSE <- apply(parameter_combinations, FUN=STLerrors, MARGIN=c(1), rotations = rotations, deviations = deviations)

  top_parameters <- parameter_combinations[order(MSE)[1:10], ]
  
  return(top_parameters)
}


STLoptimization <- function(top_parameters, rotations, deviations) { # run optimx on the best starting positions:
  allfits <- do.call("rbind",
                     apply( top_parameters,
                            MARGIN=c(1),
                            FUN=optimx::optimx,
                            fn=STLerrors,
                            method = "L-BFGS-B", 
                            lower = c(0.1, 1), 
                            upper = c(1, 60),
                            rotations = rotations,
                            deviations = deviations
                     ))
  
  win <- allfits[order(allfits$value)[1], ]
  return(unlist(win[1:2]))
}


STL <- function(processed_data) {
  # 
  # Preprocess data for the given identifier
  # processed_data <- preprocess_file(identifier)
  
  # Extract rotations and deviations
  rotations <- processed_data$output$rotation
  deviations <- processed_data$output$response
  
  # Run grid search to find best parameters
  top_parameters <- STLgridsearch(rotations, deviations)
  
  cat("Top 10 Parameters:\n")
  print(top_parameters)
  
  # run optimx on the best starting positions:
  winning_par <- STLoptimization(top_parameters, rotations, deviations)
  
  return(winning_par)
  
}


STLIndividualFits <- function(identifiers_list) {
  # Create an empty data frame to store the best fits for all participants
  all_fits_df <- data.frame(participant = character(),
                            slope = numeric(),
                            width = numeric(),
                            stringsAsFactors = FALSE)
  
  # Loop through each participant
  for (identifier in identifiers_list) {
    # Preprocess data for the current participant
    processed_data <- preprocess_file(identifier)
    
    # Extract rotations and deviations
    rotations <- processed_data$output$rotation
    deviations <- processed_data$output$response
    
    # Run grid search to find best parameters
    top_parameters <- STLgridsearch(rotations, deviations)
    
    # Run optimization on the best starting positions
    winning_par <- STLoptimization(top_parameters, rotations, deviations)
    
    # Add the participant id, slope, and width to the data frame
    all_fits_df[nrow(all_fits_df) + 1, ] <- c(identifier, winning_par[1], winning_par[2])
  }
  
  # Set column names
  colnames(all_fits_df) <- c("participant", "slope", "width")
  
  return(all_fits_df)
}

# Function to perform STL grid search and optimization on raw data
STLfit <- function(participant, data_file) {
  # Read raw data from CSV
  data <- read.csv(data_file)
  
  # If participant is a single ID, convert it to a list
  if (!is.list(participant)) {
    participant <- list(participant)
  }
  
  # Extract data for all participants
  all_participant_data <- lapply(participant, function(id) subset(data, Participant_ID == id))
  all_data <- do.call(rbind, all_participant_data)
  
  # Extract rotations and deviations from combined data
  rotations <- all_data$Rotation
  deviations <- all_data$Response
  
  # Run grid search to find best parameters
  top_parameters <- STLgridsearch(rotations, deviations)
  
  # Run optimization on the best parameters
  winning_par <- STLoptimization(top_parameters, rotations, deviations)
  
  # Create a data frame for the results
  results_df <- data.frame(Participant_ID = "Combined",
                           Slope = winning_par[1],
                           Width = winning_par[2],
                           stringsAsFactors = FALSE)
  
  return(results_df)
}


# Graph Predicted Responses using Winning Parameter

# # Extract the optimized parameters
# opt_slope <- best_parameters$slope
# opt_width <- best_parameters$width
# 
# # Run STLpredict with optimized parameters
# predicted_values <- STLpredict(opt_slope, opt_width, rotations)
# 
# # Open a PDF
# layout(mat = matrix(c(1:6), byrow = TRUE, nrow = 3))
# pdf_filename <- paste0(identifier, "_predict.pdf")
# pdf(pdf_filename, width = 8.5, height = 11)
# 
# # Plot Observed data
# plot(rotations, deviations, type = "p", pch = 16, col = "black",
#      xlab = "Rotation in Degrees", ylab = "Average Response",
#      main = identifier)
# 
# # Plot Mean and Median Lines
# abline(h = 0, lty = 2, col = "grey")
# lines(processed_data$result_mean$rotation, processed_data$result_mean$response, type = "l", col = "red")
# lines(processed_data$result_median$rotation, processed_data$result_median$response, type = "l", col = "blue")
# 
# # Plot STL Predicted Line with optimized parameters
# stl_response <- STLpredict(opt_slope, opt_width, rotations = c(1,5,10,15,20,25,30,35,40,45))
# lines(x = c(1,5,10,15,20,25,30,35,40,45), stl_response, type = 'l', col = 'black')
# 
# dev.off()
# 
# cat("PDF file saved as:", pdf_filename, "\n")
# 
# # Print PDF file path
# pdf_full_path <- normalizePath(pdf_filename)
# cat("PDF file saved at:", pdf_full_path, "\n")
# 
# return(predicted_values)
    
    
    
    
    
    