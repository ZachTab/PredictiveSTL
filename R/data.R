library(Reach)
library(optimx)
library(dplyr)
library(ggplot2)

Download_data <- function() {
  Reach::downloadOSFdata(repository='6m24e',
                         filelist=list('data'=c('data.zip')), 
                         folder='data', 
                         unzip=TRUE, 
                         removezips=TRUE)
}

# Preprocess participant data
Preprocess_file <- function(id) {
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
RawData_CSV <- function(data_portion = "all", phase, file_name) {
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
    
    if (phase == "STL") {
      # Extract relevant data for STL phase
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
      
    } else if (phase %in% c("Short", "Short1", "Short2")) {
      # Extract relevant data for Short phase
      short_indices <- which(grepl("^short", df$label))
      
      if (length(short_indices) < 50) {
        warning("Not enough short trials found for participant: ", identifier)
        next
      }
      
      if (phase == "Short1") {
        short_start <- short_indices[1]
        short_end <- tail(which(df$label == "short-washout" & lead(df$label) == "short-acquire"), 1)
      } else if (phase == "Short2") {
        short_start <- head(which(df$label == "short-acquire" & lag(df$label) == "short-washout"), 1)
        short_end <- short_indices[length(short_indices)]
      } else {
        short_start <- short_indices[1]
        short_end <- short_indices[length(short_indices)]
      }
      
      if (is.na(short_start) || is.na(short_end)) {
        warning("Could not find valid start and end indices for participant: ", identifier)
        next
      }
      
      short_phase_indices <- short_start:short_end
      
      # Store into data frame
      output <- data.frame(Participant_ID = rep(identifier, length(short_phase_indices)),
                           Label = df$label[short_phase_indices],
                           Rotation = df$rotation[short_phase_indices],
                           Reach_Deviation = df$reachdeviation_deg[short_phase_indices])
      
    } else {
      stop("Invalid phase. Please choose 'STL', 'Short', 'Short1', or 'Short2'.")
    }
    
    # Debug print to check the number of rows before subsetting
    cat("Total rows for", identifier, "before subsetting:", nrow(output), "\n")
    
    # Subset the data based on the specified portion
    if (data_portion == "first_half") {
      output <- output[1:floor(nrow(output) / 2), ]
    } else if (data_portion == "last_half") {
      output <- output[(floor(nrow(output) / 2) + 1):nrow(output), ]
    }
    
    # Debug print to check the number of rows after subsetting
    cat("Total rows for", identifier, "after subsetting:", nrow(output), "\n")
    
    # Add data for the current participant to the all_data data frame
    all_data <- rbind(all_data, output)
  }
  
  # Write processed data to CSV file
  write.csv(all_data, file_name, row.names = FALSE)
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

Preprocess_all_files <- function(identifiers_list) {
  processed_data_list <- list()
  
  for (identifier in identifiers_list) {
    processed_data <- Preprocess_file(identifier)
    processed_data_list[[identifier]] <- processed_data
  }
  
  return(processed_data_list)
}

# STL Prediction Function
STLPredict <- function(par, rotations) {
  
  return((rotations * par['s']) * (dnorm(rotations, mean=0, sd= par['w']) / dnorm(0, mean=0, sd= par['w'])))
}

# Function to establish & return predicted adaptions for given fit data and rotations  
STLPredictions <- function(fits_data, rotations) {
  results <- data.frame(participant = character(), rotations = numeric(), predictions = numeric(), stringsAsFactors = FALSE)
  
  fits_data$slope <- as.numeric(fits_data$slope)
  fits_data$width <- as.numeric(fits_data$width)
  
  for (i in 1:nrow(fits_data)) {
    participant <- fits_data$participant[i]
    par <- c(s = fits_data$slope[i], w = fits_data$width[i])
    predictions <- STLPredict(par = par, rotations = rotations)
    
    participant_results <- data.frame(participant = rep(participant, length(rotations)),
                                      rotations = rotations,
                                      predictions = predictions)
    
    results <- rbind(results, participant_results)
  }
  
  return(results)
}    


# STL Error Function
STLErrors <- function(par, rotations, deviations) {
  
  # Generate model predictions using STLpredict function
  model_predictions <- STLPredict(par, rotations)
  
  # Calculate squared errors between model predictions and actual deviations
  squared_errors <- (deviations - model_predictions)^2
  
  # Compute the mean squared error (MSE)
  MSE <- mean(squared_errors)
  
  return(MSE)
}

# STL Grid Search Function
STLGridsearch <- function(rotations, deviations) {
  # Define parameter grids
  s_values <- seq(0.1, 1, length.out = 10)  # Example range for parameter s
  w_values <- seq(0, 60, length.out = 60)  # Example range for parameter w
  
  # Generate all combinations of parameters
  parameter_combinations <- expand.grid(s = s_values, w = w_values)
  
  MSE <- apply(parameter_combinations, FUN=STLErrors, MARGIN=c(1), rotations = rotations, deviations = deviations)

  top_parameters <- parameter_combinations[order(MSE)[1:10], ]
  
  return(top_parameters)
}


STLOptimization <- function(top_parameters, rotations, deviations) { # run optimx on the best starting positions:
  allfits <- do.call("rbind",
                     apply( top_parameters,
                            MARGIN=c(1),
                            FUN=optimx::optimx,
                            fn=STLErrors,
                            method = "L-BFGS-B", 
                            lower = c(0.1, 1), 
                            upper = c(1, 60),
                            rotations = rotations,
                            deviations = deviations
                     ))
  
  win <- allfits[order(allfits$value)[1], ]
  return(unlist(win[1:2]))
}

# Returns individual STL Fit parameters for a list of participants
STLIndividualFits <- function(identifiers_list, data_file = NULL) {
  # Create an empty data frame to store the best fits for all participants
  all_fits_df <- data.frame(participant = character(),
                            slope = numeric(),
                            width = numeric(),
                            stringsAsFactors = FALSE)
  
  # Check if a data_file argument is provided
  if (!is.null(data_file)) {
    # Read the data file
    data <- read.csv(data_file, stringsAsFactors = FALSE)
    
    # Loop through each participant
    for (identifier in identifiers_list) {
      # Extract data for the current participant from the data file
      participant_data <- data[data$Participant_ID == identifier, ]
      
      # Extract rotations and deviations
      rotations <- participant_data$Rotation
      deviations <- participant_data$Response
      
      # Run grid search to find best parameters
      top_parameters <- STLGridsearch(rotations, deviations)
      
      # Run optimization on the best starting positions
      winning_par <- STLOptimization(top_parameters, rotations, deviations)
      
      # Add the participant id, slope, and width to the data frame
      all_fits_df[nrow(all_fits_df) + 1, ] <- c(identifier, winning_par[1], winning_par[2])
    }
  } else {
    # Loop through each participant
    for (identifier in identifiers_list) {
      # Preprocess data for the current participant
      processed_data <- Preprocess_file(identifier)
      
      # Extract rotations and deviations
      rotations <- processed_data$output$rotation
      deviations <- processed_data$output$response
      
      # Run grid search to find best parameters
      top_parameters <- STLGridsearch(rotations, deviations)
      
      # Run optimization on the best starting positions
      winning_par <- STLOptimization(top_parameters, rotations, deviations)
      
      # Add the participant id, slope, and width to the data frame
      all_fits_df[nrow(all_fits_df) + 1, ] <- c(identifier, winning_par[1], winning_par[2])
    }
  }
  
  # Set column names
  colnames(all_fits_df) <- c("participant", "slope", "width")
  
  return(all_fits_df)
}


# Returns a single (or weighted) set of STL Fit parameters for a set of participant data
STLFit <- function(participant, data_file) {
  # Read raw data from CSV
  data <- read.csv(data_file)
  
  # If participant is a single ID, convert it to a list
  if (!is.vector(participant)) {
    participant <- vector(participant)
  }
  
  # Extract data for all participants
  all_data <- NA
  for (id in participant) {
    if (is.data.frame(all_data)) {
      all_data <- rbind(all_data, data[which(data$Participant_ID == id),])
    } else {
      all_data <- data[which(data$Participant_ID == id),]
    }
  }
  
  # Extract rotations and deviations from combined data
  rotations <- all_data$Rotation
  deviations <- all_data$Response
  
  # Run grid search to find best parameters
  top_parameters <- STLGridsearch(rotations, deviations)
  
  # Run optimization on the best parameters
  winning_par <- STLOptimization(top_parameters, rotations, deviations)
  
  # Create a data frame for the results
  results_df <- data.frame(Participant_ID = "Combined",
                           Slope = winning_par[1],
                           Width = winning_par[2],
                           stringsAsFactors = FALSE)
  
  return(results_df)
}

# Baselines and derives Exponential Model Fits from participant Fixed Rotation data 
STLExponentialFits <- function(ids, phase, data_file = NULL) {
  if (!is.vector(ids)) {
    ids <- c(ids)  # Ensure ids is a vector
  }
  
  if (!is.null(data_file)) {
    df <- read.csv(data_file, stringsAsFactors = FALSE)
  }
  
  results <- data.frame(participant = character(), rotation = numeric(), lambda = numeric(), N0 = numeric(), stringsAsFactors = FALSE)
  
  for (id in ids) {
    cat("Processing participant:", id, "\n")
    
    if (is.null(data_file)) {
      filename <- sprintf('data/%s_performance.csv', id)
      df <- read.csv(filename, stringsAsFactors = FALSE)
    } else {
      participant_df <- df[df$Participant_ID == id, ]
    }
    
    if (phase == "long") {
      # Extract baseline data
      baseline_indices <- which(participant_df$label == 'fixed-rotation-baseline')
      if (length(baseline_indices) < 60) {
        warning("Not enough baseline trials for participant: ", id)
        next
      }
      bias <- median(participant_df$reachdeviation_deg[baseline_indices[11:60]], na.rm = TRUE)
      learning_indices <- which(participant_df$label == 'fixed-rotation')
      learning <- participant_df$reachdeviation_deg[learning_indices] - bias
      
      # Derive Exponential Fits
      fit <- Reach::exponentialFit(signal = learning)
      
      # Initialize a data frame to store Exponential Fit data
      participant_result <- data.frame(participant = id, rotation = NA, lambda = fit['lambda'], N0 = fit['N0'])
      results <- rbind(results, participant_result)
      
    } else if (phase == "short") {
      # Identify the indices for the 'short' phases
      short_indices <- which(grepl("^short-", participant_df$Label))
      cat("Short indices for participant", id, ":", short_indices, "\n")
      
      if (length(short_indices) < 90) {
        warning("Not enough short trials found for participant: ", id)
        next
      }
      
      # Extract baseline data
      baseline_indices <- which(participant_df$Label == 'short-baseline')
      if (length(baseline_indices) < 10) {
        warning("Not enough short-baseline trials found for participant: ", id)
        next
      }
      bias <- median(participant_df$Reach_Deviation[baseline_indices], na.rm = TRUE)
      
      # Calculate learning for all short-rotated trials
      short_rotation_indices <- which(participant_df$Label == 'short-rotated')
      if (length(short_rotation_indices) == 0) {
        warning("No short-rotation trials found for participant: ", id)
        next
      }
      
      # Organize results by rotation values
      unique_rotations <- unique(participant_df$Rotation[short_rotation_indices])
      for (rotation_value in unique_rotations) {
        rotation_indices <- short_rotation_indices[participant_df$Rotation[short_rotation_indices] == rotation_value]
        learning <- (sign(rotation_indices) * -1) * (participant_df$Reach_Deviation[rotation_indices] - bias)
        
        if (any(is.na(learning))) {
          warning("NA values found in learning data for participant: ", id)
          learning <- na.omit(learning)
        }
        
        if (length(learning) == 0) {
          warning("No valid learning data found for participant: ", id)
          next
        }
        
        # Derive Exponential Fits
        fit <- Reach::exponentialFit(signal = learning)
        
        # Initialize a data frame to store Exponential Fit data
        participant_result <- data.frame(participant = id, rotation = rotation_value, lambda = fit['lambda'], N0 = fit['N0'])
        results <- rbind(results, participant_result)
      }
      
    } else {
      stop("Invalid phase. Please choose 'long' or 'short'.")
    }
  }
  
  return(results)
}
    

# Obtains graph-able Exponential Model values for participants
STLExponentialModel <- function(participants_results, mode = 'learning', setN0 = NULL, points, rotation = NULL) {
  
  participants_results$lambda <- as.numeric(participants_results$lambda)
  participants_results$N0 <- as.numeric(participants_results$N0)
  
  # Filter participants' results based on the rotation value if provided
  if (!is.null(rotation)) {
    participants_results <- participants_results[participants_results$rotation == rotation, ]
  }
  
  # Initialize an empty data frame to hold all the plot data
  all_plot_data <- data.frame()
  
  # Iterate over each participant's results
  for (i in 1:nrow(participants_results)) {
    participant_result <- participants_results[i, ]
    lambda <- participant_result$lambda
    N0 <- participant_result$N0
    participant <- participant_result$participant
    
    # Generate the time points (for example, 0 to 'points')
    time_points <- seq(0, points, by = 1)
    
    # Generate the model values using the exponentialModel function
    model_output <- Reach::exponentialModel(par = c('lambda' = lambda, 'N0' = N0), timepoints = time_points, mode = mode, setN0 = setN0)
    
    # Create a data frame for the current participant's plot data
    plot_data <- data.frame(time = model_output$trial, value = model_output$output, participant = participant)
    
    # Combine with the overall plot data
    all_plot_data <- rbind(all_plot_data, plot_data)
  }
  
  return(all_plot_data)
}
    
# Function to plot a participant's observed STL reach data along with average & predicted deviation lines
PlotReachdata <- function(identifiers_list, filename, predictions_data) {
  
  dir_path <- "~/Desktop/Masters/Thesis/PredictiveSTL/Plots/Raw Data Plots/"
  
  pdf_name <- paste0(dir_path, filename, "_Reach_Data.pdf")
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  pdf(pdf_name, width = 8.5, height = 11)  # Open PDF device
  
  layout(mat=matrix(c(1:6),byrow=TRUE,nrow=3))
  
  data_file <- read.csv("~/Desktop/Masters/Thesis/PredictiveSTL/Data Files/raw_data.csv", stringsAsFactors = FALSE)
  
  # Loop through each participant
  for (identifier in identifiers_list) {
    # Extract data for the current participant from the data file
    participant_data <- data_file[data_file$Participant_ID == identifier, ]
    
    # Extract rotations and deviations
    rotations <- participant_data$Rotation
    deviations <- participant_data$Response
    
    # Calculate Mean Deviation per rotation value
    mean_deviation <- participant_data %>%
      group_by(Rotation) %>%
      summarize(mean_deviation = mean(Response))
    
    # Extract predicted data for the current participant
    participant_predictions <- predictions_data %>%
      filter(participant == identifier)
    
    p <- ggplot(participant_data, aes(x = Rotation, y = Response)) +
      geom_point(color = "black") +
      
      # Plot mean data points with red line
      geom_line(data = mean_deviation, aes(x = Rotation, y = mean_deviation), color = "red", alpha = 0.6) +
      
      # Plot predicted data points with blue line
      geom_line(data = participant_predictions, aes(x = rotations, y = predictions), color = "blue", alpha = 0.6) +
      
      # Add a dotted grey line at y = 0
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey") +    
      
      labs(title = identifier,
           x = "Rotation in Degrees(°)", 
           y = "Deviation in Degrees(°)") +
      theme_minimal() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    # Print the plot to the PDF device
    print(p)
   
  }
  
  dev.off()  # Close PDF device

}

# Function to plot exponential model values 
STLExponentialPlot <- function(data, rot, vert_line = NULL, trials = NULL){
  
  main_title <- paste(paste0(rot,"°"), "Adaptation Over", trials, "Trials")
  
  # Allows to plot for a certain number of trials
  if (is.null(trials)){
  
  ggplot(data, aes(x = time, y = value, color = participant)) +
    geom_line(alpha = 0.7) +
      labs(title = main_title,
           x = "Trials",
           y = "Model Value") +
      theme_minimal() +
      guides(color = "none") +
      geom_vline(xintercept = vert_line, linetype = "dotted", color = "black")
    
    
  } else{
    # Add the trials column, resetting at 0 for each new participant
      data <- data %>%
        group_by(participant) %>%
        mutate(trials = row_number() - 1) %>%
        ungroup()
      
      # Filter data to include only up to the specified number of trials
      data <- data %>%
        filter(trials <= trials)
      
      ggplot(data, aes(x = trials, y = value, color = participant)) +
        geom_line(alpha = 0.7) +
        labs(title = main_title,
             x = "Trials",
             y = "Model Value") +
        theme_minimal() +
        guides(color = "none") +
        geom_vline(xintercept = vert_line, linetype = "dotted", color = "black") +
        scale_x_continuous(limits = c(0, trials))
      
  }
}

# Function to plot predicted STL adaptation values
STLPredictPlot <- function(data, vert_line = NULL){
  
main <- "Predicted STL Adaption"
  
ggplot(data, aes(x = rotations, y = predictions, group = participant)) +
  geom_line(color = "red", alpha = 0.6) +  # Lower alpha for more transparency
  labs(title = main,
       x = "Rotations in Degrees(°)",
       y = "Predicted Adaptation in Degrees(°)") +
  theme_minimal() +
  geom_vline(xintercept = vert_line, linetype = "dotted", color = "black")    
}   

# To plot XY Scatter plots
STL_XY <- function(data, rot, x, y) {
  # Filter data based on the rotation value
  filtered_data <- data[data$rotation == rot, ]
  
  # Calculate the regression slope
  regression_model <- lm(as.formula(paste(y, "~", x)), data = filtered_data)
  slope <- coef(regression_model)[2]
  
  # Determine the range for the x and y axes
  x_range <- range(filtered_data[[x]], na.rm = TRUE)
  y_range <- range(filtered_data[[y]], na.rm = TRUE)
  common_range <- range(c(x_range, y_range))
  
  # Create the plot
  p <- ggplot(filtered_data, aes_string(x = x, y = y)) +
    geom_point(color = "black") +  # Plot points as black dots
    geom_abline(intercept = 0, slope = slope, linetype = "dotted", color = "red") +  # Add red, dotted regression line from origin
    labs(title = paste("STL vs Short-Fixed Prediction", rot, "°"),
         x = "STL Deviation (°)",
         y = "Exponential Deviation (°)") +
    coord_fixed(ratio = 1, xlim = common_range, ylim = common_range) +  # Ensure identical scaling for both axes
    theme_minimal()  # Use a minimal theme
  
  # Print the plot
  print(p)
}
