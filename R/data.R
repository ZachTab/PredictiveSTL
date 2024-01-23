library(Reach)

download_data <- function() {
  Reach::downloadOSFdata(repository='j67bv',
                         filelist=list('data'=c('performance.zip')), 
                         folder='data', 
                         unzip=TRUE, 
                         removezips=TRUE)
}


preprocess_file <- function(id) {
  
  # read file for participant id
  filename <- sprintf('data/%s_performance.csv', id) 
  df <- read.csv(filename, stringsAsFactors = FALSE)
  
  # extract relevant data 
  idx <- which(df$label == 'stl-target-rotation')
  response <- df$reachdeviation_deg[idx + 1] - df$reachdeviation_deg[idx - 1]
  rotation <- df$rotation[idx]
  
  # store into data frame
  output <- data.frame(rotation, response, 'trial' = idx)
  output$participant <- id
  
  # normalize data
  output$response[which(output$rotation > 0)] <- -1 *  output$response[which(output$rotation > 0)]
  output$rotation[which(output$rotation < 0)] <- -1 *  output$rotation[which(output$rotation < 0)]
  return(output)
  
}






