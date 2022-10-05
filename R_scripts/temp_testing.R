hexp_windowed <- calc_sliding_window(pos_hexp, GENOME_LENGTH, windowSize = WINDOW_SIZE, pointSpacing = WINDOW_SPACING)

frequencies_windowed <- calc_sliding_window(pos_freq, GENOME_LENGTH, windowSize = WINDOW_SIZE, pointSpacing = WINDOW_SPACING)




posValData <- pos_hexp
totalLength <- GENOME_LENGTH
windowSize <- WINDOW_SIZE
pointSpacing <- WINDOW_SPACING


calc_sliding_window <- function(posValData, totalLength, windowSize, pointSpacing){
  centers <- seq(0, totalLength, by=pointSpacing)
  output <- data.frame(position=centers, average=0)
  for(i in 1:length(centers)){
    correspondingValues <- posValData[2][posValData[1] > centers[i]-windowSize & posValData[1] <= centers[i]+windowSize]
    output[i, 2] <- mean(correspondingValues)
  }
  return(output)
}
