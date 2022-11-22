
# read in files (values: selection coefficient, migration rate, replicate #)
files <- list.files(path=PATH, pattern="*.txt", full.names=F, recursive=FALSE)
n_files <- length(files)
# pre-calculate window centers' positions
window_centers <- seq(0, GENOME_LENGTH, by=WINDOW_SPACING)

tags_index <- data.frame(population=character(n_files), sel_coef=numeric(n_files), migration=numeric(n_files), 
                         repl=integer(n_files), stringsAsFactors=F)
correlations_3d <- array(numeric(), dim=c(N_TILES, N_TILES, n_files))



if(INVERSION_PRESENT && generation > 5000){
  n_repl <- length(unique(tags_index$repl))
  
  correlations_3d_normal <- array(numeric(), dim=c(N_TILES, N_TILES, n_repl))
  correlations_3d_inverted <- array(numeric(), dim=c(N_TILES, N_TILES, n_repl))  
  
  for(repl in 1:max(tags_index$repl)){
    # get and prepare ms_data
    filepath_p1 <- paste0(PATH, "/", files[which(tags_index$population=='p1' & tags_index$repl==repl)])
    filepath_p2 <- paste0(PATH, "/", files[which(tags_index$population=='p2' & tags_index$repl==repl)])
    ms_p1 <- get_ms_data(filepath_p1)
    ms_p2 <- get_ms_data(filepath_p2)
    pos_p1 <- get_positions(filepath_p1)
    pos_p2 <- get_positions(filepath_p2)
    pos_both <- unique(sort(c(pos_p1, pos_p2)))
    colnames(ms_p1) <- pos_p1
    colnames(ms_p2) <- pos_p2
    n_indiv <- dim(ms_p1)[1]
    
    ms_both <- matrix(0, nrow=2*n_indiv, ncol=length(pos_both))
    colnames(ms_both) <- pos_both
    # top half of new matrix is p1 data, bottom half is p2 data. All missing rows in a population 0 by default
    ms_both[1:200,as.character(pos_p1)] <- ms_p1
    ms_both[201:400,as.character(pos_p2)] <- ms_p2
    # extract rows based on the presence of both inversion markers
    ms_normal <- ms_both[ms_both[,as.character(INV_START)]==0 & ms_both[,as.character(INV_START)]==0, ]
    ms_inverted <- ms_both[ms_both[,as.character(INV_START)]==1 & ms_both[,as.character(INV_START)]==1, ]
    # unname here for compatability with reduce_to_long() function
    corr_data_normal <- unname(get_correlations(ms_normal, pos_both, numTiles = N_TILES))  
    corr_data_inverted <- unname(get_correlations(ms_inverted, pos_both, numTiles = N_TILES))
    corr_long_normal <- reduce_to_long(corr_data_normal, pos_both, numTiles = N_TILES)
    corr_long_inverted <- reduce_to_long(corr_data_inverted, pos_both, numTiles = N_TILES)
    
    correlations_3d_normal[,,repl] <- as.matrix(dcast(corr_long_normal, Var1 ~ Var2)[,-1])
    correlations_3d_inverted[,,repl] <- as.matrix(dcast(corr_long_inverted, Var1 ~ Var2)[,-1])
    
  }
}

# correlation heatmap
corr_summ_p1 <- apply(correlations_3d_normal, c(1, 2), mean, na.rm = TRUE)
corr_summ_p1_long <- melt(corr_summ_p1)
# correct group values to bin centers
bin_size <- GENOME_LENGTH / N_TILES
corr_summ_p1_long$Var1 <- corr_summ_p1_long$Var1*bin_size - bin_size/2
corr_summ_p1_long$Var2 <- corr_summ_p1_long$Var2*bin_size - bin_size/2

corr_summ_p2 <- apply(correlations_3d_inverted, c(1, 2), mean, na.rm = TRUE)
corr_summ_p2_long <- melt(corr_summ_p2)
# correct group values to bin centers
bin_size <- GENOME_LENGTH / N_TILES
corr_summ_p2_long$Var1 <- corr_summ_p2_long$Var1*bin_size - bin_size/2
corr_summ_p2_long$Var2 <- corr_summ_p2_long$Var2*bin_size - bin_size/2

corr_a <- ggplot(corr_summ_p1_long, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low='white', high='blue') +
  ggtitle('P1') + xlab('Position') + ylab('Position')

corr_b <- ggplot(corr_summ_p2_long, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low='white', high='blue') +
  ggtitle('P2') + xlab('Position') + ylab('Position')

plot_correlation <- grid.arrange(corr_a, corr_b, nrow=1)


# -----------------------------------
#  Calc nucdiv in Kim's Windows of replicate 1
# -----------------------------------

FILEPATH1 <- "Outputs/inversionLAA_2pop_s0.01_m0.001_mu1e-6/15000/linkage2P_p1_0.01_0.001_1_.txt"
FILEPATH2 <- "Outputs/inversionLAA_2pop_s0.01_m0.001_mu1e-6/15000/linkage2P_p2_0.01_0.001_1_.txt"

ms_p1 <- get_ms_data(FILEPATH1)
ms_p2 <- get_ms_data(FILEPATH2)
pos_p1 <- get_positions(FILEPATH1)
pos_p2 <- get_positions(FILEPATH2)

ms_p1_inv <- ms_p1[ms_p1[,which(pos_p1==INV_START)]==1 & ms_p1[,which(pos_p1==INV_END-1)]==1,]
ms_p1_nor <- ms_p1[ms_p1[,which(pos_p1==INV_START)]!=1 & ms_p1[,which(pos_p1==INV_END-1)]!=1,]
ms_p2_inv <- ms_p2[ms_p2[,which(pos_p2==INV_START)]==1 & ms_p2[,which(pos_p2==INV_END-1)]==1,]
ms_p2_nor <- ms_p2[ms_p2[,which(pos_p2==INV_START)]!=1 & ms_p2[,which(pos_p2==INV_END-1)]!=1,]


ms_p1_inv_win1 <- ms_p1_inv[,pos_p1 <= 2500 & pos_p1 >= 1]
ms_p1_inv_win2 <- ms_p1_inv[,pos_p1 <= 12500 & pos_p1 >= 10001]
ms_p2_inv_win1 <- ms_p2_inv[,pos_p2 <= 2500 & pos_p2 >= 1]
ms_p2_inv_win2 <- ms_p2_inv[,pos_p2 <= 12500 & pos_p2 >= 10001]
ms_p1_nor_win1 <- ms_p1_nor[,pos_p1 <= 2500 & pos_p1 >= 1]
ms_p1_nor_win2 <- ms_p1_nor[,pos_p1 <= 12500 & pos_p1 >= 10001]
ms_p2_nor_win1 <- ms_p2_nor[,pos_p2 <= 2500 & pos_p2 >= 1]
ms_p2_nor_win2 <- ms_p2_nor[,pos_p2 <= 12500 & pos_p2 >= 10001]



nucdiv_short <- function(ms_in_seq, seqLen){
  num_of_seq <- dim(ms_in_seq)[1]
  # dist calculates distance between every combination of rows in a matrix. Manhattan method avoids "diagonal" distance
  distances_all <- dist(ms_in_seq, method='manhattan')
  # use this to adjust the sequence length when the window exceeds the range of the genome
  output <- (2*sum(distances_all) / (num_of_seq^2)) / seqLen
  return(output)
}


nucdiv_p1_inv_win1 <- nucdiv_short(ms_p1_inv_win1, 2500)
nucdiv_p1_inv_win2 <- nucdiv_short(ms_p1_inv_win2, 2500)
nucdiv_p2_inv_win1 <- nucdiv_short(ms_p2_inv_win1, 2500)
nucdiv_p2_inv_win2 <- nucdiv_short(ms_p2_inv_win2, 2500)
nucdiv_p1_nor_win1 <- nucdiv_short(ms_p1_nor_win1, 2500)
nucdiv_p1_nor_win2 <- nucdiv_short(ms_p1_nor_win2, 2500)
nucdiv_p2_nor_win1 <- nucdiv_short(ms_p2_nor_win1, 2500)
nucdiv_p2_nor_win2 <- nucdiv_short(ms_p2_nor_win2, 2500)
















# -----------------------------------
#  Gaussian sliding window
# -----------------------------------

testdata <- c(1,5,1,6,8,2,3,6,2,7,5,7,5,7,2,4,6,2)
testpositions <- c(1,12,35,41,45,58,65,67,74,88,99,101,109,110,121,131,136,140)

testall <- data.frame(testpositions, testdata)

center <- 8

dnorm(testpositions, mean=testpositions[center], sd=20)

# takes a dataframe where first column is position and second column is the value
calc_sliding_window_gaussian <- function(posValData, totalLength, windowSize, pointSpacing, stdev){
  centers <- seq(0, totalLength, by=pointSpacing)
  output <- data.frame(position=centers, average=NA)
  for(i in 1:length(centers)){
    correspondingIndeces <- posValData[1] > centers[i]-windowSize & posValData[1] <= centers[i]+windowSize
    correspondingPositions <- posValData[1][correspondingIndeces]
    correspondingValues <- posValData[2][correspondingIndeces]
    gaussianValues <- dnorm(correspondingPositions, mean=centers[i], sd=stdev)
    gaussianNormalized <- gaussianValues / sum(gaussianValues)
    output[i, 2] <- mean(correspondingValues)
  }
  return(output)
}

gaussian_test <- calc_sliding_window_gaussian(testall, totalLength=150, windowSize=20, pointSpacing=20)


plot(testpositions, testdata, type='l')

plot(gaussian_test$position, gaussian_test$average, type='l')








