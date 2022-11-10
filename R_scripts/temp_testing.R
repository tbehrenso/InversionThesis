
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
#     Comparing original and SFS nucdiv functions
# -----------------------------------

positions_test <- c(1,5,10,15,20,25)

ms_test <- matrix(c(0,0,0,0,1,0,
                    0,1,0,0,1,1,
                    0,1,0,1,1,0,
                    1,0,0,1,1,1,
                    0,1,1,0,1,0,
                    0,0,0,1,0,1), nrow = 6, byrow = T)

ms_test_long <- rbind(ms_test, ms_test, ms_test)

calc_nuc_div <- function(msdata, positions, totalLength, seqLen=200, centerSpacing=100){
  centers <- seq(0, totalLength, by=centerSpacing)
  # prepare storage for nucleotide diversity at each center position
  output <- data.frame(position=centers, nuc_div=NA)
  num_of_seq <- dim(msdata)[1]
  for(seqCenter in centers){
    positions_in_sequence <- which(   # select positions within window
      positions > seqCenter-seqLen & positions <= seqCenter+seqLen & !(positions %in% c(INV_START, INV_END-1))
    )  
    ms_in_seq <- as.matrix(msdata[,positions_in_sequence])
    # only do the calculations if more than one position
    if(dim(ms_in_seq)[2] > 0){
      # dist calculates distance between every combination of rows in a matrix. Manhattan method avoids "diagonal" distance
      distances_all <- dist(ms_in_seq, method='manhattan')
      # use this to adjust the sequence length when the window exceeds the range of the genome
      adjusted_seq_len <- sum((c((seqCenter-seqLen):(seqCenter+seqLen)))>=0 & (c((seqCenter-seqLen):(seqCenter+seqLen)))<=totalLength)
      # at a given position (center), nucdiv is the average number of differences divided by the length of the sequence window
      # Note: average number of differences is the total number of pairwise differences / the number of pairwise difference nChoosek
    
      # output[output$position==seqCenter,2] <- (sum(distances_all) / choose(num_of_seq,2)) / adjusted_seq_len
      # here using 2x sum of distances_all to account for both parts of the pairwise comparison matrix (above and below the diagonal)
      output[output$position==seqCenter,2] <- (2*sum(distances_all) / (num_of_seq^2)) / adjusted_seq_len
    }
  }
  return(output)
}

calc_nuc_div_sfs <- function(msdata, positions, totalLength, seqLen=200, centerSpacing=100){
  centers <- seq(0, totalLength, by=centerSpacing)
  # prepare storage for nucleotide diversity at each center position
  output <- data.frame(position=centers, nuc_div=NA)
  num_of_seq <- dim(msdata)[1]
  for(seqCenter in centers){
    positions_in_sequence <- which(   # select positions within window
      positions > seqCenter-seqLen & positions <= seqCenter+seqLen & !(positions %in% c(INV_START, INV_END-1))
    )  
    ms_in_seq <- as.matrix(msdata[,positions_in_sequence])
    # only do the calculations if more than one position
    if(dim(ms_in_seq)[2] > 0){
      # use this to adjust the sequence length when the window exceeds the range of the genome
      adjusted_seq_len <- sum((c((seqCenter-seqLen):(seqCenter+seqLen)))>=0 & (c((seqCenter-seqLen):(seqCenter+seqLen)))<=totalLength)
      
      sfs.raw <- colSums(ms_in_seq)
      #sfs.inverse <- num_of_seq - sfs.raw
      # adjust so it always considers in respect to the less frequent allele
      #sfs.total <- pmin(sfs.raw, sfs.inverse)
      # just some line above to get the counts of alleles per site, across all individuals you sample. 
      sfs.total <- sfs.raw
      
      p.all <- sfs.total/num_of_seq # your sample size is 200
      q.all <- 1-p.all
      numerator.all <- 2*p.all*q.all*sfs.total
      pi.all <- sum(numerator.all)/adjusted_seq_len     # want to divide by all possible sites, not just where SNPs are
      output[output$position==seqCenter,2] <- pi.all
    }
  }
  return(output)
}

calc_nuc_div(ms_test_long, positions_test, totalLength=26, seqLen=6,centerSpacing=10)

calc_nuc_div_sfs(ms_test_long, positions_test, totalLength=26, seqLen=6,centerSpacing=10)
