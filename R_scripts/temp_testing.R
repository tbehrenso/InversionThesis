
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

positions_test <- c(1,6,11,16,21,26)

ms_test <- matrix(c(0,0,0,0,1,0,
                    0,1,0,0,1,1,
                    0,1,0,1,1,0,
                    1,0,0,1,1,1,
                    0,1,1,0,1,0,
                    0,0,0,1,0,1), nrow = 6, byrow = T)

ms_test_long <- rbind(ms_test, ms_test, ms_test, ms_test, ms_test, ms_test)


positions_test_addedzeros <- c(1,6,11,16,21,25,26)
ms_test_addedzeros <- matrix(c(0,0,0,0,1,0,0,
                               0,1,0,0,1,0,1,
                               0,1,0,1,1,0,0,
                               1,0,0,1,1,0,1,
                               0,1,1,0,1,0,0,
                               0,0,0,1,0,0,1), nrow = 6, byrow = T)

calc_nuc_div <- function(msdata, positions, totalLength, seqLen=200, centerSpacing=100){
  centers <- seq(0, totalLength, by=centerSpacing)
  # prepare storage for nucleotide diversity at each center position
  output <- data.frame(position=centers, nuc_div=NA)
  num_of_seq <- dim(msdata)[1]
  for(seqCenter in centers){
    positions_in_sequence <- which(   # select positions within window
      positions > seqCenter-seqLen & positions <= seqCenter+seqLen
    ) 
    print(positions_in_sequence)
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
      positions > seqCenter-seqLen & positions <= seqCenter+seqLen
    )  
    ms_in_seq <- as.matrix(msdata[,positions_in_sequence])
    # only do the calculations if more than one position
    if(dim(ms_in_seq)[2] > 0){
      # use this to adjust the sequence length when the window exceeds the range of the genome
      adjusted_seq_len <- sum((c((seqCenter-seqLen):(seqCenter+seqLen)))>=0 & (c((seqCenter-seqLen):(seqCenter+seqLen)))<=totalLength)
      
      sfs.raw <- table(colSums(ms_in_seq))

      # alternate to adjust so it always considers in respect to the less frequent allele
      # sfs.inverse <- num_of_seq - sfs.raw
      # sfs.total <- pmin(sfs.raw, sfs.inverse)
      
      sfs.total <- sfs.raw
      
      counts.sfs.total <- as.numeric(names(sfs.total))
      
      p.all <- counts.sfs.total/num_of_seq # your sample size is 200
      q.all <- 1-p.all
      numerator.all <- 2*p.all*q.all*sfs.total
      pi.all <- sum(numerator.all)/adjusted_seq_len     # want to divide by all possible sites, not just where SNPs are
      output[output$position==seqCenter,2] <- pi.all
    }
  }
  return(output)
}
calc_nuc_div_popgenome <- function(msdata, positions, totalLength, seqLen=200, centerSpacing=100){
  centers <- seq(0, totalLength, by=centerSpacing)
  # prepare storage for nucleotide diversity at each center position
  output <- data.frame(position=centers, nuc_div=NA)
  num_of_seq <- dim(msdata)[1]
  for(seqCenter in centers){
    positions_in_sequence <- which(   # select positions within window
      positions > seqCenter-seqLen & positions <= seqCenter+seqLen
    )  
    ms_in_seq <- as.matrix(msdata[,positions_in_sequence])
    # only do the calculations if more than one position
    if(dim(ms_in_seq)[2] > 0){
      ones <- colSums(ms_in_seq)
      zeros <- num_of_seq - ones
      n.comparisons <- (num_of_seq * (num_of_seq-1)) / 2
      nuc_div_all <- (ones * zeros) / n.comparisons
      nuc_div_mean <- mean(nuc_div_all)
      
      output[output$position==seqCenter,2] <- nuc_div_mean
    }
  }
  return(output)
}

# testing on default version
calc_nuc_div(ms_test, positions_test, totalLength=26, seqLen=7,centerSpacing=7)
calc_nuc_div_sfs(ms_test, positions_test, totalLength=26, seqLen=7,centerSpacing=7)
calc_nuc_div_popgenome(ms_test, positions_test, totalLength=26, seqLen=7,centerSpacing=7)

# testing on longer version (more samples)
calc_nuc_div(ms_test_long, positions_test, totalLength=26, seqLen=7,centerSpacing=7)
calc_nuc_div_sfs(ms_test_long, positions_test, totalLength=26, seqLen=7,centerSpacing=7)
calc_nuc_div_popgenome(ms_test_long, positions_test, totalLength=26, seqLen=7,centerSpacing=7)

# testing with smaller window (and hence more windows)
calc_nuc_div(ms_test, positions_test, totalLength=26, seqLen=5,centerSpacing=5)
calc_nuc_div_sfs(ms_test, positions_test, totalLength=26, seqLen=5,centerSpacing=5)
calc_nuc_div_popgenome(ms_test, positions_test, totalLength=26, seqLen=5,centerSpacing=5)

# testing on version with a row of zeros in final window
calc_nuc_div(ms_test_addedzeros, positions_test_addedzeros, totalLength=26, seqLen=7,centerSpacing=7)
calc_nuc_div_sfs(ms_test_addedzeros, positions_test_addedzeros, totalLength=26, seqLen=7,centerSpacing=7)
calc_nuc_div_popgenome(ms_test_addedzeros, positions_test_addedzeros, totalLength=26, seqLen=7,centerSpacing=7)


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













