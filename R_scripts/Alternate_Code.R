# -----------------------------------------------------   (probably need to run main script first)
#    EFFECT OF WINDOW SIZE ON NUCLEOTIDE DIVERSITY
# -----------------------------------------------------

nucdiv_windowed <- calc_nuc_div(ms_binary, abs_positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
nucdiv_df[i,] <- nucdiv_windowed[[2]]

msdata <- ms_binary
positions <- abs_positions
totalLength <- GENOME_LENGTH
seqLen <- WINDOW_SIZE
centerSpacing <- WINDOW_SPACING




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
    if(dim(ms_in_seq)[2] > 1){
      # dist calculates distance between every combination of rows in a matrix. Manhattan method avoids "diagonal" distance
      distances_all <- dist(ms_in_seq, method='manhattan')
      # use this to adjust the sequence length when the window exceeds the range of the genome
      adjusted_seq_len <- sum((c((seqCenter-seqLen):(seqCenter+seqLen)))>=0 & (c((seqCenter-seqLen):(seqCenter+seqLen)))<=totalLength)
      # at a given position (center), nucdiv is the average number of differences divided by the length of the sequence window
      # Note: average number of differences is the total number of pairwise differences / the number of pairwise difference nChoosek
      
      #output[output$position==seqCenter,2] <- (sum(distances_all) / choose(num_of_seq,2)) / adjusted_seq_len
      output[output$position==seqCenter,2] <- (2*sum(distances_all) / num_of_seq^2) / adjusted_seq_len
    }
  }
  return(output)
}


m1 <- round(matrix(runif(200*5), 200, 5))
dist1 <- dist(m1, method='manhattan')
sum(dist1)/choose(200,2)/101

m2 <- round(matrix(runif(200*10), 200, 10))
dist2 <- dist(m2, method='manhattan')
sum(dist2)/choose(200,2)/200

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
}



# TAKEN FROM SCRIPT FOR QUICKER TESTING
window_centers <- seq(0, GENOME_LENGTH, by=WINDOW_SPACING)
nucdiv_df <- matrix(0, nrow=n_files, ncol=length(window_centers))

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
  # extract metadata from filename
  tags <- strsplit(files[i], split='_')[[1]]
  tags_index[i,] <- list(tags[2], as.numeric(tags[3]), as.numeric(tags[4]), as.integer(tags[5]))
  
  # calc nucleotide diversity (over sliding window by default)
  nucdiv_windowed <- calc_nuc_div(ms_binary, abs_positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
  nucdiv_df[i,] <- nucdiv_windowed[[2]]
}

nucdiv_summ_p1 <- data.frame(center = window_centers, nucdiv = colMeans(nucdiv_df[which(tags_index$population=='p1'),], na.rm=T),
                             stdev=apply(nucdiv_df[which(tags_index$population=='p1'),], 2, sd, na.rm=T))
nucdiv_summ_p2 <- data.frame(center = window_centers, nucdiv = colMeans(nucdiv_df[which(tags_index$population=='p2'),], na.rm=T),
                             stdev=apply(nucdiv_df[which(tags_index$population=='p2'),], 2, sd, na.rm=T))

nucdiv_a <- ggplot(nucdiv_summ_p1, aes(x=center, y=nucdiv)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('c) Nucleotide Diversity - P1') +
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red')} +
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue', alpha=0.4)} +
  xlab('Position') + ylab('\u03c0')
#ylim(c(0.0014, 0.0029))

nucdiv_b <- ggplot(nucdiv_summ_p2, aes(x=center, y=nucdiv)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('d) Nucleotide Diversity - P2') +
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red')} +
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue', alpha=0.4)} +
  #geom_errorbar(aes(ymin=nucdiv-stdev, ymax=nucdiv+stdev), width=1, position=position_dodge(0.1)) +
  xlab('Position') + ylab('\u03c0')

grid.arrange(nucdiv_a, nucdiv_b, nrow=1)


##### Here comparing histograms at first two positions

pos1 <- data.frame(x=nucdiv_df[which(tags_index$population=='p1'),1])
pos2 <- data.frame(x=nucdiv_df[which(tags_index$population=='p1'),2])

pos_both <- data.frame(x=c(nucdiv_df[which(tags_index$population=='p1'),1],nucdiv_df[which(tags_index$population=='p1'),2]),
                       pos=as.factor(c(rep(1,100), rep(2, 100))))

ggplot(pos_both, aes(x=x, fill=pos)) + geom_histogram(alpha=0.5, position="identity")


#### COMPARING DIFFERENT WINDOWS SIZES
seqLen <- WINDOW_SIZE / 2
centerSpacing <- WINDOW_SPACING / 2


window_centers <- seq(0, GENOME_LENGTH, by=centerSpacing)
nucdiv_df_large <- matrix(0, nrow=n_files, ncol=length(window_centers))

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
  # extract metadata from filename
  tags <- strsplit(files[i], split='_')[[1]]
  tags_index[i,] <- list(tags[2], as.numeric(tags[3]), as.numeric(tags[4]), as.integer(tags[5]))
  
  # calc nucleotide diversity (over sliding window by default)
  nucdiv_windowed <- calc_nuc_div(ms_binary, abs_positions, GENOME_LENGTH, seqLen = seqLen, centerSpacing = centerSpacing)
  nucdiv_df_large[i,] <- nucdiv_windowed[[2]]
}

nucdiv_summ_p1_LARGE <- data.frame(center = window_centers, nucdiv = colMeans(nucdiv_df_large[which(tags_index$population=='p1'),], na.rm=T),
                                   stdev=apply(nucdiv_df_large[which(tags_index$population=='p1'),], 2, sd, na.rm=T))
nucdiv_summ_p2_LARGE <- data.frame(center = window_centers, nucdiv = colMeans(nucdiv_df_large[which(tags_index$population=='p2'),], na.rm=T),
                                   stdev=apply(nucdiv_df_large[which(tags_index$population=='p2'),], 2, sd, na.rm=T))

nucdiv_a_large <- ggplot(nucdiv_summ_p1_LARGE, aes(x=center, y=nucdiv)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('c) Nucleotide Diversity - P1') +
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red')} +
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue', alpha=0.4)} +
  xlab('Position') + ylab('\u03c0')
#ylim(c(0.0014, 0.0029))

nucdiv_b_large <- ggplot(nucdiv_summ_p2_LARGE, aes(x=center, y=nucdiv)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('d) Nucleotide Diversity - P2') +
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red')} +
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue', alpha=0.4)} +
  #geom_errorbar(aes(ymin=nucdiv-stdev, ymax=nucdiv+stdev), width=1, position=position_dodge(0.1)) +
  xlab('Position') + ylab('\u03c0')


grid.arrange(nucdiv_a_large, nucdiv_b_large, nrow=1)

#### Histogram, but comparing runs with different window size (and spacing)
# 1 is window size 200, 2 is window size 100
winpos1 <- data.frame(x=nucdiv_df[which(tags_index$population=='p1'),2])
winpos2 <- data.frame(x=nucdiv_df_large[which(tags_index$population=='p1'),2])

winpos_both <- data.frame(x=c(nucdiv_df[which(tags_index$population=='p1'),2],nucdiv_df_large[which(tags_index$population=='p1'),2]),
                          pos=as.factor(c(rep(1,100), rep(2, 100))))

ggplot(winpos_both, aes(x=x, fill=pos)) + geom_histogram(alpha=0.5, position="identity")


# -----------------------------------------------------
#    
# -----------------------------------------------------
