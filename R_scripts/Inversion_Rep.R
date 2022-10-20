setwd("C:/Users/tbehr/Desktop/Thesis")

library(ggplot2)
library(gridExtra)
library(reshape2)
library(tidyr)

#args = commandArgs(trailingOnly=TRUE)

#-----------------------------------------------------------
# PARAMETERS
#-----------------------------------------------------------

GENOME_LENGTH <- 22000
FIXED_MUTATION_POS1 <- 8000
FIXED_MUTATION_POS2 <- 12000
INV_START <-6000
INV_END <- 16000  # this value should NOT be the '-1' value that the SLiM script uses. This script does that correction later
# PATH <- "RepOutput_Inversion_Burnin_20percent/6000/"
PATH <- "Outputs/inversionLAA_2pop_s0.01_m0.001_mu1e-6/7000"
WINDOW_SPACING <- 200
WINDOW_SIZE <- 100   # NOTE: window size is added on each side (so the full size is more like twice this value)
N_TILES <- 60

# record presence or absence of inversion and locally adapted alleles
INVERSION_PRESENT <- TRUE
LAA_PRESENT <- TRUE


#-----------------------------------------------------------
# FUNCTIONS
#-----------------------------------------------------------

get_ms_data <- function(filename){
  # read in ms file (ignoring first 3 lines), and as strings
  ms_binary_conjoined <- read.table(filename, skip = 3, colClasses = 'character')
  # split each line into its elements as separate columns
  ms_binary <- apply(ms_binary_conjoined, 1, function(x) as.integer(unlist(strsplit(x,""))))
  # flip axes so it matches with original file
  ms_binary <- t(ms_binary)
  return(ms_binary)
}

get_positions <- function(filename){
  # read in third line of text file containing relative positions
  positions_row <- readLines(filename)[3]
  relative_positions <- as.numeric(unlist(strsplit(positions_row, " "))[-1])
  abs_positions <- relative_positions * (GENOME_LENGTH-1)
  #return(abs_positions)
  return(round(abs_positions))
}

calc_hexp <- function(msdata){
  num_individuals <- nrow(msdata) 
  # calculate expected heterozygosity as 2*p*(1-p)
  expected_heterozygosity <- 2 * (colSums(msdata) / num_individuals) * (1-(colSums(msdata) / num_individuals))
  return(expected_heterozygosity)
}

# calculating nucleotide diversity. Returns dataframe with nucdiv at spaced positions  --> NOTE: edited to remove marker mutations
# NOTE: seqLen is only really HALF the window size (adds seqLen in both directions of a point)
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
      
      
      # output[output$position==seqCenter,2] <- (sum(distances_all) / choose(num_of_seq,2)) / adjusted_seq_len
      output[output$position==seqCenter,2] <- (2*sum(distances_all) / (num_of_seq^2)) / adjusted_seq_len
    }
  }
  return(output)
}

# takes a dataframe where first column is position and second column is the value
calc_sliding_window <- function(posValData, totalLength, windowSize, pointSpacing){
  centers <- seq(0, totalLength, by=pointSpacing)
  output <- data.frame(position=centers, average=0)
  for(i in 1:length(centers)){
    correspondingValues <- posValData[2][posValData[1] > centers[i]-windowSize & posValData[1] <= centers[i]+windowSize]
    output[i, 2] <- mean(correspondingValues)
  }
  return(output)
}

# MODIFIED TO REMOVE MARKER MUTATIONS
get_correlations <- function(msdata, positions, numTiles=20){
  # remove inversion marker mutations
  inv_start_index_reversed <- which(positions!=INV_START)
  inv_end_index_reversed <- which(positions!=INV_END-1)
  not_inv_markers <- intersect(inv_start_index_reversed, inv_end_index_reversed)
  positions <- positions[not_inv_markers]
  msdata <- msdata[,not_inv_markers]
  
  # prep storage (with default value of 1 for the diagonal)
  num_sites <- length(positions)
  # use default method (pearson)
  corr_all <- cor(msdata, method="pearson")
  return(corr_all)
}

#MODIFIED TO REMOVE MARKER MUTATIONS
reduce_to_long <- function(corrData, positions, numTiles=20){
  # split positions into bins (using range up to full length so replicates can be combined)
  # remove inversion marker mutations
  inv_start_index <- which(positions!=INV_START)
  inv_end_index <- which(positions!=INV_END-1)
  not_inv_markers <- intersect(inv_start_index, inv_end_index)
  positions <- positions[not_inv_markers]
  max_pos <- max(positions)
  
  groups <- cut(c(0, positions, GENOME_LENGTH), breaks=numTiles, labels=F)
  # associate positions with their groups. Remove first and last group which were only included to specify range
  pos_grouping <- data.frame(position=positions, group=groups[-c(1,length(groups))]) 
  
  # convert to long, then convert position indeces to corresponding group numbers
  data_long <- melt(corrData)
  data_long$Var1 <- pos_grouping$group[data_long$Var1]
  data_long$Var2 <- pos_grouping$group[data_long$Var2]
  
  # averaging correlations within each combination (i,j) of bins
  red_long_incomplete <- aggregate(value ~ Var1 + Var2, data=data_long, FUN=mean)
  # dataframe including all possible tile coordinates, to combine with aggregated means
  red_long_empty <- data.frame(Var1=sort(rep(1:numTiles, numTiles)), Var2=rep(1:numTiles, numTiles))
  red_long <- merge(red_long_empty, red_long_incomplete, by=c('Var1', 'Var2'), all=T)
  # scale group numbers to nucleotide positions
  red_long[1] <- red_long[1] * (max(max_pos)/numTiles)
  red_long[2] <- red_long[2] * (max(max_pos)/numTiles)
  #red_long[is.na(red_long)] <- 0                           # does it make sense to fill with zero?????????
  
  return(red_long)
}

# Getting average (and SD) number of polymorphic sites. 
get_average_polymorphism_count <- function(){
  polymorphism_counts <- numeric(n_files)
  
  for(i in 1:n_files){
    filepath <- paste0(PATH, "/", files[i])
    abs_positions <- get_positions(filepath)
    
    polymorphism_counts[i] <- length(abs_positions)
    average <- mean(polymorphism_counts)
    stdev <- sd(polymorphism_counts)
  }
  print(paste('Polymorphism Count - Mean:', average, 'SD:', stdev))
}

# INVERSION FREQUENCIES      ---> can make generic to take any two (or set) of positions
get_inversion_frequencies <- function(){
  freq_all <- numeric(n_files)
  for(i in 1:n_files){
    filepath <- paste0(PATH, "/", files[i])
    ms_binary <- get_ms_data(filepath)
    abs_positions <- get_positions(filepath)
    
    if(all(c(INV_START, INV_END-1) %in% abs_positions)){
      inv_start_index <- which(abs_positions==INV_START)
      inv_start_index <- inv_start_index[1] ##### TEMPORARY: TO "FIX" MULTIPLE MUTATIONS AT A SITE
      inv_end_index <- which(abs_positions==INV_END-1)
      inv_end_index <- inv_end_index[1]  ##### TEMPORARY: TO "FIX" MULTIPLE MUTATIONS AT A SITE
      
      freq_all[i] <- mean(ms_binary[,inv_start_index])
    }
  }
  
  inv_freq_p1 <- mean(freq_all[tags_index$population=='p1'])
  inv_freq_p2 <- mean(freq_all[tags_index$population=='p2'])
  print(paste('P1 Frequency:', inv_freq_p1))
  print(paste('P2 Frequency:', inv_freq_p2))
  return(c(inv_freq_p1, inv_freq_p2))
}

# get frequency of neutral alleles (only does mean, haven't implemented standard deviation)
get_neutral_frequency <- function(){
  neutral_frequencies <- matrix(0, nrow = n_files, ncol = 2)
  for(i in 1:n_files){
    filepath <- paste0(PATH, "/", files[i])
    ms_binary <- get_ms_data(filepath)
    abs_positions <- get_positions(filepath)
    
    ms_neutrals <- ms_binary[,!(abs_positions %in% c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2, INV_START, INV_END-1))]
    neutral_frequencies[i, 1] <- mean(ms_neutrals)
    neutral_frequencies[i, 2] <- dim(ms_neutrals)[2]
  }
  overall_neutral_frequency <- sum(apply(neutral_frequencies, 1, prod) / sum(neutral_frequencies[,2]))
  return(overall_neutral_frequency)
}

# get frequency of locally adapted alleles. 
get_allele_frequency <- function(){
  freq_allele1 <- numeric(n_files)
  freq_allele2 <- numeric(n_files)
  for(i in 1:n_files){
    filepath <- paste0(PATH, "/", files[i])
    ms_binary <- get_ms_data(filepath)
    abs_positions <- get_positions(filepath)
    
    if(FIXED_MUTATION_POS1 %in% abs_positions){
      freq_allele1[i] <- mean(ms_binary[,which(abs_positions==FIXED_MUTATION_POS1)])
    }
    if(FIXED_MUTATION_POS2 %in% abs_positions){
      freq_allele2[i] <- mean(ms_binary[,which(abs_positions==FIXED_MUTATION_POS2)])
    }
  }
  allele1_freq_p1 <- mean(freq_allele1[tags_index$population=='p1'])
  allele1_freq_p2 <- mean(freq_allele1[tags_index$population=='p2'])
  allele2_freq_p1 <- mean(freq_allele2[tags_index$population=='p1'])
  allele2_freq_p2 <- mean(freq_allele2[tags_index$population=='p2'])
  
  print(paste('P1 Allele1 Frequency:', allele1_freq_p1))
  print(paste('P1 Allele2 Frequency:', allele2_freq_p1))
  print(paste('P2 Allele1 Frequency:', allele1_freq_p2))
  print(paste('P2 Allele2 Frequency:', allele2_freq_p2))
}


# -------------------------------------------------------------------------------------------
# DATA EXTRACTION
# -------------------------------------------------------------------------------------------

# read in files (values: selection coefficient, migration rate, replicate #)
files <- list.files(path=PATH, pattern="*.txt", full.names=F, recursive=FALSE)
n_files <- length(files)
# pre-calculate window centers' positions
window_centers <- seq(0, GENOME_LENGTH, by=WINDOW_SPACING)

# STORAGE DATAFRAMES
tags_index <- data.frame(population=character(n_files), sel_coef=numeric(n_files), migration=numeric(n_files), 
                         repl=integer(n_files), stringsAsFactors=F)
hexp_red_df <- matrix(0, nrow=n_files, ncol=length(window_centers))
nucdiv_df <- matrix(0, nrow=n_files, ncol=length(window_centers))
frequencies <- matrix(0, nrow=n_files, ncol=length(window_centers))
correlations_3d <- array(numeric(), dim=c(N_TILES, N_TILES, n_files))  # NOTE: end up switching a lot between long and wide here
#      maybe change so its all in long?

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
  # extract metadata from filename
  tags <- strsplit(files[i], split='_')[[1]]
  tags_index[i,] <- list(tags[2], as.numeric(tags[3]), as.numeric(tags[4]), as.integer(tags[5]))
  
  # calc hexp at each position, and then a sliding window across the genome
  exp_het <- calc_hexp(ms_binary)
  pos_hexp <- data.frame(position=abs_positions, hexp=exp_het)
  # remove marker mutations
  pos_hexp <- pos_hexp[pos_hexp$position!=INV_START & pos_hexp$position!=INV_END-1,]
  hexp_windowed <- calc_sliding_window(pos_hexp, GENOME_LENGTH, windowSize = WINDOW_SIZE, pointSpacing = WINDOW_SPACING)
  hexp_red_df[i,] <- hexp_windowed[[2]]
  
  # calc nucleotide diversity (over sliding window by default)
  nucdiv_windowed <- calc_nuc_div(ms_binary, abs_positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
  nucdiv_df[i,] <- nucdiv_windowed[[2]]
  
  # calc mutation frequency, then sliding window
  pos_freq <- data.frame(position=abs_positions, freq=colMeans(ms_binary))
  frequencies_windowed <- calc_sliding_window(pos_freq, GENOME_LENGTH, windowSize = WINDOW_SIZE, pointSpacing = WINDOW_SPACING)
  frequencies[i,] <- frequencies_windowed[[2]]
  
  # correlation matrix into a 3D array (third dimension is file index)
  corr_data <- get_correlations(ms_binary, abs_positions, numTiles = N_TILES)
  corr_long <- reduce_to_long(corr_data, abs_positions, numTiles = N_TILES)
  correlations_3d[,,i] <- as.matrix(dcast(corr_long, Var1 ~ Var2)[,-1]) # exclude first column (variable names)
}

get_inversion_frequencies()

#-----------------------------------------------------------
# DIVERSITY
#-----------------------------------------------------------
hexp_summ_p1 <- data.frame(center= window_centers, hexp=colMeans(hexp_red_df[which(tags_index$population=='p1'),], na.rm=T), 
                           stdev=apply(hexp_red_df[which(tags_index$population=='p1'),], 2, sd, na.rm=T))
hexp_summ_p2 <- data.frame(center= window_centers, hexp=colMeans(hexp_red_df[which(tags_index$population=='p2'),], na.rm=T),
                           stdev=apply(hexp_red_df[which(tags_index$population=='p2'),], 2, sd, na.rm=T))

nucdiv_summ_p1 <- data.frame(center = window_centers, nucdiv = colMeans(nucdiv_df[which(tags_index$population=='p1'),], na.rm=T),
                             stdev=apply(nucdiv_df[which(tags_index$population=='p1'),], 2, sd, na.rm=T))
nucdiv_summ_p2 <- data.frame(center = window_centers, nucdiv = colMeans(nucdiv_df[which(tags_index$population=='p2'),], na.rm=T),
                             stdev=apply(nucdiv_df[which(tags_index$population=='p2'),], 2, sd, na.rm=T))

hexp_a <- ggplot(hexp_summ_p1, aes(x=center, y=hexp)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('a) Expected Heterozygosity - P1') +
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red')} +
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue', alpha=0.4)} +
  xlab('Position') + ylab(expression(H[exp]))
#ylim(c(0.061, 0.106))

hexp_b <- ggplot(hexp_summ_p2, aes(x=center, y=hexp)) +
  geom_line() +
  scale_color_brewer(palette="Dark2") +
  ggtitle('b) Expected Heterozygosity - P2') +
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red')} +
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue', alpha=0.4)} +
  xlab('Position') + ylab(expression(H[exp]))
#ylim(c(0.061, 0.106))

#  grid.arrange(hexp_a, hexp_b, nrow=1)

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
#ylim(c(0.0014, 0.0029))

#  grid.arrange(nucdiv_a, nucdiv_b, nrow=1)

grid.arrange(hexp_a, hexp_b, nucdiv_a, nucdiv_b, nrow=2)

# nucleotide diversity - separating haplotypes
nucdiv_inverted <- matrix(0, nrow=n_files, ncol=length(window_centers))
nucdiv_normal <- matrix(0, nrow=n_files, ncol=length(window_centers))

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
  # if at least one sample individual has the inversion
  if(all(c(INV_START, INV_END-1) %in% abs_positions)){
    inv_start_index <- which(abs_positions==INV_START)
    inv_start_index <- inv_start_index[1] ##### TEMPORARY: TO "FIX" MULTIPLE MUTATIONS AT A SITE
    inv_end_index <- which(abs_positions==INV_END-1)
    inv_end_index <- inv_end_index[1]  ##### TEMPORARY: TO "FIX" MULTIPLE MUTATIONS AT A SITE
    
    ms_inverted <- ms_binary[ms_binary[,inv_start_index]==1 & ms_binary[,inv_end_index]==1,]
    ms_normal <- ms_binary[ms_binary[,inv_start_index]!=1 & ms_binary[,inv_end_index]!=1,]
    
    # account for if only one individual is inverted/normal (in which case cannot calculate nucleotide diversity)
    if(is.null(dim(ms_inverted))){    # NOTE: This part is a mess right now. 
      nucdiv_inverted[i,] <- NA
      nucdiv_normal_windowed <- calc_nuc_div(ms_normal, abs_positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
      nucdiv_normal[i,] <- nucdiv_normal_windowed[[2]]
    } else if(is.null(dim(ms_normal))){
      nucdiv_normal[i,] <- NA
      nucdiv_inverted_windowed <- calc_nuc_div(ms_inverted, abs_positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
      nucdiv_inverted[i,] <- nucdiv_inverted_windowed[[2]]
    } else if (dim(ms_inverted)[1] == 0 ){
      nucdiv_inverted[i,] <- NA
      nucdiv_normal_windowed <- calc_nuc_div(ms_normal, abs_positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
      nucdiv_normal[i,] <- nucdiv_normal_windowed[[2]]
    } else if (dim(ms_inverted)[1] == 0){
      nucdiv_normal[i,] <- NA
      nucdiv_inverted_windowed <- calc_nuc_div(ms_inverted, abs_positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
      nucdiv_inverted[i,] <- nucdiv_inverted_windowed[[2]]
    }else {
      nucdiv_inverted_windowed <- calc_nuc_div(ms_inverted, abs_positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
      nucdiv_normal_windowed <- calc_nuc_div(ms_normal, abs_positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
      
      nucdiv_inverted[i,] <- nucdiv_inverted_windowed[[2]]
      nucdiv_normal[i,] <- nucdiv_normal_windowed[[2]]
    }
    # if no individual has the inversion (hence markers not in list of positions)
  } else {
    ms_normal <- ms_binary
    nucdiv_normal_windowed <- calc_nuc_div(ms_normal, abs_positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
    nucdiv_normal[i,] <- nucdiv_normal_windowed[[2]]
    
    nucdiv_inverted[i,] <- NA
  }
}

nucdiv_inverted_summ_p1 <- data.frame(center = window_centers, 
                                      nucdiv = colMeans(nucdiv_inverted[which(tags_index$population=='p1'),], na.rm=T),
                                      stdev = apply(nucdiv_inverted[which(tags_index$population=='p1'),], 2, sd, na.rm=T))
nucdiv_inverted_summ_p2 <- data.frame(center = window_centers, 
                                      nucdiv = colMeans(nucdiv_inverted[which(tags_index$population=='p2'),], na.rm=T),
                                      stdev = apply(nucdiv_inverted[which(tags_index$population=='p2'),], 2, sd, na.rm=T))
nucdiv_normal_summ_p1 <- data.frame(center = window_centers, 
                                    nucdiv = colMeans(nucdiv_normal[which(tags_index$population=='p1'),], na.rm=T),
                                    stdev = apply(nucdiv_normal[which(tags_index$population=='p1'),], 2, sd, na.rm=T))
nucdiv_normal_summ_p2 <- data.frame(center = window_centers, 
                                    nucdiv = colMeans(nucdiv_normal[which(tags_index$population=='p2'),], na.rm=T),
                                    stdev = apply(nucdiv_normal[which(tags_index$population=='p2'),], 2, sd, na.rm=T))
# compile into one data frame
if(INVERSION_PRESENT){
  nucdiv_all <- cbind(nucdiv_inverted_summ_p1$center, nucdiv_inverted_summ_p1$nucdiv, nucdiv_inverted_summ_p2$nucdiv, nucdiv_normal_summ_p1$nucdiv, nucdiv_normal_summ_p2$nucdiv)
  nucdiv_all <- as.data.frame(nucdiv_all)
  names(nucdiv_all) <- c('center', 'inv_P1', 'inv_P2', 'nor_P1', 'nor_P2')
  # exclude empty inverted columns if no inversion present
}else{
  nucdiv_all <- cbind(nucdiv_inverted_summ_p1$center, nucdiv_normal_summ_p1$nucdiv, nucdiv_normal_summ_p2$nucdiv)
  nucdiv_all <- as.data.frame(nucdiv_all)
  names(nucdiv_all) <- c('center', 'nor_P1', 'nor_P2')
}

nucdiv_all_long <- melt(nucdiv_all, id='center')

ggplot(nucdiv_all_long, aes(x=center, y=value, col=variable)) +
  geom_line() +
  ggtitle('Nucleotide Diversity for Different Haplotypes') +
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red')} +
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue', alpha=0.3)}

#-----------------------------------------------------------
# CORRELATION
#-----------------------------------------------------------

# correlation heatmap
corr_summ_p1 <- apply(correlations_3d[, , which(tags_index$population == "p1")], c(1, 2), mean, na.rm = TRUE)
corr_summ_p1_long <- melt(corr_summ_p1)
# correct group values to bin centers
bin_size <- GENOME_LENGTH / N_TILES
corr_summ_p1_long$Var1 <- corr_summ_p1_long$Var1*bin_size - bin_size/2
corr_summ_p1_long$Var2 <- corr_summ_p1_long$Var2*bin_size - bin_size/2

corr_summ_p2 <- apply(correlations_3d[, , which(tags_index$population == "p2")], c(1, 2), mean, na.rm = TRUE)
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

grid.arrange(corr_a, corr_b, nrow=1)

#-----------------------------------------------------------
# DIFFERENTIATION --> F_ST
#-----------------------------------------------------------

# F_ST statistic
# storage of Fst across sliding windows
n_repl <- length(unique(tags_index$repl))
fst_windowed_all <- matrix(0, nrow=n_repl, ncol=length(window_centers))

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
  
  #### IGNORING IF IT IS MISSING IN ONE POP OR THE OTHER ####
  #pos_both <- intersect(pos_p1, pos_p2)
  
  ms_both <- matrix(0, nrow=2*n_indiv, ncol=length(pos_both))
  colnames(ms_both) <- pos_both
  # top half of new matrix is p1 data, bottom half is p2 data. All missing rows in a population 0 by default
  ms_both[1:200,as.character(pos_p1)] <- ms_p1
  ms_both[201:400,as.character(pos_p2)] <- ms_p2
  # extract rows again (so ms_p1 and ms_p2 now have all positions, with missing positions filled with zeros)
  ms_p1 <- ms_both[1:200,]
  ms_p2 <- ms_both[201:400,]
  
  # storage
  hexp_df <- data.frame(pos=pos_both, p1=numeric(length(pos_both)), p2=numeric(length(pos_both)), total=numeric(length(pos_both)))
  fst_all <- data.frame(pos=pos_both, fst=numeric(length(pos_both)))
  
  for(i in 1:length(pos_both)){
    # extract columns at current positions
    ms_vect_p1 <- ms_p1[ ,as.character(pos_both[i])]
    ms_vect_p2 <- ms_p2[ ,as.character(pos_both[i])]
    # calc hexp as 2pq
    hexp_df$p1[i] <- 2 * mean(ms_vect_p1) * (1 - mean(ms_vect_p1))
    hexp_df$p2[i] <- 2 * mean(ms_vect_p2) * (1 - mean(ms_vect_p2))
    av_hexp <- mean(hexp_df$p1[i], hexp_df$p2[i])
    # combine ms data to get hexp of total metapopulation
    ms_vect_both <- c(ms_vect_p1, ms_vect_p2)
    hexp_total <- 2 * mean(ms_vect_both) * (1 - mean(ms_vect_both))
    hexp_df$total[i] <- hexp_total
    # calculate F_ST as (Ht - Hs) / Ht
    fst_all$fst[i] <- (hexp_total - av_hexp) / hexp_total
  }
  fst_all <- fst_all[fst_all$pos!=INV_START & fst_all$pos!=INV_END-1,]
  
  fst_windowed <- calc_sliding_window(fst_all, GENOME_LENGTH, windowSize = WINDOW_SIZE, pointSpacing = WINDOW_SPACING)
  fst_windowed_all[repl,] <- fst_windowed[,2]
}

fst_average <- data.frame(pos=window_centers, av_fst=colMeans(fst_windowed_all, na.rm = T), 
                          stdev=apply(fst_windowed_all, 2, sd, na.rm=T))

ggplot(fst_average, aes(x=pos, y=av_fst)) +
  geom_line() +
  scale_fill_gradient(low='white', high='blue') +
  ggtitle('F_ST between Populations') + 
  xlab('Position') + ylab(expression(F[ST])) +
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red')} +
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue')}

###### FST -- SEPARATING HaPLOTYPES

fst_windowed_haplotypes_all <- matrix(0, nrow=n_repl, ncol=length(window_centers))

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
  
  # prep storage
  hexp_df <- data.frame(pos=pos_both, normal=numeric(length(pos_both)), inverted=numeric(length(pos_both)), total=numeric(length(pos_both)))
  fst_all <- data.frame(pos=pos_both, fst=numeric(length(pos_both)))
  
  for(i in 1:length(pos_both)){
    # extract columns at current positions
    ms_vect_normal <- ms_normal[ ,as.character(pos_both[i])]
    ms_vect_inverted <- ms_inverted[ ,as.character(pos_both[i])]
    # calc hexp as 2pq
    hexp_df$normal[i] <- 2 * mean(ms_vect_normal) * (1 - mean(ms_vect_normal))
    hexp_df$inverted[i] <- 2 * mean(ms_vect_inverted) * (1 - mean(ms_vect_inverted))
    av_hexp <- mean(hexp_df$normal[i], hexp_df$inverted[i])
    # combine ms data to get hexp of total metapopulation
    ms_vect_both <- c(ms_vect_normal, ms_vect_inverted)
    hexp_total <- 2 * mean(ms_vect_both) * (1 - mean(ms_vect_both))
    hexp_df$total[i] <- hexp_total
    # calculate F_ST as (Ht - Hs) / Ht
    fst_all$fst[i] <- (hexp_total - av_hexp) / hexp_total
  }
  fst_all <- fst_all[fst_all$pos!=INV_START & fst_all$pos!=INV_END-1,]
  
  fst_windowed <- calc_sliding_window(fst_all, GENOME_LENGTH, windowSize = WINDOW_SIZE, pointSpacing = WINDOW_SPACING)
  fst_windowed_haplotypes_all[repl,] <- fst_windowed[,2]
}

fst_haplotypes_average <- data.frame(pos=window_centers, av_fst=colMeans(fst_windowed_haplotypes_all, na.rm = T), 
                          stdev=apply(fst_windowed_haplotypes_all, 2, sd, na.rm=T))

ggplot(fst_haplotypes_average, aes(x=pos, y=av_fst)) +
  geom_line() +
  scale_fill_gradient(low='white', high='blue') +
  ggtitle('F_ST between Haplotypes') + 
  xlab('Position') + ylab(expression(F[ST])) +
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red')} +
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue')}


# for testing different window sizes
WINDOW_SPACING <- 100
WINDOW_SIZE <- 50
window_centers <- seq(0, GENOME_LENGTH, by=WINDOW_SPACING)

