# determine whether on cluster or local
on_cluster <- Sys.info()['sysname']=='Linux'

if(on_cluster){
  # first argument is directory name, second is generation time point
  args = commandArgs(trailingOnly=TRUE)
}else{
  setwd("C:/Users/tbehr/Desktop/Thesis")
  }

library(ggplot2)
library(gridExtra)
library(reshape2)
library(tidyr)

#-----------------------------------------------------------
# PARAMETERS
#-----------------------------------------------------------

GENOME_LENGTH <- 120000
FIXED_MUTATION_POS1 <- 30000
FIXED_MUTATION_POS2 <- 70000
INV_START <- 10000
INV_END <- 110000  # this value should NOT be the '-1' value that the SLiM script uses. This script does that correction later
WINDOW_SPACING <- 100
WINDOW_SIZE <- 100   # NOTE: window size is added on each side (so the full size is more like twice this value)
N_TILES <- 600   # number of tiles along each axis of the correlation heatmap
FIRST_GEN <- 5000  # first generation where inversion/locally adapted alleles are introduced

if(on_cluster){
  PATH <- paste("Outputs", args[1], args[2], sep="/")
  simtype <- strsplit(args[1], split='_')[[1]][1]
  generation <- as.integer(args[2])
}else{
  PATH <- "Outputs/inversionLAA_2pop_s0.01_m0.001_mu1e-6/15000"
  simtype <- strsplit(strsplit(PATH, split='/')[[1]][2], split='_')[[1]][1]
  generation <- as.integer(strsplit(PATH, split='/')[[1]][3])
}

# record presence or absence of inversion and locally adapted alleles
INVERSION_PRESENT <- ifelse(simtype=='adaptiveInversion' || simtype=='inversionLAA' ,TRUE, FALSE)
LAA_PRESENT <- ifelse(simtype=='locallyAdapted' || simtype=='inversionLAA' ,TRUE, FALSE)

# reusable layer for ggplot to include marker lines for inversion bounds (blue) and locally adapted alleles (red)
gglayer_markers <- list(
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red', alpha=0.3)},
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue', alpha=0.4)}
)

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
  # if only one polymorphic site it becomes a vector. If just one individual, cannot calculate nucdiv (or rather, would be all zero)
  if(is.null(dim(msdata))){
    # msdata <- t(as.matrix(msdata))
    return(output)
  }
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
      # here using 2x sum of distances_all to account for both parts of the pairwise comparison matrix (above and below the diagonal)
      output[output$position==seqCenter,2] <- (2*sum(distances_all) / (num_of_seq^2)) / adjusted_seq_len
      
      # slightly different calculation which ignores comparisons along diagonal
      # output[output$position==seqCenter,2] <- (sum(distances_all) / choose(num_of_seq,2)) / adjusted_seq_len
    }
  }
  return(output)
}

# alternative function for calculating nucleotide diversity, using site frequency spectrum (SFS)
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
      
      sfs.raw <- table(colSums(ms_in_seq))
      
      # alternate to adjust so it always considers in respect to the less frequent allele
      #sfs.inverse <- num_of_seq - sfs.raw
      #sfs.total <- pmin(sfs.raw, sfs.inverse)

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

# calculating nucleotide diversity using the method in PopGenome package
# calculates nucdiv PER SITE first, and then averages those values across each window
calc_nuc_div_popgenome <- function(msdata, positions, totalLength, seqLen=200, centerSpacing=100){
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

# takes a dataframe where first column is position and second column is the value
calc_sliding_window <- function(posValData, totalLength, windowSize, pointSpacing){
  centers <- seq(0, totalLength, by=pointSpacing)
  output <- data.frame(position=centers, average=NA)
  for(i in 1:length(centers)){
    correspondingValues <- posValData[2][posValData[1] > centers[i]-windowSize & posValData[1] <= centers[i]+windowSize]
    output[i, 2] <- mean(correspondingValues, na.rm=T)
  }
  return(output)
}

# takes a dataframe where first column is position and second column is the value
calc_sliding_window_gaussian <- function(posValData, totalLength, windowSize, pointSpacing, stdev=100){
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

get_correlations <- function(msdata, positions, numTiles=20){
  # remove inversion marker mutations
  #positions <- positions[! positions %in% c(INV_START, INV_END-1)]
  #msdata <- msdata[, ! positions %in% c(INV_START, INV_END-1)]
  
  num_sites <- length(positions)
  # use default method (pearson)
  corr_all <- cor(msdata, method="pearson")
  # taking absolute value of correlation
  corr_all <- abs(corr_all)
  return(corr_all)
}

reduce_to_long <- function(corrData, positions, numTiles=20){
  # split positions into bins (using range up to full length so replicates can be combined)
  # remove inversion marker mutations
  #positions_reduced <- positions[! positions %in% c(INV_START, INV_END-1)]
  positions_reduced <- positions
  #corrData <- corrData[! positions %in% c(INV_START, INV_END-1), ! positions %in% c(INV_START, INV_END-1)]
  
  groups <- cut(c(0, positions_reduced, GENOME_LENGTH), breaks=numTiles, labels=F)
  # associate positions with their groups. Remove first and last group which were only included to specify range
  pos_grouping <- data.frame(position=positions_reduced, group=groups[-c(1,length(groups))]) 
  
  # convert to long, then convert position indeces to corresponding group numbers
  data_long <- melt(corrData)
  
  data_long$Var1 <- pos_grouping$group[data_long$Var1]
  data_long$Var2 <- pos_grouping$group[data_long$Var2]
  
  #data_long$Var1 <- pos_grouping$group[match(data_long$Var1, pos_grouping$position)]
  #data_long$Var2 <- pos_grouping$group[match(data_long$Var2, pos_grouping$position)]
  
  # averaging correlations within each combination (i,j) of bins
  red_long_incomplete <- aggregate(value ~ Var1 + Var2, data=data_long, FUN=mean, drop=F, na.rm=T)
  # By chance, some bins may be empty, so dataframe including all possible tile coordinates is combined with aggregated means
  red_long_empty <- data.frame(Var1=rep(1:numTiles, numTiles), Var2=sort(rep(1:numTiles, numTiles)))
  red_long <- merge(red_long_empty, red_long_incomplete, by=c('Var1', 'Var2'), all=T)
  # scale group numbers to nucleotide positions
  red_long[1] <- red_long[1] * (max(GENOME_LENGTH)/numTiles)
  red_long[2] <- red_long[2] * (max(GENOME_LENGTH)/numTiles)

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

# get frequency of locally adapted alleles. --> doesn't work if multiple mutations at a site
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

# calculate differentiation (Fst) between two groups in ms format, with column names of their positions
# Note: removes inversion markers
calc_fst_between <- function(msGroup1, msGroup2){
  # extract positions from column names
  allPositions <- unique(sort(c(as.integer(colnames(msGroup1)), as.integer(colnames(msGroup2)))))
  
  # storage
  hexp_df <- data.frame(pos=allPositions, group1=numeric(length(allPositions)), group2=numeric(length(allPositions)), 
                        total=numeric(length(allPositions)))
  fst_all <- data.frame(pos=allPositions, fst=numeric(length(allPositions)))
  
  # convert to matrix of one row if the msdata has only one sample (and hence was converted to a vector)
  if(is.null(dim(msGroup1))){
    msGroup1 <- t(as.matrix(msGroup1))
  }
  if(is.null(dim(msGroup2))){
    msGroup2 <- t(as.matrix(msGroup2))
  }
  
  for(i in 1:length(allPositions)){
    # extract columns at current positions
    ms_vect_1 <- msGroup1[ ,i]
    ms_vect_2 <- msGroup2[ ,i]
    # calc hexp as 2pq
    hexp_df$group1[i] <- 2 * mean(ms_vect_1) * (1 - mean(ms_vect_1))
    hexp_df$group2[i] <- 2 * mean(ms_vect_2) * (1 - mean(ms_vect_2))
    av_hexp <- mean(c(hexp_df$group1[i], hexp_df$group2[i]), na.rm = TRUE)
    # combine ms data to get hexp of total metapopulation
    ms_vect_both <- c(ms_vect_1, ms_vect_2)
    hexp_total <- 2 * mean(ms_vect_both) * (1 - mean(ms_vect_both))
    hexp_df$total[i] <- hexp_total
    # calculate F_ST as (Ht - Hs) / Ht
    fst_all$fst[i] <- (hexp_total - av_hexp) / hexp_total
  }
  # remove inversion markers
  fst_all <- fst_all[fst_all$pos!=INV_START & fst_all$pos!=INV_END-1,]
  return(fst_all)
}

# quick way to count # of windows and see the ranges
get_windows_ranges <- function(){
  centers <- seq(0, GENOME_LENGTH, by=WINDOW_SPACING)
  
  range_list <- character(length=length(centers))
  for(i in 1:length(centers)){
    range_list[i] <- paste(centers[i]-WINDOW_SIZE, "-", centers[i]+WINDOW_SIZE)
  }
  print(paste('Number of Windows:', length(centers)))
  return(range_list)
}

# get the indeces of the breakpoints in a vector of absolute positions
# requires the ms data to disentangle when multiple polymorphisms exist at one or both breakpoints
get_breakpoint_indeces <- function(msdata, positions, breakpoints){
  inv_start <- breakpoints[1]
  inv_end <- breakpoints[2]
  
  inv_start_index <- which(positions==inv_start)
  inv_end_index <- which(positions==inv_end-1)
  
  # Fix for multiple mutations at a site
  # this first statement is if there are MORE THAN 2 mutations at a breakpoint, which is kinda weird. In code, added warnings for when this happens
  if(length(inv_start_index)>2 | length(inv_end_index)>2){
    return(c(NA, NA))
  } else if(length(inv_start_index)==2 & length(inv_end_index)==2){
    # if both indeces are duplicated, need to find the pair of columns that are identical
    comparison_indeces <- which(colSums(msdata[,inv_start_index]!=msdata[,inv_end_index])==0)
    if(length(comparison_indeces)==0){
      # if no columns are the same, then flip one of the matrices for the other two comparisons
      inv_start_index <- inv_start_index[c(2,1)]
    }
    inv_start_index <- inv_start_index[comparison_indeces[1]]  # this last index is in case all columns are identical
    inv_end_index <- inv_end_index[comparison_indeces[1]]
    
  } else if(length(inv_start_index)>1){
    # if only one index is duplicated, pick the index that is identical to the ms of the single index
    # works by taking the index of the comparison matrix with sum of zero (ie. no differences). 
    # Index at end is needed if its identical to both, in which case in doesn't matter which to take
    inv_start_index <- inv_start_index[which(colSums(ms_binary[,inv_end_index]!=ms_binary[,inv_start_index])==0)[1]]
  } else if(length(inv_end_index)>1){
    # same as previous else if, but if the end breakpoint is duplicated
    inv_end_index <- inv_end_index[which(colSums(ms_binary[,inv_start_index]!=ms_binary[,inv_end_index])==0)[1]]
  }
  
  return(c(inv_start_index, inv_end_index))
}

# separate ms file into two separate ms (outputted as list along with new positions) by haplotype
# also removes marker mutations
split_ms_by_haplotype <- function(msdata, positions, breakpoints, indeces){
  inv_start_index <- indeces[1]
  inv_end_index <- indeces[2]
  inv_start <- breakpoints[1]
  inv_end <- breakpoints[2] - 1
  
  ms_normal <- ms_binary[ms_binary[,inv_start_index]==0 & ms_binary[,inv_end_index]==0, ]
  ms_inverted <- ms_binary[ms_binary[,inv_start_index]==1 & ms_binary[,inv_end_index]==1, ]
  
  # convert to matrix of one row if the msdata has only one sample (and hence was converted to a vector)
  if(is.null(dim(ms_normal))){
    ms_normal <- t(as.matrix(ms_normal))
  }
  if(is.null(dim(ms_inverted))){
    ms_inverted <- t(as.matrix(ms_inverted))
  }
  
  # remove marker mutations
  ms_normal <- ms_normal[,! abs_positions %in% c(inv_start, inv_end)]
  ms_inverted <- ms_inverted[,! abs_positions %in% c(inv_start, inv_end)]
  positions_reduced <- abs_positions[! abs_positions %in% c(inv_start, inv_end)]
  
  if(is.null(dim(ms_normal))){
    ms_normal <- t(as.matrix(ms_normal))
  }
  if(is.null(dim(ms_inverted))){
    ms_inverted <- t(as.matrix(ms_inverted))
  }
  
  colnames(ms_normal) <- positions_reduced
  colnames(ms_inverted) <- positions_reduced
  
  return(list(ms_normal, ms_inverted, positions_reduced))
}

# function for easily loading and renaming a Rds file. Taken from stackoverflow
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#-----------------------------------------------------------
# DATA EXTRACTION
#-----------------------------------------------------------

# read in files (values: selection coefficient, migration rate, replicate #)
files <- list.files(path=PATH, pattern="*.txt", full.names=F, recursive=FALSE)
n_files <- length(files)
# pre-calculate window centers' positions
# Ideally is perfectly divisble so first and last windows lie on the first and last positions of the genome
if(GENOME_LENGTH%%WINDOW_SPACING){warning("Genome length not divisible by window spacing: Window centers don't fit cleanly")}
window_centers <- seq(0, GENOME_LENGTH, by=WINDOW_SPACING)

# STORAGE DATAFRAMES
tags_index <- data.frame(population=character(n_files), sel_coef=numeric(n_files), migration=numeric(n_files), 
                         repl=integer(n_files), stringsAsFactors=F)
pos_frequency <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("pop", "position", "frequency"))
hexp_red_df <- matrix(0, nrow=n_files, ncol=length(window_centers))
nucdiv_df <- matrix(0, nrow=n_files, ncol=length(window_centers))
frequencies <- matrix(0, nrow=n_files, ncol=length(window_centers))
# correlations_3d <- array(numeric(), dim=c(N_TILES, N_TILES, n_files))  # NOTE: end up switching a lot between long and wide here
#      maybe change so its all in long?

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
  # extract metadata from filename
  tags <- strsplit(files[i], split='_')[[1]]
  tags_index[i,] <- list(tags[2], as.numeric(tags[3]), as.numeric(tags[4]), as.integer(tags[5]))
  
  # get allele frequencies at all positions
  pos_frequency_subset <- data.frame(pop=tags[2], position=abs_positions, frequency=colMeans(ms_binary))
  pos_frequency <- rbind(pos_frequency, pos_frequency_subset)
  
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
  
  # # correlation matrix into a 3D array (third dimension is file index)
  # corr_data <- get_correlations(ms_binary, abs_positions, numTiles = N_TILES)
  # corr_long <- reduce_to_long(corr_data, abs_positions, numTiles = N_TILES)
  # correlations_3d[,,i] <- as.matrix(dcast(corr_long, Var1 ~ Var2)[,-1]) # exclude first column (variable names)
}

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

#(FOR JUST ONE REPLCIATE)
#nucdiv_summ_p1 <- data.frame(center = window_centers, nucdiv = nucdiv_df[which(tags_index$population=='p1'),], stdev=NA)

# find max and min values between both populations to get shared axis
min_hexp <- min(hexp_summ_p1$hexp, hexp_summ_p2$hexp)
max_hexp <- max(hexp_summ_p1$hexp, hexp_summ_p2$hexp)
min_nucdiv <- min(nucdiv_summ_p1$nucdiv, nucdiv_summ_p2$nucdiv)
max_nucdiv <- max(nucdiv_summ_p1$nucdiv, nucdiv_summ_p2$nucdiv)

hexp_a <- ggplot(hexp_summ_p1, aes(x=center, y=hexp)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('a) Expected Heterozygosity - P1') +
  gglayer_markers +
  xlab('Position') + ylab(expression(H[exp])) +
  ylim(c(min_hexp, max_hexp))

hexp_b <- ggplot(hexp_summ_p2, aes(x=center, y=hexp)) +
  geom_line() +
  scale_color_brewer(palette="Dark2") +
  ggtitle('b) Expected Heterozygosity - P2') +
  gglayer_markers +
  xlab('Position') + ylab(expression(H[exp])) +
  ylim(c(min_hexp, max_hexp))

nucdiv_a <- ggplot(nucdiv_summ_p1, aes(x=center, y=nucdiv)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('c) Nucleotide Diversity - P1') +
  gglayer_markers +
  xlab('Position') + ylab('\u03c0') +
  ylim(c(min_nucdiv, max_nucdiv))

nucdiv_b <- ggplot(nucdiv_summ_p2, aes(x=center, y=nucdiv)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('d) Nucleotide Diversity - P2') +
  gglayer_markers +
  #geom_errorbar(aes(ymin=nucdiv-stdev, ymax=nucdiv+stdev), width=1, position=position_dodge(0.1)) +
  xlab('Position') + ylab('\u03c0') +
  ylim(c(min_nucdiv, max_nucdiv))

plot_nucdiv_hexp <- grid.arrange(hexp_a, hexp_b, nucdiv_a, nucdiv_b, nrow=2)

######  nucleotide diversity - separating haplotypes #######

nucdiv_inverted <- matrix(0, nrow=n_files, ncol=length(window_centers))
nucdiv_normal <- matrix(0, nrow=n_files, ncol=length(window_centers))

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
  # if at least one sample individual has the inversion
  if(all(c(INV_START, INV_END-1) %in% abs_positions)){
    
    breakpoint_indeces <- get_breakpoint_indeces(ms_binary, abs_positions, c(INV_START, INV_END))
    if(any(is.na(breakpoint_indeces))){
      warning(paste("More than two mutations at a breakpoint. Did not determine true breakpoint, setting to NA -", filepath))}
    
    inv_start_index <- breakpoint_indeces[1]
    inv_end_index <- breakpoint_indeces[2]
    
    # if both NA, just skip the loop
    if(all(is.na(breakpoint_indeces))){
      nucdiv_inverted[i,] <- NA
      nucdiv_normal[i,] <- NA
      next
    }
    
    # extract ms rows based on presence of inversion markers
    ms_split <- split_ms_by_haplotype(ms_binary, abs_positions, c(INV_START, INV_END), breakpoint_indeces)
    
    ms_normal <- ms_split[[1]]
    ms_inverted <- ms_split[[2]]
    abs_positions <- ms_split[[3]]

    nucdiv_normal_windowed <- calc_nuc_div(ms_normal, abs_positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
    nucdiv_inverted_windowed <- calc_nuc_div(ms_inverted, abs_positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
    
    nucdiv_normal[i,] <- nucdiv_normal_windowed[[2]]
    nucdiv_inverted[i,] <- nucdiv_inverted_windowed[[2]]
    
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

plot_nucdiv_haplotypes <- ggplot(nucdiv_all_long, aes(x=center, y=value, col=variable)) +
  geom_line() +
  ggtitle('Nucleotide Diversity for Different Haplotypes') +
  gglayer_markers

#-----------------------------------------------------------
# CORRELATION
#-----------------------------------------------------------

# # correlation heatmap
# corr_summ_p1 <- apply(correlations_3d[, , which(tags_index$population == "p1")], c(1, 2), mean, na.rm = TRUE)
# corr_summ_p1_long <- melt(corr_summ_p1)
# # correct group values to bin centers
# bin_size <- GENOME_LENGTH / N_TILES
# corr_summ_p1_long$Var1 <- corr_summ_p1_long$Var1*bin_size - bin_size/2
# corr_summ_p1_long$Var2 <- corr_summ_p1_long$Var2*bin_size - bin_size/2
# 
# corr_summ_p2 <- apply(correlations_3d[, , which(tags_index$population == "p2")], c(1, 2), mean, na.rm = TRUE)
# corr_summ_p2_long <- melt(corr_summ_p2)
# # correct group values to bin centers
# bin_size <- GENOME_LENGTH / N_TILES
# corr_summ_p2_long$Var1 <- corr_summ_p2_long$Var1*bin_size - bin_size/2
# corr_summ_p2_long$Var2 <- corr_summ_p2_long$Var2*bin_size - bin_size/2
# 
# corr_a <- ggplot(corr_summ_p1_long, aes(x=Var1, y=Var2, fill=value)) +
#   geom_tile() +
#   scale_fill_gradient(low='white', high='blue') +
#   ggtitle('P1') + xlab('Position') + ylab('Position')
# 
# corr_b <- ggplot(corr_summ_p2_long, aes(x=Var1, y=Var2, fill=value)) +
#   geom_tile() +
#   scale_fill_gradient(low='white', high='blue') +
#   ggtitle('P2') + xlab('Position') + ylab('Position')
# 
# plot_correlation <- grid.arrange(corr_a, corr_b, nrow=1)

#-----------------------------------------------------------
# Breakpoint Linkage 
#-----------------------------------------------------------
if(INVERSION_PRESENT && generation > FIRST_GEN){
  breakpoint_corr_windowed_all <- matrix(0, nrow=n_files, ncol=length(window_centers))
  breakpoint_corr_windowed_all_filt <- matrix(0, nrow=n_files, ncol=length(window_centers))
  
  for(i in 1:n_files){
    filepath <- paste0(PATH, "/", files[i])
    ms_binary <- get_ms_data(filepath)
    abs_positions <- get_positions(filepath)
    
    breakpoint_indeces <- get_breakpoint_indeces(ms_binary, abs_positions, c(INV_START, INV_END))
    if(any(is.na(breakpoint_indeces))){
      warning(paste("More than two mutations at a breakpoint. Did not determine true breakpoint, setting to NA -", filepath))}
    
    inv_start_index <- breakpoint_indeces[1]
    inv_end_index <- breakpoint_indeces[2]
    
    # if both NA, just skip the loop
    if(all(is.na(breakpoint_indeces))){
      breakpoint_corr_windowed_all[i,] <- NA
      breakpoint_corr_windowed_all_filt[i,] <- NA
      next
    }
    
    breakpoint_vector <- ms_binary[,inv_start_index]
    
    ##### filter out low freq. alleles (as they may be perfectly correlated with breakpoint)
    allele_frequencies <- colMeans(ms_binary)
    
    ms_binary_filtered <- ms_binary[, allele_frequencies > 0.1 & allele_frequencies < 0.9]
    abs_positions_filtered <- abs_positions[allele_frequencies > 0.1 & allele_frequencies < 0.9]
    #####
    
    breakpoint_corr <- cor(breakpoint_vector, ms_binary)
    breakpoints_corr_abs <- abs(breakpoint_corr)
    breakpoints_corr_df <- data.frame(pos=abs_positions, corr=as.vector(breakpoints_corr_abs))
    breakpoint_corr_windowed <- calc_sliding_window(breakpoints_corr_df, GENOME_LENGTH, WINDOW_SIZE, WINDOW_SPACING)
    breakpoint_corr_windowed_all[i,] <- breakpoint_corr_windowed$average
    
    breakpoint_corr_filt <- cor(breakpoint_vector, ms_binary_filtered)
    breakpoints_corr_abs_filt <- abs(breakpoint_corr_filt)
    breakpoints_corr_df_filt <- data.frame(pos=abs_positions_filtered, corr=as.vector(breakpoints_corr_abs_filt))
    breakpoint_corr_windowed_filt <- calc_sliding_window(breakpoints_corr_df_filt, GENOME_LENGTH, WINDOW_SIZE, WINDOW_SPACING)
    breakpoint_corr_windowed_all_filt[i,] <- breakpoint_corr_windowed_filt$average
  }
  
  breakpoints_corr_mean <- data.frame(pos=window_centers, corr_mean=colMeans(breakpoint_corr_windowed_all, na.rm=T))
  plot_corr_breakpoints <- ggplot(breakpoints_corr_mean, aes(x=pos, y=corr_mean)) + geom_line() + gglayer_markers
  
  breakpoints_corr_mean_filt <- data.frame(pos=window_centers, corr_mean=colMeans(breakpoint_corr_windowed_all_filt, na.rm=T))
  plot_corr_breakpoints_filt <- ggplot(breakpoints_corr_mean_filt, aes(x=pos, y=corr_mean)) + geom_line() + gglayer_markers
  
  if(on_cluster){
    ggsave('corr_breakpoint_win100.png', plot_corr_breakpoints, device="png", path=paste("Plots", args[1], args[2], sep="/"), width=12, height=8)
    ggsave('corr_breakpoint_filtered_win100.png', plot_corr_breakpoints_filt, device="png", path=paste("Plots", args[1], args[2], sep="/"), width=12, height=8)
  }else{
    print(plot_corr_breakpoints)
    print(plot_corr_breakpoints_filt)
  }
}

#-----------------------------------------------------------
# DIFFERENTIATION --> F_ST
#-----------------------------------------------------------

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
  ms_both[1:n_indiv,as.character(pos_p1)] <- ms_p1
  ms_both[(n_indiv+1):(n_indiv*2),as.character(pos_p2)] <- ms_p2
  # extract rows again (so ms_p1 and ms_p2 now have all positions, with missing positions filled with zeros)
  ms_p1 <- ms_both[1:n_indiv,]
  ms_p2 <- ms_both[(n_indiv+1):(n_indiv*2),]
  
  fst_all <- calc_fst_between(ms_p1, ms_p2)
  
  fst_windowed <- calc_sliding_window(fst_all, GENOME_LENGTH, windowSize = WINDOW_SIZE, pointSpacing = WINDOW_SPACING)
  fst_windowed_all[repl,] <- fst_windowed[,2]
}

fst_average <- data.frame(pos=window_centers, av_fst=colMeans(fst_windowed_all, na.rm = T), 
                          stdev=apply(fst_windowed_all, 2, sd, na.rm=T))

plot_fst <- ggplot(fst_average, aes(x=pos, y=av_fst)) +
  geom_line() +
  scale_fill_gradient(low='white', high='blue') +
  ggtitle('F_ST between Populations') + 
  xlab('Position') + ylab(expression(F[ST])) +
  gglayer_markers

##### separating by both haplotype and population
if(INVERSION_PRESENT && generation > FIRST_GEN){
  fst_hudson_windowed_all <- matrix(0, nrow=n_files, ncol=length(window_centers))
  
  for(i in 1:n_files){
    filepath <- paste0(PATH, "/", files[i])
    ms_binary <- get_ms_data(filepath)
    abs_positions <- get_positions(filepath)
    # extract metadata from filename
    tags <- strsplit(files[i], split='_')[[1]]
    tags_index[i,] <- list(tags[2], as.numeric(tags[3]), as.numeric(tags[4]), as.integer(tags[5]))
    
    if(INVERSION_PRESENT && generation > FIRST_GEN){
      breakpoint_indeces <- get_breakpoint_indeces(ms_binary, abs_positions, c(INV_START, INV_END))
      if(any(is.na(breakpoint_indeces))){
        warning(paste("More than two mutations at a breakpoint. Did not determine true breakpoint, setting to NA -", filepath))}
      
      inv_start_index <- breakpoint_indeces[1]
      inv_end_index <- breakpoint_indeces[2]
      
      # if both NA, just skip the loop
      if(all(is.na(breakpoint_indeces))){
        fst_hudson_windowed_all[i,] <- NA
        next
      }
      
      # extract ms rows based on presence of inversion markers
      ms_split <- split_ms_by_haplotype(ms_binary, abs_positions, c(INV_START, INV_END), breakpoint_indeces)
      abs_positions <- abs_positions[! abs_positions %in% c(INV_START, INV_END-1)]
      
      ms_normal <- ms_split[[1]]
      ms_inverted <- ms_split[[2]]
      positions_reduced <- ms_split[[3]]
      
      normal_is_valid <- TRUE
      inverted_is_valid <- TRUE
      
      if(dim(ms_normal)[1]<=1){
        normal_is_valid <- FALSE
      }
      if(dim(ms_inverted)[1]<=1){
        inverted_is_valid <- FALSE
      }
      if(normal_is_valid & inverted_is_valid){
        freq_nor <- colMeans(ms_normal)
        freq_inv <- colMeans(ms_inverted)
        # n is sample size
        n_nor <- dim(ms_normal)[[1]]
        n_inv <- dim(ms_inverted)[[1]]
        # Hudson estimator for Fst
        fst_hudson_numerator <- ((freq_nor - freq_inv)^2) - ((freq_nor*(1-freq_nor))/(n_nor-1)) - ((freq_inv*(1-freq_inv))/(n_inv-1))
        fst_hudson_demoninator <- (freq_nor*(1-freq_inv)) + (freq_inv*(1-freq_nor))
        fst_hudson <- data.frame(pos=positions_reduced, fst=fst_hudson_numerator / fst_hudson_demoninator)
        
        fst_hudson_windowed <- calc_sliding_window(fst_hudson, GENOME_LENGTH, WINDOW_SIZE, WINDOW_SPACING)
        
        fst_hudson_windowed_all[i,] <- fst_hudson_windowed[[2]]
      }
    }
  }
  
  fst_windowed_p1 <- fst_hudson_windowed_all[tags_index$population=='p1',]
  fst_windowed_p2 <- fst_hudson_windowed_all[tags_index$population=='p2',]
  
  fst_p1_average <- data.frame(pos=window_centers, av_fst=colMeans(fst_windowed_p1, na.rm = T),
                               stdev=apply(fst_windowed_p1, 2, sd, na.rm=T))
  fst_p2_average <- data.frame(pos=window_centers, av_fst=colMeans(fst_windowed_p2, na.rm = T),
                               stdev=apply(fst_windowed_p2, 2, sd, na.rm=T))
  
  plot_fst_p1 <- ggplot(fst_p1_average, aes(x=pos, y=av_fst)) +
    geom_line()+
    scale_fill_gradient(low='white', high='blue') +
    ggtitle('P1 - Between Haplotypes') +
    xlab('Position') + ylab(expression(F[ST])) +
    gglayer_markers
  
  plot_fst_p2 <- ggplot(fst_p2_average, aes(x=pos, y=av_fst)) +
    geom_line() +
    scale_fill_gradient(low='white', high='blue') +
    ggtitle('P2 - Between Haplotypes') +
    xlab('Position') + ylab(expression(F[ST])) +
    gglayer_markers
  
  plot_fst_hudson <- grid.arrange(plot_fst_p1, plot_fst_p2, nrow=1)
  
  if(on_cluster){
    ggsave('fst_hudson_win100.png', plot_fst_hudson, device="png", path=paste("Plots", args[1], args[2], sep="/"), width=20, height=8)
  }else{
    print(plot_fst_hudson)
  }
}

#-----------------------------------------------------------
# FINALIZATION - PLOTTING
#-----------------------------------------------------------

if(on_cluster){
  ggsave('nucdiv_hexp_win100.png', plot_nucdiv_hexp, device="png", path=paste("Plots", args[1], args[2], sep="/"), width=16, height=12)
  ggsave('nucdiv_haplotypes_win100.png', plot_nucdiv_haplotypes, device="png", path=paste("Plots", args[1], args[2], sep="/"), width=18, height=12)
  #ggsave('correlation.png', plot_correlation, path=paste("Plots", args[1], args[2], sep="/"), width=12, height=5.5)
  ggsave('fst_pops.png', plot_fst, device="png", path=paste("Plots", args[1], args[2], sep="/"), width=12, height=9)
}else{
  # view plots (the ones created with grid.arrange are displayed automatically)
  print(plot_nucdiv_haplotypes)
  print(plot_fst)
}

if(on_cluster){
  save(hexp_summ_p1, file=paste("data_summary", args[1], args[2],"hexp_summ_p1.Rds", sep="/"))
  save(hexp_summ_p2, file=paste("data_summary", args[1], args[2],"hexp_summ_p2.Rds", sep="/"))
  save(nucdiv_summ_p1, file=paste("data_summary", args[1], args[2],"nucdiv_summ_p1.Rds", sep="/"))
  save(nucdiv_summ_p2, file=paste("data_summary", args[1], args[2],"nucdiv_summ_p2.Rds", sep="/"))
  save(nucdiv_all_long, file=paste("data_summary", args[1], args[2],"nucdiv_all_long.Rds", sep="/"))
  
  if(INVERSION_PRESENT && generation > FIRST_GEN){
    save(breakpoints_corr_mean_filt, file=paste("data_summary", args[1], args[2],"breakpoints_corr_mean_filt.Rds", sep="/"))
    save(fst_p1_average, file=paste("data_summary", args[1], args[2],"fst_p1_average.Rds", sep="/"))
    save(fst_p2_average, file=paste("data_summary", args[1], args[2],"fst_p2_average.Rds", sep="/"))
  }
  
  save(fst_average, file=paste("data_summary", args[1], args[2],"fst_average.Rds", sep="/"))

}