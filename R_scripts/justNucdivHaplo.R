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
WINDOW_SPACING <- 400
WINDOW_SIZE <- 400   # NOTE: window size is added on each side (so the full size is more like twice this value)
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

# calculating nucleotide diversity. Returns dataframe with nucdiv at spaced positions  --> NOTE: edited to remove marker mutations
# NOTE: seqLen is only really HALF the window size (adds seqLen in both directions of a point)
calc_nuc_div <- function(msdata, positions, totalLength, seqLen=200, centerSpacing=100){
  centers <- seq(0, totalLength, by=centerSpacing)
  # prepare storage for nucleotide diversity at each center position
  output <- data.frame(position=centers, nuc_div=NA)
  # if only one individual it becomes a vector. If just one individual, cannot calculate nucdiv (or rather, would be all zero)
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

get_breakpoint_indeces <- function(msdata, positions, breakpoints){
  inv_start <- breakpoints[1]
  inv_end <- breakpoints[2]
  
  inv_start_index <- which(positions==inv_start)
  inv_end_index <- which(positions==inv_end-1)
  
  # Fix for multiple mutations at a site
  # this first statement is if there are MORE THAN 2 mutations at a breakpoint, which is really weird. Dunno why thats happening
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

#-----------------------------------------------------------
# DATA EXTRACTION
#-----------------------------------------------------------

# read in files (values: selection coefficient, migration rate, replicate #)
files <- list.files(path=PATH, pattern="*.txt", full.names=F, recursive=FALSE)
n_files <- length(files)
# pre-calculate window centers' positions
window_centers <- seq(0, GENOME_LENGTH, by=WINDOW_SPACING)

# STORAGE DATAFRAMES
tags_index <- data.frame(population=character(n_files), sel_coef=numeric(n_files), migration=numeric(n_files), 
                         repl=integer(n_files), stringsAsFactors=F)
pos_frequency <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("pop", "position", "frequency"))
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

#---------------------------------------------------------------
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



if(on_cluster){
  ggsave('nucdiv_haplotypes_win400.png', plot_nucdiv_haplotypes, path=paste("Plots", args[1], args[2], sep="/"), width=8, height=6)
  save(nucdiv_all_long, file=paste("data_summary", args[1], args[2],"nucdiv_all_long_win200.Rds", sep="/"))
}else{
  # view plots (the ones created with grid.arrange are displayed automatically)
  print(plot_nucdiv_haplotypes)
}

