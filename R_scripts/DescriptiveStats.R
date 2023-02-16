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
  PATH <- "Outputs/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000"
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
frequencies <- matrix(0, nrow=n_files, ncol=length(window_centers))
neutral_frequencies <- matrix(0, nrow = n_files, ncol = 2)
freq_allele1 <- numeric(n_files)
freq_allele2 <- numeric(n_files)
inv_freq_all <- numeric(n_files)
polymorphism_counts <- numeric(n_files)
# correlations_3d <- array(numeric(), dim=c(N_TILES, N_TILES, n_files))  # NOTE: end up switching a lot between long and wide here
#      maybe change so its all in long?

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
  # extract metadata from filename
  tags <- strsplit(files[i], split='_')[[1]]
  tags_index[i,] <- list(tags[2], as.numeric(tags[3]), as.numeric(tags[4]), as.integer(tags[5]))
  
  breakpoint_indeces <- get_breakpoint_indeces(ms_binary, abs_positions, c(INV_START, INV_END))
  
  # get allele frequencies at all positions
  pos_frequency_subset <- data.frame(pop=tags[2], position=abs_positions, frequency=colMeans(ms_binary))
  pos_frequency <- rbind(pos_frequency, pos_frequency_subset)
  
  # calc mutation frequency, then sliding window
  pos_freq <- data.frame(position=abs_positions, freq=colMeans(ms_binary))
  frequencies_windowed <- calc_sliding_window(pos_freq, GENOME_LENGTH, windowSize = WINDOW_SIZE, pointSpacing = WINDOW_SPACING)
  frequencies[i,] <- frequencies_windowed[[2]]
  
  ms_neutrals <- ms_binary[,!(abs_positions %in% c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2, INV_START, INV_END-1))]
  neutral_frequencies[i, 1] <- mean(ms_neutrals)
  neutral_frequencies[i, 2] <- dim(ms_neutrals)[2]
  
  if(FIXED_MUTATION_POS1 %in% abs_positions){
    freq_allele1[i] <- mean(ms_binary[,tail(which(abs_positions==FIXED_MUTATION_POS1), n=1)])  # tail() is in case there are multiple mutations at a site
  }
  if(FIXED_MUTATION_POS2 %in% abs_positions){
    freq_allele2[i] <- mean(ms_binary[,tail(which(abs_positions==FIXED_MUTATION_POS2), n=1)])
  }
  
  polymorphism_counts[i] <- length(abs_positions)

  inv_freq_all[i] <- mean(ms_binary[,breakpoint_indeces[1]])
}


overall_neutral_frequency <- sum(apply(neutral_frequencies, 1, prod) / sum(neutral_frequencies[,2]))

# mean allele frequency 
allele1_freq_p1 <- mean(freq_allele1[tags_index$population=='p1'], na.rm=T)
allele1_freq_p2 <- mean(freq_allele1[tags_index$population=='p2'], na.rm=T)
allele2_freq_p1 <- mean(freq_allele2[tags_index$population=='p1'], na.rm=T)
allele2_freq_p2 <- mean(freq_allele2[tags_index$population=='p2'], na.rm=T)

# polymorphism count
average_polymorphism_count <- mean(polymorphism_counts)
stdev_polymorphism_count <- sd(polymorphism_counts)

# inversion frequency
inv_freq_p1 <- mean(inv_freq_all[tags_index$population=='p1'], na.rm=T)
inv_freq_p2 <- mean(inv_freq_all[tags_index$population=='p2'], na.rm=T)
inv_stdev_p1 <- sd(inv_freq_all[tags_index$population=='p1'], na.rm=T)
inv_stdev_p2 <- sd(inv_freq_all[tags_index$population=='p2'], na.rm=T)

print(paste("Overall neutral frequency:", overall_neutral_frequency))
print(paste("Average polymorphism count:", average_polymorphism_count))
print(paste("STdev polymorphism count:", stdev_polymorphism_count))
print(paste("LAA 1 Freq, P1:", allele1_freq_p1))
print(paste("LAA 2 Freq, P1:", allele2_freq_p1))
print(paste("LAA 1 Freq, P2:", allele1_freq_p2))
print(paste("LAA 2 Freq, P2:", allele2_freq_p2))
print(paste("Inversion Freq P1:", inv_freq_p1, "with SD:", inv_stdev_p1))
print(paste("Inversion Freq P2:", inv_freq_p2, "with SD:", inv_stdev_p2))


