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

GENOME_LENGTH <- 22000
FIXED_MUTATION_POS1 <- 8000
FIXED_MUTATION_POS2 <- 12000
INV_START <- 6000
INV_END <- 16000  # this value should NOT be the '-1' value that the SLiM script uses. This script does that correction later
WINDOW_SPACING <- 100
WINDOW_SIZE <- 100   # NOTE: window size is added on each side (so the full size is more like twice this value)
N_TILES <- 200    # number of tiles along each axis of the correlation heatmap

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
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red')},
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
    output[i, 2] <- mean(correspondingValues)
  }
  return(output)
}

# MODIFIED TO REMOVE MARKER MUTATIONS
get_correlations <- function(msdata, positions, numTiles=20){
  # remove inversion marker mutations
  positions <- positions[! positions %in% c(INV_START, INV_END-1)]
  msdata <- msdata[, ! positions %in% c(INV_START, INV_END-1)]
  
  num_sites <- length(positions)
  # use default method (pearson)
  corr_all <- cor(msdata, method="pearson")
  # taking absolute value of correlation
  corr_all <- abs(corr_all)
  return(corr_all)
}

#MODIFIED TO REMOVE MARKER MUTATIONS
reduce_to_long <- function(corrData, positions, numTiles=20){
  # split positions into bins (using range up to full length so replicates can be combined)
  # remove inversion marker mutations
  positions_reduced <- positions[! positions %in% c(INV_START, INV_END-1)]
  
  corrData <- corrData[! positions %in% c(INV_START, INV_END-1), ! positions %in% c(INV_START, INV_END-1)]
  
  groups <- cut(c(0, positions_reduced, GENOME_LENGTH), breaks=numTiles, labels=F)
  # associate positions with their groups. Remove first and last group which were only included to specify range
  pos_grouping <- data.frame(position=positions_reduced, group=groups[-c(1,length(groups))]) 
  
  # convert to long, then convert position indeces to corresponding group numbers
  data_long <- melt(corrData)
  
  #data_long$Var1 <- pos_grouping$group[data_long$Var1]
  #data_long$Var2 <- pos_grouping$group[data_long$Var2]
  
  data_long$Var1 <- pos_grouping$group[match(data_long$Var1, pos_grouping$position)]
  data_long$Var2 <- pos_grouping$group[match(data_long$Var2, pos_grouping$position)]
  
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
corr_outliers_all <- data.frame(pos1=numeric(), pos2=numeric(), value=numeric())

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
  # extract metadata from filename
  tags <- strsplit(files[i], split='_')[[1]]
  tags_index[i,] <- list(tags[2], as.numeric(tags[3]), as.numeric(tags[4]), as.integer(tags[5]))
  
  corr_data <- get_correlations(ms_binary, abs_positions, numTiles = N_TILES)
  data_long <- melt(corr_data)
  
  mean_corr <- mean(data_long$value, na.rm=TRUE)
  sd_corr <- sd(data_long$value, na.rm=TRUE)
  
  corr_outliers <- subset(data_long, value > (mean_corr + 2*sd_corr))
  
  corr_outliers_all <- rbind(corr_outliers_all, corr_outliers)
}

corr_outliers_all <- subset(corr_outliers_all, value<1)

#outlier_scatterplot <- ggplot(corr_outliers_all, aes(x=Var1, y=Var2)) + geom_point(alpha=0.01, pch='.')

outlier_hexplot <- ggplot(corr_outliers_all, aes(x=Var1, y=Var2)) +
  stat_binhex()


if(on_cluster){
  #ggsave('corr_outlier_scatter.png', outlier_scatterplot, path=paste("Plots", args[1], args[2], sep="/"), width=8, height=6)
  ggsave('corr_outlier_hex.png', outlier_hexplot, path=paste("Plots", args[1], args[2], sep="/"), width=8, height=6)
}


