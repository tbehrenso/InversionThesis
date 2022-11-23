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
  
  #data_long$Var1 <- pos_grouping$group[pos_grouping$position==data_long$Var1]
  #data_long$Var2 <- pos_grouping$group[pos_grouping$position==data_long$Var2]
  
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
correlations_3d <- array(numeric(), dim=c(N_TILES, N_TILES, n_files))

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
  # extract metadata from filename
  tags <- strsplit(files[i], split='_')[[1]]
  tags_index[i,] <- list(tags[2], as.numeric(tags[3]), as.numeric(tags[4]), as.integer(tags[5]))
  
  inv_start_index <- which(abs_positions==INV_START)
  inv_end_index <- which(abs_positions==INV_END-1)
  
  # Fix for multiple mutations at a site
  # this first statement is if there are MORE THAN 2 mutations at a breakpoint, which is really weird. Dunno why thats happening
  if(length(inv_start_index)>2 | length(inv_end_index)>2){
    nucdiv_inverted[i,] <- NA
    nucdiv_normal[i,] <- NA
    next
  } else if(length(inv_start_index)==2 & length(inv_end_index)==2){
    # if both indeces are duplicated, need to find the pair of columns that are identical
    comparison_indeces <- which(colSums(ms_binary[,inv_start_index]!=ms_binary[,inv_end_index])==0)
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
  
  
  
  
  
}







