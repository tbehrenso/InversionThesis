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

calc_sliding_window <- function(posValData, totalLength, windowSize, pointSpacing){
  centers <- seq(0, totalLength, by=pointSpacing)
  output <- data.frame(position=centers, average=NA)
  for(i in 1:length(centers)){
    correspondingValues <- posValData[2][posValData[1] > centers[i]-windowSize & posValData[1] <= centers[i]+windowSize]
    output[i, 2] <- mean(correspondingValues)
  }
  return(output)
}

calc_fst_between <- function(msGroup1, msGroup2){
  # extract positions from column names
  allPositions <- unique(sort(c(as.integer(colnames(msGroup1)), as.integer(colnames(msGroup2)))))
  
  # storage
  hexp_df <- data.frame(pos=allPositions, group1=numeric(length(allPositions)), group2=numeric(length(allPositions)), 
                        total=numeric(length(allPositions)))
  fst_all <- data.frame(pos=allPositions, fst=numeric(length(allPositions)))
  
  for(i in 1:length(allPositions)){
    # convert to matrix of one row if the msdata has only one sample (and hence was converted to a vector)
    if(is.null(dim(msGroup1))){
      msGroup1 <- t(as.matrix(msGroup1))
    }
    if(is.null(dim(msGroup2))){
      msGroup2 <- t(as.matrix(msGroup2))
    }
    # extract columns at current positions
    ms_vect_1 <- msGroup1[ ,as.character(allPositions[i])]
    ms_vect_2 <- msGroup2[ ,as.character(allPositions[i])]
    # calc hexp as 2pq
    hexp_df$group1[i] <- 2 * mean(ms_vect_1) * (1 - mean(ms_vect_1))
    hexp_df$group2[i] <- 2 * mean(ms_vect_2) * (1 - mean(ms_vect_2))
    av_hexp <- mean(hexp_df$group1[i], hexp_df$group2[i], na.rm = TRUE)
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

calc_kst_between <- function(msGroup1, msGroup2, positions){
  # extract positions from column names
  allPositions <- unique(sort(c(as.integer(colnames(msGroup1)), as.integer(colnames(msGroup2)))))
  
  # storage
  num_of_windows <- length(window_centers)
  nucdiv_df <- data.frame(pos=window_centers, group1=numeric(num_of_windows), group2=numeric(num_of_windows), 
                        total=numeric(num_of_windows))
  kst_all <- data.frame(pos=window_centers, fst=numeric(num_of_windows))
  
  # calc nucleotide diversity (over sliding window by default)
  nucdiv_windowed_1 <- calc_nuc_div(msGroup1, positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
  nucdiv_windowed_2 <- calc_nuc_div(msGroup2, positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
  nucdiv_df$group1 <- nucdiv_windowed_1[[2]]
  nucdiv_df$group2 <- nucdiv_windowed_2[[2]]
  
  ms_both <- rbind(msGroup1, msGroup2)
  nucdiv_windowed_both <- calc_nuc_div(msGroup1, positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
  nucdiv_df$total <- nucdiv_windowed_both[[2]]
  
  av_nucdiv <- rowMeans(cbind(nucdiv_df$group1, nucdiv_df$group2))
  kst_all$fst <- (nucdiv_df$total - av_nucdiv) / nucdiv_df$total

  return(kst_all)
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
kst_windowed_all <- matrix(0, nrow=n_files, ncol=length(window_centers))
nucdiv_df <- matrix(0, nrow=n_files, ncol=length(window_centers))

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
  # extract metadata from filename
  tags <- strsplit(files[i], split='_')[[1]]
  tags_index[i,] <- list(tags[2], as.numeric(tags[3]), as.numeric(tags[4]), as.integer(tags[5]))
  
  if(INVERSION_PRESENT && generation > 5000){
    inv_start_index <- which(abs_positions==INV_START)
    inv_end_index <- which(abs_positions==INV_END-1)
    
    # Fix for multiple mutations at a site
    # this first statement is if there are MORE THAN 2 mutations at a breakpoint, which is really weird. Dunno why thats happening
    if(length(inv_start_index)>2 | length(inv_end_index)>2){
      kst_windowed_all[i,] <- NA
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
    
    ms_normal <- ms_binary[ms_binary[,inv_start_index]==0 & ms_binary[,inv_end_index]==0, ]
    ms_inverted <- ms_binary[ms_binary[,inv_start_index]==1 & ms_binary[,inv_end_index]==1, ]
    
    # convert to matrix of one row if the msdata has only one sample (and hence was converted to a vector)
    if(is.null(dim(ms_normal))){
      ms_normal <- t(as.matrix(ms_normal))
    }
    if(is.null(dim(ms_inverted))){
      ms_inverted <- t(as.matrix(ms_inverted))
    }
    
    colnames(ms_normal) <- abs_positions
    colnames(ms_inverted) <- abs_positions
    
    normal_is_valid <- TRUE
    inverted_is_valid <- TRUE
    
    if(dim(ms_normal)[1]<=1){
      normal_is_valid <- FALSE
    }
    if(dim(ms_inverted)[1]<=1){
      inverted_is_valid <- FALSE
    }
    if(normal_is_valid & inverted_is_valid){
      kst_windowed <- calc_kst_between(ms_normal, ms_inverted, abs_positions)
      kst_windowed_all[i,] <- kst_windowed[[2]]
    }
  }
}

kst_windowed_all[is.infinite(kst_windowed_all)] <- NA

kst_windowed_p1 <- kst_windowed_all[tags_index$population=='p1',]
kst_windowed_p2 <- kst_windowed_all[tags_index$population=='p2',]

kst_p1_average <- data.frame(pos=window_centers, av_kst=colMeans(kst_windowed_p1, na.rm = T), 
                             stdev=apply(kst_windowed_all, 2, sd, na.rm=T))
kst_p2_average <- data.frame(pos=window_centers, av_kst=colMeans(kst_windowed_p2, na.rm = T), 
                             stdev=apply(kst_windowed_all, 2, sd, na.rm=T))

plot_kst_p1 <- ggplot(kst_p1_average, aes(x=pos, y=av_kst)) +
  geom_line() +
  scale_fill_gradient(low='white', high='blue') +
  ggtitle('P1 - Between Haplotypes') + 
  xlab('Position') + ylab(expression(F[ST])) +
  gglayer_markers

plot_kst_p2 <- ggplot(kst_p2_average, aes(x=pos, y=av_kst)) +
  geom_line() +
  scale_fill_gradient(low='white', high='blue') +
  ggtitle('P2 - Between Haplotypes') + 
  xlab('Position') + ylab(expression(F[ST])) +
  gglayer_markers

plot_kst_haps_pops <- grid.arrange(plot_kst_p1, plot_kst_p2, nrow=1)

ggsave('kst_haps_pops.png', plot_kst_haps_pops, path=paste("Plots", args[1], args[2], sep="/"), width=12, height=5.5)




