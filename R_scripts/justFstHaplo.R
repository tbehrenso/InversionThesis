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
WINDOW_SPACING <- 25
WINDOW_SIZE <- 25   # NOTE: window size is added on each side (so the full size is more like twice this value)
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

# adjusting for IMBALANCED SAMPLE SIZE
calc_fst_between <- function(msGroup1, msGroup2){
  # here can just use abs_positions (but generally can use colnames)
  allPositions <- abs_positions
  
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
    ms_vect_1 <- msGroup1[ ,as.character(allPositions[i])]
    ms_vect_2 <- msGroup2[ ,as.character(allPositions[i])]
    # calc hexp as 2pq
    hexp_df$group1[i] <- 2 * mean(ms_vect_1) * (1 - mean(ms_vect_1))
    hexp_df$group2[i] <- 2 * mean(ms_vect_2) * (1 - mean(ms_vect_2))
    # use weighted average to calculate Hs
    #av_hexp <- mean(c(hexp_df$group1[i], hexp_df$group2[i]), na.rm = TRUE)
    av_hexp <- weighted.mean(c(hexp_df$group1[i], hexp_df$group2[i]), w=c(length(ms_vect_1), length(ms_vect_2)), na.rm = TRUE)
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
fst_windowed_all <- matrix(0, nrow=n_files, ncol=length(window_centers))

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
      fst_windowed_all[i,] <- NA
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
      fst_all <- calc_fst_between(ms_normal, ms_inverted)
      fst_windowed <- calc_sliding_window(fst_all, GENOME_LENGTH, windowSize = WINDOW_SIZE, pointSpacing = WINDOW_SPACING)
      fst_windowed_all[i,] <- fst_windowed[,2]
    }
    
  }
}

fst_windowed_p1 <- fst_windowed_all[tags_index$population=='p1',]
fst_windowed_p2 <- fst_windowed_all[tags_index$population=='p2',]

fst_p1_average <- data.frame(pos=window_centers, av_fst=colMeans(fst_windowed_p1, na.rm = T), 
                          stdev=apply(fst_windowed_all, 2, sd, na.rm=T))
fst_p2_average <- data.frame(pos=window_centers, av_fst=colMeans(fst_windowed_p2, na.rm = T), 
                             stdev=apply(fst_windowed_all, 2, sd, na.rm=T))

plot_fst_p1 <- ggplot(fst_p1_average, aes(x=pos, y=av_fst)) +
  geom_line() +
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

plot_fst_haps_pops <- grid.arrange(plot_fst_p1, plot_fst_p2, nrow=1)

ggsave('fst_haps_pops_weighted_win25.png', plot_fst_haps_pops, path=paste("Plots", args[1], args[2], sep="/"), width=12, height=5.5)



