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


if(on_cluster){
  ggsave('nucdiv_hexp_win400.png', plot_nucdiv_hexp, device="png", path=paste("Plots", args[1], args[2], sep="/"), width=16, height=12)
  save(hexp_summ_p1, file=paste("data_summary", args[1], args[2],"hexp_summ_p1_400.Rds", sep="/"))
  save(hexp_summ_p2, file=paste("data_summary", args[1], args[2],"hexp_summ_p2_400.Rds", sep="/"))
  save(nucdiv_summ_p1, file=paste("data_summary", args[1], args[2],"nucdiv_summ_p1_400.Rds", sep="/"))
  save(nucdiv_summ_p2, file=paste("data_summary", args[1], args[2],"nucdiv_summ_p2_400.Rds", sep="/"))
}






