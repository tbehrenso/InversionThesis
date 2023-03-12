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
N_TILES <- 200   # number of tiles along each axis of the correlation heatmap
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

get_sfs <- function(msdata, positions, totalLength, seqLen=100, centerSpacing=100){
  centers <- seq(0, totalLength, by=centerSpacing)
  # prepare storage for nucleotide diversity at each center position
  output <- data.frame(position=centers, nuc_div=NA)
  num_of_seq <- dim(msdata)[1]

  sfs.raw <- table(colSums(msdata))

  return(OUTHERE)
}

#-----------------------------------------------------------
# DATA EXTRACTION
#-----------------------------------------------------------

# read in files (values: selection coefficient, migration rate, replicate #)
files <- list.files(path=PATH, pattern="*.txt", full.names=F, recursive=FALSE)
n_files <- length(files)

# STORAGE DATAFRAMES
tags_index <- data.frame(population=character(n_files), sel_coef=numeric(n_files), migration=numeric(n_files), 
                         repl=integer(n_files), stringsAsFactors=F)
# correlations_3d <- array(numeric(), dim=c(N_TILES, N_TILES, n_files))  # NOTE: end up switching a lot between long and wide here
#      maybe change so its all in long?

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
  # extract metadata from filename
  tags <- strsplit(files[i], split='_')[[1]]
  tags_index[i,] <- list(tags[2], as.numeric(tags[3]), as.numeric(tags[4]), as.integer(tags[5]))
}

n_indiv <- dim(ms_binary)[1]
n_repl <- length(unique(tags_index$repl))
sfs_2D <- matrix(0, nrow=n_indiv, ncol=n_indiv)

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
  ms_both[1:n_indiv,as.character(pos_p1)] <- ms_p1
  ms_both[(n_indiv+1):(n_indiv*2),as.character(pos_p2)] <- ms_p2
  # extract rows again (so ms_p1 and ms_p2 now have all positions, with missing positions filled with zeros)
  ms_p1 <- ms_both[1:n_indiv,]
  ms_p2 <- ms_both[(n_indiv+1):(n_indiv*2),]
  
  colsums_p1 <- colSums(ms_p1)
  colsums_p2 <- colSums(ms_p2)
  
  for(i in 1:length(pos_both)){
    freq1 <- colsums_p1[[i]]
    freq2 <- colsums_p2[[i]]
    
    sfs_2D[freq1, freq2] <- sfs_2D[freq1, freq2] + 1
  }
}

max_freq <- dim(sfs_2D)

sfs_2D[max_freq,max_freq] <- 0
sfs_2D[max_freq,max_freq] <- max(sfs_2D)

sfs_2D_long <- melt(sfs_2D)
# taking LOG scale (+1 to avoid NA values)
sfs_2D_long$value <- log10(sfs_2D_long$value + 1)

sfs_2D_plot <- ggplot(sfs_2D_long, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low='white', high='blue', limits=c(0,4.3)) +
  #ggtitle('2D SFS') + 
  xlab('Allele Frequency - P1') + ylab('Allele Frequency - P2') +
  scale_x_continuous(n.breaks=10, expand=c(0,0)) +
  scale_y_continuous(n.breaks=10, expand=c(0,0)) +
  theme(axis.line = element_line(color="black", size = 0.5)) +
  coord_fixed()

if(on_cluster){
  ggsave('sfs_2D.png', sfs_2D_plot, path=paste("Plots", args[1], args[2], sep="/"), width=20, height=18)
  save(sfs_2D_long, file=paste("data_summary", args[1], args[2],"sfs_2D_long_LOG.Rds", sep="/"))
}else{
  print(sfs_2D_plot)
}


# -------------------------------------------------------------
# PLOTTING TWO NEXT TO EACH OTHER



# (LOAD DATA FIRST)
sfs_2D_long_invLAA <- loadRData('C:/Users/tbehr/Desktop/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/sfs_2D_long_LOG.Rds')
sfs_2D_long_adapInv <- loadRData('C:/Users/tbehr/Desktop/data_summary/adaptiveInversion_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/sfs_2D_long_LOG.Rds')

sfs_2D_long_invLAA$Var1 <- sfs_2D_long_invLAA$Var1/200
sfs_2D_long_invLAA$Var2 <- sfs_2D_long_invLAA$Var2/200
sfs_2D_long_adapInv$Var1 <- sfs_2D_long_adapInv$Var1/200
sfs_2D_long_adapInv$Var2 <- sfs_2D_long_adapInv$Var2/200


max_value <- max(c(sfs_2D_long_invLAA$value, sfs_2D_long_adapInv$value))

sfs_2D_plot_invLAA <- ggplot(sfs_2D_long_invLAA, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  ggtitle('a) invLAA') +
  scale_fill_gradient(low='white', high='blue', limits=c(0,max_value), name='Allele \nFrequency\n(log10)') +
  #ggtitle('2D SFS') + 
  xlab('Allele Frequency - P1') + ylab('Allele Frequency - P2') +
  scale_x_continuous(n.breaks=10, expand=c(0,0)) +
  scale_y_continuous(n.breaks=10, expand=c(0,0)) +
  theme(axis.line = element_line(color="black", size = 0.5)) +
  coord_fixed()

sfs_2D_plot_adapInv <- ggplot(sfs_2D_long_adapInv, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  ggtitle('b) adapInv') +
  scale_fill_gradient(low='white', high='blue', limits=c(0,max_value), name='Allele \nFrequency\n(log10)') +
  #ggtitle('2D SFS') + 
  xlab('Allele Frequency - P1') + ylab('Allele Frequency - P2') +
  scale_x_continuous(n.breaks=10, expand=c(0,0)) +
  scale_y_continuous(n.breaks=10, expand=c(0,0)) +
  theme(axis.line = element_line(color="black", size = 0.5)) +
  coord_fixed()

sfs_2D_plot_both <- grid.arrange(sfs_2D_plot_invLAA, sfs_2D_plot_adapInv, nrow=1)


