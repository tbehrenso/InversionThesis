# -----------------------------------------------------   (probably need to run main script first)
#    EFFECT OF WINDOW SIZE ON NUCLEOTIDE DIVERSITY
# -----------------------------------------------------

nucdiv_windowed <- calc_nuc_div(ms_binary, abs_positions, GENOME_LENGTH, seqLen = WINDOW_SIZE, centerSpacing = WINDOW_SPACING)
nucdiv_df[i,] <- nucdiv_windowed[[2]]

msdata <- ms_binary
positions <- abs_positions
totalLength <- GENOME_LENGTH
seqLen <- WINDOW_SIZE
centerSpacing <- WINDOW_SPACING




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
      
      #output[output$position==seqCenter,2] <- (sum(distances_all) / choose(num_of_seq,2)) / adjusted_seq_len
      output[output$position==seqCenter,2] <- (2*sum(distances_all) / num_of_seq^2) / adjusted_seq_len
    }
  }
  return(output)
}


m1 <- round(matrix(runif(200*5), 200, 5))
dist1 <- dist(m1, method='manhattan')
sum(dist1)/choose(200,2)/101

m2 <- round(matrix(runif(200*10), 200, 10))
dist2 <- dist(m2, method='manhattan')
sum(dist2)/choose(200,2)/200

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
}



# TAKEN FROM SCRIPT FOR QUICKER TESTING
window_centers <- seq(0, GENOME_LENGTH, by=WINDOW_SPACING)
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

nucdiv_summ_p1 <- data.frame(center = window_centers, nucdiv = colMeans(nucdiv_df[which(tags_index$population=='p1'),], na.rm=T),
                             stdev=apply(nucdiv_df[which(tags_index$population=='p1'),], 2, sd, na.rm=T))
nucdiv_summ_p2 <- data.frame(center = window_centers, nucdiv = colMeans(nucdiv_df[which(tags_index$population=='p2'),], na.rm=T),
                             stdev=apply(nucdiv_df[which(tags_index$population=='p2'),], 2, sd, na.rm=T))

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

grid.arrange(nucdiv_a, nucdiv_b, nrow=1)


##### Here comparing histograms at first two positions

pos1 <- data.frame(x=nucdiv_df[which(tags_index$population=='p1'),1])
pos2 <- data.frame(x=nucdiv_df[which(tags_index$population=='p1'),2])

pos_both <- data.frame(x=c(nucdiv_df[which(tags_index$population=='p1'),1],nucdiv_df[which(tags_index$population=='p1'),2]),
                       pos=as.factor(c(rep(1,100), rep(2, 100))))

ggplot(pos_both, aes(x=x, fill=pos)) + geom_histogram(alpha=0.5, position="identity")


#### COMPARING DIFFERENT WINDOWS SIZES
seqLen <- WINDOW_SIZE / 2
centerSpacing <- WINDOW_SPACING / 2


window_centers <- seq(0, GENOME_LENGTH, by=centerSpacing)
nucdiv_df_large <- matrix(0, nrow=n_files, ncol=length(window_centers))

for(i in 1:n_files){
  filepath <- paste0(PATH, "/", files[i])
  ms_binary <- get_ms_data(filepath)
  abs_positions <- get_positions(filepath)
  # extract metadata from filename
  tags <- strsplit(files[i], split='_')[[1]]
  tags_index[i,] <- list(tags[2], as.numeric(tags[3]), as.numeric(tags[4]), as.integer(tags[5]))
  
  # calc nucleotide diversity (over sliding window by default)
  nucdiv_windowed <- calc_nuc_div(ms_binary, abs_positions, GENOME_LENGTH, seqLen = seqLen, centerSpacing = centerSpacing)
  nucdiv_df_large[i,] <- nucdiv_windowed[[2]]
}

nucdiv_summ_p1_LARGE <- data.frame(center = window_centers, nucdiv = colMeans(nucdiv_df_large[which(tags_index$population=='p1'),], na.rm=T),
                                   stdev=apply(nucdiv_df_large[which(tags_index$population=='p1'),], 2, sd, na.rm=T))
nucdiv_summ_p2_LARGE <- data.frame(center = window_centers, nucdiv = colMeans(nucdiv_df_large[which(tags_index$population=='p2'),], na.rm=T),
                                   stdev=apply(nucdiv_df_large[which(tags_index$population=='p2'),], 2, sd, na.rm=T))

nucdiv_a_large <- ggplot(nucdiv_summ_p1_LARGE, aes(x=center, y=nucdiv)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('c) Nucleotide Diversity - P1') +
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red')} +
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue', alpha=0.4)} +
  xlab('Position') + ylab('\u03c0')
#ylim(c(0.0014, 0.0029))

nucdiv_b_large <- ggplot(nucdiv_summ_p2_LARGE, aes(x=center, y=nucdiv)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('d) Nucleotide Diversity - P2') +
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red')} +
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue', alpha=0.4)} +
  #geom_errorbar(aes(ymin=nucdiv-stdev, ymax=nucdiv+stdev), width=1, position=position_dodge(0.1)) +
  xlab('Position') + ylab('\u03c0')


grid.arrange(nucdiv_a_large, nucdiv_b_large, nrow=1)

#### Histogram, but comparing runs with different window size (and spacing)
# 1 is window size 200, 2 is window size 100
winpos1 <- data.frame(x=nucdiv_df[which(tags_index$population=='p1'),2])
winpos2 <- data.frame(x=nucdiv_df_large[which(tags_index$population=='p1'),2])

winpos_both <- data.frame(x=c(nucdiv_df[which(tags_index$population=='p1'),2],nucdiv_df_large[which(tags_index$population=='p1'),2]),
                          pos=as.factor(c(rep(1,100), rep(2, 100))))

ggplot(winpos_both, aes(x=x, fill=pos)) + geom_histogram(alpha=0.5, position="identity")


# -----------------------------------------------------
#     Testing allele frequency plot
# -----------------------------------------------------
PATH <- "Outputs/inversionLAA_2pop_s0.1_m0.001_mu1e-5_r1e-6/15000"
simtype <- strsplit(strsplit(PATH, split='/')[[1]][2], split='_')[[1]][1]

INVERSION_PRESENT <- ifelse(simtype=='adaptiveInversion' || simtype=='inversionLAA' ,TRUE, FALSE)
LAA_PRESENT <- ifelse(simtype=='locallyAdapted' || simtype=='inversionLAA' ,TRUE, FALSE)


# read in files (values: selection coefficient, migration rate, replicate #)
files <- list.files(path=PATH, pattern="*.txt", full.names=F, recursive=FALSE)
n_files <- length(files)
# pre-calculate window centers' positions
window_centers <- seq(0, GENOME_LENGTH, by=WINDOW_SPACING)

# STORAGE DATAFRAMES
tags_index <- data.frame(population=character(n_files), sel_coef=numeric(n_files), migration=numeric(n_files), 
                         repl=integer(n_files), stringsAsFactors=F)
pos_frequency <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("pop", "position", "frequency"))

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
}

ggplot(pos_frequency, aes(x=position, y=frequency)) + geom_point(size=0.05) + gglayer_markers +  
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", fx = TRUE, k = 100))

# -----------------------------------------------------
#     Comparing original and SFS nucdiv functions
# -----------------------------------------------------
positions_test <- c(1,6,11,16,21,26)

ms_test <- matrix(c(0,0,0,0,1,0,
                    0,1,0,0,1,1,
                    0,1,0,1,1,0,
                    1,0,0,1,1,1,
                    0,1,1,0,1,0,
                    0,0,0,1,0,1), nrow = 6, byrow = T)

ms_test_long <- rbind(ms_test, ms_test, ms_test, ms_test, ms_test, ms_test)


positions_test_addedzeros <- c(1,6,11,16,21,25,26)
ms_test_addedzeros <- matrix(c(0,0,0,0,1,0,0,
                               0,1,0,0,1,0,1,
                               0,1,0,1,1,0,0,
                               1,0,0,1,1,0,1,
                               0,1,1,0,1,0,0,
                               0,0,0,1,0,0,1), nrow = 6, byrow = T)

calc_nuc_div <- function(msdata, positions, totalLength, seqLen=200, centerSpacing=100){
  centers <- seq(0, totalLength, by=centerSpacing)
  # prepare storage for nucleotide diversity at each center position
  output <- data.frame(position=centers, nuc_div=NA)
  num_of_seq <- dim(msdata)[1]
  for(seqCenter in centers){
    positions_in_sequence <- which(   # select positions within window
      positions > seqCenter-seqLen & positions <= seqCenter+seqLen
    ) 
    print(positions_in_sequence)
    ms_in_seq <- as.matrix(msdata[,positions_in_sequence])
    # only do the calculations if more than one position
    if(dim(ms_in_seq)[2] > 0){
      # dist calculates distance between every combination of rows in a matrix. Manhattan method avoids "diagonal" distance
      distances_all <- dist(ms_in_seq, method='manhattan')
      # use this to adjust the sequence length when the window exceeds the range of the genome
      adjusted_seq_len <- sum((c((seqCenter-seqLen):(seqCenter+seqLen)))>=0 & (c((seqCenter-seqLen):(seqCenter+seqLen)))<=totalLength)
      # at a given position (center), nucdiv is the average number of differences divided by the length of the sequence window
      # Note: average number of differences is the total number of pairwise differences / the number of pairwise difference nChoosek
      
      # output[output$position==seqCenter,2] <- (sum(distances_all) / choose(num_of_seq,2)) / adjusted_seq_len
      # here using 2x sum of distances_all to account for both parts of the pairwise comparison matrix (above and below the diagonal)
      output[output$position==seqCenter,2] <- (2*sum(distances_all) / (num_of_seq^2)) / adjusted_seq_len
    }
  }
  return(output)
}
calc_nuc_div_sfs <- function(msdata, positions, totalLength, seqLen=200, centerSpacing=100){
  centers <- seq(0, totalLength, by=centerSpacing)
  # prepare storage for nucleotide diversity at each center position
  output <- data.frame(position=centers, nuc_div=NA)
  num_of_seq <- dim(msdata)[1]
  for(seqCenter in centers){
    positions_in_sequence <- which(   # select positions within window
      positions > seqCenter-seqLen & positions <= seqCenter+seqLen
    )  
    ms_in_seq <- as.matrix(msdata[,positions_in_sequence])
    # only do the calculations if more than one position
    if(dim(ms_in_seq)[2] > 0){
      # use this to adjust the sequence length when the window exceeds the range of the genome
      adjusted_seq_len <- sum((c((seqCenter-seqLen):(seqCenter+seqLen)))>=0 & (c((seqCenter-seqLen):(seqCenter+seqLen)))<=totalLength)
      
      sfs.raw <- table(colSums(ms_in_seq))
      
      # alternate to adjust so it always considers in respect to the less frequent allele
      # sfs.inverse <- num_of_seq - sfs.raw
      # sfs.total <- pmin(sfs.raw, sfs.inverse)
      
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
calc_nuc_div_popgenome <- function(msdata, positions, totalLength, seqLen=200, centerSpacing=100){
  centers <- seq(0, totalLength, by=centerSpacing)
  # prepare storage for nucleotide diversity at each center position
  output <- data.frame(position=centers, nuc_div=NA)
  num_of_seq <- dim(msdata)[1]
  for(seqCenter in centers){
    positions_in_sequence <- which(   # select positions within window
      positions > seqCenter-seqLen & positions <= seqCenter+seqLen
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

# testing on default version
calc_nuc_div(ms_test, positions_test, totalLength=26, seqLen=7,centerSpacing=7)
calc_nuc_div_sfs(ms_test, positions_test, totalLength=26, seqLen=7,centerSpacing=7)
calc_nuc_div_popgenome(ms_test, positions_test, totalLength=26, seqLen=7,centerSpacing=7)

# testing on longer version (more samples)
calc_nuc_div(ms_test_long, positions_test, totalLength=26, seqLen=7,centerSpacing=7)
calc_nuc_div_sfs(ms_test_long, positions_test, totalLength=26, seqLen=7,centerSpacing=7)
calc_nuc_div_popgenome(ms_test_long, positions_test, totalLength=26, seqLen=7,centerSpacing=7)

# testing with smaller window (and hence more windows)
calc_nuc_div(ms_test, positions_test, totalLength=26, seqLen=5,centerSpacing=5)
calc_nuc_div_sfs(ms_test, positions_test, totalLength=26, seqLen=5,centerSpacing=5)
calc_nuc_div_popgenome(ms_test, positions_test, totalLength=26, seqLen=5,centerSpacing=5)

# testing on version with a row of zeros in final window
calc_nuc_div(ms_test_addedzeros, positions_test_addedzeros, totalLength=26, seqLen=7,centerSpacing=7)
calc_nuc_div_sfs(ms_test_addedzeros, positions_test_addedzeros, totalLength=26, seqLen=7,centerSpacing=7)
calc_nuc_div_popgenome(ms_test_addedzeros, positions_test_addedzeros, totalLength=26, seqLen=7,centerSpacing=7)
# -----------------------------------------------------
#     Original Fst code that separates by JUST haplotype (merges populations)
# -----------------------------------------------------


###### FST -- SEPARATING HaPLOTYPES

if(INVERSION_PRESENT && generation > 5000){
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
    
    # if in a sample the inversion markers are not present, set to NA and go to next replicate
    if(!all(c(INV_START, INV_END) %in% pos_both)){
      fst_all <- NA
      next
    }
    
    ms_both <- matrix(0, nrow=2*n_indiv, ncol=length(pos_both))
    colnames(ms_both) <- pos_both
    # top half of new matrix is p1 data, bottom half is p2 data. All missing rows in a population 0 by default
    ms_both[1:n_indiv,as.character(pos_p1)] <- ms_p1
    ms_both[(n_indiv+1):(n_indiv*2),as.character(pos_p2)] <- ms_p2
    # extract rows based on the presence of both inversion markers
    ms_normal <- ms_both[ms_both[,as.character(INV_START)]==0 & ms_both[,as.character(INV_END)]==0, ]
    ms_inverted <- ms_both[ms_both[,as.character(INV_START)]==1 & ms_both[,as.character(INV_END)]==1, ]
    
    # convert to matrix of one row if the msdata has only one sample (and hence was converted to a vector)
    if(is.null(dim(ms_normal))){
      ms_normal <- t(as.matrix(ms_normal))
    }
    if(is.null(dim(ms_inverted))){
      ms_inverted <- t(as.matrix(ms_inverted))
    }
    
    if(dim(ms_normal)[1]==0 | dim(ms_inverted)[1]==0){
      fst_all <- NA
    } else {
      fst_all <- calc_fst_between(ms_normal, ms_inverted)
    }
    
    fst_windowed <- calc_sliding_window(fst_all, GENOME_LENGTH, windowSize = WINDOW_SIZE, pointSpacing = WINDOW_SPACING)
    fst_windowed_haplotypes_all[repl,] <- fst_windowed[,2]
  }
  
  fst_haplotypes_average <- data.frame(pos=window_centers, av_fst=colMeans(fst_windowed_haplotypes_all, na.rm = T), 
                                       stdev=apply(fst_windowed_haplotypes_all, 2, sd, na.rm=T))
  
  plot_fst_haplotypes <- ggplot(fst_haplotypes_average, aes(x=pos, y=av_fst)) +
    geom_line() +
    scale_fill_gradient(low='white', high='blue') +
    ggtitle('F_ST between Haplotypes') + 
    xlab('Position') + ylab(expression(F[ST])) +
    gglayer_markers
  
  
  if(on_cluster){
    ggsave('fst_haps.png', plot_fst_haplotypes, path=paste("Plots", args[1], args[2], sep="/"), width=8, height=6)
  }else{
    print(plot_fst_haplotypes)
  }
}

# -----------------------------------------------------
#     Overlap some Nucdiv data (original vs added deleterious mutations)
# -----------------------------------------------------

GENOME_LENGTH <- 120000
FIXED_MUTATION_POS1 <- 30000
FIXED_MUTATION_POS2 <- 70000
INV_START <- 10000
INV_END <- 110000  # this value should NOT be the '-1' value that the SLiM script uses. This script does that correction later
WINDOW_SPACING <- 100
WINDOW_SIZE <- 100   # NOTE: window size is added on each side (so the full size is more like twice this value)
N_TILES <- 600   # number of tiles along each axis of the correlation heatmap
FIRST_GEN <- 5000  # first generation where inversion/locally adapted alleles are introduced

# record presence or absence of inversion and locally adapted alleles
INVERSION_PRESENT <- TRUE
LAA_PRESENT <- TRUE

# reusable layer for ggplot to include marker lines for inversion bounds (blue) and locally adapted alleles (red)
gglayer_markers <- list(
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red', alpha=0.3)},
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='dotted', colour='blue', alpha=0.3)}
)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

setwd("C:/Users/tbehr/Desktop")

nucdiv_original <- loadRData("nucdiv_summ_p1_original.Rds")
nucdiv_deleterious <- loadRData("nucdiv_summ_p1_deleterious.Rds")

min_nucdiv <- 0
max_nucdiv <- max(c(nucdiv_original$nucdiv, nucdiv_deleterious$nucdiv))

fst_average_all <- data.frame(pos=nucdiv_original$center, original=nucdiv_original$nucdiv, deleterious=nucdiv_deleterious$nucdiv)

fst_average_melt <- melt(fst_average_all, id='pos')

ggplot(fst_average_melt, aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('Nucleotide Diversity') + 
  xlab('Position') + ylab('\u03c0') +
  gglayer_markers





