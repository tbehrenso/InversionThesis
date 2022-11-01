
# read in files (values: selection coefficient, migration rate, replicate #)
files <- list.files(path=PATH, pattern="*.txt", full.names=F, recursive=FALSE)
n_files <- length(files)
# pre-calculate window centers' positions
window_centers <- seq(0, GENOME_LENGTH, by=WINDOW_SPACING)

tags_index <- data.frame(population=character(n_files), sel_coef=numeric(n_files), migration=numeric(n_files), 
                         repl=integer(n_files), stringsAsFactors=F)
correlations_3d <- array(numeric(), dim=c(N_TILES, N_TILES, n_files))



# correlation heatmap
corr_summ_p1 <- apply(correlations_3d[, , which(tags_index$population == "p1")], c(1, 2), mean, na.rm = TRUE)
corr_summ_p1_long <- melt(corr_summ_p1)
# correct group values to bin centers
bin_size <- GENOME_LENGTH / N_TILES
corr_summ_p1_long$Var1 <- corr_summ_p1_long$Var1*bin_size - bin_size/2
corr_summ_p1_long$Var2 <- corr_summ_p1_long$Var2*bin_size - bin_size/2

corr_summ_p2 <- apply(correlations_3d[, , which(tags_index$population == "p2")], c(1, 2), mean, na.rm = TRUE)
corr_summ_p2_long <- melt(corr_summ_p2)
# correct group values to bin centers
bin_size <- GENOME_LENGTH / N_TILES
corr_summ_p2_long$Var1 <- corr_summ_p2_long$Var1*bin_size - bin_size/2
corr_summ_p2_long$Var2 <- corr_summ_p2_long$Var2*bin_size - bin_size/2

corr_a <- ggplot(corr_summ_p1_long, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low='white', high='blue') +
  ggtitle('P1') + xlab('Position') + ylab('Position')

corr_b <- ggplot(corr_summ_p2_long, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low='white', high='blue') +
  ggtitle('P2') + xlab('Position') + ylab('Position')

plot_correlation <- grid.arrange(corr_a, corr_b, nrow=1)




if(INVERSION_PRESENT && generation > 5000){
  n_repl <- length(unique(tags_index$repl))
  
  correlations_3d_normal <- array(numeric(), dim=c(N_TILES, N_TILES, n_repl))
  correlations_3d_inverted <- array(numeric(), dim=c(N_TILES, N_TILES, n_repl))  
  
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
    ms_both[1:200,as.character(pos_p1)] <- ms_p1
    ms_both[201:400,as.character(pos_p2)] <- ms_p2
    # extract rows based on the presence of both inversion markers
    ms_normal <- ms_both[ms_both[,as.character(INV_START)]==0 & ms_both[,as.character(INV_START)]==0, ]
    ms_inverted <- ms_both[ms_both[,as.character(INV_START)]==1 & ms_both[,as.character(INV_START)]==1, ]
    
    corr_data_normal <- get_correlations(ms_normal, pos_both, numTiles = N_TILES)
    corr_data_inverted <- get_correlations(ms_inverted, pos_both, numTiles = N_TILES)
    corr_long_normal <- reduce_to_long(corr_data_normal, pos_both, numTiles = N_TILES)
    corr_long_inverted <- reduce_to_long(corr_data_inverted, pos_both, numTiles = N_TILES)
    
    correlations_3d_normal[,,repl] <- as.matrix(dcast(corr_long_normal, Var1 ~ Var2)[,-1])
    correlations_3d_inverted[,,repl] <- as.matrix(dcast(corr_long_inverted, Var1 ~ Var2)[,-1])
    
    
    
  }
}



