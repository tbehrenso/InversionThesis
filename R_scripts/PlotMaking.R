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
library(patchwork)

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

gglayer_markers_invonly <- list(
  {if(F)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red', alpha=0.3)},
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue', alpha=0.4)}
)

gglayer_markers_short <- list(
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1/10, FIXED_MUTATION_POS2/10), linetype='dashed', colour='red', alpha=0.3)},
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START/10, INV_END/10), linetype='solid', colour='blue', alpha=0.4)}
)

gglayer_markers_short_invonly <- list(
  {if(F)geom_vline(xintercept = c(FIXED_MUTATION_POS1/10, FIXED_MUTATION_POS2/10), linetype='dashed', colour='red', alpha=0.3)},
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START/10, INV_END/10), linetype='solid', colour='blue', alpha=0.4)}
)

# function for easily loading and renaming a Rds file. Taken from stackoverflow
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# ------------------------------------------------------------------------------------------
# Hudson FST

hudson_p1 <- loadRData("C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_p1_average_FIXED.Rds")
hudson_p2 <- loadRData("C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_p2_average_FIXED.Rds")

max_y <- max(c(hudson_p1$av_fst, hudson_p2$av_fst))

plot_fst_p1 <- ggplot(hudson_p1, aes(x=pos, y=av_fst)) +
  geom_line()+
  ggtitle('a) P1') +
  xlab('Position') + ylab(bquote(F["ST"]*" between Arrangements")) +
  ylim(0, max_y+max_y*0.05) +
  gglayer_markers + theme_bw() + coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))

plot_fst_p2 <- ggplot(hudson_p2, aes(x=pos, y=av_fst)) +
  geom_line() +
  ggtitle('b) P2') +
  xlab('Position') + ylab('') +
  ylim(0, max_y+max_y*0.05) +
  gglayer_markers + theme_bw() + coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))

plot_fst_hudson <- grid.arrange(plot_fst_p1, plot_fst_p2, nrow=1)

plot_fst_p1 + plot_fst_p2 + plot_layout(ncol=2)



# ------------------------------------------------------------------------------------------
# NucDiv and Hexp


invLAA_nucdiv_p1 <- loadRData("C:/Users/tbehr/Desktop/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/nucdiv_summ_p1.Rds")
invLAA_nucdiv_p2 <- loadRData("C:/Users/tbehr/Desktop/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/nucdiv_summ_p2.Rds")
adapInv_nucdiv_p1 <- loadRData("C:/Users/tbehr/Desktop/data_summary/adaptiveInversion_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/nucdiv_summ_p1.Rds")
adapInv_nucdiv_p2 <- loadRData("C:/Users/tbehr/Desktop/data_summary/adaptiveInversion_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/nucdiv_summ_p2.Rds")

invLAA_hexp_p1 <- loadRData("C:/Users/tbehr/Desktop/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/hexp_summ_p1.Rds")
invLAA_hexp_p2 <- loadRData("C:/Users/tbehr/Desktop/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/hexp_summ_p2.Rds")
adapInv_hexp_p1 <- loadRData("C:/Users/tbehr/Desktop/data_summary/adaptiveInversion_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/hexp_summ_p1.Rds")
adapInv_hexp_p2 <- loadRData("C:/Users/tbehr/Desktop/data_summary/adaptiveInversion_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/hexp_summ_p2.Rds")

min_hexp_invLAA <- min(invLAA_hexp_p1$hexp, invLAA_hexp_p2$hexp)
max_hexp_invLAA <- max(invLAA_hexp_p1$hexp, invLAA_hexp_p2$hexp)
min_nucdiv_invLAA <- min(invLAA_nucdiv_p1$nucdiv, invLAA_nucdiv_p2$nucdiv)
max_nucdiv_invLAA <- max(invLAA_nucdiv_p1$nucdiv, invLAA_nucdiv_p2$nucdiv)

min_hexp_adapInv <- min(adapInv_hexp_p1$hexp, adapInv_hexp_p2$hexp)
max_hexp_adapInv <- max(adapInv_hexp_p1$hexp, adapInv_hexp_p2$hexp)
min_nucdiv_adapInv <- min(adapInv_nucdiv_p1$nucdiv, adapInv_nucdiv_p2$nucdiv)
max_nucdiv_adapInv <- max(adapInv_nucdiv_p1$nucdiv, adapInv_nucdiv_p2$nucdiv)

invLAA_hexp_a <- ggplot(invLAA_hexp_p1, aes(x=center, y=hexp)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('a) Expected Heterozygosity - P1') +
  gglayer_markers +
  xlab('Position') + ylab(expression(H[exp])) +
  ylim(c(min_hexp_invLAA, max_hexp_invLAA)) +
  coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))

invLAA_hexp_b <- ggplot(invLAA_hexp_p2, aes(x=center, y=hexp)) +
  geom_line() +
  scale_color_brewer(palette="Dark2") +
  ggtitle('b) Expected Heterozygosity - P2') +
  gglayer_markers +
  xlab('Position') + ylab(expression(H[exp])) +
  ylim(c(min_hexp_invLAA, max_hexp_invLAA)) +
  coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))

invLAA_nucdiv_a <- ggplot(invLAA_nucdiv_p1, aes(x=center, y=nucdiv)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('c) Nucleotide Diversity - P1') +
  gglayer_markers +
  xlab('Position') + ylab('\u03c0') +
  ylim(c(min_nucdiv_invLAA, max_nucdiv_invLAA)) +
  coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))

invLAA_nucdiv_b <- ggplot(invLAA_nucdiv_p2, aes(x=center, y=nucdiv)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('d) Nucleotide Diversity - P2') +
  gglayer_markers +
  xlab('Position') + ylab('\u03c0') +
  ylim(c(min_nucdiv_invLAA, max_nucdiv_invLAA)) +
  coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))

adapInv_hexp_a <- ggplot(adapInv_hexp_p1, aes(x=center, y=hexp)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('a) Expected Heterozygosity - P1') +
  gglayer_markers_invonly +
  xlab('Position') + ylab(expression(H[exp])) +
  ylim(c(min_hexp_adapInv, max_hexp_adapInv)) +
  coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))

adapInv_hexp_b <- ggplot(adapInv_hexp_p2, aes(x=center, y=hexp)) +
  geom_line() +
  scale_color_brewer(palette="Dark2") +
  ggtitle('b) Expected Heterozygosity - P2') +
  gglayer_markers_invonly +
  xlab('Position') + ylab(expression(H[exp])) +
  ylim(c(min_hexp_adapInv, max_hexp_adapInv)) +
  coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))

adapInv_nucdiv_a <- ggplot(adapInv_nucdiv_p1, aes(x=center, y=nucdiv)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('c) Nucleotide Diversity - P1') +
  gglayer_markers_invonly +
  xlab('Position') + ylab('\u03c0') +
  ylim(c(min_nucdiv_adapInv, max_nucdiv_adapInv)) +
  coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))

adapInv_nucdiv_b <- ggplot(adapInv_nucdiv_p2, aes(x=center, y=nucdiv)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('d) Nucleotide Diversity - P2') +
  gglayer_markers_invonly +
  xlab('Position') + ylab('\u03c0') +
  ylim(c(min_nucdiv_adapInv, max_nucdiv_adapInv)) +
  coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))


plot_nucdiv_hexp_invLAA <- grid.arrange(invLAA_hexp_a, invLAA_hexp_b, invLAA_nucdiv_a, invLAA_nucdiv_b, nrow=2)

plot_nucdiv_hexp_adapInv <- grid.arrange(adapInv_hexp_a, adapInv_hexp_b, adapInv_nucdiv_a, adapInv_nucdiv_b, nrow=2)


# ------------------------------------------------------------------------------------------
# nucdiv by haplotype

nucdiv_haps_invLAA <- loadRData('C:/Users/tbehr/Desktop/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/nucdiv_all_long.Rds')
nucdiv_haps_adapInv <- loadRData('C:/Users/tbehr/Desktop/data_summary/adaptiveInversion_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/nucdiv_all_long.Rds')

nucdiv_haps_invLAA <- loadRData('C:/Users/tbehr/Desktop/data_summary/inversionLAA_2pop_DELETERIOUS_s0.1_m0.01_mu1e-5_r1e-6/15000/nucdiv_all_long.Rds')
nucdiv_haps_adapInv <- loadRData('C:/Users/tbehr/Desktop/data_summary/adaptiveInversion_2pop_DELETERIOUS_s0.1_m0.01_mu1e-5_r1e-6/15000/nucdiv_all_long.Rds')

ylow <- min(c(nucdiv_haps_adapInv$value, nucdiv_haps_invLAA$value))
yhigh <- max(c(nucdiv_haps_adapInv$value, nucdiv_haps_invLAA$value))

plot_nucdiv_haplotypes_invLAA <- ggplot(nucdiv_haps_invLAA, aes(x=center, y=value, col=variable)) +
  geom_line() +
  #ggtitle('a) invLAA') +
  gglayer_markers +
  xlab('Position') + ylab('\u03c0') +
  ylim(c(ylow, yhigh)) +
  coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(guide = guide_axis(angle = 45)) +
  labs(colour = 'Arrangement') +
  scale_color_hue(labels = c("Inv. P1", "Inv. P2","Col. P1","Col. P2"))


plot_nucdiv_haplotypes_adapInv <- ggplot(nucdiv_haps_adapInv, aes(x=center, y=value, col=variable)) +
  geom_line() +
  ggtitle('b) adapInv') +
  gglayer_markers_invonly +
  xlab('Position') + ylab('\u03c0') +
  ylim(c(ylow, yhigh)) +
  coord_cartesian(expand=F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(guide = guide_axis(angle = 45)) +
  labs(colour = 'Arrangement') +
  scale_color_hue(labels = c("Inv. P1", "Inv. P2","Col. P1","Col. P2"))

plot_nucdiv_haplotypes_both <- grid.arrange(plot_nucdiv_haplotypes_invLAA, plot_nucdiv_haplotypes_adapInv, nrow=1)


# ------------------------------------------------------------------------------------------
# correlation heatmap 

corr_summ_p1 <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/corr_summ_p1.Rds')
corr_summ_p2 <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/corr_summ_p2.Rds')

corr_summ_p1_long <- melt(corr_summ_p1)
bin_size <- GENOME_LENGTH / N_TILES
corr_summ_p1_long$Var1 <- corr_summ_p1_long$Var1*bin_size - bin_size/2
corr_summ_p1_long$Var2 <- corr_summ_p1_long$Var2*bin_size - bin_size/2

corr_summ_p2_long <- melt(corr_summ_p2)
bin_size <- GENOME_LENGTH / N_TILES
corr_summ_p2_long$Var1 <- corr_summ_p2_long$Var1*bin_size - bin_size/2
corr_summ_p2_long$Var2 <- corr_summ_p2_long$Var2*bin_size - bin_size/2


###### Extracting just the inversion to plot in the heatmap
corr_summ_p1_long <- subset(corr_summ_p1_long, Var1 >= INV_START & Var1 <= INV_END)
corr_summ_p1_long <- subset(corr_summ_p1_long, Var2 >= INV_START & Var2 <= INV_END)

corr_summ_p2_long <- subset(corr_summ_p2_long, Var1 >= INV_START & Var1 <= INV_END)
corr_summ_p2_long <- subset(corr_summ_p2_long, Var2 >= INV_START & Var2 <= INV_END)

corr_a <- ggplot(corr_summ_p1_long, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low='white', high='blue', name='Correlation') +
  ggtitle('a)') + xlab('Position') + ylab('Position') +
  coord_fixed() +
  scale_x_continuous(n.breaks=5, expand=c(0,0), guide = guide_axis(angle = 45)) +
  scale_y_continuous(n.breaks=5, expand=c(0,0)) +
  theme(axis.line = element_line(color="black", size = 0.5))

corr_b <- ggplot(corr_summ_p2_long, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low='white', high='blue', name='Correlation') +
  ggtitle('P2') + xlab('Position') + ylab('Position') +
  coord_fixed() +
  scale_x_continuous(n.breaks=5, expand=c(0,0), guide = guide_axis(angle = 45)) +
  scale_y_continuous(n.breaks=5, expand=c(0,0)) +
  theme(axis.line = element_line(color="black", size = 0.5)) +
  labs(colour = 'Correlation')

plot_correlation <- grid.arrange(corr_a, corr_b, nrow=1)

# ------------------------------------------------------------------------------------------
# correlation with breakpoints

breakpoint_corr_invLAA <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/breakpoints_corr_mean_filt.Rds')

ymax <- max(breakpoint_corr_invLAA$corr_mean)

plot_corr_breakpoints_filt <- ggplot(breakpoint_corr_invLAA, aes(x=pos, y=corr_mean)) + 
  geom_line() + gglayer_markers +
  coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(guide = guide_axis(angle = 45)) +
  ylim(c(0, ymax+ymax*0.05)) +
  xlab('Position') + ylab('Correlation with\nInversion Breakpoint') +
  ggtitle('b)')


corr_a + plot_corr_breakpoints_filt + plot_layout(nrow=2, heights=c(3,1))

# ------------------------------------------------------------------------------------------
# SHORT - Fst


fst_short_invLAA <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_SHORTINV_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds')
fst_orig_invLAA <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds')

fst_short_adapInv <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/adaptiveInversion_2pop_SHORTINV_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds')
fst_orig_adapInv <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/adaptiveInversion_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds')

fst_orig_invLAA_reduced = fst_orig_invLAA[seq(1, nrow(fst_orig_invLAA), 10), ]
fst_invLAA_both <- data.frame(pos=fst_short_invLAA$pos, fst_orig=fst_orig_invLAA_reduced$av_fst, fst_short=fst_short_invLAA$av_fst)

fst_average_melt_invLAA <- melt(fst_invLAA_both, id='pos')

ymin_invLAA <- min(fst_average_melt_invLAA$value)
ymax_invLAA <- max(fst_average_melt_invLAA$value)


fst_invLAA_both <- ggplot(fst_average_melt_invLAA, aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('a)') +
  xlab('Position') + ylab(bquote(F["ST"]*" between P1 and P2")) +
  gglayer_markers_short +
  scale_color_manual(values=c('black', '#FF4900'), labels = c("Original\n(Scaled)", "Short")) +
  labs(colour = 'Simulation\nType') +
  theme_bw() + coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45)) +
  ylim(c(ymin_invLAA-ymin_invLAA*0.05, ymax_invLAA+ymax_invLAA*0.05))


fst_orig_adapInv_reduced = fst_orig_adapInv[seq(1, nrow(fst_orig_adapInv), 10), ]
fst_adapInv_both <- data.frame(pos=fst_short_adapInv$pos, fst_orig=fst_orig_adapInv_reduced$av_fst, fst_short=fst_short_adapInv$av_fst)

fst_average_melt_adapInv <- melt(fst_adapInv_both, id='pos')

ymin_adapInv <- min(fst_average_melt_adapInv$value)
ymax_adapInv <- max(fst_average_melt_adapInv$value)


fst_adapInv_both <- ggplot(fst_average_melt_adapInv, aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('b)') +
  xlab('Position') + ylab(bquote(F["ST"]*" between P1 and P2")) +
  gglayer_markers_short_invonly +
  scale_color_manual(values=c('black', '#0AD5B6'), labels = c("Original\n(Scaled)", "Short")) +
  labs(colour = 'Simulation\nType') +
  theme_bw() + coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45)) +
  ylim(c(ymin_adapInv-ymin_adapInv*0.05, ymax_adapInv+ymax_adapInv*0.05))


grid.arrange(fst_invLAA_both, fst_adapInv_both, nrow=1)


# ------------------------------------------------------------------------------------------
# DELETERIOUS - Fst

gglayer_markers <- list(
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red', alpha=0.3)},
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue', alpha=0.4)}
)

gglayer_markers_invonly <- list(
  {if(F)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red', alpha=0.3)},
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='solid', colour='blue', alpha=0.4)}
)

fst_short_invLAA <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_DELETERIOUS_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds')
fst_orig_invLAA <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds')

fst_short_adapInv <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/adaptiveInversion_2pop_DELETERIOUS_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds')
fst_orig_adapInv <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/adaptiveInversion_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds')

fst_invLAA_both <- data.frame(pos=fst_short_invLAA$pos, fst_orig=fst_orig_invLAA$av_fst, fst_short=fst_short_invLAA$av_fst)

fst_average_melt_invLAA <- melt(fst_invLAA_both, id='pos')

ymin_invLAA <- min(fst_average_melt_invLAA$value)
ymax_invLAA <- max(fst_average_melt_invLAA$value)


fst_invLAA_both <- ggplot(fst_average_melt_invLAA, aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('a)') +
  xlab('Position') + ylab(bquote(F["ST"]*" between P1 and P2")) +
  gglayer_markers +
  scale_color_manual(values=c('black', '#FF4900'), labels = c("Original\n(Scaled)", "Deleterious")) +
  labs(colour = 'Simulation\nType') +
  theme_bw() + coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45)) +
  ylim(c(ymin_invLAA-ymin_invLAA*0.05, ymax_invLAA+ymax_invLAA*0.05))

fst_adapInv_both <- data.frame(pos=fst_short_adapInv$pos, fst_orig=fst_orig_adapInv$av_fst, fst_short=fst_short_adapInv$av_fst)

fst_average_melt_adapInv <- melt(fst_adapInv_both, id='pos')

ymin_adapInv <- min(fst_average_melt_adapInv$value)
ymax_adapInv <- max(fst_average_melt_adapInv$value)


fst_adapInv_both <- ggplot(fst_average_melt_adapInv, aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('b)') +
  xlab('Position') + ylab(bquote(F["ST"]*" between P1 and P2")) +
  gglayer_markers_invonly +
  scale_color_manual(values=c('black', '#0AD5B6'), labels = c("Original\n(Scaled)", "Deleterious")) +
  labs(colour = 'Simulation\nType') +
  theme_bw() + coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45)) +
  ylim(c(ymin_adapInv-ymin_adapInv*0.05, ymax_adapInv+ymax_adapInv*0.05))


grid.arrange(fst_invLAA_both, fst_adapInv_both, nrow=1)


# ------------------------------------------------------------------------------------------
# DELETERIOUS - Nucdiv

nucdiv_orig_invLAA <- loadRData('C:/Users/tbehr/Desktop/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/nucdiv_summ_p1.Rds')
nucdiv_deleterious_invLAA <- loadRData('C:/Users/tbehr/Desktop/data_summary/inversionLAA_2pop_DELETERIOUS_s0.1_m0.01_mu1e-5_r1e-6/15000/nucdiv_summ_p1.Rds')

nucdiv_invLAA_both <- data.frame(pos=nucdiv_deleterious_invLAA$center, nucdiv_orig=nucdiv_orig_invLAA$nucdiv, nucdiv_short=nucdiv_deleterious_invLAA$nucdiv)

nucdiv_average_melt_invLAA <- melt(nucdiv_invLAA_both, id='pos')
#Remove one, can comment out
nucdiv_melt_reduced <- melt(subset(nucdiv_invLAA_both, select=-c(nucdiv_short)), id='pos')


ymin_invLAA <- min(nucdiv_average_melt_invLAA$value)
ymax_invLAA <- max(nucdiv_average_melt_invLAA$value)

fst_invLAA_both <- ggplot(nucdiv_average_melt_invLAA, aes(x=pos, y=value, col=variable)) +
  geom_line() +
  xlab('Position') + ylab('\u03c0') +
  gglayer_markers +
  scale_color_manual(values=c('black', '#FF4900'), labels = c("Original\n(Scaled)", "Deleterious")) +
  labs(colour = 'Simulation\nType') +
  theme_bw() + coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45)) +
  ylim(c(ymin_invLAA-ymin_invLAA*0.05, ymax_invLAA+ymax_invLAA*0.05))

fst_invLAA_justOne <- ggplot(nucdiv_melt_reduced, aes(x=pos, y=value, col=variable)) +
  geom_line() +
  xlab('Position') + ylab('\u03c0') +
  gglayer_markers +
  scale_color_manual(values=c('black', '#FF4900'), labels = c("Original\n(Scaled)", "Deleterious")) +
  labs(colour = 'Simulation\nType') +
  theme_bw() + coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45)) +
  ylim(c(ymin_invLAA-ymin_invLAA*0.05, ymax_invLAA+ymax_invLAA*0.05))

# ------------------------------------------------------------------------------------------
# FST of last two generations of SHORT

fst_shortinv_14000 <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_SHORTINV_s0.1_m0.01_mu1e-5_r1e-6/14000/fst_average.Rds')
fst_shortinv_13000 <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_SHORTINV_s0.1_m0.01_mu1e-5_r1e-6/13000/fst_average.Rds')

max_y <- max(c(fst_shortinv_14000$av_fst, fst_shortinv_13000$av_fst))

# correcting errors from running script. Ideally would just rerun with correction genome size parameters
fst_shortinv_13000_reduced = fst_shortinv_13000[seq(1, nrow(fst_shortinv_13000), 10), ]
fst_shortinv_14000_reduced = fst_shortinv_14000[seq(1, nrow(fst_shortinv_14000), 10), ]

fst_shortinv_13000_reduced$pos <- fst_shortinv_13000_reduced$pos / 10
fst_shortinv_14000_reduced$pos <- fst_shortinv_14000_reduced$pos / 10

plot_fst13000 <- ggplot(fst_shortinv_13000_reduced, aes(x=pos, y=av_fst)) +
  geom_line()+
  ggtitle('a) Generation 13000') +
  xlab('Position') + ylab(bquote(F["ST"]*" between P1 and P2")) +
  ylim(0.09, max_y+max_y*0.05) +
  gglayer_markers_short + theme_bw() + coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))

plot_fst14000 <- ggplot(fst_shortinv_14000_reduced, aes(x=pos, y=av_fst)) +
  geom_line()+
  ggtitle('b) Generation 14000') +
  xlab('Position') + ylab(bquote(F["ST"]*" between P1 and P2")) +
  ylim(0.09, max_y+max_y*0.05) +
  gglayer_markers_short + theme_bw() + coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))

grid.arrange(plot_fst13000, plot_fst14000, nrow=1)



# ------------------------------------------------------------------------------------------
# SHORT v SHORTINV

fst_shortinv <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_SHORTINV_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds')
fst_short <- loadRData('C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_SHORT_s0.1_m0.01_mu1e-5_r1e-6/1500/fst_average.Rds')

max_y <- max(c(fst_shortinv$av_fst, fst_short$av_fst))

fst_shortinv_13000_reduced$pos <- fst_shortinv_13000_reduced$pos / 10
fst_shortinv_14000_reduced$pos <- fst_shortinv_14000_reduced$pos / 10

plot_fst_shortinv <- ggplot(fst_short, aes(x=pos, y=av_fst)) +
  geom_line()+
  ggtitle('a)') +
  xlab('Position') + ylab(bquote(F["ST"]*" between P1 and P2")) +
  ylim(0.09, max_y+max_y*0.05) +
  gglayer_markers_short + theme_bw() + coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))

plot_fst_short <- ggplot(fst_shortinv, aes(x=pos, y=av_fst)) +
  geom_line()+
  ggtitle('b)') +
  xlab('Position') + ylab(bquote(F["ST"]*" between P1 and P2")) +
  ylim(0.09, max_y+max_y*0.05) +
  gglayer_markers_short + theme_bw() + coord_cartesian(expand=F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(guide = guide_axis(angle = 45))

grid.arrange(plot_fst_shortinv, plot_fst_short, nrow=1)

























