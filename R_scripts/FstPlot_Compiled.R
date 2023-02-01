library(ggplot2)
library(gridExtra)
library(reshape2)
library(tidyr)

if(on_cluster){
  # first argument is directory name, second is generation time point
  PATH_NEUT <- "data_summary/neutral_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds"
  PATH_LAA <- "data_summary/locallyAdapted_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds"
  PATH_ADAP <- "data_summary/adaptiveInversion_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds"
  PATH_INVLAA <- "data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds"
}else{
  setwd("C:/Users/tbehr/Desktop/Thesis/data_summary")
}

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

fst_average_neutral <- loadRData("neutral_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds")
fst_average_laa <- loadRData("locallyAdapted_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds")
fst_average_adap <- loadRData("adaptiveInversion_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds")
fst_average_invlaa <- loadRData("inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds")

min_fst <- 0
max_fst <- max(c(fst_average_neutral$av_fst, fst_average_laa$av_fst, fst_average_adap$av_fst, fst_average_invlaa$av_fst))

fst_average_all <- data.frame(pos=fst_average_neutral$pos, neutral=fst_average_neutral$av_fst, LAA=fst_average_laa$av_fst,
                              inv=fst_average_adap$av_fst, invLAA=fst_average_invlaa$av_fst)

fst_average_melt <- melt(fst_average_all, id='pos')

# ggplot(fst_average_all, aes(x=pos)) +
#   geom_line(aes(y=neut), color="violet") +
#   geom_line(aes(y=laa), color="darkgreen") +
#   geom_line(aes(y=adap), color="blue") +
#   geom_line(aes(y=invlaa), color="red") +
#   ggtitle('F_ST between Populations') + 
#   xlab('Position') + ylab(expression(F[ST])) +
#   gglayer_markers +
#   ylim(c(min_fst, max_fst))

ggplot(fst_average_melt, aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('F_ST between Populations') + 
  xlab('Position') + ylab(expression(F[ST])) +
  gglayer_markers

# individual plots but with scaled axis
library(scales)
hex <- hue_pal()(4)

ggplot(melt(subset(fst_average_all, select=-c(LAA, inv, invLAA)), id='pos'), aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('F_ST between Populations') + 
  xlab('Position') + ylab(expression(F[ST])) +
  gglayer_markers +
  ylim(NA, max_fst) +
  scale_color_manual(values=c(hex[1]))

ggplot(melt(subset(fst_average_all, select=-c(inv, invLAA)), id='pos'), aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('F_ST between Populations') + 
  xlab('Position') + ylab(expression(F[ST])) +
  gglayer_markers +
  ylim(NA, max_fst) +
  scale_color_manual(values=c(hex[1], hex[2]))

ggplot(melt(subset(fst_average_all, select=-c(invLAA)), id='pos'), aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('F_ST between Populations') + 
  xlab('Position') + ylab(expression(F[ST])) +
  gglayer_markers +
  ylim(NA, max_fst) +
  scale_color_manual(values=c(hex[1], hex[2], hex[3]))

ggplot(melt(fst_average_all, id='pos'), aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('F_ST between Populations') + 
  xlab('Position') + ylab(expression(F[ST])) +
  gglayer_markers +
  ylim(NA, max_fst) +
  scale_color_manual(values=c(hex[1], hex[2], hex[3], hex[4]))


# -------------------------------------------------------------------------------------------


setwd("C:/Users/tbehr/Desktop/Thesis/data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/")

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

fst_average_mean <- loadRData("fst_average.Rds")
fst_average_1 <- loadRData("fst_average_1.Rds")
fst_average_2 <- loadRData("fst_average_2.Rds")
fst_average_3 <- loadRData("fst_average_3.Rds")

min_fst <- 0
max_fst <- max(c(fst_average_mean$av_fst, fst_average_1$av_fst, fst_average_2$av_fst, fst_average_3$av_fst))

fst_average_all <- data.frame(pos=fst_average_mean$pos, overall_mean=fst_average_mean$av_fst, samp1=fst_average_1$av_fst,
                              samp2=fst_average_2$av_fst, samp3=fst_average_3$av_fst)

fst_average_melt <- melt(fst_average_all, id='pos')

ggplot(melt(subset(fst_average_all, select=-c(samp2,samp3)), id='pos'), aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('F_ST between Populations') + 
  xlab('Position') + ylab(expression(F[ST])) +
  gglayer_markers +
  ylim(NA, max_fst) +
  scale_color_manual(values=c('black', hex[1]))

ggplot(melt(subset(fst_average_all, select=-c(samp1,samp3)), id='pos'), aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('F_ST between Populations') + 
  xlab('Position') + ylab(expression(F[ST])) +
  gglayer_markers +
  ylim(NA, max_fst) +
  scale_color_manual(values=c('black', hex[2]))

ggplot(melt(subset(fst_average_all, select=-c(samp2,samp1)), id='pos'), aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('F_ST between Populations') + 
  xlab('Position') + ylab(expression(F[ST])) +
  gglayer_markers +
  ylim(NA, max_fst) +
  scale_color_manual(values=c('black', hex[4]))



# -------------------------------------------------------------------------------------------
# ORIGINAL FST VS ADDED DELETERIOUs   
gglayer_markers <- list(
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1, FIXED_MUTATION_POS2), linetype='dashed', colour='red', alpha=0.3)},
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START, INV_END), linetype='dotted', colour='blue', alpha=0.3)}
)



fst_original_invlaa <- loadRData("adaptiveInversion_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds")
fst_deleterious_invlaa <- loadRData("adaptiveInversion_2pop_DELETERIOUS_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds")

fst_all <- data.frame(pos=fst_original_invlaa$pos, original=fst_original_invlaa$av_fst, deleterious=fst_deleterious_invlaa$av_fst)

ggplot(melt(fst_all, id='pos'), aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('F_ST between Populations') + 
  xlab('Position') + ylab(expression(F[ST])) +
  gglayer_markers +
  #ylim(NA, max_fst) +
  scale_color_manual(values=c('black', hex[3]))

# -------------------------------------------------------------------------------------------
# ORIGINAL FST VS SHORT 

gglayer_markers <- list(
  {if(LAA_PRESENT)geom_vline(xintercept = c(FIXED_MUTATION_POS1/10, FIXED_MUTATION_POS2/10), linetype='dashed', colour='red', alpha=0.3)},
  {if(INVERSION_PRESENT)geom_vline(xintercept = c(INV_START/10, INV_END/10), linetype='dotted', colour='blue', alpha=0.3)}
)

fst_original_invlaa <- loadRData("adaptiveInversion_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/fst_average.Rds")
fst_short_invlaa <- loadRData("adaptiveInversion_2pop_SHORT_s0.1_m0.01_mu1e-5_r1e-6/1500/fst_average.Rds")

fst_all <- data.frame(pos=fst_short_invlaa$pos, original=fst_original_invlaa$av_fst[seq(1,length(fst_original_invlaa$av_fst),by=10)], 
                      short=fst_short_invlaa$av_fst)

ggplot(melt(fst_all, id='pos'), aes(x=pos, y=value, col=variable)) +
  geom_line() +
  ggtitle('F_ST between Populations') + 
  xlab('Position') + ylab(expression(F[ST])) +
  gglayer_markers +
  #ylim(NA, max_fst) +
  scale_color_manual(values=c('black', hex[3]))



# -------------------------------------------------------------------------------------------
# ORIGINAL NUCDIV VS INDIVIDUAL RUNS

nucdiv_all <- loadRData("data_summary/inversionLAA_2pop_s0.1_m0.01_mu1e-5_r1e-6/15000/nucdiv_summ_p1.Rds")

theframe <- data.frame(pos=nucdiv_all$center, mean=nucdiv_all$nucdiv, single=nucdiv_summ_p1$nucdiv)


ggplot(melt(theframe, id='pos'), aes(x=pos, y=value, col=variable)) +
  geom_line() +
  scale_color_brewer(palette="Dark2")+
  ggtitle('c) Nucleotide Diversity - P1') +
  gglayer_markers +
  xlab('Position') + ylab('\u03c0')





