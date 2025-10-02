# This script performs clonal diversity analysis via Alakazam and mutation analysis via Shazam

install.packages("ggsignif")
install.packages("ggpubr")
install.packages("ggplot2")

library("alakazam")
library("shazam")
library("dplyr")
library("ggplot2")
library("ggsignif")
library("ggpubr")
library("readr")

# Import dataset including duplicate count
HCtsv_DUPCOUNT <- read_tsv("Data/ChangeO HC_LC tsv_deduplicated_filtered/Reprocessed_HC_Combined_filtered_dist_ham_clone-pass_germ-pass_DUPCOUNT.tsv")

# Partition data based on subset, weigh clone sizes by the DUPCOUNT column 
clone_count <- countClones(HCtsv_DUPCOUNT, group="subset", copy = "DUPCOUNT", clone="clone_id")

# Partition data based on subset, and calculate a 95% confidence interval via 100 bootstrap realizations 
curve <- estimateAbundance(HCtsv_DUPCOUNT, group="subset", ci=0.95, nboot=100, clone="clone_id")

# Plot rank abundance curve of the relative clonal abundances 
## set color palette
sample_colors <- c("day0"="black", "ASC"="#69BED8", "FOB"="#E7899E", "CD21lo"="#CF76D6", "GC"="#F8AC6D")
## plot abundance curve 
plot(curve, colors = sample_colors, legend_title="Subset", 
     axis.title = element_text(size = 14),
     plot.title = element_text(size = 16),
     legend.text = element_text(size = 12), 
     panel.border = element_blank(),
     axis.line = element_line(color = "black"))

# Generate diversity curve
## Compare diversity curve across values in subset column
## q ranges from 0 (min_q=0) to 4 (max_q=4) in 0.05 increments (step_q=0.05)
## A 95% confidence interval will be calculated (ci=0.95)
## 100 resampling realizations are performed (nboot=100)
subset_diversity_curve <- alphaDiversity(HCtsv_DUPCOUNT, group = "subset", clone="clone_id",
                                         min_q=0, max_q=, step_q=0.1, 
                                         ci=0.95, nboot=100)
## Plot a log-log (log_q=TRUE, log_d=TRUE) plot of sample diversity
subset_main <- paste0("Subset diversity")
p <- plot(subset_diversity_curve, 
          axis= element_text(size=20),
          colors=sample_colors, 
          main_title=subset_main,
          legend_title="Subset",
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16),
          legend.text = element_text(size = 12), 
          panel.border = element_blank(),
          axis.line = element_line(color = "black"),
          )

p + geom_vline(xintercept = c(0,1,2), color = "grey50", linetype = "dashed") +
  geom_text(data = data.frame(q = c(0,1,2), y = round(max(p$data$d_upper)/2),
                              label = c("Richeness", "Shannon", "Simpson")),
            aes(x = q, y = y, label = label), size = 5, angle = 90, vjust = -0.4, inherit.aes = F, color = "grey50") 

#### --------------------------Mutation analysis-----------------------------------

# Filter to non-chimeric sequences 
is_chimeric <- slideWindowDb(
  HCtsv_DUPCOUNT, 
  sequenceColumn = "sequence_alignment",
  germlineColumn = "germline_alignment_d_mask", 
  mutThres = 6,
  windowSize = 10
) # extract chimeric sequence indexes 

is_chimeric <- data.frame(is_chimeric) # change chimeric indexes to dataframe

HCtsv_DUPCOUNT_non_chimeric <- HCtsv_DUPCOUNT[!is_chimeric,]

# Calculate total R and S mutation counts
## Add R and S mutation counts as separate columns 
HCtsv_DUPCOUNT_obs_count <- observedMutations(HCtsv_DUPCOUNT_non_chimeric, sequenceColumn = "sequence_alignment", 
                                        germlineColumn = "germline_alignment_d_mask",
                                        regionDefinition = NULL,
                                        frequency=FALSE, #new column added as mutation counts
                                        nproc=1)

## Add total mutation counts as a single column 
HCtsv_DUPCOUNT_obs_count_combined <- observedMutations(HCtsv_DUPCOUNT_non_chimeric, sequenceColumn = "sequence_alignment", 
                                                 germlineColumn = "germline_alignment_d_mask",
                                                 regionDefinition = NULL,
                                                 combine = TRUE, #combine R and S mutations into a single column
                                                 frequency=FALSE,#new column added as mutation counts  
                                                 nproc=1)

## Add R and S mutation frequencies as separate columns 
HCtsv_DUPCOUNT_obs <- observedMutations(HCtsv_DUPCOUNT_non_chimeric, sequenceColumn = "sequence_alignment", 
                                        germlineColumn = "germline_alignment_d_mask",
                                        regionDefinition = NULL,
                                        frequency=TRUE, #new column added as mutation frequencies
                                        nproc=1)

## Add total mutation frequencies as a single column 
HCtsv_DUPCOUNT_obs_combined <- observedMutations(HCtsv_DUPCOUNT_non_chimeric, sequenceColumn = "sequence_alignment", 
                                        germlineColumn = "germline_alignment_d_mask",
                                        regionDefinition = NULL,
                                        combine = TRUE, #combine R and S mutations into a single column
                                        frequency=TRUE, #new column added as mutation frequencies
                                        nproc=1)

# Plot mutation frequencies (overall)
MutFreq_plot <- ggplot(HCtsv_DUPCOUNT_obs_combined, aes(x=subset, y=mu_freq, fill=subset)) +
  geom_violin(width=1.4, trim = TRUE) +
  theme_classic() + ggtitle("Total mutations") +
  xlab("Subset") + ylab("Mutation frequency") +
  scale_fill_manual(name="subset", values = sample_colors) +
  coord_cartesian(ylim=c(0,0.06))
plot(MutFreq_plot)

# Plot mutation count (overall)

sample_colors_mutation <- c("day0"="grey", "ASC"="#69BED8", "FOB"="#E7899E", "CD21lo"="#CF76D6", "GC"="#F8AC6D")

MutCount_plot <- ggplot(HCtsv_DUPCOUNT_obs_count_combined, aes(x=subset, y=mu_count, fill=subset)) +
  geom_violin(width=1.4, trim = TRUE, adjust = 6, linewidth = 0.5) +
  geom_boxplot(outlier.shape = NA, width=0.08, color="black", alpha=0.2) +
  theme_classic() + ggtitle("Total mutations") +
  xlab("Subset") + ylab("Mutation count") +
  scale_fill_manual(name="subset", values = sample_colors_mutation) +
  coord_cartesian(ylim=c(0,20)) +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14))
plot(MutCount_plot)

# Plot replacement mutation count

MutCount_plot_R <- ggplot(HCtsv_DUPCOUNT_obs_count, aes(x=subset, y=mu_count_seq_r, fill=subset)) +
  geom_violin(width=1, trim = TRUE, adjust = 5) +
  geom_boxplot(outlier.shape = NA, width=0.08, color="black", alpha=0.2) +
  theme_classic() + ggtitle("replacement mutations") +
  xlab("Subset") + ylab("Replacement mutation count") +
  scale_fill_manual(name="subset", values = sample_colors_mutation) +
  coord_cartesian(ylim=c(0,20)) +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14))

plot(MutCount_plot_R)




