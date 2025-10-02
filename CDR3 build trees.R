# Build trees for expanded clones using CDR3 after filtering 

library(dplyr)
library(readr)
library(ggplot2)
library(alakazam)
library(RColorBrewer)  # Load the RColorBrewer package
library(shazam)
library(ggtree)
library(dowser)

# read in file 
HCtsv_DUPCOUNT <- read_tsv("Data/Reprocessed_HC_Combined_filtered_dist_ham_clone-pass_germ-pass_DUPCOUNT.tsv")

# subset clones 
clone_237071 <- subset(HCtsv_DUPCOUNT,  clone_id == "237071")
clone_212983 <- subset(HCtsv_DUPCOUNT,  clone_id == "212983")

# Filter to non-chimeric sequences 
is_chimeric_237071 <- slideWindowDb(
  clone_237071, 
  sequenceColumn = "sequence_alignment",
  germlineColumn = "germline_alignment_d_mask", 
  mutThres = 6,
  windowSize = 10
) # extract chimeric sequence indexes 

is_chimeric_237071 <- data.frame(is_chimeric_237071) # change chimeric indexes to dataframe

is_chimeric_212983 <- slideWindowDb(
  clone_212983, 
  sequenceColumn = "sequence_alignment",
  germlineColumn = "germline_alignment_d_mask", 
  mutThres = 6,
  windowSize = 10
) # extract chimeric sequence indexes 

is_chimeric_212983 <- data.frame(is_chimeric_212983) # change chimeric indexes to dataframe

clone_237071_non_chimeric <- clone_237071[!is_chimeric_237071,]
clone_212983_non_chimeric <- clone_212983[!is_chimeric_212983,]

# Filter to primary V gene assignment 
clone_237071_non_chimeric <- clone_237071_non_chimeric %>% 
  filter(v_call == "Musmus IGHV1-72*01 F")

clone_212983_non_chimeric <- clone_212983_non_chimeric %>% 
  filter(v_call == "Musmus IGHV1-42*01 F")

# Extract CDR3 germline sequences
## 237071
clone_237071_non_chimeric <- clone_237071_non_chimeric %>% 
  mutate(cdr3_germline = substr(germline_alignment_d_mask, 313, 333))

# build trees
formatclones_237071 <- formatClones(clone_237071_non_chimeric, 
                                    seq = "cdr3",
                                    germ = "cdr3_germline",
                                    traits=c("subset"), 
                                    num_fields=c("DUPCOUNT"), 
                                    columns=c("subset", "mouse"), 
                                    minseq = 20, 
                                    nproc = 8)
cloneTree_mp_237071 <- getTrees(formatclones_237071, nproc=8, collapse = TRUE)
mp_plot_237071 <- plotTrees(cloneTree_mp_237071, tips="subset", tipsize="DUPCOUNT", palette = "Set1", layout = "circular", scale = 0.05)

cloneTree_ml_237071 <- getTrees(formatclones_237071, build="pml", collapse = TRUE)
ml_plots_237071 <- plotTrees(cloneTree_ml_237071, tips="subset", tipsize="DUPCOUNT", palette = "Set1", layout = "circular", scale = 0.01)


## 212983
clone_212983_non_chimeric <- clone_212983_non_chimeric %>% 
  mutate(cdr3_germline = substr(germline_alignment_d_mask, 313, 333))

# build trees
formatclones_212983 <- formatClones(clone_212983_non_chimeric, 
                                    seq = "cdr3",
                                    germ = "cdr3_germline",
                                    traits=c("subset"), 
                                    num_fields=c("DUPCOUNT"), 
                                    columns=c("subset", "mouse"), 
                                    minseq = 20, 
                                    nproc = 8)
cloneTree_mp_212983 <- getTrees(formatclones_212983, nproc=8, collapse = TRUE)

mp_plots_212983 <- plotTrees(cloneTree_mp_212983, tips="subset", tipsize="DUPCOUNT", palette = "Set2", layout = "rectangular", scale = 0.01)

cloneTree_ml_212983 <- getTrees(formatclones_212983, build="pml", collapse = TRUE)
ml_plots_212983 <- plotTrees(cloneTree_ml_212983, tips="subset", tipsize="DUPCOUNT", palette = "Set2", layout = "circular", scale = 0.01)

