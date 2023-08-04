# Create Alluvial plots to illustrate clonal expansion 

install.packages("ggalluvial")
install.packages("cowplot")

library(dplyr)
library(ggalluvial)
library(readr)
library(ggplot2)
library(alakazam)
library(RColorBrewer)  # Load the RColorBrewer package
library(cowplot)

# Read in files 
HCtsv_DUPCOUNT <- read_tsv("Data/Reprocessed_HC_Combined_filtered_dist_ham_clone-pass_germ-pass_DUPCOUNT.tsv")

# Count and filter clones 
clones <- countClones(HCtsv_DUPCOUNT, group="subset") # generate clone table
clones <- clones[clones$seq_count > 1,] # subset to clones with more than 1 sequence

write.table(clones, file="Total clone count by subset.tsv", quote=FALSE, sep="\t", row.names = FALSE)

# Generate subset clone tables 
day0_clones <- subset(clones, subset == "day0")
ASC_clones <- subset(clones, subset == "ASC")
CD21lo_clones <- subset(clones, subset == "CD21lo")
FOB_clones <- subset(clones, subset == "FOB")
GC_clones <- subset(clones, subset == "GC")
ASC_FOB_CD21lo_GC <- subset(clones, subset != "day0") # subset all clones found in ASC, FOB, CD21lo, GC

# Find shared clones
common_clones <- intersect(intersect(intersect(day0_clones$clone_id, ASC_FOB_CD21lo_GC$clone_id), FOB_clones$clone_id), CD21lo_clones$clone_id) # identify clones shared by day0, ASC, FOB, and CD21lo

# Generate alluvial plot table
alluvial_common <- clones[clones$clone_id %in% common_clones, ] # subset common clones and count frequencies 
write.table(alluvial_common, file="clones shared by day0, ASC, FOB and CD21lo.tsv", quote=FALSE, sep="\t", row.names = FALSE) # export clone files for external analysis 
alluvial_common <- subset(alluvial_common, subset != "GC") # remove GC clones 
alluvial_common <- alluvial_common[c("clone_id", "subset", "seq_count", "seq_freq")] %>%
  subset(select = -seq_freq) %>% # rearrange dataframe and drop seq_freq column 
  mutate(clone_id=as.character(clone_id),
       seq_count=as.numeric(seq_count),
       subset = factor(subset, levels= (c("day0", "FOB", "CD21lo", "ASC")))
       )

# generate colors vector
values = c("237071"="#BC5C7C", "212983" = "#2DBDCB", "258497" = "#E781A3", "258551" = "#ACCB99", 
           "226086" = "#b3c4d2", "135564" = "#d9e0e6", "221617" = "#99CBBA", "254088" = "#80e0e0", 
           "258509" = "#8061e6", "228298" = "#FD7373", "232052" = "#FDA873" , "239718" = "#FDD873", 
           "243126" = "#B1FFD7", "247070" = "#B1F9FF", "234414" = "#B1E6FF" , "220708" = "#B1D1FF", 
           "230537" = "#B5B1FF", "217635" = "#E8B1FF" , "243277" = "#CCB1FF", "230530" = "#FFB1F8", 
           "146397" = "#FFB1D9", "234430" = "#FFB1B8" , "226086" = "#FFCEB1", "230530" = "#B5E4A1", 
           "254062" = "#A1E4C1" , "243277" = "#A1DEE4", "224647" = "#EF818D", "196673" = "#81EFBB",
           "258551" = "#EFD281", "228298" = "#81EFCE", "258551" = "#81EFE5", "206144" = "#A295CB",
           "230530" = "#B2E1C1", "243126" = "#E0A1A1", "223089" = "#A1D4E0")

# generate alluvial plot
ggplot(alluvial_common,
       aes(y = seq_count, x = subset)) +
  geom_flow(aes(alluvium = clone_id, fill = clone_id), 
            alpha= 0.5, 
            lty = 1, color = "black",
            curve_type = "cubic", 
            width = 0.3,
            linewidth = 0) +
  geom_col(aes(fill = clone_id), width = .3, color = "black", alpha = 1, linewidth = 0.3) +
  scale_fill_manual(values = values, na.value = "white") +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() + 
  ylab("Sequence Count") +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 30)) + 
  theme(axis.title = element_text(size = 30))

# Alluvial plot without ASCs 
## subset alluvial data 
alluvial_common_minusASC <- subset(alluvial_common, subset != "ASC") # remove ASC clones 

# plot
ggplot(alluvial_common_minusASC,
       aes(y = seq_count, x = subset)) +
  geom_flow(aes(alluvium = clone_id, fill = clone_id), 
            alpha= 0.5, 
            lty = 1, color = "black",
            curve_type = "cubic", 
            width = 0.3,
            lwd = 0) +
  geom_col(aes(fill = clone_id), width = .3, color = "black", alpha = 1, linewidth = 0.3) +
  scale_fill_manual(values = values, na.value = "white") +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() + 
  ylab("Sequence Count") +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 20)) + 
  theme(axis.title = element_text(size = 20))


