#Immcantation - Alakazam for BCR-seq repertoire analysis

install.packages("BiocManager")
install.packages("alakazam") #install alakazam from CRAN
install.packages("shazam") #install shazam from CRAN
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "Rcpp")) #install build dependencies 
library(devtools)
install_bitbucket("kleinstein/alakazam@master")
install_bitbucket("kleinstein/shazam@master")

library("alakazam")
library("shazam")
library("dplyr")
library("scales")

# Read heavy chain tsv files generated from ChangeO 
day0_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/day0_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M1_ASC_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M1_ASC_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE) 
M1_CD21lo_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M1_CD21lo_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M1_FOB_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M1_FOB_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M1_GC_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M1_GC_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M2_ASC_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M2_ASC_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE) 
M2_CD21lo_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M2_CD21lo_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M2_FOB_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M2_FOB_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M2_GC_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M2_GC_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M3_ASC_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M3_ASC_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M3_CD21lo_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M3_CD21lo_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M3_FOB_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M3_FOB_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M3_GC_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M3_GC_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M4_ASC_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M4_ASC_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M4_CD21lo_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M4_CD21lo_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M4_FOB_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M4_FOB_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)
M4_GC_HC <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/M4_GC_reprocess_heavy_parse-select.tsv", select = NULL, drop = NULL, seq_upper = TRUE)

# Add metadata columns 
day0_HC_meta <- mutate(day0_HC, Mouse="M0", Subset="day0")
M1_ASC_HC_meta <- mutate(M1_ASC_HC, Mouse="M1", Subset="ASC")
M1_CD21lo_HC_meta <- mutate(M1_CD21lo_HC, Mouse="M1", Subset="CD21lo")
M1_FOB_HC_meta <- mutate(M1_FOB_HC, Mouse="M1", Subset="FOB")
M1_GC_HC_meta <- mutate(M1_GC_HC, Mouse="M1", Subset="GC")
M2_ASC_HC_meta <- mutate(M2_ASC_HC, Mouse="M2", Subset="ASC")
M2_CD21lo_HC_meta <- mutate(M2_CD21lo_HC, Mouse="M2", Subset="CD21lo")
M2_FOB_HC_meta <- mutate(M2_FOB_HC, Mouse="M2", Subset="FOB")
M2_GC_HC_meta <- mutate(M2_GC_HC, Mouse="M2", Subset="GC")
M3_ASC_HC_meta <- mutate(M3_ASC_HC, Mouse="M3", Subset="ASC")
M3_CD21lo_HC_meta <- mutate(M3_CD21lo_HC, Mouse="M3", Subset="CD21lo")
M3_FOB_HC_meta <- mutate(M3_FOB_HC, Mouse="M3", Subset="FOB")
M3_GC_HC_meta <- mutate(M3_GC_HC, Mouse="M3", Subset="GC")
M4_ASC_HC_meta <- mutate(M4_ASC_HC, Mouse="M4", Subset="ASC")
M4_CD21lo_HC_meta <- mutate(M4_CD21lo_HC, Mouse="M4", Subset="CD21lo")
M4_FOB_HC_meta <- mutate(M4_FOB_HC, Mouse="M4", Subset="FOB")
M4_GC_HC_meta <- mutate(M4_GC_HC, Mouse="M4", Subset="GC")

# Merge all dataframes into one mega data file
HC_combineAll <- 
  bind_rows(day0_HC_meta, M1_ASC_HC_meta, M1_CD21lo_HC_meta, M1_FOB_HC_meta, M1_GC_HC_meta,
            M2_ASC_HC_meta, M2_CD21lo_HC_meta, M2_FOB_HC_meta, M2_GC_HC_meta,
            M3_ASC_HC_meta, M3_CD21lo_HC_meta, M3_FOB_HC_meta, M3_GC_HC_meta,
            M4_ASC_HC_meta, M4_CD21lo_HC_meta, M4_FOB_HC_meta, M4_GC_HC_meta)

# Determine the number of informative sites of sequence alignment
n_informative <- str_count(HC_combineAll$sequence_alignment, "[ACTG]")
greater_than_250 <- n_informative > 250 # Create a logical vector indicating which counts are greater than 250
count_greater_than_250 <- sum(greater_than_250) # Use the sum function to count the number of TRUE values
print(count_greater_than_250) # Print the result

# Filter to sequences with 250 or more informative sites
HC_combineAll_250 <- HC_combineAll %>% 
  filter(n_informative > 250)

# Write file to table for germline assignment  
write.table(HC_combineAll_250, file="HC_combineAll_250.tsv", quote=FALSE, sep="\t", row.names = FALSE)

##################### ChangeO germline assignment #######################

# Read germline assigned tsv file generated from ChangeO 
HC_combineAll_250_germline <- readChangeoDb("Data/ChangeO HC_LC tsv_deduplicated_filtered/HC_combineAll_250_germ-pass.tsv", select = NULL, drop = NULL, seq_upper = TRUE)

# Filter to non-chimeric sequences 
is_chimeric <- slideWindowDb(
  HC_combineAll_250_germline, 
  sequenceColumn = "sequence_alignment",
  germlineColumn = "germline_alignment_d_mask", 
  mutThres = 6,
  windowSize = 10
) # extract chimeric sequence indexes 

is_chimeric <- data.frame(is_chimeric) # change chimeric indexes to dataframe

HC_combineAll_250_germline_non_chimeric <- HC_combineAll_250_germline[!is_chimeric,]

# Define clonal assignment threshold 
## Use Nucleotide Hamming distance and normalize by junction length
dist_ham <- distToNearest(HC_combineAll_250_germline_non_chimeric,
                          sequenceColumn = "junction",
                          vCallColumn = "v_call", jCallColumn = "j_call",
                          model = "ham", normalize = "len", nproc=1)

# Export dist_ham as TSV file
write.table(dist_ham, file="Reprocessed_HC_Combined_filtered_dist_ham.tsv", quote=FALSE, sep="\t", row.names = FALSE)

##Threshold determination by manual inspection 
### Generate Hamming distance histogram 
HamPlot <- ggplot(subset(dist_ham, !is.na(dist_nearest)),
                  aes(x=dist_nearest)) +
  theme_bw() +
  xlab("Hamming distance") +
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.005) +
  geom_vline(xintercept=0.045, color="firebrick", linetype=2)
plot(HamPlot)

##Generate grouped histograms
HamPlot_by_subset <- ggplot(subset(dist_ham, !is.na(dist_nearest)),
                            aes(x=dist_nearest)) +
  theme_bw() +
  xlab("Hamming distance by subset") +
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.005) +
  geom_vline(xintercept=0.045, color="firebrick", linetype=2) +
  facet_grid(subset ~., scales="free_y")
plot(HamPlot_by_subset)

##################### ChangeO clonal clustering #######################


