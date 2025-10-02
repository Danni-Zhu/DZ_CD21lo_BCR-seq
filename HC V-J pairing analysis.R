# compute Heavy chain V-J pairs for different subsets 

library(dplyr)
library(circlize)
library(readr)

# Read in file
HCtsv_DUPCOUNT <- read_tsv("Data/Reprocessed_HC_Combined_filtered_dist_ham_clone-pass_germ-pass_DUPCOUNT.tsv")

# subset data 
## subset to day0 data 
day0_HC_all <- HCtsv_DUPCOUNT[HCtsv_DUPCOUNT$subset == "day0", ]

## subset to ASC data
ASC_HC_all <- HCtsv_DUPCOUNT[HCtsv_DUPCOUNT$subset == "ASC", ]
ASC_HC_M1 <- ASC_HC_all[ASC_HC_all$mouse == "M1", ]
ASC_HC_M2 <- ASC_HC_all[ASC_HC_all$mouse == "M2", ]
ASC_HC_M3 <- ASC_HC_all[ASC_HC_all$mouse == "M3", ]
ASC_HC_M4 <- ASC_HC_all[ASC_HC_all$mouse == "M4", ]

## subset to FOB data
FOB_HC_all <- HCtsv_DUPCOUNT[HCtsv_DUPCOUNT$subset == "FOB", ]
FOB_HC_M1 <- FOB_HC_all[FOB_HC_all$mouse == "M1", ]
FOB_HC_M2 <- FOB_HC_all[FOB_HC_all$mouse == "M2", ]
FOB_HC_M3 <- FOB_HC_all[FOB_HC_all$mouse == "M3", ]
FOB_HC_M4 <- FOB_HC_all[FOB_HC_all$mouse == "M4", ]

## subset to CD21lo data 
CD21lo_HC_all <- HCtsv_DUPCOUNT[HCtsv_DUPCOUNT$subset == "CD21lo", ]
CD21lo_HC_M1 <- CD21lo_HC_all[CD21lo_HC_all$mouse == "M1", ]
CD21lo_HC_M2 <- CD21lo_HC_all[CD21lo_HC_all$mouse == "M2", ]
CD21lo_HC_M3 <- CD21lo_HC_all[CD21lo_HC_all$mouse == "M3", ]
CD21lo_HC_M4 <- CD21lo_HC_all[CD21lo_HC_all$mouse == "M4", ]

## subset to GC data 
GC_HC_all <- HCtsv_DUPCOUNT[HCtsv_DUPCOUNT$subset == "GC", ]
GC_HC_M1 <- GC_HC_all[GC_HC_all$mouse == "M1", ]
GC_HC_M2 <- GC_HC_all[GC_HC_all$mouse == "M2", ]
GC_HC_M3 <- GC_HC_all[GC_HC_all$mouse == "M3", ]
GC_HC_M4 <- GC_HC_all[GC_HC_all$mouse == "M4", ]

# Compute V-J pair counts
## day0 FoB
pair_counts_day0 <- day0_HC_all %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_day0, file="day0 VJ pair counts.tsv", quote=FALSE, sep="\t", row.names = FALSE)

## ASCs
pair_counts_ASC <- ASC_HC_all %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

### mouse 1 ### 
pair_counts_ASC_M1 <- ASC_HC_M1 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_ASC_M1, file="ASC VJ pair counts M1.tsv", quote=FALSE, sep="\t", row.names = FALSE)

### mouse 2 ### 
pair_counts_ASC_M2 <- ASC_HC_M2 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_ASC_M2, file="ASC VJ pair counts M2.tsv", quote=FALSE, sep="\t", row.names = FALSE)

### mouse 3 ### 
pair_counts_ASC_M3 <- ASC_HC_M3 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_ASC_M3, file="ASC VJ pair counts M3.tsv", quote=FALSE, sep="\t", row.names = FALSE)

### mouse 4 ### 
pair_counts_ASC_M4 <- ASC_HC_M4 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_ASC_M4, file="ASC VJ pair counts M4.tsv", quote=FALSE, sep="\t", row.names = FALSE)

## FoB
pair_counts_FoB <- FOB_HC_all %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

### mouse 1 ### 
pair_counts_FoB_M1 <- FOB_HC_M1 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_FoB_M1, file="FOB VJ pair counts M1.tsv", quote=FALSE, sep="\t", row.names = FALSE)

### mouse 2 ### 
pair_counts_FoB_M2 <- FOB_HC_M2 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_FoB_M2, file="FOB VJ pair counts M2.tsv", quote=FALSE, sep="\t", row.names = FALSE)

### mouse 3 ### 
pair_counts_FoB_M3 <- FOB_HC_M3 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_FoB_M3, file="FOB VJ pair counts M3.tsv", quote=FALSE, sep="\t", row.names = FALSE)

### mouse 4 ### 
pair_counts_FoB_M4 <- FOB_HC_M4 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_FoB_M4, file="FOB VJ pair counts M4.tsv", quote=FALSE, sep="\t", row.names = FALSE)

## CD21lo
pair_counts_CD21lo <- CD21lo_HC_all %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

### mouse 1 ###
pair_counts_CD21lo_M1 <- CD21lo_HC_M1 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_CD21lo_M1, file="CD21lo VJ pair counts M1.tsv", quote=FALSE, sep="\t", row.names = FALSE)

### mouse 2 ###
pair_counts_CD21lo_M2 <- CD21lo_HC_M2 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_CD21lo_M2, file="CD21lo VJ pair counts M2.tsv", quote=FALSE, sep="\t", row.names = FALSE)

### mouse 3 ###
pair_counts_CD21lo_M3 <- CD21lo_HC_M3 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_CD21lo_M3, file="CD21lo VJ pair counts M3.tsv", quote=FALSE, sep="\t", row.names = FALSE)

### mouse 4 ###
pair_counts_CD21lo_M4 <- CD21lo_HC_M4 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_CD21lo_M4, file="CD21lo VJ pair counts M4.tsv", quote=FALSE, sep="\t", row.names = FALSE)

## GC
pair_counts_GC <- GC_HC_all %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

### mouse 1 ###
pair_counts_GC_M1 <- GC_HC_M1 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_GC_M1, file="GC VJ pair counts M1.tsv", quote=FALSE, sep="\t", row.names = FALSE)

### mouse 2 ###
pair_counts_GC_M2 <- GC_HC_M2 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_GC_M2, file="GC VJ pair counts M2.tsv", quote=FALSE, sep="\t", row.names = FALSE)

### mouse 3 ### 
pair_counts_GC_M3 <- GC_HC_M3 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_GC_M3, file="GC VJ pair counts M3.tsv", quote=FALSE, sep="\t", row.names = FALSE)

### mouse 4 ### 
pair_counts_GC_M4 <- GC_HC_M4 %>% 
  group_by(germline_v_call, germline_j_call) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts)) %>%
  ungroup()

write.table(pair_counts_GC_M4, file="GC VJ pair counts M4.tsv", quote=FALSE, sep="\t", row.names = FALSE)

# Build edge lists
## day 0 FoB

day0_top_pairs <- pair_counts_day0 %>% 
  slice_max(counts, n = 100) 

day0_VJlist <- day0_HC_all %>% 
  select(germline_v_call, germline_j_call)

day0_VJlist <- inner_join(day0_VJlist, day0_top_pairs, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

day0_VJlist <- subset(day0_VJlist, select = -c(counts))

day0_VJ_adj <- with(day0_VJlist, table(germline_v_call, germline_j_call))

day0_gridcol <- rep("grey", nrow(day0_VJ_adj)+ncol(day0_VJ_adj)) 

day0_gridcol[c(18, 20, 22, 36, 37, 38, 39)] <- c("#FFA500", "#DB7093", "#87CEEB", "#963C78", "#B06D9A", "#CB9DBB", "#E5CEDD")

day0_linkcol <- matrix("grey", nrow = nrow(day0_VJ_adj), ncol = ncol(day0_VJ_adj))

IGHV1_82 <- match("IGHV1-82*01", rownames(day0_VJ_adj))
IGHV1_80 <- match("IGHV1-80*01", rownames(day0_VJ_adj))
IGHV1_72 <- match("IGHV1-72*01", rownames(day0_VJ_adj))

day0_linkcol[IGHV1_72, ] <- "#FFA500"
day0_linkcol[IGHV1_80, ] <- "#DB7093"
day0_linkcol[IGHV1_82, ] <- "#87CEEB"

chordDiagram(day0_VJ_adj, transparency = 0.5, annotationTrack = "grid", col = day0_linkcol, grid.col = day0_gridcol,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(day0_VJ_adj))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  label <- CELL_META$sector.index
  if (label == "IGHV1-72*01" | label == "IGHV1-80*01" | label == "IGHV1-82*01" | label == "IGHJ1*03" | label == "IGHJ2*01" | label == "IGHJ3*01" | label == "IGHJ4*01") {
    col <- "black"
  } else {
    col <- "white"
  }
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 1.2, col = col)
}, bg.border = NA)

## FOB 

FOB_top_pairs <- pair_counts_FoB %>% 
  slice_max(counts, n = 100) 

FOB_VJlist <- FOB_HC_all %>% 
  select(germline_v_call, germline_j_call)

FOB_VJlist <- inner_join(FOB_VJlist, FOB_top_pairs, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

FOB_VJlist <- subset(FOB_VJlist, select = -c(counts))

FOB_VJ_adj <- with(FOB_VJlist, table(germline_v_call, germline_j_call))

FOB_gridcol <- rep("grey", nrow(FOB_VJ_adj)+ncol(FOB_VJ_adj)) 

FOB_gridcol[c(18, 20, 22, 34, 35, 36, 37)] <- c("#FFA500", "#DB7093", "#87CEEB", "#963C78", "#B06D9A", "#CB9DBB", "#E5CEDD")

FOB_linkcol <- matrix("grey", nrow = nrow(FOB_VJ_adj), ncol = ncol(FOB_VJ_adj))

IGHV1_82 <- match("IGHV1-82*01", rownames(FOB_VJ_adj))
IGHV1_80 <- match("IGHV1-80*01", rownames(FOB_VJ_adj))
IGHV1_72 <- match("IGHV1-72*01", rownames(FOB_VJ_adj))

FOB_linkcol[IGHV1_72, ] <- "#FFA500"
FOB_linkcol[IGHV1_80, ] <- "#DB7093"
FOB_linkcol[IGHV1_82, ] <- "#87CEEB"

chordDiagram(FOB_VJ_adj, transparency = 0.5, annotationTrack = "grid", col = FOB_linkcol, grid.col = FOB_gridcol,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(FOB_VJ_adj))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  label <- CELL_META$sector.index
  if (label == "IGHV1-72*01" | label == "IGHV1-80*01" | label == "IGHV1-82*01" | label == "IGHJ1*03" | label == "IGHJ2*01" | label == "IGHJ3*01" | label == "IGHJ4*01") {
    col <- "black"
  } else {
    col <- "white"
  }
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 1.2, col = col)
}, bg.border = NA)

### mouse 1 ### 
FOB_top_pairs_M1 <- pair_counts_FoB_M1 %>% 
  slice_max(counts, n = 100) 

FOB_VJlist_M1 <- FOB_HC_M1 %>% 
  select(germline_v_call, germline_j_call)

FOB_VJlist_M1 <- inner_join(FOB_VJlist_M1, FOB_top_pairs_M1, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

FOB_VJlist_M1 <- subset(FOB_VJlist_M1, select = -c(counts))

FOB_VJ_adj_M1 <- with(FOB_VJlist_M1, table(germline_v_call, germline_j_call))

chordDiagram(FOB_VJ_adj_M1, transparency = 0.5, annotationTrack = "grid",
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(FOB_VJ_adj_M1))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.8 )
}, bg.border = NA)

### mouse 2 ### 
FOB_top_pairs_M2 <- pair_counts_FoB_M2 %>% 
  slice_max(counts, n = 100) 

FOB_VJlist_M2 <- FOB_HC_M2 %>% 
  select(germline_v_call, germline_j_call)

FOB_VJlist_M2 <- inner_join(FOB_VJlist_M2, FOB_top_pairs_M2, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

FOB_VJlist_M2 <- subset(FOB_VJlist_M2, select = -c(counts))

FOB_VJ_adj_M2 <- with(FOB_VJlist_M2, table(germline_v_call, germline_j_call))

chordDiagram(FOB_VJ_adj_M2, transparency = 0.5, annotationTrack = "grid",
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(FOB_VJ_adj_M2))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.8 )
}, bg.border = NA)

### mouse 3 ### 
FOB_top_pairs_M3 <- pair_counts_FoB_M3 %>% 
  slice_max(counts, n = 100) 

FOB_VJlist_M3 <- FOB_HC_M3 %>% 
  select(germline_v_call, germline_j_call)

FOB_VJlist_M3 <- inner_join(FOB_VJlist_M3, FOB_top_pairs_M3, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

FOB_VJlist_M3 <- subset(FOB_VJlist_M3, select = -c(counts))

FOB_VJ_adj_M3 <- with(FOB_VJlist_M3, table(germline_v_call, germline_j_call))

chordDiagram(FOB_VJ_adj_M3, transparency = 0.5, annotationTrack = "grid",
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(FOB_VJ_adj_M3))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.8 )
}, bg.border = NA)

### mouse 4 ### 
FOB_top_pairs_M4 <- pair_counts_FoB_M4 %>% 
  slice_max(counts, n = 100) 

FOB_VJlist_M4 <- FOB_HC_M4 %>% 
  select(germline_v_call, germline_j_call)

FOB_VJlist_M4 <- inner_join(FOB_VJlist_M4, FOB_top_pairs_M4, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

FOB_VJlist_M4 <- subset(FOB_VJlist_M4, select = -c(counts))

FOB_VJ_adj_M4 <- with(FOB_VJlist_M4, table(germline_v_call, germline_j_call))

chordDiagram(FOB_VJ_adj_M4, transparency = 0.5, annotationTrack = "grid",
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(FOB_VJ_adj_M4))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.8 )
}, bg.border = NA)


## ASCs 

ASC_top_pairs <- pair_counts_ASC %>% 
  slice_max(counts, n = 100) 

ASC_VJlist <- ASC_HC_all %>% 
  select(germline_v_call, germline_j_call)

ASC_VJlist <- inner_join(ASC_VJlist, ASC_top_pairs, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

ASC_VJlist <- subset(ASC_VJlist, select = -c(counts))

ASC_VJ_adj <- with(ASC_VJlist, table(germline_v_call, germline_j_call))

ASC_gridcol <- rep("grey", nrow(ASC_VJ_adj)+ncol(ASC_VJ_adj)) 

ASC_gridcol[c(20, 23, 25, 37, 38, 39, 40)] <- c("#FFA500", "#DB7093", "#87CEEB", "#963C78", "#B06D9A", "#CB9DBB", "#E5CEDD")

ASC_linkcol <- matrix("grey", nrow = nrow(ASC_VJ_adj), ncol = ncol(ASC_VJ_adj))

IGHV1_82 <- match("IGHV1-82*01", rownames(ASC_VJ_adj))
IGHV1_80 <- match("IGHV1-80*01", rownames(ASC_VJ_adj))
IGHV1_72 <- match("IGHV1-72*01", rownames(ASC_VJ_adj))

ASC_linkcol[IGHV1_72, ] <- "#FFA500"
ASC_linkcol[IGHV1_80, ] <- "#DB7093"
ASC_linkcol[IGHV1_82, ] <- "#87CEEB"

chordDiagram(ASC_VJ_adj, transparency = 0.5, annotationTrack = "grid", col = ASC_linkcol, grid.col = ASC_gridcol,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(ASC_VJ_adj))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  label <- CELL_META$sector.index
  if (label == "IGHV1-72*01" | label == "IGHV1-80*01" | label == "IGHV1-82*01" | label == "IGHJ1*03" | label == "IGHJ2*01" | label == "IGHJ3*01" | label == "IGHJ4*01") {
    col <- "black"
  } else {
    col <- "white"
  }
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 1.2, col = col)
}, bg.border = NA)

### mouse 1 ###
ASC_top_pairs_M1 <- pair_counts_ASC_M1 %>% 
  slice_max(counts, n = 100) 

ASC_VJlist_M1 <- ASC_HC_M1 %>% 
  select(germline_v_call, germline_j_call)

ASC_VJlist_M1 <- inner_join(ASC_VJlist_M1, ASC_top_pairs_M1, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

ASC_VJlist_M1 <- subset(ASC_VJlist_M1, select = -c(counts))

ASC_VJ_adj_M1 <- with(ASC_VJlist_M1, table(germline_v_call, germline_j_call))

chordDiagram(ASC_VJ_adj_M1, transparency = 0.5, annotationTrack = "grid",
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(ASC_VJ_adj_M1))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.8 )
}, bg.border = NA)

### mouse 2 ###
ASC_top_pairs_M2 <- pair_counts_ASC_M2 %>% 
  slice_max(counts, n = 100) 

ASC_VJlist_M2 <- ASC_HC_M2 %>% 
  select(germline_v_call, germline_j_call)

ASC_VJlist_M2 <- inner_join(ASC_VJlist_M2, ASC_top_pairs_M2, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

ASC_VJlist_M2 <- subset(ASC_VJlist_M2, select = -c(counts))

ASC_VJ_adj_M2 <- with(ASC_VJlist_M2, table(germline_v_call, germline_j_call))

chordDiagram(ASC_VJ_adj_M2, transparency = 0.5, annotationTrack = "grid",
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(ASC_VJ_adj_M2))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.8 )
}, bg.border = NA)

### mouse 3 ### 
ASC_top_pairs_M3 <- pair_counts_ASC_M3 %>% 
  slice_max(counts, n = 100) 

ASC_VJlist_M3 <- ASC_HC_M3 %>% 
  select(germline_v_call, germline_j_call)

ASC_VJlist_M3 <- inner_join(ASC_VJlist_M3, ASC_top_pairs_M3, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

ASC_VJlist_M3 <- subset(ASC_VJlist_M3, select = -c(counts))

ASC_VJ_adj_M3 <- with(ASC_VJlist_M3, table(germline_v_call, germline_j_call))

chordDiagram(ASC_VJ_adj_M3, transparency = 0.5, annotationTrack = "grid",
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(ASC_VJ_adj_M3))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.8 )
}, bg.border = NA)

### mouse 4 ### 
ASC_top_pairs_M4 <- pair_counts_ASC_M4 %>% 
  slice_max(counts, n = 100) 

ASC_VJlist_M4 <- ASC_HC_M4 %>% 
  select(germline_v_call, germline_j_call)

ASC_VJlist_M4 <- inner_join(ASC_VJlist_M4, ASC_top_pairs_M4, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

ASC_VJlist_M4 <- subset(ASC_VJlist_M4, select = -c(counts))

ASC_VJ_adj_M4 <- with(ASC_VJlist_M4, table(germline_v_call, germline_j_call))

chordDiagram(ASC_VJ_adj_M4, transparency = 0.5, annotationTrack = "grid",
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(ASC_VJ_adj_M4))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.8 )
}, bg.border = NA)

## CD21lo 

CD21lo_top_pairs <- pair_counts_CD21lo %>% 
  slice_max(counts, n = 100) 

CD21lo_VJlist <- CD21lo_HC_all %>% 
  select(germline_v_call, germline_j_call)

CD21lo_VJlist <- inner_join(CD21lo_VJlist, CD21lo_top_pairs, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

CD21lo_VJlist <- subset(CD21lo_VJlist, select = -c(counts))

CD21lo_VJ_adj <- with(CD21lo_VJlist, table(germline_v_call, germline_j_call))

CD21lo_gridcol <- rep("grey", nrow(CD21lo_VJ_adj)+ncol(CD21lo_VJ_adj)) 

CD21lo_gridcol[c(18, 22, 24, 45, 46, 47, 48)] <- c("#FFA500", "#DB7093", "#87CEEB", "#963C78", "#B06D9A", "#CB9DBB", "#E5CEDD")

CD21lo_textcol <- rep("white", nrow(CD21lo_VJ_adj)+ncol(CD21lo_VJ_adj)) 

CD21lo_textcol[c(18, 22, 24, 45, 46, 47, 48)] <- rep("black", length(c(18, 22, 24, 45, 46, 47, 48)))

CD21lo_linkcol <- matrix("grey", nrow = nrow(CD21lo_VJ_adj), ncol = ncol(CD21lo_VJ_adj))

IGHV1_82 <- match("IGHV1-82*01", rownames(CD21lo_VJ_adj))
IGHV1_80 <- match("IGHV1-80*01", rownames(CD21lo_VJ_adj))
IGHV1_72 <- match("IGHV1-72*01", rownames(CD21lo_VJ_adj))

CD21lo_linkcol[IGHV1_72, ] <- "#FFA500"
CD21lo_linkcol[IGHV1_80, ] <- "#DB7093"
CD21lo_linkcol[IGHV1_82, ] <- "#87CEEB"

chordDiagram(CD21lo_VJ_adj, transparency = 0.5, annotationTrack = "grid", col = CD21lo_linkcol, grid.col = CD21lo_gridcol,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(CD21lo_VJ_adj))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  label <- CELL_META$sector.index
  if (label == "IGHV1-72*01" | label == "IGHV1-80*01" | label == "IGHV1-82*01" | label == "IGHJ1*03" | label == "IGHJ2*01" | label == "IGHJ3*01" | label == "IGHJ4*01") {
    col <- "black"
  } else {
    col <- "white"
  }
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 1.2, col = col)
}, bg.border = NA)

### mouse 1 ### 
CD21lo_top_pairs_M1 <- pair_counts_CD21lo_M1%>% 
  slice_max(counts, n = 100) 

CD21lo_VJlist_M1 <- CD21lo_HC_M1 %>% 
  select(germline_v_call, germline_j_call)

CD21lo_VJlist_M1 <- inner_join(CD21lo_VJlist_M1, CD21lo_top_pairs_M1, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

CD21lo_VJlist_M1 <- subset(CD21lo_VJlist_M1, select = -c(counts))

CD21lo_VJ_adj_M1 <- with(CD21lo_VJlist_M1, table(germline_v_call, germline_j_call))

chordDiagram(CD21lo_VJ_adj_M1, transparency = 0.5, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(CD21lo_VJ_adj_M1))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.7 )
}, bg.border = NA)

### mouse 2 ### 
CD21lo_top_pairs_M2 <- pair_counts_CD21lo_M2%>% 
  slice_max(counts, n = 100) 

CD21lo_VJlist_M2 <- CD21lo_HC_M2 %>% 
  select(germline_v_call, germline_j_call)

CD21lo_VJlist_M2 <- inner_join(CD21lo_VJlist_M2, CD21lo_top_pairs_M2, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

CD21lo_VJlist_M2 <- subset(CD21lo_VJlist_M2, select = -c(counts))

CD21lo_VJ_adj_M2 <- with(CD21lo_VJlist_M2, table(germline_v_call, germline_j_call))

chordDiagram(CD21lo_VJ_adj_M2, transparency = 0.5, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(CD21lo_VJ_adj_M2))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.7 )
}, bg.border = NA)

### mouse 3 ### 
CD21lo_top_pairs_M3 <- pair_counts_CD21lo_M3%>% 
  slice_max(counts, n = 100) 

CD21lo_VJlist_M3 <- CD21lo_HC_M3 %>% 
  select(germline_v_call, germline_j_call)

CD21lo_VJlist_M3 <- inner_join(CD21lo_VJlist_M3, CD21lo_top_pairs_M3, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

CD21lo_VJlist_M3 <- subset(CD21lo_VJlist_M3, select = -c(counts))

CD21lo_VJ_adj_M3 <- with(CD21lo_VJlist_M3, table(germline_v_call, germline_j_call))

chordDiagram(CD21lo_VJ_adj_M3, transparency = 0.5, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(CD21lo_VJ_adj_M3))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.7 )
}, bg.border = NA)

### mouse 4 ### 
CD21lo_top_pairs_M4 <- pair_counts_CD21lo_M4 %>% 
  slice_max(counts, n = 100) 

CD21lo_VJlist_M4 <- CD21lo_HC_M4 %>% 
  select(germline_v_call, germline_j_call)

CD21lo_VJlist_M4 <- inner_join(CD21lo_VJlist_M4, CD21lo_top_pairs_M4, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

CD21lo_VJlist_M4 <- subset(CD21lo_VJlist_M4, select = -c(counts))

CD21lo_VJ_adj_M4 <- with(CD21lo_VJlist_M4, table(germline_v_call, germline_j_call))

chordDiagram(CD21lo_VJ_adj_M4, transparency = 0.5, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(CD21lo_VJ_adj_M4))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.7 )
}, bg.border = NA)

## GC

GC_top_pairs <- pair_counts_GC %>% 
  slice_max(counts, n = 100) 

GC_VJlist <- GC_HC_all %>% 
  select(germline_v_call, germline_j_call)

GC_VJlist <- inner_join(GC_VJlist, GC_top_pairs, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

GC_VJlist <- subset(GC_VJlist, select = -c(counts))

GC_VJ_adj <- with(GC_VJlist, table(germline_v_call, germline_j_call))

GC_gridcol <- rep("grey", nrow(GC_VJ_adj)+ncol(GC_VJ_adj)) 

GC_gridcol[c(21, 24, 26, 40, 41, 42, 43)] <- c("#FFA500", "#DB7093", "#87CEEB", "#963C78", "#B06D9A", "#CB9DBB", "#E5CEDD")

GC_linkcol <- matrix("grey", nrow = nrow(GC_VJ_adj), ncol = ncol(GC_VJ_adj))

IGHV1_82 <- match("IGHV1-82*01", rownames(GC_VJ_adj))
IGHV1_80 <- match("IGHV1-80*01", rownames(GC_VJ_adj))
IGHV1_72 <- match("IGHV1-72*01", rownames(GC_VJ_adj))

GC_linkcol[IGHV1_72, ] <- "#FFA500"
GC_linkcol[IGHV1_80, ] <- "#DB7093"
GC_linkcol[IGHV1_82, ] <- "#87CEEB"

chordDiagram(GC_VJ_adj, transparency = 0.5, annotationTrack = "grid", col = GC_linkcol, grid.col = GC_gridcol,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(GC_VJ_adj))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  label <- CELL_META$sector.index
  if (label == "IGHV1-72*01" | label == "IGHV1-80*01" | label == "IGHV1-82*01" | label == "IGHJ1*03" | label == "IGHJ2*01" | label == "IGHJ3*01" | label == "IGHJ4*01") {
    col <- "black"
  } else {
    col <- "white"
  }
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 1.2, col = col)
}, bg.border = NA)

### mouse 1 ### 
GC_top_pairs_M1 <- pair_counts_GC_M1 %>% 
  slice_max(counts, n = 100) 

GC_VJlist_M1 <- GC_HC_M1 %>% 
  select(germline_v_call, germline_j_call)

GC_VJlist_M1 <- inner_join(GC_VJlist_M1, GC_top_pairs_M1, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

GC_VJlist_M1 <- subset(GC_VJlist_M1, select = -c(counts))

GC_VJ_adj_M1 <- with(GC_VJlist_M1, table(germline_v_call, germline_j_call))

chordDiagram(GC_VJ_adj_M1, transparency = 0.5, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(GC_VJ_adj_M1))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.8 )
}, bg.border = NA)

### mouse 2 ### 
GC_top_pairs_M2 <- pair_counts_GC_M2 %>% 
  slice_max(counts, n = 100) 

GC_VJlist_M2 <- GC_HC_M2 %>% 
  select(germline_v_call, germline_j_call)

GC_VJlist_M2 <- inner_join(GC_VJlist_M2, GC_top_pairs_M2, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

GC_VJlist_M2 <- subset(GC_VJlist_M2, select = -c(counts))

GC_VJ_adj_M2 <- with(GC_VJlist_M2, table(germline_v_call, germline_j_call))

chordDiagram(GC_VJ_adj_M2, transparency = 0.5, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(GC_VJ_adj_M2))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.8 )
}, bg.border = NA)

### mouse 3 ### 
GC_top_pairs_M3 <- pair_counts_GC_M3 %>% 
  slice_max(counts, n = 100) 

GC_VJlist_M3 <- GC_HC_M3 %>% 
  select(germline_v_call, germline_j_call)

GC_VJlist_M3 <- inner_join(GC_VJlist_M3, GC_top_pairs_M3, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

GC_VJlist_M3 <- subset(GC_VJlist_M3, select = -c(counts))

GC_VJ_adj_M3 <- with(GC_VJlist_M3, table(germline_v_call, germline_j_call))

chordDiagram(GC_VJ_adj_M3, transparency = 0.5, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(GC_VJ_adj_M3))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.8 )
}, bg.border = NA)

### mouse 4 ### 
GC_top_pairs_M4 <- pair_counts_GC_M4 %>% 
  slice_max(counts, n = 100) 

GC_VJlist_M4 <- GC_HC_M4 %>% 
  select(germline_v_call, germline_j_call)

GC_VJlist_M4 <- inner_join(GC_VJlist_M4, GC_top_pairs_M4, by = c("germline_v_call" = "germline_v_call", "germline_j_call" = "germline_j_call"))

GC_VJlist_M4 <- subset(GC_VJlist_M4, select = -c(counts))

GC_VJ_adj_M4 <- with(GC_VJlist_M4, table(germline_v_call, germline_j_call))

chordDiagram(GC_VJ_adj_M4, transparency = 0.5, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(GC_VJ_adj_M4))))))

circos.track(track.index = 1, panel.fun = function(x,y){
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.8 )
}, bg.border = NA)


