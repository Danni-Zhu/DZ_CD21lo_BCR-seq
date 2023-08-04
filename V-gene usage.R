#Immcantation - count V gene usage frequency for post-hoc analyses

library("alakazam")
library("shazam")
library("dplyr")
library("scales")

# Import dataset including duplicate count
HCtsv_DUPCOUNT <- read_tsv("Data/ChangeO HC_LC tsv_deduplicated_filtered/Reprocessed_HC_Combined_filtered_dist_ham_clone-pass_germ-pass_DUPCOUNT.tsv")

# subset to HC table by mouse 
## Mouse 1
M1_HC_all <- HCtsv_DUPCOUNT[HCtsv_DUPCOUNT$mouse == "M1", ]
## Mouse 2
M2_HC_all <- HCtsv_DUPCOUNT[HCtsv_DUPCOUNT$mouse == "M2", ]
## Mouse 3
M3_HC_all <- HCtsv_DUPCOUNT[HCtsv_DUPCOUNT$mouse == "M3", ]
## Mouse 4
M4_HC_all <- HCtsv_DUPCOUNT[HCtsv_DUPCOUNT$mouse == "M4", ]
## Input follicular B cells 
Day0FOB_HC_all <- HCtsv_DUPCOUNT[HCtsv_DUPCOUNT$mouse == "M0", ]

# Count V gene frequencies for each mouse and convert to tables
## Input follicular B cells 
Day0_gene <- countGenes(Day0FOB_HC_all, gene="v_call", groups="subset", mode="gene")
write.table(Day0_gene, file="Day 0 Vh usage.tsv", quote=FALSE, sep="\t", row.names = FALSE)

## Mouse 1
M1gene <- countGenes(M1_HC_all, gene="v_call", groups="subset", mode="gene")
write.table(M1gene, file="Mouse 1 Vh usage.tsv", quote=FALSE, sep="\t", row.names = FALSE)

## Mouse 2
M2gene <- countGenes(M2_HC_all, gene="v_call", groups="subset", mode="gene")
write.table(M2gene, file="Mouse 2 Vh usage.tsv", quote=FALSE, sep="\t", row.names = FALSE)

## Mouse 3
M3gene <- countGenes(M3_HC_all, gene="v_call", groups="subset", mode="gene")
write.table(M3gene, file="Mouse 3 Vh usage.tsv", quote=FALSE, sep="\t", row.names = FALSE)

## Mouse 4
M4gene <- countGenes(M4_HC_all, gene="v_call", groups="subset", mode="gene")
write.table(M4gene, file="Mouse 4 Vh usage.tsv", quote=FALSE, sep="\t", row.names = FALSE)







