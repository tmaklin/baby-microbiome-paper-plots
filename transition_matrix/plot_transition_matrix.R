
source("data_wrangling.R") ## Functions to build the transition matrix
source("plotting.R") ## Functions to plot the transition matrix

## Read in the sample metadata and filter accessions to exclude Mother/Infancy samples
metadata <- read.table("../wgs_meta_delivery.tsv", header=TRUE)
##vaginal.accessions <- metadata$err_accession[!grepl("Mother|Infancy", metadata$Time_point) & metadata$Delivery_mode == "Vaginal"]
##caesarean.accessions <- metadata$err_accession[!grepl("Mother|Infancy", metadata$Time_point) & metadata$Delivery_mode == "Caesarean"]
##vaginal.accessions <- metadata$err_accession[!grepl("Mother", metadata$Time_point) & metadata$Delivery_mode == "Vaginal"]
##caesarean.accessions <- metadata$err_accession[!grepl("Mother", metadata$Time_point) & metadata$Delivery_mode == "Caesarean"]
##vaginal.accessions <- metadata$err_accession[!grepl("Infancy", metadata$Time_point) & metadata$Delivery_mode == "Vaginal"]
##caesarean.accessions <- metadata$err_accession[!grepl("Infancy", metadata$Time_point) & metadata$Delivery_mode == "Caesarean"]
vaginal.accessions <- metadata$err_accession[metadata$Delivery_mode == "Vaginal"]
caesarean.accessions <- metadata$err_accession[!metadata$Delivery_mode == "Caesarean"]

## Read in the results from demix_check
demix.results <- read.table("../ecoli-new-reference/E_col_demix_results_top2.tsv", sep='\t', header=FALSE, stringsAsFactors=FALSE)

## Only use samples that have a relative abundance highe than 1%
demix.results <- demix.results[demix.results[, 3] > 0.01, ]


## Extract cohorts
vaginal <- demix.results[demix.results$V1 %in% vaginal.accessions, ]
caesarean <- demix.results[demix.results$V1 %in% caesarean.accessions, ]

## Rename the clusters from E_col_Pop[0-999] to E_col_ST[0-9999]_SC[0-999]
new.clusters <- read.table(file = "../ecoli-new-reference/new_clusters.tsv", sep='\t', header = TRUE)
vaginal$V2 <- paste("E_col_", new.clusters$Cluster[match(vaginal$V2, new.clusters$PopPUNK)], sep='')
caesarean$V2 <- paste("E_col_", new.clusters$Cluster[match(caesarean$V2, new.clusters$PopPUNK)], sep='')

## Build the transition matrix for each cohort
vaginal.mat <- GetTransitionMatrix(vaginal, metadata, "E_col")
caesarean.mat <- GetTransitionMatrix(caesarean, metadata, "E_col")

gradient <- colorRampPalette(c("#feebe2", "#fbb4b9", "#f768a1", "#c51b8a", "#7a0177"))

pdf(file = "E_coli_transition_matrix_no_mother.pdf", width = 12, height = 6)
n.max <- max(max(vaginal.mat), max(caesarean.mat))
par(mar = c(1, 24, 20, 4))
layout(matrix(1:3, ncol = 3), widths = c(2, 1, 0.3), heights = c(1))
more.than.once <- which(rowSums(vaginal.mat) > 1 | colSums(vaginal.mat) > 1)
PlotTransitionMatrix(vaginal.mat[more.than.once, more.than.once], "a)", 0.15, -8, 23.3, n.max, gradient)
par(mar = c(1, 8, 20, 4))
more.than.once <- which(rowSums(caesarean.mat) > 1 | colSums(caesarean.mat) > 1)
PlotTransitionMatrix(caesarean.mat[more.than.once, more.than.once], "b)", 0.60, -8, 9.9, n.max, gradient)
par(mar = c(1, 0, 20, 4))
IntensityLegend(gradient, n.max)
dev.off()


##pdf(file = "E_coli_transition_matrix_caesarean_cohort.pdf", width = 18, height = 18)
##dev.off()

###

library("vegan")

hospital.b <- grepl("^B", rownames(subset))
hospital.c <- grepl("^C", rownames(subset))
hospital.a <- !hospital.b & !hospital.c

hospital.b <- grepl("^B", rownames(mt.subset))
hospital.c <- grepl("^C", rownames(mt.subset))
hospital.a <- !hospital.b & !hospital.c

cae.shannon <- diversity(mt.subset, index = "shannon")
cae.a.shannon <- diversity(mt.subset[hospital.a, ], index = "shannon")
cae.b.shannon <- diversity(mt.subset[hospital.b, ], index = "shannon")
cae.c.shannon <- diversity(mt.subset[hospital.c, ], index = "shannon")

wilcox.test(cae.shannon, vag.shannon, paired = FALSE, alternative = "less")
wilcox.test(cae.a.shannon, vag.a.shannon, paired = FALSE, alternative = "less")
wilcox.test(cae.c.shannon, vag.b.shannon, paired = FALSE, alternative = "less")
wilcox.test(cae.b.shannon, vag.c.shannon, paired = FALSE, alternative = "less")

wilcox.test(vag.shannon, vag.a.shannon, paired = FALSE)
wilcox.test(vag.shannon, vag.b.shannon, paired = FALSE)
wilcox.test(vag.shannon, vag.c.shannon, paired = FALSE)

wilcox.test(vag.a.shannon, vag.b.shannon, paired = FALSE)
wilcox.test(vag.a.shannon, vag.c.shannon, paired = FALSE)
wilcox.test(vag.b.shannon, vag.c.shannon, paired = FALSE)

wilcox.test(cae.shannon, cae.a.shannon, paired = FALSE)
wilcox.test(cae.shannon, cae.b.shannon, paired = FALSE)
wilcox.test(cae.shannon, cae.c.shannon, paired = FALSE)

wilcox.test(cae.a.shannon, cae.b.shannon, paired = FALSE)
wilcox.test(cae.a.shannon, cae.c.shannon, paired = FALSE)
wilcox.test(cae.b.shannon, cae.c.shannon, paired = FALSE)
