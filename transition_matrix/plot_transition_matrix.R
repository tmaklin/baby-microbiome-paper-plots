
source("data_wrangling.R") ## Functions to build the transition matrix
source("plotting.R") ## Functions to plot the transition matrix

## Read in the sample metadata and filter accessions to exclude Mother/Infancy samples
metadata <- read.table("../wgs_meta_delivery.tsv", header=TRUE)
vaginal.accessions <- metadata$err_accession[!grepl("Mother|Infancy", metadata$Time_point) & metadata$Delivery_mode == "Vaginal"]
caesarean.accessions <- metadata$err_accession[!grepl("Mother|Infancy", metadata$Time_point) & metadata$Delivery_mode == "Caesarean"]

## Read in the results from demix_check
demix.results <- read.table("../ecoli-new-reference/E_col_demix_results_top2.tsv", sep='\t', header=FALSE, stringsAsFactors=FALSE)
##demix.results <- demix.results[demix.results[, 3] > 0.01, ]

## Extract cohorts
vaginal <- demix.results[demix.results$V1 %in% vaginal.accessions, ]
caesarean <- demix.results[demix.results$V1 %in% caesarean.accessions, ]

## Build the transition matrix for each cohort
vaginal.mat <- GetTransitionMatrix(vaginal, metadata, "E_col")
caesarean.mat <- GetTransitionMatrix(caesarean, metadata, "E_col")

par(mar = c(1, 23, 23, 1))
PlotTransitionMatrix(caesarean.mat)

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
