
source("data_wrangling.R") ## Functions to build the transition matrix
source("plotting.R") ## Functions to plot the transition matrix
source("fix_cluster_labels.R") ## Functions to change the cluster labels to aesthetically pleasing format

vaginal.mats <- list("no_mother_no_infancy" = ReadTransitionMatrix("../wgs_meta_delivery.tsv",
                                                                   "../ecoli-new-reference/E_col_demix_results_top2.tsv",
                                                                   "../ecoli-new-reference/new_clusters.tsv",
                                                                   "Vaginal", "Mother|Infancy", "E_col"),
                     "no_mother" = ReadTransitionMatrix("../wgs_meta_delivery.tsv",
                                                        "../ecoli-new-reference/E_col_demix_results_top2.tsv",
                                                        "../ecoli-new-reference/new_clusters.tsv",
                                                        "Vaginal", "Mother", "E_col"))
caesarean.mats <- list("no_mother_no_infancy" = ReadTransitionMatrix("../wgs_meta_delivery.tsv",
                                                                     "../ecoli-new-reference/E_col_demix_results_top2.tsv",
                                                                     "../ecoli-new-reference/new_clusters.tsv",
                                                                     "Caesarean", "Mother|Infancy", "E_col"),
                       "no_mother" = ReadTransitionMatrix("../wgs_meta_delivery.tsv",
                                                          "../ecoli-new-reference/E_col_demix_results_top2.tsv",
                                                          "../ecoli-new-reference/new_clusters.tsv",
                                                          "Caesarean", "Mother", "E_col"))

## Change the cluster labels to what we want in the plot
vaginal.mats <- lapply(vaginal.mats, function(x) { rownames(x) <- colnames(x) <- PrettyClusterLabels(rownames(x), EcolPrettyLabels); return(x) })
caesarean.mats <- lapply(caesarean.mats, function(x) { rownames(x) <- colnames(x) <- PrettyClusterLabels(rownames(x), EcolPrettyLabels); return(x) })

## Reorder the matrices so they ascend in ST code
vaginal.mats <- lapply(vaginal.mats, function(x) { new.order <- order(as.numeric(gsub("[- ].*$", "", gsub("ST", "", rownames(x))))); return(x[new.order, new.order]) })
caesarean.mats <- lapply(caesarean.mats, function(x) { new.order <- order(as.numeric(gsub("[- ].*$", "", gsub("ST", "", rownames(x))))); return(x[new.order, new.order]) })

## Single hue color for intensity
gradient <- colorRampPalette(c("#feebe2", "#fbb4b9", "#f768a1", "#c51b8a", "#7a0177"))

## Get the maximum number of transitions across all samples
n.max <- max(unlist(lapply(vaginal.mats, max)), unlist(lapply(caesarean.mats, max)))

pdf(file = "E_coli_transition_matrix_no_mother.pdf", width = 9, height = 8)
layout(matrix(c(1, 2, 5, 3, 4, 5), ncol = 3, byrow = TRUE), widths = c(2, 1.5, 0.3), heights = c(0.65, 0.9))
par(mar = c(0, 8, 8, 0.5))
PlotTransitionMatrix(vaginal.mats$no_mother_no_infancy, "a)", 0.07, -5, 18, n.max, 0.3, gradient)
par(mar = c(7.5, 6, 8, 1))
PlotTransitionMatrix(caesarean.mats$no_mother_no_infancy, "b)", 0.565, -5, 10.6, n.max, 0.25, gradient)
par(mar = c(0, 8, 8, 0.5))
PlotTransitionMatrix(vaginal.mats$no_mother, "c)", 0.07, -30.75, 30, n.max, 0.5, gradient)
par(mar = c(11, 6, 8, 1))
PlotTransitionMatrix(caesarean.mats$no_mother, "d)", 0.565, -30.75, 16.4, n.max, 0.4, gradient)
par(mar = c(1, 0, 8, 4))
IntensityLegend(gradient, n.max)
dev.off()

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
