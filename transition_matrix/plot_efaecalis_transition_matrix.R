
source("data_wrangling.R") ## Functions to build the transition matrix
source("plotting.R") ## Functions to plot the transition matrix
source("fix_cluster_labels.R") ## Functions to change the cluster labels to aesthetically pleasing format

vaginal.mats <- list("no_mother_no_infancy" = ReadTransitionMatrix("../wgs_meta_delivery.tsv",
                                                                   "../wgs_demix_check_high_confidence.tsv",
                                                                   "../poppunk_cluster_info/E_fcs_new_clusters.tsv",
                                                                   "Vaginal", "Mother|Infancy", "E_fcs"),
                     "no_mother" = ReadTransitionMatrix("../wgs_meta_delivery.tsv",
                                                        "../wgs_demix_check_high_confidence.tsv",
                                                        "../poppunk_cluster_info/E_fcs_new_clusters.tsv",
                                                        "Vaginal", "Mother", "E_fcs"))

caesarean.mats <- list("no_mother_no_infancy" = ReadTransitionMatrix("../wgs_meta_delivery.tsv",
                                                                     "../wgs_demix_check_high_confidence.tsv",
                                                                     "../poppunk_cluster_info/E_fcs_new_clusters.tsv",
                                                                     "Caesarean", "Mother|Infancy", "E_fcs"),
                       "no_mother" = ReadTransitionMatrix("../wgs_meta_delivery.tsv",
                                                          "../wgs_demix_check_high_confidence.tsv",
                                                          "../poppunk_cluster_info/E_fcs_new_clusters.tsv",
                                                          "Caesarean", "Mother", "E_fcs"))

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

pdf(file = "E_faecalis_transition_matrix_first_month.pdf", width = 9, height = 5)
layout(matrix(c(1, 2, 3), ncol = 3, byrow = TRUE), widths = c(1.7, 2, 0.35), heights = c(0.65, 0.9))
par(mar = c(5, 8, 8, 2))
PlotTransitionMatrix(vaginal.mats$no_mother_no_infancy, "a)", 0.07, -5, nrow(vaginal.mats$no_mother_no_infancy) + 1, n.max, 0.25, gradient)
par(mar = c(0, 6, 8, 1))
PlotTransitionMatrix(caesarean.mats$no_mother_no_infancy, "b)", 0.5, -5, nrow(caesarean.mats$no_mother_no_infancy) + 1.25, n.max, 0.4, gradient)
## par(mar = c(10, 8, 8, 2))
## PlotTransitionMatrix(vaginal.mats$no_mother, "c)", 0.07, -30.75, nrow(vaginal.mats$no_mother) + 1, n.max, 0.4, gradient)
## par(mar = c(6, 6, 8, 1))
## PlotTransitionMatrix(caesarean.mats$no_mother, "d)", 0.565, -30.75, nrow(caesarean.mats$no_mother) + 2, n.max, 0.4, gradient)
par(mar = c(1, 0, 8, 4))
IntensityLegend(gradient, 18, 1.5)
## Optinally add some titles as well
## title(main = expression(italic("E. faecalis")~"competition"), line = -57, outer = TRUE, adj = 0.5)
## title(main = "Vaginal delivery cohort", line = -1, outer = TRUE, adj = 0.25)
## mtext(side = 2, "Mother removed", line = -1.2, outer = TRUE, adj = 0.25, cex = 0.8)
## mtext(side = 2, "Mother & Infancy removed", line = -1.2, outer = TRUE, adj = 0.8, cex = 0.8)
## title(main = "Caesarean delivery cohort", line = -1, outer = TRUE, adj = 0.75)
dev.off()
