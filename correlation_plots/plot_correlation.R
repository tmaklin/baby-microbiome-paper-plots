
source("plotting.R")

library("RColorBrewer")

metadata <- read.table("../wgs_meta_delivery.tsv", header=TRUE)
vaginal.accessions <- metadata$err_accession[!grepl("Mother|Infancy", metadata$Time_point) & metadata$Delivery_mode == "Vaginal"]
caesarean.accessions <- metadata$err_accession[!grepl("Mother|Infancy", metadata$Time_point) & metadata$Delivery_mode == "Caesarean"]

demix.filter <- ReadDemixResults()

OtuFilter <- function(x) { grepl("[EKS][_][a-z][a-z][a-z]", x) & x != "S_pne" }

vaginal.correlations <- ReadCorrelations("data/vaginal_correlation.tsv", "data/vaginal_pvalues.tsv", OtuFilter)
caesarean.correlations <- ReadCorrelations("data/caesarean_correlation.tsv", "data/caesarean_pvalues.tsv", OtuFilter)

vaginal.correlations <- vaginal.correlations[OtuFilter(rownames(vaginal.correlations)), OtuFilter(rownames(vaginal.correlations))]
caesarean.correlations <- caesarean.correlations[OtuFilter(rownames(caesarean.correlations)), OtuFilter(rownames(caesarean.correlations))]

vaginal.obs.counts <- DemixCheckCounts(demix.filter, vaginal.accessions, vaginal.correlations)
caesarean.obs.counts <- DemixCheckCounts(demix.filter, caesarean.accessions, caesarean.correlations)

gradient <- colorRampPalette(rev(c("#ca0020", "#f4a582", "white", "#92c5de", "#0571b0")))

pdf(file = "out/correlation_plots.pdf", width = 8, height = 4)
par(mar = c(4, 4, 4, 2))
layout(matrix(1:3, ncol = 3), widths = c(1, 1, 0.3), heights = c(1, 0.2))
CorrelationPlot(vaginal.correlations, vaginal.obs.counts, gsub("_", ". ", rownames(vaginal.correlations)), "a)", 0.01, gradient)
par(mar = c(4, 4, 4, 2))
CorrelationPlot(caesarean.correlations, caesarean.obs.counts, gsub("_", ". ", rownames(caesarean.correlations)), "b)", 0.46, gradient)
par(mar = c(4, 0, 4, 0))
CorrelationIntensityLegend(gradient)
par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
PointAreaLegend()
dev.off()
