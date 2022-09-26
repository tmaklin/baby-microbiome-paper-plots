
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

a <- match(rownames(caesarean.correlations), rownames(vaginal.correlations))
b <- match(colnames(caesarean.correlations), colnames(vaginal.correlations))

tmp <- as.data.frame(as.matrix(vaginal.correlations)[a, b])
rownames(tmp) <- rownames(caesarean.correlations)
colnames(tmp) <- colnames(caesarean.correlations)
tmp[lower.tri(tmp)] <- caesarean.correlations[lower.tri(caesarean.correlations)]
tmp[is.na(tmp)] <- 0

tmp.counts <- as.data.frame(as.matrix(vaginal.obs.counts)[a, b])
rownames(tmp.counts) <- rownames(caesarean.obs.counts)
colnames(tmp.counts) <- colnames(caesarean.obs.counts)
tmp.counts[lower.tri(tmp.counts)] <- caesarean.obs.counts[lower.tri(caesarean.obs.counts)]
tmp.counts[is.na(tmp.counts)] <- 0

pdf(file = "correlation_plot_revised.pdf", width = 8, height = 8)
par(mar = c(4, 4, 4, 2))
layout(matrix(1:2, ncol = 2), widths = c(1, 0.3), heights = c(1, 0.2))
CorrelationPlot(tmp, tmp.counts, gsub("_", ". ", rownames(tmp)), "", 0.01, gradient)
par(mar = c(4, 0, 4, 0))
CorrelationIntensityLegend(gradient, -33.5)
par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
PointAreaLegend()
dev.off()
