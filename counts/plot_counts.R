
EcolPrettyLabels <- function(cluster.labels) {
    ## Rename ST10 subclusters
    cluster.labels <- gsub("E. col ST10 SC6", "E. col ST10-1", cluster.labels)
    cluster.labels <- gsub("E. col ST10 SC351", "E. col ST10-2", cluster.labels)
    cluster.labels <- gsub("E. col ST10 SC79", "E. col ST10-3", cluster.labels)
    cluster.labels
}

RenameClusters <- function(data, new.cluster.names.path) {
    ## ## ## Rename the clusters from [Species]_Pop[0-999] to [Species]_ST[0-9999]_SC[0-999]
    new.clusters <- read.table(file = new.cluster.names.path, sep='\t', header = TRUE)
    has.new.names <- data$cluster %in% new.clusters$PopPUNK
    data$cluster[has.new.names] <- new.clusters$Cluster[match(data$cluster[has.new.names], new.clusters$PopPUNK)]
    data$cluster[has.new.names] <- sub("SC([0-9]*)_ST([0-9]*)", "ST\\2_SC\\1", data$cluster[has.new.names])
    data$cluster <- sub("^([A-Z])(_)([a-z]*)", "\\1. \\3", data$cluster)
    data$cluster <- gsub("_", " ", data$cluster)
    data
}

ReadDemixCheckData <- function(metadata.path, demix.path, new.cluster.names.path, cohort.name, exclude.time.points, species.name) {
    metadata <- read.table(metadata.path, header=TRUE)

    ## Read in the results from demix_check
    demix.results <- read.table(demix.path, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    ## Filter by species
    demix.results <- demix.results[grepl(species.name, demix.results$cluster), ]
    ## Only use samples that have a relative abundance higher than 1%
    demix.results <- demix.results[demix.results$abundance > 0.01, ]

    ## Extract cohort and remove excluded time points
    accessions <- metadata$err_accession[!grepl(exclude.time.points, metadata$Time_point) & metadata$Delivery_mode == cohort.name]
    samples <- demix.results[demix.results$accession %in% accessions, ]

    ## Rename the clusters
    samples <- RenameClusters(samples, new.cluster.names.path)
    samples$cluster <- EcolPrettyLabels(samples$cluster)
##    samples$cluster <- gsub(" SC[0-9]*$", "", samples$cluster)

    samples
}

vaginal.counts <- list("no_mother_no_infancy" = ReadDemixCheckData("../wgs_meta_delivery.tsv",
                                                                   "../wgs_demix_check_high_confidence.tsv",
                                                                   "../poppunk_cluster_info/all_new_clusters.tsv",
                                                                   "Vaginal", "Mother|Infancy", "E_col|K_[a-z]*|S_aur|E_fcs|E_fcm*"),
                       "no_mother" = ReadDemixCheckData("../wgs_meta_delivery.tsv",
                                                        "../wgs_demix_check_high_confidence.tsv",
                                                        "../poppunk_cluster_info/all_new_clusters.tsv",
                                                        "Vaginal", "Mother", "E_col|K_[a-z]*|S_aur|E_fcs|E_fcm*"))


caesarean.counts <- list("no_mother_no_infancy" = ReadDemixCheckData("../wgs_meta_delivery.tsv",
                                                                   "../wgs_demix_check_high_confidence.tsv",
                                                                   "../poppunk_cluster_info/all_new_clusters.tsv",
                                                                   "Caesarean", "Mother|Infancy", "E_col|K_[a-z]*|S_aur|E_fcs|E_fcm*"),
                       "no_mother" = ReadDemixCheckData("../wgs_meta_delivery.tsv",
                                                        "../wgs_demix_check_high_confidence.tsv",
                                                        "../poppunk_cluster_info/all_new_clusters.tsv",
                                                        "Caesarean", "Mother", "E_col|K_[a-z]*|S_aur|E_fcs|E_fcm*"))

vaginal.no_mother.counts <- table(gsub("([A-Z][.][ ][a-z][a-z][a-z]).*$", "\\1", vaginal.counts$no_mother$cluster))
caesarean.no_mother.counts <- table(gsub("([A-Z][.][ ][a-z][a-z][a-z]).*$", "\\1", caesarean.counts$no_mother$cluster))

clusters <- sort(unique(c(names(vaginal.no_mother.counts), names(caesarean.no_mother.counts))))
n.clusters <- length(clusters)

obs.mat <- matrix(0, nrow = 2, ncol = n.clusters)
colnames(obs.mat) <- clusters
rownames(obs.mat) <- c("Vaginal delivery", "Caesarean delivery")
for (i in 1:n.clusters) {
    if (clusters[i] %in% names(vaginal.no_mother.counts)) {
        obs.mat[1, clusters[i]] <- vaginal.no_mother.counts[clusters[i]]
    }
    if (clusters[i] %in% names(caesarean.no_mother.counts)) {
        obs.mat[2, clusters[i]] <- caesarean.no_mother.counts[clusters[i]]
    }
}

diffs <- obs.mat[1, ] - obs.mat[2, ]
pdf(file = "counts_differences.pdf", width = 6, height = 6)
layout(matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(0.85, 0.15))
par(mar = c(5, 4, 0.5, 1))
plot(y = n.clusters:1, x = diffs/colSums(obs.mat), xaxt = 'n', yaxt = 'n', ylab = '', xlab = "Relative difference", cex = log(abs(diffs) + 1, base = 5), col = ifelse(diffs < 0, "#92c5de", "#f4a582"), pch = 19, bty = 'n', cex.lab = 1.2)
abline(v = 0, col = "gray80", lty = "dashed")
axis(side = 2, at = n.clusters:1, labels = clusters, las = 2, font = 3, cex.axis = 1.2)
x.axis.labels <- c(paste('-', seq(100, 25, by = -25), '%', sep=''), "0%", paste('+', seq(25, 100, by = 25), '%', sep=''))
axis(side = 1, at = seq(-1.0, 1.0, by = 0.25), labels = x.axis.labels, cex.axis = 1.25)
par(mar = c(0, 0, 0, 0))
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
legend("topleft", legend = c(1, 25, 50, 100, 250), pt.cex = log(c(1, 25, 50, 100, 250) + 1, base = 5), pch = rep(19, 5), cex = 1, ncol = 5, bty = 'n')
title("Absolute difference", font.main = 1, cex.main = 1.2, adj = 0.15, line = -3.75)
legend("topright", legend = c("Caesarean section cohort", "Vaginal delivery cohort"), pch = rep(19, 2), cex = 1, ncol = 1, bty = 'n', col = c("#92c5de", "#f4a582"))
title("Enrichened in", font.main = 1, cex.main = 1.2, adj = 0.85, line = -3.75)
dev.off()
