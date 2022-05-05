
library("UpSetR")

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

SetupUpsetData <- function(data) {
    by.accession <- by(data$cluster, data$accession, unique)
    by.accession <- lapply(by.accession, function(x) sub("^K. ([a-z]*).*", "K. \\1", x))
##    by.accession <- lapply(by.accession, function(x) sub("^E. ([a-z]*).*", "E. \\1", x))

    all.identified <- unique(unlist(by.accession))

    n.samples <- length(by.accession)
    n.lineages <- length(all.identified)

    presabs <- matrix(FALSE, nrow = n.samples, ncol = n.lineages)
    rownames(presabs) <- names(by.accession)
    colnames(presabs) <- sort(all.identified) ## TODO: sort by species and ST

    for (i in 1:n.samples) {
        sample.name <- names(by.accession)[i]
        identified <- by.accession[[i]]
        presabs[sample.name, identified] <- TRUE
    }

    filtered.presabs <- presabs[, colSums(presabs) >= 5]

    by.cluster <- apply(filtered.presabs, 2, which)

    list("presabs" = filtered.presabs, "by_cluster" = by.cluster)
}


vaginal.counts <- list("no_mother_no_infancy" = ReadDemixCheckData("../wgs_meta_delivery.tsv",
                                                                   "../wgs_demix_check_high_confidence.tsv",
                                                                   "../poppunk_cluster_info/all_new_clusters.tsv",
                                                                   "Vaginal", "Mother|Infancy", "E_col|K_[a-z]*"),
                       "no_mother" = ReadDemixCheckData("../wgs_meta_delivery.tsv",
                                                        "../wgs_demix_check_high_confidence.tsv",
                                                        "../poppunk_cluster_info/all_new_clusters.tsv",
                                                        "Vaginal", "Mother", "E_col|K_[a-z]*"))


caesarean.counts <- list("no_mother_no_infancy" = ReadDemixCheckData("../wgs_meta_delivery.tsv",
                                                                   "../wgs_demix_check_high_confidence.tsv",
                                                                   "../poppunk_cluster_info/all_new_clusters.tsv",
                                                                   "Caesarean", "Mother|Infancy", "E_col|K_[a-z]*"),
                       "no_mother" = ReadDemixCheckData("../wgs_meta_delivery.tsv",
                                                        "../wgs_demix_check_high_confidence.tsv",
                                                        "../poppunk_cluster_info/all_new_clusters.tsv",
                                                        "Caesarean", "Mother", "E_col|K_[a-z]*"))



pdf(file = "upset_vaginal_no_mother.pdf", width = 8, height = 4)
vnm <- SetupUpsetData(vaginal.counts$no_mother)
upset(fromList(vnm$by_cluster), order.by = "freq", nsets = length(vnm$by_cluster), mb.ratio = c(0.3, 0.7), nintersects = 200)
dev.off()
pdf(file = "upset_caesarean_no_mother.pdf", width = 8, height = 4)
vnm <- SetupUpsetData(caesarean.counts$no_mother)
upset(fromList(vnm$by_cluster), order.by = "freq", nsets = length(vnm$by_cluster), mb.ratio = c(0.3, 0.7), nintersects = 200)
dev.off()
