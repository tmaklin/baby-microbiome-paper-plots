
ReadAbundances <- function(files, species, genome.len, coverage) {
    ## Read in species abundances
    species.files <- files[grepl(species, files)]

    abundances <- read.table(species.files[1], sep = '\t', comment.char = '#')
    n.clusters <- nrow(abundances)
    n.files <- length(species.files)

    abundances.mat <- matrix(0, nrow = n.files, ncol = n.clusters)
    alncounts <- numeric(n.files)
    identifiers <- gsub(".*[/]", "", gsub("_abundances.txt", "", species.files))
    rownames(abundances.mat) <- identifiers
    colnames(abundances.mat) <- abundances[, 1]
    names(alncounts) <- identifiers

    for (i in 1:n.files) {
        alncounts[i] <- as.numeric(read.table(species.files[i], comment.char = '@')[2, 2])
        abundances.mat[i, ] <- read.table(species.files[i], sep = '\t', comment.char = '#')[, 2]
    }

    list("abundances" = abundances.mat, "aln_counts" = alncounts)
}

FilterByCoverage <- function(species.data, strain.data, species.name, genome.len, min.coverage, read.len = 150) {
    ## Filters the output from ReadAbundances for species & strains by
    ## n.aligned.reads*read.len > genome.length*min.cov
    min.bases <- genome.len*min.coverage

    n.strain.samples <- length(strain.data$aln_counts)
    pass.or.fail <- numeric(n.strain.samples)
    names(pass.or.fail) <- rownames(strain.data$abundances)
    read.counts <- numeric(n.strain.samples)
    names(read.counts) <- rownames(strain.data$abundances)
    abus <- numeric(n.strain.samples)
    names(abus) <- rownames(strain.data$abundances)

    species.pos <- which(grepl(species.name, colnames(species.data$abundances)))
    for (i in 1:n.strain.samples) {
        accession <- sub(".*(ERR[0-9]*).*", "\\1", rownames(strain.data$abundances)[i])
        pos <- match(accession, rownames(species.data$abundances))
        read.count <- species.data$aln_counts[pos]
        species.read.count <- floor(read.count*species.data$abundances[pos, species.pos])
        species.pass.qc <- species.read.count*read.len > min.bases
        pass.or.fail[i] <- species.pass.qc
        read.counts[i] <- species.read.count
        abus[i] <- species.data$abundances[pos, species.pos]
    }

    list("filter" = pass.or.fail, "aligned_reads" = read.counts, "abundances" = abus)
}


RenameClusters <- function(old.labels, new.labels) {
    ## Find new names for the demix.check.data clusters
    order.in.labels <- match(old.labels, new.labels$PopPUNK)
    ordered.labels <- new.labels[order.in.labels, ]
    stopifnot(all(ordered.labels$PopPUNK == old.labels) && length(ordered.labels) > 0)

    ## Rename from [Species][SC]_[ST] to [Species]_[ST]_[SC]
    renamed.clusters <- sub("(SC[0-9]*)_(ST[0-9]*)", "\\2_\\1", ordered.labels$Cluster)

    renamed.clusters
}

