
## Read in the E. coli demix check data that has abundance > 0.01
demix.check.data <- read.table("../wgs_demix_check_high_confidence.tsv", header = TRUE, comment.char='@', sep = '\t')
e.colis <- demix.check.data[grepl("^E_col", demix.check.data$cluster), ]
e.colis <- e.colis[e.colis$abundance > 0.01, ]

## Read in the metadata (accession, delivery_mode, individual, timepoint) for above samples
metadata <- read.table("../wgs_meta_delivery.tsv", header=TRUE)
metadata.for.ecolis <- metadata$err_accession %in% e.colis$accession
metadata <- metadata[metadata.for.ecolis, ]

## Rename samples from the mothers as [Individual]-Mother
metadata$Individual[grepl("Mother", metadata$Time_point)] <- paste(metadata$Individual[grepl("Mother", metadata$Time_point)], "Mother", sep = '-')

## Combine the demix.check output and the metadata
order.for.demix.data <- match(metadata$err_accession, e.colis$accession)
combined.data <- cbind(e.colis[order.for.demix.data, ], metadata)

## Collapse families
combined.data$Individual <- gsub("[-_].*$", "", combined.data$Individual)

## Find identified clusters in each individuals times series
cluster.by.individual <- by(combined.data$cluster, combined.data$Individual, function(x) unique(c(x)))

## Create a presence-absence matrix for the identified clusters
n.individuals <- length(cluster.by.individual)
n.clusters <- length(unique(unlist(cluster.by.individual)))
presabs <- matrix(0, nrow = n.individuals, ncol = n.clusters)
rownames(presabs) <- names(cluster.by.individual)
colnames(presabs) <- unique(unlist(cluster.by.individual))

## Fill the presence-absence matrix
for (i in 1:n.individuals) {
    individual <- names(cluster.by.individual)[i]
    n.clusters.in.series <- length(cluster.by.individual[[i]])
    ## R lists are weird so have to handle length > 1 and length == 1 separately.
    if (n.clusters.in.series > 1) {
        for (j in 1:n.clusters.in.series) {
            cluster <- cluster.by.individual[[i]][j]
            presabs[individual, cluster] <- 1
        }
    } else {
        cluster <- cluster.by.individual[[i]]
        presabs[individual, cluster] <- 1
    }
}

## Rename the clusters to the ST SC format
new.clusters <- read.table("../poppunk_cluster_info/E_col_new_clusters.tsv", sep='\t', header=TRUE)
new.cluster.order <- match(colnames(presabs), new.clusters$PopPUNK)
new.cluster.names <- sub("E_col_(SC[0-9]*)_(ST[0-9]*)", "E_col_\\2_\\1", new.clusters[new.cluster.order, ]$Cluster)
colnames(presabs) <- new.cluster.names

## Write results
presabs <- cbind("Individual" = rownames(presabs), presabs)
write.table(presabs, file = "E_col_collapsed_colonization_matrix.tsv", sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
