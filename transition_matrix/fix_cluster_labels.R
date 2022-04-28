
PrettyClusterLabels <- function(cluster.labels, SpeciesSpecificRenames) {
    ## Change the cluster labels to the format used in the final figure
    ## Run the species-specific fixes first
    cluster.labels <- SpeciesSpecificRenames(cluster.labels)

    ## Run the rest of the fixes
    ## Remove SC designation from the leading SC in coexisting clusters
    cluster.labels <- gsub("_SC[0-9]* ", " ", cluster.labels)

    ## Remove the remaining SC designations
    cluster.labels <- gsub("_SC[0-9]*$", "", cluster.labels)

    cluster.labels
}

EcolPrettyLabels <- function(cluster.labels) {
    ## Change the special E. col cluster names
    ## ST1434 and ST548 coexistence are _probably_ false detections based on looking at the data manually
    cluster.labels <- gsub("ST1434_SC224 x ST1434_SC402 x ST1434_SC723", "ST1434", cluster.labels)
    cluster.labels <- gsub("ST543_SC305 x ST543_SC408", "ST543", cluster.labels)

    ## Rename ST10 subclusters
    cluster.labels <- gsub("ST10_SC6", "ST10-1", cluster.labels)
    cluster.labels <- gsub("ST10_SC351", "ST10-2", cluster.labels)
    cluster.labels <- gsub("ST10_SC79", "ST10-3", cluster.labels)
    cluster.labels
}

EfcsPrettyLabels <- function(cluster.labels) {
    ## Change the special E. fcs cluster names
    cluster.labels
}
