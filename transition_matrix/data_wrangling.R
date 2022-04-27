
ReadEcolAbundances <- function(cohort, full.metadata, what.to.extract = "E_col") {
    ## Reads in the abundance matrix with rows corresponding to
    ## samples from each individual at each time point and cols to the
    ## identified lineages. Nonzero values in the matrix mean that the
    ## lineage passed demix_check. The values themselves are filtered
    ## relative abundances, where the filtering has set the values
    ## that did not pass demix_check to zero.

    ## Find cohort-specific samples in the metadata
    order.in.metadata <- match(cohort$V1, full.metadata$err_accession)

    ## Rename time points so they sort correctly with sort()
    new.time.points <- full.metadata[order.in.metadata, ]$Time_point
    new.time.points <- gsub("Mother", "10", new.time.points)
    new.time.points <- gsub("^21$", "23", new.time.points)
    new.time.points <- gsub("^18$", "22", new.time.points)
    new.time.points <- gsub("^17$", "21", new.time.points)
    new.time.points <- gsub("^14$", "20", new.time.points)
    new.time.points <- gsub("^13$", "19", new.time.points)
    new.time.points <- gsub("^12$", "18", new.time.points)
    new.time.points <- gsub("^11$", "17", new.time.points)
    new.time.points <- gsub("^10$", "16", new.time.points)
    new.time.points <- gsub("^9$", "15", new.time.points)
    new.time.points <- gsub("^8$", "14", new.time.points)
    new.time.points <- gsub("^7$", "13", new.time.points)
    new.time.points <- gsub("^6$", "12", new.time.points)
    new.time.points <- gsub("^4$", "11", new.time.points)
    new.time.points <- gsub("Infancy", "30", new.time.points)

    ## Extract the column and row names
    cols <- unique(cohort$V2) ## V2 is the identified lineage's name
    rows <- cbind("Individual" = full.metadata[order.in.metadata, ]$Individual, "Time_point" = new.time.points)
    row.order <- order(rows[, 1]) ## rows[, 1] contains individual names
    vals.ordered <- cbind(cohort[row.order, ], rows[row.order, 1], rows[row.order, 2]) ## rows[, 2] contains the sampling timepoint (in days since birth)

    ## Read in the abundances
    abundances <- matrix(0, nrow = nrow(rows), ncol = length(cols))
    rownames(abundances) <- paste(vals.ordered[, 13], vals.ordered[, 14], sep='-') ## Merge individual names and timepoints
    colnames(abundances) <- cols
    for (i in 1:nrow(rows)) {
        ## vals.ordered[, 13] is again the individual name and vals.ordered[, 14] the timepoint
        ## Find the correct row in abundances by merging the names as above before the loop
        abundances[paste(vals.ordered[i, 13], vals.ordered[i, 14], sep='-'), vals.ordered[i, 2]] <- vals.ordered[i, 3]
    }

    ## Only want to extract E. colis by default
    subset <- abundances[, grepl(what.to.extract, colnames(abundances))]
    subset
}

ExtractStateNames <- function(abundances) {
    ## Extracts the names of states that were visited at least once in
    ## the abundance matrix produced by ReadEcolAbundances() States
    ## that are of the form "State_1 x State_2" mean that State_1 and
    ## State_2 were present in the same sample at the same time point.
    states <- sort(unique(unlist(lapply(apply(abundances, 1, function(x) colnames(abundances)[x > 0]), function(x) paste(x, collapse=' x ')))))
    states <- c(states, "No E. coli")
    states
}

ExtractTimeSeriesSamples <- function(abundances) {
    ## Extracts samples that have at least 2 observations from each
    ## individual aka a time series.

    ## Merge timepoints that have multiple samples taken from them
    ## T1 T2 T3 etc should be merged by hour, they're _probably_ the same sample but just duplicate
    abundances <- do.call("rbind", by(abundances, factor(gsub("_T[0-9]", "", rownames(abundances))), colSums))

    ## Turn the sample names to individual labels by removing the time point designation
    individual.identifier <- gsub("-.*$", "", rownames(abundances))

    ## Extract individual labels that appear at least 2 times
    transition.samples <- names(which(table(individual.identifier) > 1))

    ## Subset the values and return the time series samples
    time.series.subset <- abundances[individual.identifier %in% transition.samples, ]
    time.series.individuals <- gsub("-.*$", "", rownames(time.series.subset))

    ## Return a list containing the time series abundances and the individuals that have a time series
    list("abundances" = time.series.subset, "individuals" = time.series.individuals)
}

GetTransitions <- function(time.series) {
    ## Given an abundance matrix of observations from the same
    ## individual in temporal order (first observation == first row
    ## etc.), return the set of transitions observed in the time
    ## seires.
    n.steps <- nrow(time.series)
    transitions <- matrix("", nrow = n.steps - 1, ncol = 2)
    oldstate <- paste(gsub("E_col_", "", colnames(time.series)[time.series[1, ] > 0]), collapse=' x ')
    for (i in 2:n.steps) {
        oldstate[oldstate == ""] <- "No E. coli"
        transitions[i - 1, 1] <- oldstate
        newstate <- paste(gsub("E_col_", "", colnames(time.series)[time.series[i, ] > 0]), collapse=' x ')
        newstate[newstate == ""] <- "No E. coli"
        transitions[i - 1, 2] <- newstate
        oldstate <- newstate
    }
    transitions
}

BuildTransitionMatrix <- function(time.series, states) {
    ## Builds a transition matrix from the individual transitions
    ## returned by GetTransitions() and the states returned by
    ## ExtractStateNames()
    n.obs <- length(time.series)
    n.states <- length(states)
    transition.mat <- matrix(0, n.states, n.states)
    rownames(transition.mat) <- states
    colnames(transition.mat) <- states
    for (i in 1:n.obs) {
        transitions <- time.series[[i]]
        transition.mat[transitions] <- transition.mat[transitions] + 1
    }
    transition.mat
}

GetTransitionMatrix <- function(cohort, full.metadata, species.name) {
    ## Runs the above functions in correct order to return the
    ## transition matrix.

    ## Read in abundances
    subset <- ReadEcolAbundances(cohort, full.metadata, species.name)

    ## Extract samples that have at least 2 observations for the individual
    time.series <- ExtractTimeSeriesSamples(subset)

    ## Extract transition matrix states
    states <- ExtractStateNames(time.series$abundances)
    ## Remove the leading "E_col_" part
    states <- gsub("E_col_", "", states)

    ## Extract transitions within each time series observation
    transitions.within.individual <- by(time.series$abundances, time.series$individuals, GetTransitions)

    ## Build and return the transition matrix
    transition.mat <- BuildTransitionMatrix(transitions.within.individual, states)
    transition.mat
}


ReadTransitionMatrix <- function(metadata.path, demix.path, new.cluster.names.path, cohort.name, exclude.time.points, species.name) {
    metadata <- read.table(metadata.path, header=TRUE)

    ## Read in the results from demix_check
    demix.results <- read.table(demix.path, sep='\t', header=FALSE, stringsAsFactors=FALSE)
    ## Filter by species
    demix.results <- demix.results[grepl(species.name, demix.results$V2), ]
    ## Only use samples that have a relative abundance higher than 1%
    demix.results <- demix.results[demix.results$V3 > 0.01, ]

    ## Extract cohort and remove excluded time points
    accessions <- metadata$err_accession[!grepl(exclude.time.points, metadata$Time_point) & metadata$Delivery_mode == cohort.name]
    samples <- demix.results[demix.results$V1 %in% accessions, ]

    ## ## Rename the clusters from [Species]_Pop[0-999] to [Species]_ST[0-9999]_SC[0-999]
    new.clusters <- read.table(file = new.cluster.names.path, sep='\t', header = TRUE)
    samples$V2 <- paste(species.name, "_", new.clusters$Cluster[match(samples$V2, new.clusters$PopPUNK)], sep = '')

    ## ## Build the transition matrix
    mat <- GetTransitionMatrix(samples, metadata, paste(species.name, "_", sep = ''))

    ## Only return lineages which had more than 2 transitions in either the column _or_ the row
    more.than.once <- which(rowSums(mat) > 1 | colSums(mat) > 1)
    mat[more.than.once, more.than.once]
}
