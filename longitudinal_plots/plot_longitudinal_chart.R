RenameTimePointsForSorting <- function(time.points) {
    ## Rename the time points so they sort correctly as characters
    new.time.points <- gsub("Infancy", "30", time.points)
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
    new.time.points <- gsub("Mother", "10", new.time.points)
    new.time.points
}

ExtractCohortBySpecies <- function(demix.check.data, metadata, cohort.name, species.abbrv) {
    ## Extracts the cohort-specific data for some species, with a
    ## possibility to exclude some time points using a regular
    ## expression.
    ##
    demix.check.data <- demix.check.data[demix.check.data[, 3] > 0.01, ]
    demix.check.data <- demix.check.data[grepl(species.abbrv, demix.check.data[, 2]), ]

    ## Extract the cohort-specific metadata
    cohort.metadata <- metadata[metadata$Delivery_mode == cohort.name, ]


    ## Exclude time points
    ## excluded.points <- grepl(exclude.time.points, cohort.metadata$Time_point)
    ## filtered.metadata <- cohort.metadata[!excluded.points, ]
    filtered.metadata <- cohort.metadata
    stopifnot(nrow(filtered.metadata) > 0)
    stopifnot(ncol(filtered.metadata) > 0)

    ## Extract cohort specific, time point filtered data
    cohort.accessions <- demix.check.data$accession %in% filtered.metadata$err_accession
    cohort.data <- demix.check.data[cohort.accessions, ]

    ## Extract the asked species
    cohort.species <- grepl(species.abbrv, cohort.data$cluster)
    species.data <- cohort.data[cohort.species, ]
    stopifnot(nrow(species.data) > 0)
    stopifnot(ncol(species.data) > 0)

    ## Sort the metadata to same order as species.data by accession number
    order.in.metadata <- match(species.data$accession, filtered.metadata$err_accession)
    ordered.metadata <- filtered.metadata[order.in.metadata, ]

    ## Check that the output is sorted correctly
    stopifnot(all(filtered.metadata$err_accession[order.in.metadata] == species.data$accession))

    ## Return a list of the species data and the sorted metadata
    list("demix_data" = species.data, "metadata" = ordered.metadata)
}

ExtractTimeSeries <- function(cohort.data) {
    ## Process the output from ExtractCohortBySpecies: removes time
    ## points that have only one observation and renames the time
    ## points so that they sort in ascending order when using sort()
    
    ## Remove time points that have only one observation
    individual.obs.counts <- table(cohort.data$metadata$Individual)
    more.than.once <- which(individual.obs.counts > 1)
    more.than.once <- cohort.data$metadata$Individual %in% names(individual.obs.counts[more.than.once])
    cohort.data$metadata <- cohort.data$metadata[more.than.once, ]
    cohort.data$demix_data <- cohort.data$demix_data[more.than.once, ]

    ## Rename the [time]-T1 [time]-T2 [time]-T3 observations to just [time]
    cohort.data$metadata$Time_point <- gsub("-T[0-9]*$", "", cohort.data$metadata$Time_point)

    ## Rename the [time] format points so they sort correctly
    cohort.data$metadata$Time_point <- RenameTimePointsForSorting(cohort.data$metadata$Time_point)

    ## Sort the data according to the time point within each individual
    sorted.data <- do.call("rbind", by(cbind(cohort.data$demix_data, cohort.data$metadata), cohort.data$metadata$Individual, function(x) x[order(x$Time_point), ]))

    sorted.data
}

RenameClusters <- function(demix.check.data, cluster.labels) {
    ## Find new names for the demix.check.data clusters
    order.in.labels <- match(demix.check.data$cluster, cluster.labels$PopPUNK)
    ordered.labels <- cluster.labels[order.in.labels, ]
    stopifnot(all(ordered.labels$PopPUNK == demix.check.data$cluster) && length(ordered.labels) > 0)
    demix.check.data$cluster <- ordered.labels$Cluster

    ## Rename from [Species][SC]_[ST] to [Species]_[ST]_[SC]
    demix.check.data$cluster <- sub("(SC[0-9]*)_(ST[0-9]*)", "\\2_\\1", demix.check.data$cluster)
    demix.check.data
}

PersistencePlot <- function(input.data) {
    ## Sort the data within individual according to time point and merge to a single data frame
    plotting.data <- ExtractTimeSeries(input.data)

    cluster.counts <- table(sub("(ST[0-9]*).*$", "\\1", plotting.data$cluster))
    most.common <- names(cluster.counts)[cluster.counts > 4]
    clusters.to.plot <- most.common
    clusters.to.plot <- rev(clusters.to.plot[order(as.numeric(sub("E_col_ST([0-9]*).*$", "\\1", clusters.to.plot)))])

    times.to.plot <- sort(unique(plotting.data$Time_point))
    individuals.to.plot <- unique(plotting.data$Individual)

    n.times.to.plot <- length(times.to.plot)
    n.clusters.to.plot <- length(clusters.to.plot)
    n.individuals.to.plot <- length(individuals.to.plot)

    ## plot(0, type = 'n', xlim = c(0, n.times.to.plot + 1), ylim = c(0, n.clusters.to.plot + 1), xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", bty = 'n')
    ## abline(h = 0:(n.clusters.to.plot - 1) + 0.5, lty = "dashed", col = "gray90")
    ## axis(side = 1, at = 1:n.times.to.plot, labels = times.to.plot)
    ## axis(side = 2, at = 1:n.clusters.to.plot, labels = gsub("E_col_", "", clusters.to.plot), las = 2)
    ## for (i in 1:n.individuals.to.plot) {
    ##     plotting.subset <- plotting.data[plotting.data$Individual %in% individuals.to.plot[i], ]
    ##     lineages.identified <- unique(plotting.subset$cluster)
    ##     n.lineages.identified <- length(lineages.identified)
    ##     for (j in 1:n.lineages.identified) {
    ##         lineage.data <- plotting.subset[plotting.subset$cluster %in% lineages.identified[j], ]
    ##         x.coords <- match(lineage.data$Time_point, times.to.plot)
    ##         y.coords <- match(lineage.data$cluster, clusters.to.plot) + runif(length(lineage.data$cluster), -0.45, 0.45)
    ##         lines(x = x.coords, y = y.coords, type = 'p')
    ##         lines(x = x.coords, y = y.coords)
    ##     }
    ## }
    list("clusters" = clusters.to.plot, "individuals" = individuals.to.plot, "timepoints" = times.to.plot)
}

SwitchingPlot <- function(input.data, clusters.to.plot, individuals.to.plot, times.to.plot) {
    ## Sort the data within individual according to time point and merge to a single data frame
    plotting.data <- ExtractTimeSeries(input.data)

    clusters.at.timepoint <- by(plotting.data, plotting.data$Individual, function(x) x[, c("cluster", "Time_point")])
    co.colonized.at.time <- lapply(clusters.at.timepoint, function(x) by(x, x$Time_point, function(x) length(x$cluster) > 1))

    plotting.data <- plotting.data[plotting.data$cluster %in% clusters.to.plot, ]

##    has.a.switch <- by(plotting.data$cluster, plotting.data$Individual, function(x) length(unique(x)) > 1)
##    plotting.data <- plotting.data[plotting.data$Individual %in% names(has.a.switch)[unlist(has.a.switch)], ]

    times.to.plot <- c(10:23, 30)
    n.times.to.plot <- length(times.to.plot)
    n.clusters.to.plot <- length(clusters.to.plot)
    n.individuals.to.plot <- length(individuals.to.plot)

    plot(0, type = 'n', xlim = c(0, n.times.to.plot + 1), ylim = c(0, n.clusters.to.plot + 1), xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", bty = 'n')
    abline(h = 0:(n.clusters.to.plot - 1) + 0.5, lty = "solid", col = "gray90")

    xaxlabs <- c("Mother", paste(c(4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18, 21), 'd', sep=''), "Infancy")
    axis(side = 1, at = 1:n.times.to.plot, labels = xaxlabs)
    axis(side = 2, at = 1:n.clusters.to.plot, labels = gsub("E_col_", "", clusters.to.plot), las = 2, cex.axis = 1.5)

    s <- 0
    for (i in 1:n.individuals.to.plot) {
        individual.name <- individuals.to.plot[i]
        plotting.subset <- plotting.data[plotting.data$Individual %in% individual.name, ]
        timepoints.to.plot <- sort(unique(plotting.subset$Time_point))
        n.timepoints <- length(timepoints.to.plot)
        if (n.timepoints > 1) {
            for (j in 2:n.timepoints) {
                prev.time <- timepoints.to.plot[j - 1]
                prev.timepoint <- plotting.subset[plotting.subset$Time_point == prev.time, ]
                current.time <- timepoints.to.plot[j]
                current.timepoint <- plotting.subset[plotting.subset$Time_point == current.time, ]
                y.coords.prev <- match(prev.timepoint$cluster, clusters.to.plot)
                y.coords.current <- match(current.timepoint$cluster, clusters.to.plot)
                x.coords.prev <- match(prev.timepoint$Time_point, times.to.plot)
                x.coords.current <- match(current.timepoint$Time_point, times.to.plot)
                for (k in 1:length(y.coords.prev)) {
                    for (l in 1:length(y.coords.current)) {
                        noise <- runif(1, -0.45, 0.45)
                        y.coords <- c(y.coords.prev[k], y.coords.current[l])
                        x.coords <- c(x.coords.prev[k], x.coords.current[l])
                        start.or.end <- x.coords[1] == 1 || x.coords[2] == n.times.to.plot
                        switch.at.end <- (x.coords[2] == n.times.to.plot)## && (y.coords[1] != y.coords[2])
                        switch.at.start <- (x.coords[1] == 1)## && (y.coords[1] != y.coords[2])
                        if (switch.at.end && y.coords[1] != y.coords[2]) {
                            s <- s + 1
                            print(s)
                        }
                        y.coords <- y.coords + noise
                        prev.co.colonized <- co.colonized.at.time[[individual.name]][[as.character(times.to.plot[x.coords[1]])]]
                        current.co.colonized <- co.colonized.at.time[[individual.name]][[as.character(times.to.plot[x.coords[2]])]]
                        if (switch.at.start || switch.at.end) {
                            lines(x = x.coords, y = y.coords, col = ifelse(switch.at.end, "#fc8d62", "#66c2a5"), lty = ifelse(y.coords[1] != y.coords[2], "dashed", "solid"), lwd = 2)
                        } else {
                            lines(x = x.coords, y = y.coords, col = ifelse(switch.at.end, "#fc8d59", "#8da0cb"), lty = ifelse(switch.at.end, "dashed", "solid"), lwd = 1.5)
                        }
                        lines(x = x.coords[1], y = y.coords[1], type = 'p', ifelse(prev.co.colonized, 4, 1), col = ifelse(prev.co.colonized, "gray40", "gray40"), cex = 1.5)
                        lines(x = x.coords[2], y = y.coords[2], type = 'p', ifelse(current.co.colonized, 4, 1), col = ifelse(current.co.colonized, "gray40", "gray40"), cex = 1.5)
                    }
                }
            }
        }
    }
}

library("stringr")
library("graphics")

source("longtitudinal_plotting.R")
source("fix_cluster_labels.R")

## Read in the data from demix_check
demix.check.data <- read.table("../wgs_demix_check_high_confidence.tsv", sep='\t', header=TRUE)

## Falsely id'd clusters
demix.check.data <- demix.check.data[!grepl("E_col_Pop408$", demix.check.data$cluster), ]
demix.check.data <- demix.check.data[!grepl("E_col_Pop402$", demix.check.data$cluster), ]
demix.check.data <- demix.check.data[!grepl("E_col_Pop723$", demix.check.data$cluster), ]

## Read in the metadata (accession, delivery_mode, individual, timepoint)
full.metadata <- read.table("../wgs_meta_delivery.tsv", sep = '\t', header=TRUE)

## Read in the cohort specific data
caesarean <- ExtractCohortBySpecies(demix.check.data, full.metadata, "Caesarean", "E_col")
vaginal <- ExtractCohortBySpecies(demix.check.data, full.metadata, "Vaginal", "E_col")

## Rename the clusters
new.clusters <- read.table("../poppunk_cluster_info/E_col_new_clusters.tsv", sep='\t', header=TRUE)
caesarean$demix_data <- RenameClusters(caesarean$demix_data, new.clusters)
caesarean$demix_data$cluster <- PrettyClusterLabels(caesarean$demix_data$cluster, EcolPrettyLabels)
vaginal$demix_data <- RenameClusters(vaginal$demix_data, new.clusters)
vaginal$demix_data$cluster <- PrettyClusterLabels(vaginal$demix_data$cluster, EcolPrettyLabels)

pdf(file = "E_col_persistence_2.pdf", width = 18, height = 8)
params <- PersistencePlot(vaginal)
layout(matrix(c(1, 2, 3, 3), byrow = TRUE, nrow = 2, ncol = 2), heights = c(0.9, 0.1), widths = c(0.5, 0.5))
par(mar = c(2, 5, 0, 2))
SwitchingPlot(vaginal, params$clusters, params$individuals, params$timepoints)
title("a) Vaginal delivery cohort", line = -1, adj = 0.01, outer = TRUE, cex.main = 1.5)
par(mar = c(2, 5, 0, 2))
params <- PersistencePlot(caesarean)
SwitchingPlot(caesarean, params$clusters, params$individuals, params$timepoints)
title("b) Caesarean delivery cohort", line = -1, adj = 0.52, outer = TRUE, cex.main = 1.5)
par(mar = c(0, 0, 1, 0))
plot.new()
legend("top", legend = c("Co-colonized", "Single lineage", "Persistence", "Displacement", "First 21 days", "After 4-12 months", "Maternal transmission/displacement"), ncol = 4, bty = 'n', col = c("gray40", "gray40", "#8da0cb", "gray40", "gray40", "#fc8d62", "#66c2a5"), pch = c(4, 1, NA, NA, NA, NA, NA), lty = c(NA, NA, "solid", "dashed", "solid", "solid", "solid"), lwd = c(NA, NA, 2, 2, 1, 2, 2), cex = 1.5)
dev.off()


                        ##     lines(x = x.coords, y = y.coords, col = ifelse(switch.at.end, "#fc8d62", "#66c2a5"), lty = ifelse(y.coords[1] != y.coords[2], "dashed", "solid"))
                        ## } else {
                        ##     lines(x = x.coords, y = y.coords, col = ifelse(switch.at.end, "#fc8d59", "#8da0cb"), lty
## vals.to.plot <- cbind(
##    c("E_col_", "E_fcs_", "K_pne_", "K_mic_", "K_var_", "K_oxy_", "K_gri_"),
##    c("Escherichia coli", "Enterococcus faecalis", "Klebsiella pneumoniae",
##      "Klebislella michiganensis", "Klebsiella variicola", "Klebsiella oxytoca",
##      "Klebsiella grimontii"))

vals.to.plot <- cbind(c("E_col_"), c("Escherichia coli"))

for (i in 1:nrow(vals.to.plot)) {
    plot.name <- paste(vals.to.plot[i, 1], "caesarean_delivery_longitudinal.pdf", sep='')
    pdf(file = plot.name, width = ifelse(i < 3, 40, 10), height = ifelse(i < 3, 75, 20))
    par(mar = c(1, 9, 18, 2))
    PlotByHospital(caesarean, dat, "Caesarean", metadata,
                   vals.to.plot[i, 1], vals.to.plot[i, 2])
    dev.off()
    plot.name <- paste(vals.to.plot[i, 1], "vaginal_delivery_longitudinal.pdf", sep='')
    pdf(file = plot.name, width = ifelse(i < 3, 40, 10), height = ifelse(i < 3, 75, 20))
    par(mar = c(1, 9, 18, 2))
    PlotByHospital(vaginal$demix_data, dat, "Vaginal", vaginal$metadata,
                   vals.to.plot[i, 1], vals.to.plot[i, 2])
    dev.off()
}



length(unique(vaginal$metadata$Individual[vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"]))

length(unique(caesarean$metadata$Individual[caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"]))

length(unique(vaginal$metadata$Individual[vaginal$metadata$Time_point == "4" & vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"]))

length(unique(vaginal$metadata$Individual[vaginal$metadata$Time_point %in% c("6", "7") & vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"]))


length(unique(vaginal$metadata$Individual[vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"]))

in.first.time.point <- paste(vaginal$metadata$Individual[vaginal$metadata$Time_point == "4" & vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"], vaginal$demix_data$cluster[vaginal$metadata$Time_point == "4" & vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"], sep='-')

in.second.time.point <- paste(vaginal$metadata$Individual[vaginal$metadata$Time_point %in% c("6", "7") & vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"], vaginal$demix_data$cluster[vaginal$metadata$Time_point %in% c("6", "7") & vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"], sep='-')

length(unique(gsub("-.*$", "", in.first.time.point[in.first.time.point %in% in.second.time.point])))


######

length(unique(caesarean$metadata$Individual[caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"]))

length(unique(caesarean$metadata$Individual[caesarean$metadata$Time_point == "4" & caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"]))

length(unique(caesarean$metadata$Individual[caesarean$metadata$Time_point %in% c("6", "7") & caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"]))



in.first.time.point <- paste(caesarean$metadata$Individual[caesarean$metadata$Time_point == "4" & caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"], caesarean$demix_data$cluster[caesarean$metadata$Time_point == "4" & caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"], sep='-')

in.second.time.point <- paste(caesarean$metadata$Individual[caesarean$metadata$Time_point %in% c("6", "7") & caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"], caesarean$demix_data$cluster[caesarean$metadata$Time_point %in% c("6", "7") & caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"], sep='-')

length(unique(gsub("-.*$", "", in.first.time.point[in.first.time.point %in% in.second.time.point])))


###

in.first.time.point1 <- paste(vaginal$metadata$Individual[vaginal$metadata$Time_point %in% c("4", "5", "6", "7", "8") & vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"], vaginal$demix_data$cluster[vaginal$metadata$Time_point %in% c("4", "5", "6", "7", "8") & vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"], sep='-')
in.mother1 <- paste(vaginal$metadata$Individual[vaginal$metadata$Time_point == "Mother"], vaginal$demix_data$cluster[vaginal$metadata$Time_point == "Mother"], sep='-')
sum(in.mother1 %in% in.first.time.point1)
length(in.mother1)

in.first.time.point2 <- paste(caesarean$metadata$Individual[caesarean$metadata$Time_point %in% c("4", "5", "6", "7", "8") & caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"], caesarean$demix_data$cluster[caesarean$metadata$Time_point %in% c("4", "5", "6", "7", "8") & caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"], sep='-')
in.mother2 <- paste(caesarean$metadata$Individual[caesarean$metadata$Time_point == "Mother"], caesarean$demix_data$cluster[caesarean$metadata$Time_point == "Mother"], sep='-')
sum(in.mother2 %in% in.first.time.point2)
length(in.mother2)

####


in.first.time.point1 <- paste(vaginal$metadata$Individual[vaginal$metadata$Time_point == c("Infancy") & vaginal$metadata$Time_point != "Mother"], vaginal$demix_data$cluster[vaginal$metadata$Time_point == "Infancy" & vaginal$metadata$Time_point != "Mother"], sep='-')
in.mother1 <- paste(vaginal$metadata$Individual[vaginal$metadata$Time_point == "Mother"], vaginal$demix_data$cluster[vaginal$metadata$Time_point == "Mother"], sep='-')
sum(in.mother1 %in% in.first.time.point1)
length(in.mother1)

in.first.time.point2 <- paste(caesarean$metadata$Individual[caesarean$metadata$Time_point == "Infancy" & caesarean$metadata$Time_point != "Mother"], caesarean$demix_data$cluster[caesarean$metadata$Time_point %in% c("Infancy") & caesarean$metadata$Time_point != "Mother"], sep='-')
in.mother2 <- paste(caesarean$metadata$Individual[caesarean$metadata$Time_point == "Mother"], caesarean$demix_data$cluster[caesarean$metadata$Time_point == "Mother"], sep='-')
sum(in.mother2 %in% in.first.time.point2)
length(in.mother2)


########

co.colonized <- paste(vaginal$metadata$Individual[vaginal$metadata$Time_point != "Infancy" & vaginal$metadata$Time_point != "Mother"], vaginal$metadata$Time_point[vaginal$metadata$Time_point != "Infancy" & vaginal$metadata$Time_point != "Mother"], sep='-')
length(unique(gsub("-.$", "", co.colonized <- co.colonized[duplicated(co.colonized)])))
co.colonized <- paste(caesarean$metadata$Individual[caesarean$metadata$Time_point != "Infancy" & caesarean$metadata$Time_point != "Mother"], caesarean$metadata$Time_point[caesarean$metadata$Time_point != "Infancy" & caesarean$metadata$Time_point != "Mother"], sep='-')
length(unique(gsub("-.$", "", co.colonized <- co.colonized[duplicated(co.colonized)])))

co.colonized <- paste(vaginal$metadata$Individual[vaginal$metadata$Time_point == "Infancy" & vaginal$metadata$Time_point != "Mother"], vaginal$metadata$Time_point[vaginal$metadata$Time_point == "Infancy" & vaginal$metadata$Time_point != "Mother"], sep='-')
length(unique(gsub("-.$", "", co.colonized <- co.colonized[duplicated(co.colonized)])))
co.colonized <- paste(caesarean$metadata$Individual[caesarean$metadata$Time_point == "Infancy" & caesarean$metadata$Time_point != "Mother"], caesarean$metadata$Time_point[caesarean$metadata$Time_point == "Infancy" & caesarean$metadata$Time_point != "Mother"], sep='-')
length(unique(gsub("-.$", "", co.colonized <- co.colonized[duplicated(co.colonized)])))

co.colonized <- paste(vaginal$metadata$Individual[vaginal$metadata$Time_point == "Mother"], vaginal$metadata$Time_point[vaginal$metadata$Time_point == "Mother"], sep='-')
length(unique(gsub("-.$", "", co.colonized <- co.colonized[duplicated(co.colonized)])))
co.colonized <- paste(caesarean$metadata$Individual[caesarean$metadata$Time_point == "Mother"], caesarean$metadata$Time_point[caesarean$metadata$Time_point == "Mother"], sep='-')
length(unique(gsub("-.$", "", co.colonized <- co.colonized[duplicated(co.colonized)])))


in.first.time.point1 <- paste(vaginal$metadata$Individual[vaginal$metadata$Time_point == c("Infancy") & vaginal$metadata$Time_point != "Mother"], vaginal$demix_data$cluster[vaginal$metadata$Time_point == "Infancy" & vaginal$metadata$Time_point != "Mother"], sep='-')
in.mother1 <- paste(vaginal$metadata$Individual[vaginal$metadata$Time_point == "Mother"], vaginal$demix_data$cluster[vaginal$metadata$Time_point == "Mother"], sep='-')
sum(in.mother1 %in% in.first.time.point1)
length(in.mother1)

in.first.time.point2 <- paste(caesarean$metadata$Individual[caesarean$metadata$Time_point == "Infancy" & caesarean$metadata$Time_point != "Mother"], caesarean$demix_data$cluster[caesarean$metadata$Time_point %in% c("Infancy") & caesarean$metadata$Time_point != "Mother"], sep='-')
in.mother2 <- paste(caesarean$metadata$Individual[caesarean$metadata$Time_point == "Mother"], caesarean$demix_data$cluster[caesarean$metadata$Time_point == "Mother"], sep='-')
sum(in.mother2 %in% in.first.time.point2)
length(in.mother2)



## Reviewer 1 comment 14

in.first.time.point1 <- paste(vaginal$metadata$Individual[vaginal$metadata$Time_point %in% c("4", "5", "6", "7", "8") & vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"], vaginal$demix_data$cluster[vaginal$metadata$Time_point %in% c("4", "5", "6", "7", "8") & vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"], sep='-')
in.last.time.point1 <- paste(vaginal$metadata$Individual[vaginal$metadata$Time_point %in% c("21") & vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"], vaginal$demix_data$cluster[vaginal$metadata$Time_point %in% c("21") & vaginal$metadata$Time_point != "Mother" & vaginal$metadata$Time_point != "Infancy"], sep='-')
sum(in.last.time.point1 %in% in.first.time.point1)
length(in.last.time.point1)

in.first.time.point1 <- paste(caesarean$metadata$Individual[caesarean$metadata$Time_point %in% c("4", "5", "6", "7", "8") & caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"], caesarean$demix_data$cluster[caesarean$metadata$Time_point %in% c("4", "5", "6", "7", "8") & caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"], sep='-')
in.last.time.point1 <- paste(caesarean$metadata$Individual[caesarean$metadata$Time_point %in% c("21") & caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"], caesarean$demix_data$cluster[caesarean$metadata$Time_point %in% c("21") & caesarean$metadata$Time_point != "Mother" & caesarean$metadata$Time_point != "Infancy"], sep='-')
sum(in.last.time.point1 %in% in.first.time.point1)
length(in.last.time.point1)

## Reviewer 1 comment 16
sum(table(paste(vaginal$metadata$Individual, vaginal$metadata$Time_point)) > 1)
length(table(paste(vaginal$metadata$Individual, vaginal$metadata$Time_point)))

sum(table(paste(caesarean$metadata$Individual, caesarean$metadata$Time_point)) > 1)
length(table(paste(caesarean$metadata$Individual, caesarean$metadata$Time_point)))
