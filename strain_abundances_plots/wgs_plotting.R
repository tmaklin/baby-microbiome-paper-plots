## ## ##
## This file contains the functions for plotting the strain abundances
## split by delivery mode and timepoint.
## ##
##
library("stats")
library("graphics")

## https://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x
roundUpNice <- function(x, nice=1:20) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

PlotSubplot <- function(all.values, plotting.order.all, legend.order, subset, colors, ylabel, xmax) {
    ## Plot a single barplot containing the samples defined by `subset`.

    if (length(subset) > 0) {
        subset.abundances <- all.values[, subset]
    }

    ## R is a dumbass language that converts the matrix to a vector if it has just 1 column - fix here
    if (sum(subset) > 1) {
        plotting.data <- subset.abundances[legend.order, ]
    } else if (length(subset) > 0) {
        plotting.data <- matrix(subset.abundances[legend.order], ncol = 1)
    } else {
        plotting.data <- c(NA)
    }

    xmax <- max(xmax + 2, 10)
    ## Plot the actual barplot
    barplot(plotting.data, col=colors, xaxt = 'n', ylab = ylabel, cex.lab = 3, xlim = c(0, xmax), cex.axis = 1.3, width = 1, ylim = c(0, 1))
}

OrderLegendByDominant <- function(strain_abundances, plotting.order.all, frac.dominant) {
    ## Find the dominant strain in each sample and create an ordering
    ## for the legend based on their order of appearance.

    ## Find the dominant strains
    if (ncol(strain_abundances) > 1) {
        dominant.strains <- apply(strain_abundances[, plotting.order.all], 2, function(x) ifelse(any(x > frac.dominant), which(x > frac.dominant), NA))
    } else {
        dominant.strains <- which(strain_abundances > 0)
    }
    dominant.strains <- dominant.strains[!is.na(dominant.strains)]

    ## Order the legend
    legend.order <- unique(dominant.strains)

    ## Include lineages that are not dominant in any sample
    missing.legends <- setdiff(1:nrow(strain_abundances), legend.order)
    legend.order <- c(legend.order, missing.legends)
    return(legend.order)
}

PlotColumn <- function(sampleset, strain_abundances, metadata, delivery_mode, timepoints,
                       plotting.order.all, legend.order, color, n.obs.max, column.title) {
    ## Plot all plots in one column (delivery mode)

    n.timepoints <-  length(timepoints)
    for (i in 1:n.timepoints) {
        ## Leave some space at the top in the first plot
        par(mar = c(1, 6, ifelse(i > 1, 0, 3), 1))

        ## Find which samples are assigned to this timepoint
        subset <- colnames(strain_abundances) %in% metadata[metadata[, 2] == delivery_mode & metadata[, 4] == timepoints[i], 1]

        ## Plot the timepoint
        plot.label <- ifelse(grepl("Infancy|Mother", timepoints[i]), timepoints[i], paste(timepoints[i], 'd', sep=''))
        PlotSubplot(strain_abundances, plotting.order.all, legend.order, subset, color, plot.label, n.obs.max)

        ## Column title in first row
        if (i == 1)
            title(column.title, cex.main = 3.5)
    }
}

PlotBin <- function(bin.abundances, color, main.text, legend.remove, timepoints, metadata, legend.size, legend.cols) {
    ## Plots the samples in the bin divided into Vaginal/Caesarean
    ## delivery and by timepoint.
    rownames(bin.abundances) <- sub(".*(ERR[0-9]*).*", "\\1", rownames(bin.abundances))
    bin.abundances <- t(bin.abundances)

    ## Extract the number of timepoints from the data
    n.timepoints <- length(timepoints)
    ## Layout the plots - there are n.timepoints rows in 2 columns.
    layout(matrix(c(1:n.timepoints, 2*n.timepoints + 1, (n.timepoints + 1):(2*n.timepoints), 2*n.timepoints + 1),
                  nrow = n.timepoints + 1, ncol = 2), widths = 1, heights = c(rep(1, n.timepoints), 2, rep(1, n.timepoints), 2))

    ## Use hclust to group similar samples together
    ## hclust returns an error if there is only one sample, use tryCatch to avoid it
    clusts <- tryCatch(hclust(dist(t(bin.abundances), method = "euclidean"), method = "complete"), error=function(e) c(1))
    plotting.order.all <- tryCatch(clusts$order, error=function(e) c(1))

    ## Order legend by approximate order of appearance in the hclust ordering
    legend.order <- OrderLegendByDominant(bin.abundances, plotting.order.all, 0.50)
    color <- color[legend.order]

    ## Split samples by delivery mode
    strain_abundances <- bin.abundances[, plotting.order.all]
    if (!is.matrix(strain_abundances)) {
        out <- matrix(strain_abundances, ncol = 1)
        rownames(out) <- names(strain_abundances)
        colnames(out) <- colnames(bin.abundances)
        strain_abundances <- out
    }

    vaginal <- colnames(strain_abundances) %in% metadata[metadata[, 2] == "Vaginal", 1]
    caesarean <- colnames(strain_abundances) %in% metadata[metadata[, 2] == "Caesarean", 1]

    ## Find out how many samples are in the largest group, needed for
    ## fixing the column widths across all plots.
    n.obs.max <- 0
    for (i in 1:n.timepoints) {
        n.obs.max <- max(n.obs.max, sum(colnames(strain_abundances)%in% metadata[metadata[, 2] == "Vaginal" & metadata[, 4] == timepoints[i], 1]))
        n.obs.max <- max(n.obs.max, sum(colnames(strain_abundances) %in% metadata[metadata[, 2] == "Caesarean" & metadata[, 4] == timepoints[i], 1]))
    }

    ## Add some space after the last value
    n.obs.max <- n.obs.max + roundUpNice(max(sum(vaginal), sum(caesarean))/20)

    ## ## ##
    ## Plot the values
    ## Vaginal delivery

    PlotColumn(vaginal, strain_abundances, metadata, "Vaginal", timepoints,
               plotting.order.all, legend.order, color, n.obs.max, "Vaginal delivery")
    ## Caesarean delivery
    PlotColumn(caesarean, strain_abundances, metadata, "Caesarean", timepoints,
               plotting.order.all, legend.order, color, n.obs.max, "Caesarean delivery")

    ## Legend
    par(mar = c(2, 6, 8, 2))
    plot.new()
    ## Relabel the lineages as ST[0-9]* SC[0-9]*
    if (is.null(rownames(strain_abundances))) {
        ## R and its """""useful""""" type inferring
        legend.text <- names(strain_abundances)
    } else {
        legend.text <- rownames(strain_abundances)
    }
    legend.text <- gsub(legend.remove, "", legend.text[legend.order])
    legend.text <- gsub("[_]*NA", "", legend.text)
    STs <- gsub("SC[0-9]*_", "", legend.text)
    STs <- gsub("ST[ ]*$", "", STs)
    SCs <- gsub("_ST[0-9]*", "", legend.text)
    if (all(STs != SCs)) {
        no.st <- STs == ''
        SCs[!no.st] <- paste(' ', SCs[!no.st], sep='')
        legend.text <- paste(STs, SCs, sep='')
    }
    legend("topleft", legend = head(legend.text, 120), fill = head(color, 120), horiz = FALSE, ncol = legend.cols, bty = 'n', xpd = TRUE, cex=legend.size)
    ## Title for the legend
    title(main.text, adj=0, cex.main=legend.size + 0.5)
    return(cbind(legend.text, color))
}
