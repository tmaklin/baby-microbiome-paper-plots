
PointAreaLegend <- function() {
    plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
    legend("bottom", legend = c("0-1", 5, 10, 20, 100), fill = "white", pch = 19,
           border = "white", bty = 'n', ncol = 5,
           pt.cex = log(c(1, 5, 10, 20, 100),base = 6) + 1, inset = c(0, 0.07))
    title(main = "# of times a lineage was reliably\nidentified from both species", line = -30,
          font.main = 1, cex.main = 1.1)
}

CorrelationIntensityLegend <- function(ColorFunc) {
    par(xpd = TRUE)
    legend_image <- as.raster(matrix(rev(ColorFunc(41)), ncol=1))
    plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
    title(main = "Estimated\ncorrelation", font.main = 1, cex.main = 1.1,
          line = -26, adj = 0.1)
    text(x=0.9, y = seq(0,1,l=5), labels = seq(-0.20,0.20,l=5), adj = 0, cex = 1.1)
    rasterImage(legend_image, 0, 0, 0.5,1)
}

CorrelationPlot <- function(correlations, obs.counts, otu.names, title.main, title.adj, ColorFunc) {
    n.otus <- length(otu.names)
    plot(0, type = 'n', xlim = c(0, n.otus), ylim = c(0, n.otus),
         xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
    ## top axis
    axis(side = 3, at = 1:n.otus, labels = otu.names, las = 2, font = 3)
    abline(v = 1:n.otus, col = "gray90")
    ## left axis, labels in reverse order to match top axis
    axis(side = 2, at = 1:n.otus, labels = rev(otu.names), las = 2, font = 3)
    abline(h = 1:n.otus, col = "gray90")
    for (i in 1:n.otus) {
        for (j in 1:n.otus) {
            if (correlations[i, j] != 0) {
                coord.x <- i
                coord.y <- n.otus - j + 1 ## filled from top->down
                ## assume correlations range in [-0.20, 0.20] and counts in [0, 100]
                color.id <- round(correlations[i, j], digits = 2)*100 + 21
                size <- log(obs.counts[i, j] + 1, base = 6) + 1
                lines(x = coord.x, y = coord.y, type = 'p', pch = 19,
                      cex = size, col = ColorFunc(41)[color.id])
            }
        }
    }
    title(main = title.main, line = -1, adj = title.adj, outer = TRUE)
}

ReadDemixResults <- function() {
    demix.results <- read.table("../wgs_demix_check_high_confidence.tsv", sep='\t', header=TRUE)
    something.passes <- by(demix.results$cluster, demix.results$accession, function(x) unique(gsub("_NA|_SC[0-9]*_ST[0-9]*|_Pop[0-9]*|_SC[0-9]*", "", c(x))))
    demix.results.ecoli <- read.table("../ecoli-new-reference/E_col_demix_results_top2.tsv", sep='\t', header=FALSE)
    something.passes2 <- by(demix.results.ecoli$V2, demix.results.ecoli$V1, function(x) unique(gsub("_NA|_SC[0-9]*_ST[0-9]*|_Pop[0-9]*", "", c(x))))
    l <- list(something.passes, something.passes2)
    keys <- unique(unlist(lapply(l, names)))
    tst <- do.call(mapply, c(FUN=c, lapply(l, `[`, keys)))
    demix.filter <- lapply(tst, function(x) x[!is.na(x)])
    demix.filter
}

ReadCorrelations <- function(correlation.path, pvalues.path, OtuFilter) {
    correlations <- read.table(correlation.path, sep='\t')[, -1]
    pvalues <- read.table(pvalues.path, sep='\t')[, -1]
    correlations[pvalues >= 0.05] <- 0
    rownames(correlations) <- read.table(correlation.path, sep='\t')[, 1]
    colnames(correlations) <- rownames(correlations)
    diag(correlations) <- 0
    priority.pathogens <- OtuFilter(rownames(correlations))
    correlations
}

DemixCheckCounts <- function(demix.filter, accessions, correlations) {
    otu.names <- rownames(correlations)
    indices <- lapply(demix.filter[names(demix.filter) %in% accessions], function(x) t(combn(x, min(length(x), 2))))
    counts <- correlations
    counts[counts != 0] <- 0
    for (i in 1:length(indices)) {
        if (length(indices[[i]]) > 1) {
            for (j in 1:nrow(indices[[i]])) {
                if (indices[[i]][j, 1] %in% otu.names && indices[[i]][j, 2] %in% otu.names) {
                    counts[indices[[i]][j, 1], indices[[i]][j, 2]] <- counts[indices[[i]][j, 1], indices[[i]][j, 2]] + 1
                    counts[indices[[i]][j, 2], indices[[i]][j, 1]] <- counts[indices[[i]][j, 2], indices[[i]][j, 1]] + 1
                }
            }
        }
    }
    counts
}
