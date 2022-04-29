PlotLongitudinal <- function(subset, name.prefix, main.title) {
    is.nonzero <- which(by(subset, gsub("-[0-9]*$", "", rownames(subset)), sum) > 0)
    nonzeros <- match(gsub("-[0-9]*$", "", rownames(subset)), names(is.nonzero))
    subset <- subset[!is.na(nonzeros), ]

    samples.order <- by(as.numeric(gsub("Infancy", "99", gsub("Mother", "0", gsub(".*-", "", rownames(subset))))), gsub("-Infancy", "", gsub("-Mother", "", gsub("-[0-9]*$", "", rownames(subset)))), function(x) order(x, decreasing = TRUE))
    order.by.timepoint <- match(names(samples.order), gsub("-.*$", "", rownames(subset)))[match(gsub("-.*$", "", rownames(subset)), names(samples.order))] + c(unlist(samples.order)) - 1

    sample.changes <- which(diff(rev(match(names(samples.order), gsub("-.*$", "", rownames(subset)))[match(gsub("-.*$", "", rownames(subset)), names(samples.order))])) != 0)

    new.order <- order.by.timepoint

    ##colors <- read.table("K_mic_colors.tsv", comment.char='@')
    ##colors <- colors[match(gsub("K_mic_", "", colnames(subset)), paste(colors[, 2], colors[, 1], sep='_')), 3]

    if (length(new.order) > 0) {
        subset <- subset[new.order, ]
        image(1:nrow(subset), 1:ncol(subset), subset[nrow(subset):1, ], col = rev(terrain.colors(60)), axes = FALSE, xlab ='', ylab = '', zlim = c(0, 1))
        axis(3, 1:nrow(subset), rev(rownames(subset)), las = 2, cex.axis = 2  )
        axis(2, 1:ncol(subset), colnames(subset), cex.axis = 1.5, las = 2)
        title(main = main.title, adj = 0.01, outer = TRUE, line = -2.5, cex.main = 2.5)
        abline(v = rev(sample.changes) + 0.5, col = "white", lty = 'solid', lwd=4)
        ##abline(h = 1:ncol(subset), col = paste(colors, "66", sep=''), lty = 'dashed')
    } else {
        plot.new()
    }
}

PlotByHospital <- function(plotting.data, full.data, cohort.name, metadata,
                           species.abbrv, species.name) {
    cols <- unique(plotting.data[, 2])
    individuals <- full.data[full.data[, 1] %in% metadata[metadata[, 2] == cohort.name, 1], ]
    order.in.meta <- match(plotting.data[, 1], metadata[, 1])
    rows <- cbind(metadata[order.in.meta, 3], metadata[order.in.meta, 4])
    row.order <- order(rows[, 1])
    vals.ordered <- cbind(plotting.data[row.order, ], rows[row.order, 1], rows[row.order, 2])
    abundances <- matrix(0, nrow = nrow(rows), ncol = length(cols))
    rownames(abundances) <- paste(vals.ordered[, 13], vals.ordered[, 14], sep='-')
    colnames(abundances) <- cols
    for (i in 1:nrow(rows)) {
        abundances[paste(vals.ordered[i, 13], vals.ordered[i, 14], sep='-'), vals.ordered[i, 2]] <- vals.ordered[i, 3]
}
    subset <- abundances[, grepl(species.abbrv, colnames(abundances))]

    hospital.b <- grepl("^B", rownames(subset))
    hospital.c <- grepl("^C", rownames(subset))
    hospital.a <- !hospital.b & !hospital.c

    lineages <- gsub(species.abbrv, "", colnames(subset))
    if (any(grepl("ST", lineages))) {
        sts <- gsub("SC[0-9]*_", "", lineages)
        sts <- gsub("-.*$", "", sts)
        colnames(subset) <- sts
        sts <- gsub("NA", "ST9999999", sts)
        st.order <- order(as.numeric(gsub("^ST([0-9]*).*$", "\\1", sts)))
    } else {
        scs <- gsub("ST[0-9]*", "", lineages)
        scs <- gsub("_NA", "", scs)
        colnames(subset) <- scs
        st.order <- order(as.numeric(gsub("SC", "", scs)))
    }        

    subset <- subset[, rev(st.order)]

    layout(matrix(1:3, nrow = 3))
    PlotLongitudinal(subset[hospital.a, ], species.abbrv, paste(species.name, " - ", cohort.name, " delivery cohort", sep=''))
    title(main = "Hospital A", adj = 0.003, font.main = 3, cex.main = 1.8, outer = TRUE, line = -12)
    PlotLongitudinal(subset[hospital.b, ], species.abbrv, paste(species.name, " - ", cohort.name, " delivery cohort", sep=''))
    title(main = "Hospital B", adj = 0.003, font.main = 3, cex.main = 1.8, outer = TRUE, line = -75)
    PlotLongitudinal(subset[hospital.c, ], species.abbrv, paste(species.name, " - ", cohort.name, " delivery cohort", sep=''))
    title(main = "Hospital C", adj = 0.003, font.main = 3, cex.main = 1.8, outer = TRUE, line = -138)
}
