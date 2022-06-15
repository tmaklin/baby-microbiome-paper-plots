
PlotInvasiveness <- function(plotting.data) {
    plot(0, 0, pch = '', xlim = c(-1.0, 7.5), ylim = c(1, nrow(plotting.data)), xaxt = 'n', yaxt = 'n', xlab = 'Odds ratio', ylab = '', bty = 'n')
    axis(side = 1, at = seq(0.0, 7.0, by = 1.0), labels = c("0.0", "1.0", "2.0", "3.0", "4.0", "5.0", "6.0", "7.0"))
    axis(side = 2, at = 1:nrow(plotting.data), labels = plotting.data$Label, las = 2)
    abline(v = 1.0, lty = "dotted", col = "gray40")
    for (i in 1:nrow(plotting.data)) {
        linecol <- "gray80"
        if (plotting.data$p.adj.BH[i] < 0.05) {
            linecol <- ifelse(plotting.data$lower[i] > 1.0, "#d73027", "#92c5de")
        }
        x.coords <- c(plotting.data$lower[i], plotting.data$upper[i])
        y.coords <- c(i, i)
        lines(x = x.coords, y = y.coords, col = linecol, lwd = 2)
        lines(x = plotting.data$OR[i], y = i, type = 'p', pch = 19, col = linecol)
    }
}
