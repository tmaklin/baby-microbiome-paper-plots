
library("plot.matrix")
library("psych")
library("dichromat")

PlotTransitionMatrix <- function(transition.matrix, title.main, title.adj, title.line, y.lab.pos, n.max, ColorFunc) {##correlations, obs.counts, otu.names, title.main, title.adj, ColorFunc) {
    n.states <- nrow(transition.matrix)
    state.names <- colnames(transition.matrix)
    state.names <- gsub("_", " ", state.names)

    plot(0, type = 'n', xlim = c(0, n.states), ylim = c(0, n.states),
         xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')

    ## top axis, labels at 45 degree angle
    axis(side = 3, at = 1:n.states, labels = rep("", n.states))##, labels = state.names, las = 2)
    text(x = 1:n.states -0.2, y = y.lab.pos, labels = state.names, srt = 65, pos = 4, xpd = TRUE, adj = 0)
    abline(v = 1:n.states, col = "gray90")

    ## left axis, labels in reverse order to match top axis
    axis(side = 2, at = 1:n.states, labels = rev(state.names), las = 2)
    abline(h = 1:n.states, col = "gray90")

    ## Plot the values
    for (i in 1:n.states) {
        for (j in 1:n.states) {
            if (transition.matrix[i, j] > 0) {
                coord.x <- i
                coord.y <- n.states - j + 1 ## filled from top->down

                ## assume transition counts range in [1, n.max]
                color.id <- transition.matrix[i, j]

                ##size <- log(obs.counts[i, j] + 1, base = 6) + 1
                size <- 2
                lines(x = coord.x, y = coord.y, type = 'p', pch = 19,
                      cex = size, col = ColorFunc(n.max)[color.id])
            }
        }
    }
    title(main = title.main, line = title.line, adj = title.adj, outer = TRUE)
}

IntensityLegend <- function(ColorFunc, n.max) {
    par(xpd = TRUE)
    plot(c(0,2),c(0,10),type = 'n', axes = F,xlab = '', ylab = '')
##    title(main = "Transitions", font.main = 1, cex.main = 2,
##          line = -25.5, adj = 0.01)
    legend("left", legend = seq(1, n.max, l = n.max), col = ColorFunc(n.max), bty = 'n', pch = 19, cex = 3.5)
    par(xpd = FALSE)
}
