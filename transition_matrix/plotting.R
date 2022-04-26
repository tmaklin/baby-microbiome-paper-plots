
library("plot.matrix")
library("psych")
library("dichromat")

PlotTransitionMatrix <- function(transition.matrix) {
    mat.colors <- c("white", colorRampPalette(c("blue", "red"))(max(transition.matrix)))

    ## Hide states that were never visited
    nonzeros <- rowSums(transition.matrix) > 0 | colSums(transition.matrix) > 0
    no.emptys <- transition.matrix[nonzeros, nonzeros]

    ## Suppress plotting the x-axis and y-axis by passing "las = 4",
    ## which errors the axis plotting after the matrix has been
    ## plotted, since plot.matrix doesn't allow passing the usual
    ## xaxt='n' and yaxt='n' for some reason.
    err <- tryCatch(plot(no.emptys, breaks = seq(0, length(mat.colors), by=1), col = mat.colors, main="", xlab="", ylab="", las = 4), error = function(e) {})

    ## Plot the axes with labels adjusted to horizontal at the left
    ## side for y-axis and vertical at the top for x-axis
    axis(3, at = 1:nrow(no.emptys), labels = colnames(no.emptys), las = 2)
    axis(2, at = 1:nrow(no.emptys), labels = rev(rownames(no.emptys)), las = 1)
}

