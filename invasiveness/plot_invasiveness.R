
odds.ratios <- read.table("NORM_OR_E_coli_carriage_disease.csv", sep = ',', header = TRUE)

in.plot <- odds.ratios$Labe != ""

plotting.data <- odds.ratios[in.plot, ]

labels.order <- order(as.numeric(gsub("[*]", "", gsub("ST", "", plotting.data$Label))), decreasing = TRUE)

plotting.data <- plotting.data[labels.order, ]

pdf(file = "invasiveness_vs_carriage.pdf", width = 5, height = 5)
par(mar = c(4.5, 4, 0.5, 0.5))
plot(0, 0, pch = '', xlim = c(-1.0, 7.5), ylim = c(1, nrow(plotting.data)), xaxt = 'n', yaxt = 'n', xlab = 'Odds ratio', ylab = '', bty = 'n')
axis(side = 1, at = seq(0.0, 7.0, by = 1.0), labels = c("0.0", "1.0", "2.0", "3.0", "4.0", "5.0", "6.0", "7.0"))
axis(side = 2, at = 1:nrow(plotting.data), labels = plotting.data$Label, las = 2)
abline(v = 1.0, lty = "dotted", col = "gray40")
for (i in 1:nrow(plotting.data)) {
    x.coords <- c(plotting.data$lower[i], plotting.data$upper[i])
    y.coords <- c(i, i)
    linecol <- "gray80"
    if (plotting.data$p.adj.BH[i] < 0.05) {
        linecol <- ifelse(plotting.data$lower[i] > 1.0, "#d73027", "#4575b4")
    }
    lines(x = x.coords, y = y.coords, col = linecol, lwd = 2)
}
lines(x = plotting.data$OR, y = 1:nrow(plotting.data), type = 'p', pch = 19)
legend("bottomright", legend = c("Odds ratio", "95% confidence interval", "Invasive", "Opportunistic", "Commensal"), col = c("black", "black", "#d73027", "gray80", "#4575b4"), bty = 'n', lty = c(NA, 1, 1, 1, 1), lwd = c(NA, 2, 2, 2, 2), pch = c(19, NA, NA, NA, NA))
dev.off()
