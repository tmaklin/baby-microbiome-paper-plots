
source("plot_subplot.R")

odds.ratios <- read.table("NORM_OR_E_coli_carriage_disease.csv", sep = ',', header = TRUE)
in.plot <- odds.ratios$Labe != ""
plotting.data <- odds.ratios[in.plot, ]
plotting.data <- plotting.data[grepl("NORM", plotting.data$Label), ]
plotting.data$Label <- gsub("_SC[0-9]*[[:space:]]*\\(NORM\\)", "", plotting.data$Label)
labels.order <- order(as.numeric(gsub("[*]", "", gsub("ST", "", plotting.data$Label))), decreasing = TRUE)
plotting.data.all <- plotting.data[labels.order, ]
odds.ratios <- read.table("ST131_clades_OR_E_coli_carriage_disease_collapsed.csv", sep = ',', header = TRUE)
in.plot <- odds.ratios$Labe != ""
plotting.data <- odds.ratios[in.plot, ]
labels.order <- order(as.numeric(gsub("[*]", "", gsub("ST", "", plotting.data$Label))), decreasing = TRUE)
plotting.data <- plotting.data[labels.order, ]
plotting.data.st131 <- plotting.data[c(12:1, 14, 13, 15), ]

pdf(file = "invasiveness_vs_carriage.pdf", width = 8, height = 5)
layout(matrix(c(1, 2, 3, 3), byrow = TRUE, nrow = 2, ncol = 2), height = c(6, 1))
par(mar = c(4.5, 4, 0.5, 2.5))
PlotInvasiveness(plotting.data.all)
par(mar = c(4.5, 6, 0.5, 0.5))
PlotInvasiveness(plotting.data.st131)
par(mar = c(0, 0, 0, 0))
plot(0, bty ='n', xaxt='n', yaxt='n', xlab='', ylab='', pch = '')
legend("top", legend = c("Odds ratio", "95% confidence interval", "", "Invasive", "Intermediate", "Commensal"), col = c("black", "black", "black", "#d73027", "gray80", "#92c5de"), bty = 'n', lty = c(NA, 1, NA, 1, 1, 1), lwd = c(NA, 2, NA, 2, 2, 2), pch = c(19, NA, NA, NA, NA, NA), ncol = 2)
dev.off()
