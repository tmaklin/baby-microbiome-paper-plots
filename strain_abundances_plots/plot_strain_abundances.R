
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

source("wgs_read_bin.R")
source("wgs_plotting.R")

## Read in the metadata (accession, delivery_mode, individual, timepoint)
metadata <- read.table("../wgs_meta_delivery.tsv", header=TRUE)

## Read genome lengths from another file
source("genome_lengths.R")

## Extract unique timepoints and sort them in (Mother, Infancy), ascending in time. 
timepoints <- unique(metadata$Time_point)
timepoints.order <- RenameTimePointsForSorting(timepoints)
timepoints.order <- order(timepoints.order)
timepoints <- timepoints[timepoints.order]

abundances.files <- paste("../data/strain_abundances", list.files("../data/strain_abundances"), sep = '/')

abundances.files <- abundances.files[grepl("_abundances.txt", abundances.files)]

species.abundances.files <- paste("../data/species_abundances", list.files("../data/species_abundances"), sep ='/')
species.abundances <- ReadAbundances(species.abundances.files, "ERR")

species.to.plot <- c("E_col","E_fcm","E_fcs","K_aer","K_gri","K_hua","K_mic","K_orn","K_oxy","K_pas","K_pla","K_pne","K_qps","K_var","P_aer","S_aur","S_pne")

for (i in species.to.plot) {
    ecol.abundances <- ReadAbundances(abundances.files, i, 1, 1)
    ecol.qc.filter <- FilterByCoverage(species.abundances, ecol.abundances, i, genome.lens[["E_col"]], 1)
    has.new.names <- any(grepl(i, list.files("../poppunk_cluster_info/")))
    print(i)
    if (has.new.names) {
        ecol.new.names <- read.table(paste("../poppunk_cluster_info/", i, "_new_clusters.tsv", sep = ''), sep = '\t', header = TRUE)
        colnames(ecol.abundances$abundances) <- RenameClusters(colnames(ecol.abundances$abundances), ecol.new.names)
    }
    ## widths <- c(20, 20, 24, 24, 20, 20, 20, 20, 24, 20, 20, 20, 30, 20, 20, 20, 30, 24, 20)
    ## heights <- c(16, 16, 24, 24, 12, 12, 12, 12, 12, 12, 12, 12, 18, 12, 12, 12, 18, 16, 12)
    ## params <- cbind(c(2.5, 2.5, 1.5, 2.3, 2.3, 2.3, 2.5, 2.3, 1.9, 2.4, 2.5, 2.5, 1.2, 2, 2, 2, 1.9, 1.2, 2.5)
    ##               , c(8, 8, 14, 10, 17, 15, 4, 15, 19, 15, 15, 15, 15, 6, 6, 6, 10, 17, 8))
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    species.colors <- sample(color, ncol(ecol.abundances$abundances), replace = TRUE)
    source("wgs_plotting.R")
    pdf(file = paste(i, "_strain_abundances.pdf", sep = ''), width = 20, height = 20)
    dump <- PlotBin(ecol.abundances$abundances, species.colors, "abcdef", paste(i, "_", sep = ''), c("Mother", "4", "7", "21", "Infancy"), metadata, 1.6, 10)
    dev.off()
}
