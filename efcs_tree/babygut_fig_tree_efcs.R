#R script for creating E.faecalis tree figure for the baby gut microbiome project
#R v.4.2.0
#AKP 05/2022

#input 
##metadata (data_efcs.csv)
##phylogeny (poppunk_visualise_core_NJ.nwk) 

setwd("/Users/annapontinen/Desktop/baby_gut/efcs/figure_tree/")

#install packages if needed
install.packages("ggplot2")
install.packages("tidyverse")
#install.packages("ggtree")
#btw install.packages("ggtree") not available for all the latest R versions - if so, use instead e.g.
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")

#add libraries
library(ggplot2)
library(tidyverse)
library(ggtree)


#import metadata
metadata <- read.csv("data_efcs.csv", header = T, stringsAsFactors = F, sep = ",", row.names = NULL)
colnames(metadata)[1] <- 'Strain'

delivery_mode <- data.frame("delivery_mode" = metadata[,c("Delivery_mode")])
rownames(delivery_mode) <- metadata$Strain


#import phylogeny
tree <- read.tree('poppunk_visualise_core_NJ.nwk')


#nodes for 10 largest STs
#ST6
#MRCA(tree, tip=c('26975_2#149', '26975_1#113'))
MRCA(as_tibble(tree), '26975_2#149', '26975_1#113')
#2057

#ST9
#MRCA(tree, tip=c("28832_2#265","26975_1#220"))
MRCA(as_tibble(tree),"28832_2#265","26975_1#220")
#2609

#ST16
#MRCA(tree, tip=c("28099_2#365","28157_4#159"))
MRCA(as_tibble(tree),"28099_2#365","28157_4#159")
#3260

#ST21
#MRCA(tree, tip=c("26975_1#14","28157_4#384"))
MRCA(as_tibble(tree),"26975_1#14","28157_4#384")
#2805
#ST21_2
MRCA(as_tibble(tree),"ERR3404727_E_fcs_Pop2","ERR3405388_E_fcs_Pop2")
#2751
#2746 (node for ST21 cladelabel)

#ST25
#MRCA(tree, tip=c("ERR3405451_E_fcs_Pop5","28832_2#27"))
MRCA(as_tibble(tree),"ERR3405451_E_fcs_Pop5","28832_2#27")
#3010
#ST25_2
#MRCA(tree, tip=c("27688_1#68","28832_2#155"))
MRCA(as_tibble(tree),"27688_1#68","28832_2#155")
#3122
#3005 (node for ST25 cladelabel)

#ST28
#MRCA(tree, tip=c("27688_1#92","ERR3405607_E_fcs_Pop9"))
MRCA(as_tibble(tree),"27688_1#92","ERR3405607_E_fcs_Pop9")
#3175

#ST40
#MRCA(tree, tip=c("26975_2#20","ERR3405645_E_fcs_Pop4"))
MRCA(as_tibble(tree),"26975_2#20","ERR3405645_E_fcs_Pop4")
#2343

#ST55
#MRCA(tree, tip=c("27725_1#299","26975_1#34"))
MRCA(as_tibble(tree),"27725_1#299","26975_1#34")
#2509

#ST179
#MRCA(tree, tip=c("ERR3405142_E_fcs_Pop1","ERR3404659_E_fcs_Pop1"))
MRCA(as_tibble(tree),"ERR3405142_E_fcs_Pop1","ERR3404659_E_fcs_Pop1")
#3432

#ST191
#MRCA(tree, tip=c("ERR3404889_E_fcs_Pop21","ERR3406145_E_fcs_Pop21"))
MRCA(as_tibble(tree),"ERR3404889_E_fcs_Pop21","ERR3406145_E_fcs_Pop21")
#2696

#if want to check a clade, e.g.
#png("clade.png", units="in", width=8, height=25, res=300)
#viewClade(clades_tree, node = 2696)
#dev.off() 


#highlight major STs - hospital-associated STs from https://doi.org/10.1038/s41467-021-21749-5 are coloured in red
highlight <- ggtree(tree, layout = 'circular') +
geom_hilight(node=2057, fill="red", alpha=0.5) +
geom_hilight(node=2609, fill="red", alpha=0.5) +
geom_hilight(node=3260, fill="blue", alpha=0.5) +
geom_hilight(node=2805, fill="blue", alpha=0.5) +
geom_hilight(node=2751, fill="blue", alpha=0.5) +
geom_hilight(node=3010 , fill="blue", alpha=0.5) +
geom_hilight(node=3122 , fill="blue", alpha=0.5) +
geom_hilight(node=3175, fill="red", alpha=0.5) +
geom_hilight(node=2343, fill="blue", alpha=0.5) +
geom_hilight(node=2509, fill="blue", alpha=0.5) +
geom_hilight(node=3432, fill="blue", alpha=0.5) +
geom_hilight(node=2696, fill="blue", alpha=0.5) +
geom_cladelabel(node=2057, label="ST6",
                  color='black', fontsize=3.5,  offset.text=0.001) +
geom_cladelabel(node=2609, label="ST9",  
                  color='black', fontsize=3.5,  offset.text=0.003) +
geom_cladelabel(node=3260, label="ST16", 
                  color='black', fontsize=3.5, offset=0.00001, offset.text=0.0045) +
geom_cladelabel(node=2746, label="ST21", 
                  color='black', fontsize=3.5,  offset.text=0.003) +
geom_cladelabel(node=3005, label="ST25", 
                  color='black', fontsize=3.5,  offset.text=0.0035) +
geom_cladelabel(node=3175, label="ST28", 
                  color='black', fontsize=3.5,  offset.text=0.004) +
geom_cladelabel(node=2343, label="ST40",
                  color='black', fontsize=3.5,  offset.text=0.001) +
geom_cladelabel(node=2509, label="ST55", 
                  color='black', fontsize=3.5,  offset.text=0.001) +
geom_cladelabel(node=3432, label="ST179", 
                  color='black', fontsize=3.5,  offset.text=0.0045) +
geom_cladelabel(node=2696, label="ST191", 
                  color='black', fontsize=3.5,  offset.text=0.0025)

#add heatmap with delivery mode
heatmap <-  gheatmap(highlight, delivery_mode,                               
                width = 0.2, colnames = FALSE, offset = 0.004) +                               
  scale_fill_manual(name = "Delivery mode",                       
                    values = c("#92c5de", "#f4a582"),
                    breaks = c("Caesarean", "Vaginal"),
                    labels = c("Caesarean", "Vaginal")) + theme(legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.box = "horizontal", legend.margin = margin())

#save as png
png("babygut_efcs_draft_STs2.png", units="in", width=8, height=6, res=300)
heatmap
dev.off() 