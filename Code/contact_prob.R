#module load R/4.2.2-foss-2022b

#pae_matrix.R "plot title" "plot name"
# Rscript pae_matrix.R 5fvm-TORdimmer.csv "TOR-LST8 tetramer PAE matrix" pae_TOR.svg

#install.packages("cli", repos='http://cran.us.r-project.org')
my_packages <- c("cli","ggplot2","pheatmap","pak") 
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])] 
if(length(not_installed)) install.packages(not_installed, repos='http://cran.us.r-project.org')

#library(pak)
#pak::pak("r-lib/gtable")

lapply(my_packages, library, character.only = TRUE)

args<-commandArgs(TRUE)


PAEmatrix <- as.matrix(read.csv(args[1], header = FALSE))

colnames(PAEmatrix) <-  gsub("V", "",colnames(PAEmatrix))
 num_residues <- nrow(PAEmatrix)
 labels <- rep("", num_residues)  # Empty labels by default
  labels[seq(1, num_residues, by = 10)] <- seq(1, num_residues, by = 10)  # Show only every 10th residue


 Matrix_plot <- pheatmap(PAEmatrix,
           display_numbers = FALSE,
#             color = colorRampPalette(c("#006400", "#a1ffa1","white"))(100),
color = colorRampPalette(c("white","gray","black"))(100),
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
	   labels_row = labels,
           labels_col = labels,
           main = args[2] )

ggsave(filename = args[3], plot = Matrix_plot,device = "png", width=10, height=10)
