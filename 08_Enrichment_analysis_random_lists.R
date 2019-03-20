#########################
#   ENRICHR ANALYSIS    #
#########################

rm(list = ls())

# Install package
#library(devtools)
#install_github("wjawaid/enrichR")
# Or:
# install.packages("enrichR")

# Load library
library(enrichR)

# Define dirs
path <- "/Users/alejandromendozaalvarez/Desktop"
setwd(path)

# Retrieve the list of ddbb
dbs <- listEnrichrDbs()
dbs
head(dbs)

# Define ddbb to use
dbs <- "Jensen_DISEASES"

# Load and check data
data <- read.csv(file="mortality_138_genes_enrichr_input.txt", header = FALSE)
head(data)

# Transpose and convert into a vector
genes <- as.vector(t(data))
str(genes)

# Define global parameters
num.iterations <- 1000

# Define the range of genes in lists
inicio <- 30
fin <- 100

# Loop using list with 30 to 100 genes per list
# Per each number of genes (30, 31... 99, 100), extract the random list of genes, query EnrichR and save results
for (i in inicio:fin) {
  # Repeat a number of iterations
  # Create a dir to output results
  outdir <- paste(path,"/",i,"genes",sep="")
  dir.create(outdir)
  
  # Show messages
  message (paste("Outdir: ", outdir, sep=""))
  message (paste("Number of genes in list: ", i, sep=""))
  print ("*****************************************************************************************")
  
  for (j in 1:num.iterations) {
    
    # Create a random list of genes
    message (paste("Iteration: ", j, sep=""))
    message ("Creating list of genes...")
    list <- sample(genes, size = i, replace = FALSE, prob = NULL)
    #listfile <- paste(outdir,"/","list_", i, "genes","_iteration_", j, sep="")
    #write.table(list, listfile, sep="\t")
    
    # Query genes into EnrichR
    message ("Quering EnrichR with this list...")
    query <- enrichr(list, dbs)
    # str(quer)
    
    # Inspect results
    # head(query)
    
    # Save results using tab as separator (select certain cols if needed)
    message ("Saving EnrichR query results...")
    file <- paste(outdir,"/","enrichR_list_", i, "genes","_iteration_", j, sep="")
    #printEnrich(query, file, sep = "\t", columns = c(1, 2, 3, 4, 5, 6, 7, 8, 9))
    write.table(query, file, sep="\t")
    print ("-----------------------------------------------------------------------------------------")
    # End of j iteration
  }
  # End of a number of genes in the list
}
# End Of Script
