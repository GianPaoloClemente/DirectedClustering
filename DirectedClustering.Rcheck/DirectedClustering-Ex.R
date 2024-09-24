pkgname <- "DirectedClustering"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "DirectedClustering-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('DirectedClustering')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ClustBCG")
### * ClustBCG

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ClustBCG
### Title: Clustering Coefficient for Directed/Undirected and Weighted
###   Networks
### Aliases: ClustBCG

### ** Examples

if (requireNamespace("igraph", quietly = TRUE)) {
  library(igraph)
  # Generate a weighted and undirected graph
  gsim <- sample_gnp(50, 0.5, directed = FALSE, loops = FALSE)
  PESI <- runif(length(E(gsim)), 0, 1)
  E(gsim)$weight <- PESI
  A <- as_adjacency_matrix(gsim, sparse = FALSE, attr = "weight")
  BarratClust <- ClustBCG(A, "undirected")
  check <- sum(BarratClust$LocalCC - transitivity(gsim, "weighted"))

  # Generate a weighted and directed graph
  gsim <- sample_gnp(50, 0.5, directed = TRUE, loops = FALSE)
  PESI <- runif(length(E(gsim)), 0, 1)
  E(gsim)$weight <- PESI
  A <- as_adjacency_matrix(gsim, sparse = FALSE, attr = "weight")
  CGClust <- ClustBCG(A, "directed")
} else {
  cat("Please install the 'igraph' package to run this example.\n")
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ClustBCG", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ClustF")
### * ClustF

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ClustF
### Title: Clustering Coefficients for Directed/Undirected and Weighted
###   Networks
### Aliases: ClustF

### ** Examples

if (requireNamespace("igraph", quietly = TRUE)) {
  library(igraph)
  # Generate a weighted and undirected graph
  gsim <- sample_gnp(50, 0.5, directed = FALSE, loops = FALSE)
  PESI <- runif(length(E(gsim)), 0, 1)
  E(gsim)$weight <- PESI
  A <- as_adjacency_matrix(gsim, sparse = FALSE, attr = "weight")
  OnnelaClust <- ClustF(A, "undirected")

  # Generate a weighted and directed graph
  gsim <- sample_gnp(50, 0.5, directed = TRUE, loops = FALSE)
  PESI <- runif(length(E(gsim)), 0, 1)
  E(gsim)$weight <- PESI
  A <- as_adjacency_matrix(gsim, sparse = FALSE, attr = "weight")
  FagioloClust <- ClustF(A, "directed")
} else {
  cat("Please install the 'igraph' package to run this example.\n")
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ClustF", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
