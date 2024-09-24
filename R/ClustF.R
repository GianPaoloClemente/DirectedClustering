#' Clustering Coefficients for Directed/Undirected and Weighted Networks
#'
#' This function computes both Local and Global (average) Clustering Coefficients for either Directed/Undirected and Unweighted/Weighted Networks.
#' The formulas are based on Onnela et al. (2005) for undirected networks,
#' and on Fagiolo (2007) for directed networks.
#'
#' In the directed case, different components of the directed clustering coefficient are also considered.
#'
#' @param mat A weighted adjacency matrix. If weights are greater than one,
#' a normalization is provided by dividing each weight by the maximum weight observed.
#' @param type The type of clustering coefficient to calculate.
#' Possible values are: \code{"undirected"} (default) or \code{"directed"}.
#' @param isolates Character scalar, defines how to treat vertices with degree zero and one.
#' If \code{"NaN"}, their local transitivity is reported as NaN and they are not included in the averaging.
#' If \code{"zero"}, their transitivity is reported as 0 and they are included in the averaging. Default is \code{"zero"}.
#' @param norm If it is 1 (default), the link's weights are normalized by dividing by the maximum observed weight (as proposed by Fagiolo).
#' If it is 0, weights are not normalized. Weights are always normalized when the maximum weight is greater than zero, ensuring that the clustering coefficient ranges between 0 and 1.
#'
#' @details The function computes Onnela et al.'s (2005) formula for weighted and undirected networks.
#' For directed networks, Fagiolo's (2007) formula is used.
#' In the case of unweighted and undirected graphs, the classical local clustering coefficient (Watts and Strogatz) is provided.
#' Local coefficients are computed for each node, and the global coefficient is the average of these local coefficients.
#' These coefficients do not work for graphs with multiple or loop edges, so loops are removed.
#'
#' @return A list with the following components:
#' \item{LocalCC}{Local clustering coefficients for undirected networks}
#' \item{GlobalCC}{Global clustering coefficient for undirected networks}
#' \item{cycleCC}{Local Cycle clustering coefficients for directed networks}
#' \item{middlemanCC}{Local Middleman clustering coefficients for directed networks}
#' \item{inCC}{Local In clustering coefficients for directed networks}
#' \item{outCC}{Local Out clustering coefficients for directed networks}
#' \item{totalCC}{Local Total clustering coefficients for directed networks}
#' \item{GlobalcycleCC}{Global Cycle clustering coefficient for directed networks}
#' \item{GlobalmiddlemanCC}{Global Middleman clustering coefficient for directed networks}
#' \item{GlobalinCC}{Global In clustering coefficient for directed networks}
#' \item{GlobaloutCC}{Global Out clustering coefficient for directed networks}
#' \item{GlobaltotalCC}{Global Total clustering coefficient for directed networks}
#'
#' @references
#' Fagiolo, G. (2007). Clustering in complex directed networks. \emph{Physical Review E}, 76(2).
#'
#' Onnela, J.P., Saramaki, J., Kertsz, J., & Kaski, K. (2005). Intensity and coherence of motifs in weighted complex networks. \emph{Physical Review E}, 71(6).
#'
#' Watts, D.J., & Strogatz, S.H. (1998). Collective dynamics of 'small-world' networks. \emph{Nature}, 393, 440-442.
#'
#' @author Gian Paolo Clemente, \email{gianpaolo.clemente@unicatt.it}
#'
#' @examples
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   library(igraph)
#'   # Generate a weighted and undirected graph
#'   gsim <- sample_gnp(50, 0.5, directed = FALSE, loops = FALSE)
#'   PESI <- runif(length(E(gsim)), 0, 1)
#'   E(gsim)$weight <- PESI
#'   A <- as_adjacency_matrix(gsim, sparse = FALSE, attr = "weight")
#'   OnnelaClust <- ClustF(A, "undirected")
#'
#'   # Generate a weighted and directed graph
#'   gsim <- sample_gnp(50, 0.5, directed = TRUE, loops = FALSE)
#'   PESI <- runif(length(E(gsim)), 0, 1)
#'   E(gsim)$weight <- PESI
#'   A <- as_adjacency_matrix(gsim, sparse = FALSE, attr = "weight")
#'   FagioloClust <- ClustF(A, "directed")
#' } else {
#'   cat("Please install the 'igraph' package to run this example.\n")
#' }
#'
#' @export

ClustF <- function(mat, type= "undirected", isolates = "zero", norm = 1) {
    if (!is.matrix(mat)) {
        stop("Not a valid matrix object")
    }
    if (type == "undirected" & !isSymmetric(mat)) {
        stop("An asymmetric matrix has been given for an undirected graph")
    }
    if (type == "directed" & isSymmetric(mat)) {
        warning("A symmetric matrix has been given for a directed graph")
    }
    if (isolates != "zero" & isolates != "NaN") {
        isolates == "zero"
    }

    if (any(mat != 0 & mat != 1)) {
        mat1 = mat
        mat1[mat1 != 0] = 1  #Adjacency Matrix for weighted graphs
    } else {
        mat1 = mat
    }

    if (any(diag(mat) != 0)) {
        diag(mat) = 0
        diag(mat1) = 0
        warning("Loops have been removed")
    }
    if(norm==1){
      mat = mat/max(mat)
    }else{
      if (max(mat) > 1) {
        mat = mat/max(mat)
        warning("Weights have been normalized in order to assure that the clustering coefficient is between 0 and 1")
        }
    }
    mat = mat^(1/3)  #Needed for computing coefficient of weighted graphs

    if (type == "undirected") {
        c <- (diag(mat %*% mat %*% mat))/(rowSums(mat1) * (rowSums(mat1) - 1))
        if (isolates == "zero") {
            c[is.na(c) == TRUE] = 0
        }
        list(LocalCC = c, GlobalCC = mean(c, na.rm = TRUE))
    } else {

        degin = t(mat1) %*% matrix(1, dim(mat)[1], 1)
        degout = mat1 %*% matrix(1, dim(mat)[1], 1)
        dtot = (t(mat1) + mat1) %*% matrix(1, dim(mat)[1], 1)
        dbil = diag(mat1 %*% mat1)
        c_cyc = diag(mat %*% mat %*% mat)/(degin * degout - dbil)
        c_mid = diag(mat %*% t(mat) %*% mat)/(degin * degout - dbil)
        c_in = diag(t(mat) %*% (mat %*% mat))/(degin * (degin - 1))
        c_out = diag((mat %*% mat) %*% t(mat))/(degout * (degout - 1))
        c_tot = diag((mat + t(mat)) %*% (mat + t(mat)) %*% (mat + t(mat)))/(2 * (dtot * (dtot - 1) - 2 * dbil))

        if (isolates == "zero") {
            c_cyc[is.na(c_cyc) == TRUE] = 0
            c_mid[is.na(c_mid) == TRUE] = 0
            c_in[is.na(c_in) == TRUE] = 0
            c_out[is.na(c_out) == TRUE] = 0
            c_tot[is.na(c_tot) == TRUE] = 0
        }

        list(cycleCC = c_cyc, middlemanCC = c_mid, inCC = c_in, outCC = c_out, totalCC = c_tot, GlobalcycleCC = mean(c_cyc,
            na.rm = TRUE), GlobalmiddlemanCC = mean(c_mid, na.rm = TRUE), GlobalinCC = mean(c_in, na.rm = TRUE),
            GlobaloutCC = mean(c_out, na.rm = TRUE), GlobaltotalCC = mean(c_tot, na.rm = TRUE))
    }
}
