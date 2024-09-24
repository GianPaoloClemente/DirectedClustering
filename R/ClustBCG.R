#' Clustering Coefficient for Directed/Undirected and Weighted Networks
#'
#' Compute Local and Global (average) Clustering Coefficients
#' for Directed/Undirected and Unweighted/Weighted Networks.
#'
#' Formulas are based on Barrat et al. (2004) for undirected networks,
#' and on Clemente and Grassi (2018) for directed networks.
#'
#' In the directed case, different components of the directed clustering coefficient are also provided.
#'
#' @param mat A weighted adjacency matrix.
#' @param type The type of clustering coefficient to calculate.
#' Possible values are: \code{"undirected"} (default) or \code{"directed"}.
#' @param isolates Character scalar, defines how to treat vertices with degree zero and one.
#' If \code{"NaN"}, their local transitivity is reported as NaN and they are not included in the averaging. If \code{"zero"}, their transitivity is reported as 0 and they are included in the averaging. Default is \code{"zero"}.
#'
#' @details The function computes the Barrat et al. (2004) coefficient for a weighted and undirected network.
#' For a directed network, the Clemente and Grassi (2018) formula is used.
#' In case of unweighted and undirected graphs, the classical local clustering coefficient (Watts and Strogatz) is provided.
#' Local clustering coefficients are computed for each node,
#' and the global coefficient is the average of these local coefficients.
#' These coefficients do not work for graphs with multiple or loop edges, hence loops are removed.
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
#' Barrat, A., Barthelemy, M., Pastor-Satorras, R., & Vespignani, A. (2004). The architecture of complex weighted networks.
#' \emph{Proceedings of the National Academy of Sciences}, USA, 101, 3747.
#'
#' Clemente, G.P., & Grassi, R. (2018). Directed clustering in weighted networks: a new perspective.
#' \emph{Chaos, Solitons and Fractals}, 107, 26–38.
#'
#' Watts, D.J., & Strogatz, S.H. (1998). Collective dynamics of 'small-world' networks.
#' \emph{Nature}, 393, 440-442.
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
#'   BarratClust <- ClustBCG(A, "undirected")
#'   check <- sum(BarratClust$LocalCC - transitivity(gsim, "weighted"))
#'
#'   # Generate a weighted and directed graph
#'   gsim <- sample_gnp(50, 0.5, directed = TRUE, loops = FALSE)
#'   PESI <- runif(length(E(gsim)), 0, 1)
#'   E(gsim)$weight <- PESI
#'   A <- as_adjacency_matrix(gsim, sparse = FALSE, attr = "weight")
#'   CGClust <- ClustBCG(A, "directed")
#' } else {
#'   cat("Please install the 'igraph' package to run this example.\n")
#' }
#'
#' @export


ClustBCG <- function(mat, type = "undirected", isolates = "zero") {
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
        mat1[mat1 != 0] = 1  #Per ottenere matrice senza pesi per reti pesate
    } else {
        mat1 = mat
    }

    if (any(diag(mat) != 0)) {
        diag(mat) = 0
        diag(mat1) = 0
        warning("Loops are not considered")
    }

    if (type == "undirected") {
        c <- (diag(mat %*% mat1 %*% mat1))/(rowSums(mat) * (rowSums(mat1) - 1))
        if (isolates == "zero") {
            c[is.na(c) == TRUE] = 0
        }
        list(LocalCC = c, GlobalCC = mean(c, na.rm = TRUE))
    } else {
        degin = t(mat1) %*% matrix(1, dim(mat1)[1], 1)
        degout = mat1 %*% matrix(1, dim(mat1)[1], 1)
        dtot = (t(mat1) + mat1) %*% matrix(1, dim(mat1)[1], 1)
        dbil = diag(mat1 %*% mat1)
        sin = diag(t(mat1) %*% mat)
        sout = diag(mat1 %*% t(mat))
        stot = sin + sout
        sbil = diag(mat %*% mat1 + mat1 %*% mat)/2

        cyc = (diag((mat %*% mat1 %*% mat1) + (t(mat) %*% t(mat1) %*% t(mat1)))/2)/(0.5 * (sin * degout + sout *
            degin) - sbil)
        mid = (0.5 * diag((t(mat) %*% mat1 %*% t(mat1)) + (mat %*% t(mat1) %*% mat1)))/(0.5 * (sin * degout +
            sout * degin) - sbil)
        incl = (0.5 * diag(t(mat) %*% (mat1 + t(mat1)) %*% mat1))/(sin * (degin - 1))
        out = (0.5 * diag(mat %*% (mat1 + t(mat1)) %*% t(mat1)))/(sout * (degout - 1))
        tot = (0.5 * diag((mat + t(mat)) %*% (mat1 + t(mat1)) %*% (mat1 + t(mat1))))/(stot * (dtot - 1) - 2 *
            sbil)
        if (isolates == "zero") {
            cyc[is.na(cyc) == TRUE] = 0
            mid[is.na(mid) == TRUE] = 0
            incl[is.na(incl) == TRUE] = 0
            out[is.na(out) == TRUE] = 0
            tot[is.na(tot) == TRUE] = 0
        }
        list(cycleCC = cyc, middlemanCC = mid, inCC = incl, outCC = out, totalCC = tot, GlobalcycleCC = mean(cyc,
            na.rm = TRUE), GlobalmiddlemanCC = mean(mid, na.rm = TRUE), GlobalinCC = mean(incl, na.rm = TRUE),
            GlobaloutCC = mean(out, na.rm = TRUE), GlobaltotalCC = mean(tot, na.rm = TRUE))
    }
}
