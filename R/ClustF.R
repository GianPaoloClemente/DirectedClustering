ClustF <- function(mat, type, isolates = "zero", norm = 1) {
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
