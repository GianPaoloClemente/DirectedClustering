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
