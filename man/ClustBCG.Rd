\name{ClustBCG}
\alias{ClustBCG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clustering Coefficient for Directed/Undirected and Weighted Networks (Barrat et al. (2004) and Clemente, Grassi (2018) coefficients).
}
\description{
Compute Local and Global (average) Clustering Coefficients for Directed/Undirected and Unweighted/Weighted Networks. \cr Formulas are based on Barrat et al. (2004) coefficient when the network is undirected, while it is based on Clemente and Grassi (2018) proposal when the network is directed. \cr
In the directed case, different components of directed clustering coefficient are also provided.
}
\usage{
ClustBCG(mat, type = "undirected", isolates = "zero")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{A Weighted Adjacency Matrix}

  \item{type}{The type of the clustering coefficient to calculate. Possible values are: "undirected" (default) or "undirected".}

  \item{isolates}{Character scalar, defines how to treat vertices with degree zero and one. If it is 'NaN' then their local transitivity is reported as NaN and they are not included in the averaging. If it is 'zero', then we report 0 transitivity for them, and they are included in the averaging. Default value is 'zero'.}
}
\details{
The function \code{ClustBCG} computes Barrat et al. coefficient when a weighted and undirected network is considered. \cr
For directed network Clemente and Grassi formula is computed. \cr
In case of unweighted and undirected graphs, it provides classical local clustering coefficient (Watts and Strogatz).

In all cases, local clustering coefficients are obtained for each node, the global coefficient is the average of local coefficients. \cr
These clustering coefficients do not work for graphs with multiple and/or loop edges. Hence, loops are removed. \cr
In the directed case, different components of directed clustering coefficient are also considered.
}
\value{
For 'undirected case':
\item{LocalCC}{Local clustering coefficients}
\item{GlobalCC}{Global clustering coefficient}
For 'directed case':
\item{cycleCC}{Local Cycle clustering coefficients}
\item{middlemanCC}{Local Middleman clustering coefficients}
\item{inCC}{Local In clustering coefficients}
\item{outCC}{Local Out clustering coefficients}
\item{totalCC}{Local Total clustering coefficients}
\item{GlobalcycleCC}{Global Cycle clustering coefficient}
\item{GlobalmiddlemanCC}{Global Middleman clustering coefficient}
\item{GlobalinCC}{Global In clustering coefficient}
\item{GlobaloutCC}{Global Out clustering coefficient}
\item{GlobaltotalCC}{Global Total clustering coefficient}

}
\references{
Barrat, A. and Barthelemy, M. and Pastor-Satorras, R. and Vespignani, A. (2004) The architecture of complex weighted networks, \emph{Proc. Natl. Acad. Sci.}, USA \bold{101}, 3747.

Clemente, G.P. and Grassi, R. (2018) Directed clustering in weighted networks: a new perspective, \emph{Chaos, Solitons and Fractals}, \strong{107},26--38.

Watts, D.J. and Strogatz, S.H. (1998) Collective dynamics of 'small-world 'networks. \emph{Nature}, \bold{393}, 440-442.
}
\author{
Gian Paolo Clemente, gianpaolo.clemente@unicatt.it
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
if (requireNamespace("igraph", quietly = TRUE)) {

library(igraph)
#Generate a Weighted and Undirected graph with Erdos Renyi Model
gsim<-sample_gnp(50, 0.5, directed = FALSE, loops = FALSE)
PESI<-runif(length(E(gsim)), 0, 1)
E(gsim)$weight<-PESI

#Get Adjacency
A<-as_adjacency_matrix(gsim, sparse=FALSE, attr="weight")
#Compute Barrat et al. (2004) coefficient
BarratClust<-ClustBCG(A, "undirected")

#The same results can be obtained with igraph::transitivity
check<-sum(BarratClust$LocalCC-transitivity(gsim, "weighted"))

#Generate a Weighted and Directed Graph with Erdos Renyi Model
gsim<-sample_gnp(50, 0.5, directed = TRUE, loops = FALSE)
PESI<-runif(length(E(gsim)), 0, 1)
E(gsim)$weight<-PESI

#Get Adjacency
A<-as_adjacency_matrix(gsim, sparse=FALSE, attr="weight")

#Compute Clemente, Grassi (2018) coefficient
CGClust<-ClustBCG(A, "directed")
} else {
    cat("Please install the 'igraph' package to run this example.\n")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ graphs}% use one of  RShowDoc("KEYWORDS")
