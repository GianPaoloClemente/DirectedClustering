\name{ClustF}
\alias{ClustF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clustering Coefficients for Directed/Undirected and Weighted Networks (Onnela et al. (2005) and Fagiolo (2007) coefficients)
}
\description{
This function computes both Local and Global (average) Clustering Coefficients for either Directed/Undirected and Unweighted/Weighted Networks. \cr Formulas are based on Onnela et al. (2005) coefficient when the network is undirected, while it is based on Fagiolo (2007) coefficient when the network is directed. \cr
In the directed case, different components of directed clustering coefficient are also considered.
}
\usage{
ClustF(mat, type = "undirected", isolates = "zero", norm=1)
}
\arguments{
  \item{mat}{A Weighted Adjacency Matrix. If weights are greater than one, a normalization is provided by dividing each weight by the maximum weight observed.}
  \item{type}{The type of the clustering coefficient to calculate. Possible values are:
  "undirected" (default) or "undirected".}
  \item{isolates}{Character scalar, defines how to treat vertices with degree zero and one. If it is 'NaN' then their local transitivity is reported as NaN and they are not included in the averaging. If it is 'zero', then we report 0 transitivity for them, and they are included in the averaging. Default value is 'zero'.}
  \item{norm}{If it is 1 (default), link's weights are normalized by dividing by the maximum observed weight (as proposed in Fagiolo). Otherwise, weights are not normalized. It is worth pointing out that weights are always normalized when the maximum weight is greater than zero. This normalization assures that clustering coefficient ranges between 0 and 1.}
}
\details{
The function \code{ClustF} computes Onnela et al. formula when weighted and undirected networks are considered. \cr
For directed networks, Fagiolo formula is computed. In case of unweighted and undirected graphs, it provides classical local clustering coefficient (Watts and Strogatz).

Local coefficients are obtained for each node, the global coefficient is the average of local coefficients. \cr
These clustering coefficients do not work for graphs with multiple and/or loop edges. Hence, loops are removed. \cr
In the directed case, different components of directed clustering coefficient are also provided.
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
Fagiolo, G. (2007) Clustering in complex directed networks. \emph{Physical Review E}, \bold{76(2)}.

Onnela, J.P. and Saramaki, J. and Kertsz, J. and Kaski, K. (2005) Intensity and coherence of motifs in weighted complex networks. \emph{Physical Review E}, \bold{71(6)}.

Watts, D. J. and Strogatz, S. H. (1998), Collective dynamics of 'small-world' networks \emph{Nature}, \bold{393}, 440-442.
}
\author{
Gian Paolo Clemente, gianpaolo.clemente@unicatt.it
}
%\note{
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

#Compute Onnela et al. (2005) coefficient
OnnelaClust<-ClustF(A, "undirected")

#Generate a Weighted and Directed Graph with Erdos Renyi Model
gsim<-sample_gnp(50, 0.5, directed = TRUE, loops = FALSE)
PESI<-runif(length(E(gsim)), 0, 1)
E(gsim)$weight<-PESI

#Get Adjacency
A<-as_adjacency_matrix(gsim, sparse=FALSE, attr="weight")

#Compute Fagiolo (2007) coefficient
FagioloClust<-ClustF(A, "directed")
} else {
    cat("Please install the 'igraph' package to run this example.\n")
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs}% use one of  RShowDoc("KEYWORDS")
