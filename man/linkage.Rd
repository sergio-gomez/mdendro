\name{linkage}
\alias{linkage}
\alias{descplot}
\title{Extended Agglomerative Hierarchical Clustering}
\description{
  Agglomerative hierarchical clustering on a dataset of distances or
  similarities, returning a multifurcated dendrogram or \dfn{multidendrogram}.
  Descriptive measures to analyze the resulting dendrogram are additionally
  provided.
}
\usage{
linkage(prox, type.prox = "distance", digits = NULL,
        method = "arithmetic", par.method = 0, weighted = FALSE,
        group = "variable")

descplot(prox, ..., type.prox = "distance", digits = NULL,
         method = "versatile", par.method = c(-1,0,+1), weighted = FALSE,
         group = "variable", measure = "cor", slope = 10)
}
\arguments{
  \item{prox}{A structure of class \code{"dist"} containing non-negative
    proximity data (distances or similarities).}
  \item{type.prox}{A character string to indicate whether the proximity data
    represent \code{"distance"} (default) or \code{"similarity"} between
    objects.}
  \item{digits}{An integer value specifying the precision, i.e. the number of
    significant decimal digits to be used for the comparisons between proximity
    data. This is an important parameter, since equal proximity data at a
    certain precision may become different by increasing its value. Thus, it may
    be responsible of the existence of tied proximity data. If the value of this
    parameter is negative or \code{NULL} (default), then the precision is
    automatically set to that of the input proximity value with the largest
    number of significant decimal digits.}
  \item{method}{A character string specifying the linkage method to be used. For
    \code{linkage()}, this should be one of: \code{"single"}, \code{"complete"},
    \code{"arithmetic"}, \code{"geometric"}, \code{"harmonic"},
    \code{"versatile"}, \code{"ward"}, \code{"centroid"} or \code{"flexible"}.
    \code{"ward"} and \code{"centroid"} methods cannot be used with similarity
    data. \code{"versatile"} and \code{"flexible"} are the only two methods that
    can be used in \code{descplot()}. See the \emph{Details} section.}
  \item{par.method}{A real value, in the case of \code{linkage()}, or a vector
    of real values, in the case of \code{descplot()}, required as parameter for
    the methods \code{"versatile"} and \code{"flexible"}. The range of possible
    values is \code{[-Inf, +Inf]} for \code{"versatile"}, and \code{[-1, +1]}
    for \code{"flexible"}. See the \emph{Details} section.}
  \item{weighted}{A logical value to choose between the weighted and the
    unweighted (default) versions of some linkage methods. Weighted linkage
    gives merging branches in a dendrogram equal weight regardless of the number
    of objects carried on each branch. Such a procedure weights objects
    unequally, contrasting with unweighted linkage that gives equal weight to
    each object in the clusters. This parameter has no effect on the
    \code{"single"} and \code{"complete"} linkages.}
  \item{group}{A character string to choose a grouping criterion between the
    \code{"variable"}-group approach (default) that returns a multifurcated
    dendrogram (m-ary tree), and the \code{"pair"}-group approach that returns a
    bifurcated dendrogram (binary tree). See the \emph{Details} section.}
  \item{measure}{A character string specifying the descriptive measure to be
    plotted. This should be one of: \code{"cor"}, for cophenetic correlation
    coefficient; \code{"sdr"}, for space distortion ratio; \code{"ac"}, for
    agglomerative coefficient; \code{"cc"}, for chaining coefficient; or
    \code{"tb"}, for tree balance.}
  \item{slope}{A real value representing the slope of a sigmoid function to
    map the \code{"versatile"} linkage unbounded interval \code{(-Inf, +Inf)}
    onto the bounded interval \code{(-1, +1)}. It can be used to improve the
    distribution of points along the \emph{x} axis.}
  \item{\dots}{Graphical parameters (see \code{\link{par}}) may also be supplied
    and are passed to \code{\link{plot.default}}.}
}
\encoding{latin1}
\details{
  Starting from a matrix of proximity data (distances or similarities),
  \code{linkage()} calculates its dendrogram with the most commonly used
  agglomerative hierarchical clustering methods, i.e. single linkage, complete
  linkage, arithmetic linkage (also known as average linkage) and Ward's method.
  Importantly, it contains a new parameterized method named versatile linkage
  (\enc{Fernández}{Fernandez} and \enc{Gómez}{Gomez}, 2020), which includes
  single linkage, complete linkage and average linkage as particular cases, and
  which naturally defines two new methods, geometric linkage and harmonic
  linkage.

  The difference between the available hierarchical clustering methods rests in
  the way the proximity between two clusters is defined from the proximity
  between their constituent objects:
  \itemize{
    \item \code{"single"}: the proximity between clusters equals the minimum
      distance or the maximum similarity between objects.
    \item \code{"complete"}: the proximity between clusters equals the maximum
      distance or the minimum similarity between objects.
    \item \code{"arithmetic"}: the proximity between clusters equals the
      arithmetic mean proximity between objects. Also known as average linkage,
      WPGMA (weighted version) or UPGMA (unweighted version).
    \item \code{"geometric"}: the proximity between clusters equals the
      geometric mean proximity between objects.
    \item \code{"harmonic"}: the proximity between clusters equals the harmonic
      mean proximity between objects.
    \item \code{"versatile"}: the proximity between clusters equals the
      generalized power mean proximity between objects. It depends on the value
      of \code{par.method}, with the following linkage methods as particular
      cases: \code{"complete"} (\code{par.method=+Inf}), \code{"arithmetic"}
      (\code{par.method=+1}), \code{"geometric"} (\code{par.method=0}),
      \code{"harmonic"} (\code{par.method=-1}) and \code{"single"}
      (\code{par.method=-Inf}).
    \item \code{"ward"}: the distance between clusters is a weighted squared
      Euclidean distance between the centroids of each cluster. This method is
      available only for distance data.
    \item \code{"centroid"}: the distance between clusters equals the square of
      the Euclidean distance between the centroids of each cluster. Also known
      as WPGMC (weighted version) or UPGMC (unweighted version). This method is
      available only for distance data.
    \item \code{"flexible"}: the proximity between clusters is a weighted sum of
      the proximity between clusters in the previous iteration. It depends on
      the value of \code{par.method}, in the range \code{[-1, +1]}, and it is
      equivalent to \code{"arithmetic"} linkage when \code{par.method=0}.
  }

  With the argument \code{group}, users can choose between a variable-group
  approach (default) that returns a multifurcated dendrogram or multidendrogram,
  and a pair-group approach that returns a bifurcated dendrogram.
  Multidendrograms were introduced (\enc{Fernández}{Fernandez} and 
  \enc{Gómez}{Gomez}, 2008) to solve the non-uniqueness problem that arises when
  two or more minimum proximity values between different clusters are equal
  during the agglomerative process. Multidendrograms group more than two
  clusters when tied proximity values occur, what produces a uniquely determined
  solution that does not depend on the order of the input data. When there are
  no tied proximity values, the variable-group approach gives the same result as
  the pair-group one.

  \code{descplot()} can be used with methods \code{"versatile"} and
  \code{"flexible"} to analyze graphically the variation of any descriptive
  measure as a function of the corresponding method parameter.
}
\value{
  An object of class \code{"linkage"} that describes the multifurcated
  dendrogram obtained. The object is a list with the following components:
  \item{call}{The call that produced the result.}
  \item{digits}{Number of significant decimal digits used as precision. It is
    given by the user or automatically set to that of the input proximity value
    with the largest number of significant decimal digits.}
  \item{merger}{A list of vectors of integer that describes the merging of
    clusters at each step of the clustering. If a number \eqn{j} in a vector is
    negative, then singleton cluster \eqn{-j} was merged at this stage. If
    \eqn{j} is positive, then the merge was with the cluster formed at stage
    \eqn{j} of the algorithm.}
  \item{height}{A vector with the proximity values between merging clusters (for
    the particular agglomeration) at the successive stages.}
  \item{range}{A vector with the range (the maximum minus the minimum) of
    proximity values between merging clusters. It is equal to 0 for binary
    clusters.}
  \item{order}{A vector giving a permutation of the original observations to
    allow for plotting, in the sense that the branches of a clustering tree
    will not cross.}
  \item{coph}{Object of class \code{"dist"} containing the cophenetic
    (or ultrametric) proximity data in the output dendrogram, sorted in the same
    order as the input proximity data in \code{prox}.}
  \item{binary}{A logical value indicating whether the output dendrogram is a
    binary tree or, on the contrary, it contains an agglomeration of more than
    two clusters due to the existence of tied proximity data. Its value is
    always \code{TRUE} when the \code{"pair"} grouping criterion is used.}
  \item{cor}{Cophenetic correlation coefficient (Sokal and Rohlf, 1962),
    defined as the Pearson correlation coefficient between the output cophenetic
    proximity data and the input proximity data. It is a measure of how
    faithfully the dendrogram preserves the pairwise proximity between objects.}
  \item{sdr}{Space distortion ratio (\enc{Fernández}{Fernandez} and
    \enc{Gómez}{Gomez}, 2020), calculated as the difference between the maximum
    and minimum cophenetic proximity data, divided by the difference between the
    maximum and minimum initial proximity data. Space dilation occurs when the
    space distortion ratio is greater than 1.}
  \item{ac}{Agglomerative coefficient (Rousseeuw, 1986), a number between 0
    and 1 measuring the strength of the clustering structure obtained.}
  \item{cc}{Chaining coefficient (Williams \emph{et al.}, 1966), a number
    between 0 and 1 measuring the tendency for clusters to grow by the addition
    of clusters much smaller rather than by fusion with other clusters of
    comparable size.}
  \item{tb}{Tree balance (\enc{Fernández}{Fernandez} and
    \enc{Gómez}{Gomez}, 2020), a number between 0 and 1 measuring the equality
    in the number of leaves in the branches concerned at each fusion in the
    hierarchical tree.}

  Class \code{"linkage"} has methods for the following generic functions:
  \code{\link{summary}}, \code{plot} (see \code{\link{plot.linkage}}),
  \code{\link{as.dendrogram}}, \code{\link{as.hclust}} and
  \code{\link{cophenetic}}.
}
\note{
  Except for the cases containing tied proximity data, the following
  equivalences hold between function \code{linkage()} in package \emph{mdendro},
  function \code{\link{hclust}()} in package \emph{stats}, and function
  \code{\link{agnes}()} in package \emph{cluster}. When relevant, weighted
  (\code{W}) or unweighted (\code{U}) versions of the linkage methods and the
  value for \code{par.method} (\eqn{\beta}) are indicated:
  \tabular{lll}{
    \code{linkage()} \tab \code{hclust()} \tab \code{agnes()} \cr
    \code{==================} \tab \code{============} \tab
      \code{===================} \cr
    \code{"single"} \tab \code{"single"} \tab \code{"single"} \cr
    \code{"complete"} \tab \code{"complete"} \tab \code{"complete"} \cr
    \code{"arithmetic", U} \tab \code{"average"} \tab \code{"average"} \cr
    \code{"arithmetic", W} \tab \code{"mcquitty"} \tab \code{"weighted"} \cr
    \code{"ward"} \tab \code{"ward.D2"} \tab \code{"ward"} \cr
    \code{"centroid", U} \tab \code{"centroid"} \tab \code{--------} \cr
    \code{"centroid", W} \tab \code{"median"} \tab \code{--------} \cr
    \code{"flexible", U, } \eqn{\beta} \tab \code{--------} \tab
      \code{"gaverage", } \eqn{\beta} \cr
    \code{"flexible", W, } \eqn{\beta} \tab \code{--------} \tab
      \code{"flexible", } \eqn{(1-\beta)/2}
  }
}
\author{
  Alberto \enc{Fernández}{Fernandez} \email{alberto.fernandez@urv.cat}
  and Sergio \enc{Gómez}{Gomez} \email{sergio.gomez@urv.cat}.
}
\references{
  \enc{Fernández}{Fernandez}, A.; \enc{Gómez}{Gomez}, S. (2008). Solving
  non-uniqueness in agglomerative hierarchical clustering using
  multidendrograms. \emph{Journal of Classification}, \bold{25}, 43--65.

  \enc{Fernández}{Fernandez}, A.; \enc{Gómez}{Gomez}, S. (2020). Versatile
  linkage: a family of space-conserving strategies for agglomerative
  hierarchical clustering. \emph{Journal of Classification}, \bold{37},
  584--597.

  Rousseeuw, P.J. (1986). A visual display for hierarchical classification. In
  E. Diday \emph{et al.} (eds.) \emph{Data Analysis and Informatics 4}, pp.
  743--748. Amsterdam: North-Holland.

  Sokal, R.R.; Rohlf, F.J. (1962). The comparison of dendrograms by objective
  methods. \emph{Taxon}, \bold{11}, 33--40.

  Williams, W.T.; Lambert, J.M.; Lance, G.N. (1966). Multivariate methods in
  plant ecology: V. Similarity analyses and information-analysis.
  \emph{Journal of Ecology}, \bold{54}, 427--445.
}
\seealso{
  \code{\link{plot.linkage}}, \code{\link{dist}}, \code{\link{dendrogram}},
  \code{\link{hclust}}, \code{\link{agnes}}.
}
\examples{
## Plot and summary of unweighted arithmetic linkage (UPGMA) dendrogram
lnk1 <- linkage(UScitiesD)
plot(lnk1)
summary(lnk1)

## Linkage of similarity data (non-negative correlations)
sim <- as.dist(cor(EuStockMarkets))
lnk2 <- linkage(sim, type.prox = "similarity")
plot(lnk2)

## Use function as.dendrogram to plot with package dendextend
d <- dist(scale(mtcars))  # distances of standardized data
lnk <- linkage(d, digits = 1, method = "complete")
lnk.dend <- as.dendrogram(lnk)
plot(dendextend::set(lnk.dend, "branches_k_color", k = 4),
     nodePar = list(cex = 0.4, lab.cex = 0.5))

## Plot heatmap containing multidendrograms
heatmap(scale(mtcars), hclustfun = linkage)

## Plot of different versatile linkages as we increase the method parameter
d = as.dist(matrix(c( 0,  7, 16, 12,
                      7,  0,  9, 19,
                     16,  9,  0, 12,
                     12, 19, 12, 0), nrow = 4))
par(mfrow = c(2, 3))
vals <- c(-Inf, -1, 0, +1, +Inf)
names <- c("single", "harmonic", "geometric", "arithmetic", "complete")
for (i in 1:length(vals)) {
  lnk <- linkage(d, digits = 1, method = "versatile", par.method = vals[i])
  plot(lnk, main = paste0("versatile (", vals[i], ") = ", names[i]),
       ylim = c(0, 20), cex = 0.6)
}

## Analyze how descriptive measures depend on versatile linkage parameter
par(mfrow = c(2, 3))
measures <- c("cor", "sdr", "ac", "cc", "tb")
vals <- c(-Inf, (-20:+20), +Inf)
for (measure in measures) {
  descplot(UScitiesD, method = "versatile", par.method = vals,
           measure = measure,  main = measure, type = "o", col = "blue")
}
}
\keyword{cluster}
\concept{cophenetic correlation coefficient}
\concept{space distortion ratio}
\concept{agglomerative coefficient}
\concept{chaining coefficient}
\concept{tree balance}
