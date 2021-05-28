\name{plot.linkage}
\alias{plot.linkage}
\title{Plots for Extended Agglomerative Hierarchical Clustering}
\description{
  Creates plots for visualizing an object of class \code{"linkage"}.
}
\usage{
\method{plot}{linkage}(x, type = c("rectangle", "triangle"),
     center = FALSE, edge.root = FALSE,
     nodePar = NULL, edgePar = list(),
     leaflab = c("perpendicular", "textlike", "none"),
     dLeaf = NULL, xlab = "", ylab = "", xaxt = "n", yaxt = "s",
     horiz = FALSE, frame.plot = FALSE, xlim, ylim,
     col.rng = "lightgray", \dots)
}
\arguments{
  \item{x}{An object of class \code{"linkage"}, typically created by
    \code{\link{linkage}()}.}
  \item{type}{Type of plot.}
  \item{center}{Logical; if \code{TRUE}, nodes are plotted centered with respect
    to the leaves in the branch. Otherwise (default), plot them in the middle of
    all direct child nodes.}
  \item{edge.root}{Logical; if true, draw an edge to the root node.}
  \item{nodePar}{A \code{list} of plotting parameters to use for the nodes (see
    \code{\link{points}}) or \code{NULL} by default which does not draw symbols
    at the nodes. The list may contain components named \code{pch}, \code{cex},
    \code{col}, \code{xpd}, and/or \code{bg} each of which can have length two
    for specifying separate attributes for inner nodes and leaves.
    Note that the default of \code{pch} is \code{1:2}, so you may want to use
    \code{pch = NA} if you specify \code{nodePar}.}
  \item{edgePar}{A \code{list} of plotting parameters to use for the edge
    \code{\link{segments}}. The list may contain components named \code{col},
    \code{lty} and \code{lwd}. As with \code{nodePar}, each can have length two
    for differentiating leaves and inner nodes.}
  \item{leaflab}{A string specifying how leaves are labeled. The default
    \code{"perpendicular"} writes text vertically (by default),
    \code{"textlike"} writes text horizontally (in a rectangle), and
    \code{"none"} suppresses leaf labels.}
  \item{dLeaf}{A number specifying the distance in user coordinates between the
    tip of a leaf and its label. If \code{NULL} as per default, 3/4 of a letter
    width or height is used.}
  \item{xlab, ylab}{A label for the axis.}
  \item{xaxt, yaxt}{A character which specifies the axis type. Specifying "n"
    suppresses plotting, while any value other than "n" implies plotting.}
  \item{horiz}{Logical indicating if the dendrogram should be drawn
    horizontally or not.}
  \item{frame.plot}{Logical indicating if a box around the plot should be drawn,
    see \code{\link{plot.default}}.}
  \item{xlim, ylim}{Optional x- and y-limits of the plot, passed to
    \code{\link{plot.default}}. The defaults for these show the full
    dendrogram.}
  \item{col.rng}{Color ("lightgray" by default) to be used for plotting range
    rectangles in case of tied heights. If \code{NULL}, range rectangles are not
    plotted.}
  \item{\dots}{Graphical parameters (see \code{\link{par}}) may also be supplied
    and are passed to \code{\link{plot.default}}.}
}
\details{
  Based on the plot function for objects of class \code{"dendrogram"} (see
  \code{\link{plot.dendrogram}}), the plot function for objects of class
  \code{"linkage"} adds the possibility of visualizing the existence of tied
  heights in a dendrogram thanks to the \code{col.rng} parameter.
}
\seealso{
  \code{\link{linkage}}, \code{\link{dendrogram}}.
}
\examples{
## Plot complete linkage of mtcars distances, showing and hiding ranges
mtcars.dist <- dist(scale(mtcars))  # distances of standardized data
lnk <- linkage(mtcars.dist, digits = 1, method = "complete")
par(mfrow = c(1, 2))
nodePar <- list(cex = 0, lab.cex = 0.4)
plot(lnk, col.rng = "bisque", main = "show ranges", nodePar = nodePar)
plot(lnk, col.rng = NULL, main = "hide ranges", nodePar = nodePar)
}
\keyword{cluster}
\keyword{hplot}
