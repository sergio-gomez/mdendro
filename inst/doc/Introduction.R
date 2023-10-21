## ----eval = FALSE-------------------------------------------------------------
#  install.packages("mdendro")

## ----eval = FALSE-------------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("sergio-gomez/mdendro")

## -----------------------------------------------------------------------------
library(mdendro)
lnk <- linkage(UScitiesD, method = "complete")

## -----------------------------------------------------------------------------
plot(lnk)

## -----------------------------------------------------------------------------
summary(lnk)

## -----------------------------------------------------------------------------
sim <- as.dist(Harman23.cor$cov)
lnk <- linkage(sim, type.prox = "sim")
plot(lnk, main = "Harman23")

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
cars <- round(dist(scale(mtcars)), digits = 3)
nodePar <- list(cex = 0, lab.cex = 0.4)

lnk1 <- linkage(cars, method = "arithmetic")
plot(lnk1, main = "unweighted", nodePar = nodePar)

lnk2 <- linkage(cars, method = "arithmetic", weighted = TRUE)
plot(lnk2, main = "weighted", nodePar = nodePar)

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
cars <- round(dist(scale(mtcars)), digits = 1)
nodePar <- list(cex = 0, lab.cex = 0.4)

lnk1 <- linkage(cars, method = "complete")
plot(lnk1, main = "multidendrogram", nodePar = nodePar)

lnk2 <- linkage(cars, method = "complete", group = "pair")
plot(lnk2, main = "pair-group", nodePar = nodePar)

## ----collapse = TRUE----------------------------------------------------------
cars <- round(dist(scale(mtcars)), digits = 1)
lnk1 <- linkage(cars, method = "complete")
lnk2 <- linkage(cars, method = "complete", group = "pair")

# apply random permutation to data
set.seed(666)
ord <- sample(attr(cars, "Size"))
cars_p <- as.dist(as.matrix(cars)[ord, ord])
lnk1p <- linkage(cars_p, method = "complete")
lnk2p <- linkage(cars_p, method = "complete", group = "pair")

# compare original and permuted cophenetic correlation
c(lnk1$cor, lnk1p$cor)
c(lnk2$cor, lnk2p$cor)

# compare original and permuted tree balance
c(lnk1$tb, lnk1p$tb)
c(lnk2$tb, lnk2p$tb)

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
cars <- round(dist(scale(mtcars)), digits = 1)
nodePar <- list(cex = 0, lab.cex = 0.4)
lnk <- linkage(cars, method = "complete")
plot(lnk, col.rng = "bisque", main = "with ranges", nodePar = nodePar)
plot(lnk, col.rng = NULL, main = "without ranges", nodePar = nodePar)

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
cars <- round(dist(scale(mtcars)), digits = 1)
lnk <- linkage(cars, method = "complete")

lnk.dend <- as.dendrogram(lnk)
plot(dendextend::set(lnk.dend, "branches_k_color", k = 4),
     main = "dendextend package",
     nodePar = list(cex = 0.4, lab.cex = 0.5))

lnk.hcl <- as.hclust(lnk)
pal4 <- c("red", "forestgreen", "purple", "orange")
clu4 <- cutree(lnk.hcl, 4)
plot(ape::as.phylo(lnk.hcl),
     type = "fan",
     main = "ape package",
     tip.color = pal4[clu4],
     cex = 0.5)

## -----------------------------------------------------------------------------
heatmap(scale(mtcars), hclustfun = linkage)

## -----------------------------------------------------------------------------
par(mfrow = c(3, 3))
methods <- c("single", "complete", "arithmetic",
             "geometric", "harmonic", "versatile",
             "ward", "centroid", "flexible")

for (m in methods) {
  lnk <- linkage(UScitiesD, method = m)
  plot(lnk, cex = 0.6, main = m)
}

## -----------------------------------------------------------------------------
par(mfrow = c(2, 3))

vals <- c(-10.0, 0.0, 10.0)
for (v in vals) {
  lnk <- linkage(UScitiesD, method = "versatile", par.method = v)
  plot(lnk, cex = 0.6, main = sprintf("versatile (%.1f)", v))
}

vals <- c(-0.8, 0.0, 0.8)
for (v in vals) {
  lnk <- linkage(UScitiesD, method = "flexible", par.method = v)
  plot(lnk, cex = 0.6, main = sprintf("flexible (%.1f)", v))
}

## -----------------------------------------------------------------------------
par(mfrow = c(2, 3))
measures <- c("cor", "sdr", "ac", "cc", "tb")
vals <- c(-Inf, (-20:+20), +Inf)
for (m in measures)
  descplot(UScitiesD, method = "versatile",
           measure = m, par.method = vals,
           type = "o",  main = m, col = "blue")

## -----------------------------------------------------------------------------
par(mfrow = c(2, 3))
measures <- c("cor", "sdr", "ac", "cc", "tb")
vals <- seq(from = -1, to = +1, by = 0.1)
for (m in measures)
  descplot(UScitiesD, method = "flexible",
           measure = m, par.method = vals,
           type = "o",  main = m, col = "blue")

## -----------------------------------------------------------------------------
library(mdendro)
lnk <- linkage(UScitiesD, method = "complete")

library(cluster)
agn <- agnes(UScitiesD, method = "complete")

# library(stats)   # unneeded, stats included by default
hcl <- hclust(UScitiesD, method = "complete")

par(mfrow = c(1, 3))
plot(lnk)
plot(agn, which.plots = 2)
plot(hcl)

## ----collapse = TRUE----------------------------------------------------------
lnk.dend <- as.dendrogram(lnk)
agn.dend <- as.dendrogram(agn)
hcl.dend <- as.dendrogram(hcl)
identical(lnk.dend, agn.dend)

par(mfrow = c(1, 2))
plot(lnk.dend, main = "lnk.dend = agn.dend", cex = 0.7)
plot(hcl.dend, main = "hcl.dend", cex = 0.7)

## ----collapse = TRUE----------------------------------------------------------
hcl.coph <- cophenetic(hcl)
agn.coph <- cophenetic(agn)

all(lnk$coph == hcl.coph)
all(lnk$coph == agn.coph)

## ----collapse = TRUE----------------------------------------------------------
hcl.cor <- cor(UScitiesD, hcl.coph)
all.equal(lnk$cor, hcl.cor)

all.equal(lnk$ac, agn$ac)

## ----echo = FALSE, warning = FALSE, results = "hide"--------------------------
data <- "Name	VVMD5.1	VVMD5.2	VVMD7.1	VVMD7.2	VVMD27.1	VVMD27.2	VrZag62.1	VrZag62.2	VrZag79.1	VrZag79.2	VVS2-A1	VVS2.2
Alfrocheiro	226	238	249	253	179	189	188	200	251	251	145	153
Alvarinho	222	232	235	235	189	189	186	204	247	251	137	153
Antao Vaz	234	236	245	259	181	183	204	204	247	247	147	153
Aragonez	236	236	235	249	183	183	196	200	247	251	145	147
Arinto	226	238	239	247	181	185	186	188	247	251	145	153
Avesso	222	240	235	235	181	189	186	186	243	247	137	153
Azal	226	232	235	243	181	185	194	204	247	251	153	159
Baga	232	240	235	235	179	189	188	204	247	251	145	157
Bastardo	238	238	235	253	175	189	188	188	245	247	145	153
Bical	226	240	235	259	179	185	188	194	251	251	135	147
Camarate	234	236	239	249	181	189	188	200	247	251	147	153
Castelao	236	238	239	253	179	181	188	188	247	251	145	147
Cerceal Branco	226	236	245	253	179	181	188	204	247	251	145	159
Encruzado	226	232	235	253	183	189	194	194	247	251	151	153
Espadeiro	222	226	243	259	183	189	196	204	251	251	135	153
Fernao Pires	226	240	235	235	183	194	188	194	247	247	147	153
Folgasao	232	240	239	243	185	189	194	204	245	251	135	153
Galego Dourado	228	240	235	239	185	189	188	194	245	251	135	135
Gouveio	226	238	235	239	185	189	186	188	251	251	153	159
Jaen	226	236	245	253	181	189	188	194	247	251	147	153
Loureiro	232	232	247	259	181	185	186	196	251	251	145	153
Malvasia Fina	226	240	235	253	179	194	188	188	247	251	145	147
Malvasia Rei	228	240	235	245	185	194	188	194	251	257	135	147
Moreto	226	236	239	253	181	189	188	188	247	251	147	153
Moscatel Graudo	228	232	245	247	179	194	186	204	247	255	135	151
Moscatel Galego	226	228	235	245	179	189	186	188	251	255	135	153
Negra Mole	222	240	235	235	181	181	188	196	247	259	145	159
Perrum	236	240	235	239	181	185	188	188	243	247	135	147
Rabo de Ovelha	222	236	235	239	181	181	188	194	247	247	139	153
Rabigato	222	232	235	235	185	189	186	196	243	251	135	135
Rufete	226	236	235	253	181	189	188	194	245	247	135	159
Ramisco	226	238	235	259	181	185	188	196	247	251	135	159
Sercial	226	238	235	249	181	185	188	194	247	259	135	153
Siria	222	234	235	245	181	181	186	204	247	247	139	153
Talia	226	232	245	249	179	183	194	200	245	251	135	145
Terrantez	226	238	243	259	185	189	194	196	251	251	145	159
Tinta Barroca	228	236	235	239	181	183	188	192	245	247	145	153
Tinta Caiada	222	238	235	235	179	189	186	188	251	261	135	135
Tinto Cao	232	234	235	259	181	185	186	194	247	251	135	135
Tinta Miuda	226	238	235	235	179	183	186	188	251	259	141	153
Touriga Franca	226	228	235	239	181	183	192	194	245	247	145	153
Touriga Nacional	226	236	235	235	181	189	188	194	245	245	145	153
Trajadura	226	236	235	247	181	185	186	186	247	247	145	153
Trincadeira	234	238	235	245	181	185	188	204	247	251	135	153
Trincadeira das Pratas	238	240	235	253	189	194	188	188	251	257	143	145
Verdelho	222	232	235	253	181	189	194	196	247	251	135	153
Viosinho	232	232	235	239	185	189	186	188	243	245	135	153
Vinhao	222	226	235	259	189	189	188	196	245	251	135	137
Vital	222	240	235	235	181	194	188	188	247	253	147	153
Borracal	232	238	235	235	181	185	194	194	247	247	135	137
Fonte Cal	226	234	235	235	183	185	186	186	247	251	135	153"

dt <- read.table(textConnection(data), header = TRUE, sep = "\t", dec = ".")
# knitr::kable(dt)

n <- nrow(dt)
cols <- 2:ncol(dt)
nc <- length(cols)

m <- matrix(0.0, nrow = n, ncol = n,
            dimnames = list(dt$Name, dt$Name))

for (i in 1:(n-1)) {
    for (j in (i+1):n) {
        m[i, j] <- 1 - sum(dt[i, cols] == dt[j, cols]) / nc
        m[j, i] <- m[i, j]
    }
}

d <- as.dist(m)

## ----collapse = TRUE----------------------------------------------------------
length(unique(d))

## -----------------------------------------------------------------------------
lnk <- linkage(d, method = "arithmetic", digits = 3)
nodePar <- list(cex = 0, lab.cex = 0.6)
plot(lnk, col.rng = NULL, nodePar = nodePar, main = "multidendrogram")

## -----------------------------------------------------------------------------
nreps <- 1000
cors <- vector(length = nreps)
lnks <- list()
for (r in 1:nreps) {
    ord <- sample(attr(d, "Size"))
    d2 <- as.dist(as.matrix(d)[ord, ord])
    lnks[[r]] <- linkage(d2, group = "pair", digits = 3)
    cors[r] <- lnks[[r]]$cor
}
plot(sort(cors), main = "pair-group cophenetic correlations (same data)")

## -----------------------------------------------------------------------------
dend1 <- as.dendrogram(lnks[[1]])
dend2 <- as.dendrogram(lnks[[2]])
diff_12 <- dendextend::highlight_distinct_edges(dend1, dend2)
diff_21 <- dendextend::highlight_distinct_edges(dend2, dend1)

par(mfrow = c(2, 1))
nodePar <- list(cex = 0, lab.cex = 0.6)
plot(diff_12, nodePar = nodePar, main = "pair-group (original)")
plot(diff_21, nodePar = nodePar, main = "pair-group (permuted)")

## -----------------------------------------------------------------------------
lnk3 <- linkage(d, method = "arithmetic", digits = 3)
lnk1 <- linkage(d, method = "arithmetic", digits = 1)

par(mfrow = c(2, 1))
nodePar <- list(cex = 0, lab.cex = 0.6)
plot(lnk3, col.rng = NULL, nodePar = nodePar,
     main = "multidendrogram (digits = 3)")
plot(lnk1, col.rng = NULL, nodePar = nodePar,
     main = "multidendrogram (digits = 1)")

## ----echo = FALSE-------------------------------------------------------------
vl_tab <- data.frame(
  linkage = c("complete", "arithmetic", "geometric", "harmonic", "single"),
  par.method = c(+Inf, +1, 0, -1, -Inf)
)
knitr::kable(vl_tab, col.names = c("linkage", "versatile linkage (`par.method`)"))

## -----------------------------------------------------------------------------
d = as.dist(matrix(c( 0,  7, 16, 12,
                      7,  0,  9, 19,
                     16,  9,  0, 12,
                     12, 19, 12, 0), nrow = 4))

par(mfrow = c(2, 3))

vals <- c(-Inf, -1, 0, 1, Inf)
names <- c("single", "harmonic", "geometric", "arithmetic", "complete")
titles <- sprintf("versatile (%.1f) = %s", vals, names)

for (i in 1:length(vals)) {
  lnk <- linkage(d, method = "versatile", par.method = vals[i], digits = 2)
  plot(lnk, ylim = c(0, 20), cex = 0.6, main = titles[i])
}

## ----eval = FALSE-------------------------------------------------------------
#  linkage(prox, type.prox = "distance", digits = NULL,
#          method = "arithmetic", par.method = 0, weighted = FALSE,
#          group = "variable")
#  
#  descval(prox, type.prox = "distance", digits = NULL,
#          method = "versatile", par.method = c(-1,0,+1), weighted = FALSE,
#          group = "variable", measure = "cor")
#  
#  descplot(prox, ..., type.prox = "distance", digits = NULL,
#           method = "versatile", par.method = c(-1, 0, +1), weighted = FALSE,
#           group = "variable", measure = "cor", slope = 10)

## ----echo = FALSE-------------------------------------------------------------
arguments <- data.frame(
  arg = c("`prox`",
          "`type.prox`",
          "`digits`",
          "`method`",
          "`par.method`",
          "`weighted`",
          "`group`",
          "`measure`",
          "`slope`",
          "`...`"),
  desc = c("A structure of class `dist` containing non-negative proximity data (distances or similarities). All the linkage methods are used with non-squared proximity data as input, except for method `\"centroid\"`, which is meant to be used with squared Euclidean distances."
           ,
            "A character string to indicate whether the proximity data represent `\"distance\"` (default) or `\"similarity\"` between objects. Methods `\"ward\"` and `\"centroid\"` cannot be used with similarity data as input, while the rest of the linkage methods can be used with both distances and similarities."
           ,
            "An integer value specifying the precision, i.e., the number of significant decimal digits to be used for the comparisons between proximity data. This is an important parameter, since equal proximity data at a certain precision may become different by increasing its value. Thus, it may be responsible of the existence of tied proximity data. If the value of this parameter is negative or `NULL` (default), then the precision is automatically set to that of the input proximity value with the largest number of significant decimal digits."
           ,
            "A character string specifying the linkage method to be used. For `linkage()`, this should be one of: `\"single\"`, `\"complete\"`, `\"arithmetic\"`, `\"geometric\"`, `\"harmonic\"`, `\"versatile\"`, `\"ward\"`, `\"centroid\"` or `\"flexible\"`. Methods `\"versatile\"` and `\"flexible\"` are the only two methods that can be used in `descval()` and `descplot()`."
           ,
            "A real value, in the case of `linkage()`, or a vector of real values, in the case of `descval()` and `descplot()`, required as parameter for the methods `\"versatile\"` and `\"flexible\"`. The range of possible values is [-Inf, +Inf] for `\"versatile\"`, and [-1, +1] for `\"flexible\"`."
           ,
            "A logical value to choose between the weighted and the unweighted (default) versions of some linkage methods. Weighted linkage gives merging branches in a dendrogram equal weight regardless of the number of objects carried on each branch. Such a procedure weights objects unequally, contrasting with unweighted linkage that gives equal weight to each object in the clusters. This parameter has no effect on the `\"single\"` and `\"complete\"` linkages."
           ,
            "A character string to choose a grouping criterion between the `\"variable\"`-group approach (default) that returns a multidendrogram, i.e., a multifurcated dendrogram (m-ary tree), and the `\"pair\"`-group approach that returns a bifurcated dendrogram (binary tree)."
           ,
            "A character string specifying the descriptive measure to be plotted. This should be one of: `\"cor\"`, for cophenetic correlation coefficient; `\"sdr\"`, for space distortion ratio; `\"ac\"`, for agglomerative coefficient; `\"cc\"`, for chaining coefficient; or `\"tb\"`, for tree balance."
           ,
            "A real value representing the slope of a sigmoid function to map the `\"versatile\"` linkage unbounded interval (-Inf, +Inf) onto the bounded interval (-1, +1). It can be used to improve the distribution of points along the x axis."
           ,
            "Graphical parameters may also be supplied (see `par`) and are passed to `plot.default`."
            )
)
knitr::kable(arguments, col.names = c("", ""))

## ----eval = FALSE-------------------------------------------------------------
#  ## S3 method for class 'linkage'
#  plot(x, type = c("rectangle", "triangle"),
#       center = FALSE, edge.root = FALSE,
#       nodePar = NULL, edgePar = list(),
#       leaflab = c("perpendicular", "textlike", "none"),
#       dLeaf = NULL, xlab = "", ylab = "", xaxt = "n", yaxt = "s",
#       horiz = FALSE, frame.plot = FALSE, xlim, ylim,
#       col.rng = "lightgray", ...)

## ----echo = FALSE-------------------------------------------------------------
arguments <- data.frame(
  arg = c("`x`",
          "`type`",
          "`center`",
          "`edge.root`",
          "`nodePar`",
          "`edgePar`",
          "`leaflab`",
          "`dLeaf`",
          "`xlab`,`ylab`",
          "`xaxt`,`yaxt`",
          "`horiz`",
          "`frame.plot`",
          "`xlim`,`ylim`",
          "`col.rng`",
          "`...`"),
  desc = c("An object of class `linkage`, typically created by `linkage()`."
           ,
            "Type of plot"
           ,
            "Logical; if `TRUE`, nodes are plotted centered with respect to the leaves in the branch. Otherwise (default), plot them in the middle of all direct child nodes."
           ,
            "Logical; if `TRUE`, draw an edge to the root node."
           ,
            "A `list` of plotting parameters to use for the nodes (see `points`), or `NULL` by default which does not draw symbols at the nodes. The list may contain components named `pch`, `cex`, `col`, `xpd`, and/or `bg` each of which can have length two for specifying separate attributes for inner nodes and leaves. Note that the default of `pch` is `1:2`, so you may want to use `pch = NA` if you specify `nodePar`."
           ,
            "A `list` of plotting parameters to use for the edge `segments`. The list may contain components named `co`l, `lty` and `lwd`. As with `nodePar`, each can have length two for differentiating leaves and inner nodes."
           ,
            "A string specifying how leaves are labeled. The default `\"perpendicular\"` writes text vertically (by default), `\"textlike\"` writes text horizontally (in a rectangle), and `\"none\"` suppresses leaf labels."
           ,
            "A number specifying the distance in user coordinates between the tip of a leaf and its label. If `NULL` as per default, 3/4 of a letter width or height is used."
           ,
            "A label for the axis."
           ,
            "A character which specifies the axis type. Specifying `\"n\"` suppresses plotting, while any value other than `\"n\"` implies plotting."
           ,
            "Logical indicating if the dendrogram should be drawn horizontally or not."
           ,
            "Logical indicating if a box around the plot should be drawn, see `plot.default`."
           ,
            "Optional x- and y-limits of the plot, passed to `plot.default`. The defaults for these show the full dendrogram."
           ,
            "Color (`\"lightgray\"` by default) to be used for plotting range rectangles in case of tied heights. If `NULL`, range rectangles are not plotted."
           ,
            "Graphical parameters (see `par`) may also be supplied and are passed to `plot.default`."
            )
)
knitr::kable(arguments, col.names = c("", ""))

## ----echo = FALSE-------------------------------------------------------------
components <- data.frame(
  val = c("`call`",
          "`digits`",
          "`merger`",
          "`height`",
          "`range`",
          "`order`",
          "`coph`",
          "`binary`",
          "`cor`",
          "`sdr`",
          "`ac`",
          "`cc`",
          "`tb`"),
  desc = c("The call that produced the result."
           ,
            "Number of significant decimal digits used as precision. It is given by the user or automatically set to that of the input proximity value with the largest number of significant decimal digits."
           ,
            "A list of vectors of integer that describes the merging of clusters at each step of the clustering. If a number _j_ in a vector is negative, then singleton cluster _-j_ was merged at this stage. If _j_ is positive, then the merge was with the cluster formed at stage _j_ of the algorithm."
           ,
            "A vector with the proximity values between merging clusters (for the particular agglomeration) at the successive stages."
           ,
            "A vector with the range (the maximum minus the minimum) of proximity values between merging clusters. It is equal to 0 for binary clusters."
           ,
            "A vector giving a permutation of the original observations to allow for plotting, in the sense that the branches of a clustering tree will not cross."
           ,
            "Object of class `dist` containing the cophenetic (or ultrametric) proximity data in the output dendrogram, sorted in the same order as the input proximity data in `prox`."
           ,
            "A logical value indicating whether the output dendrogram is a binary tree or, on the contrary, it contains an agglomeration of more than two clusters due to the existence of tied proximity data. Its value is always `TRUE` when the `pair` grouping criterion is used."
           ,
            "Cophenetic correlation coefficient [[@Sokal1962]](#bibliography), defined as the Pearson correlation coefficient between the output cophenetic proximity data and the input proximity data. It is a measure of how faithfully the dendrogram preserves the pairwise proximity between objects."
           ,
            "Space distortion ratio [[@Fernandez2020]](#bibliography), calculated as the difference between the maximum and minimum cophenetic proximity data, divided by the difference between the maximum and minimum initial proximity data. Space dilation occurs when the space distortion ratio is greater than 1."
           ,
            "Agglomerative coefficient [[@Rousseeuw1986]](#bibliography), a number between 0 and 1 measuring the strength of the clustering structure obtained."
           ,
            "Chaining coefficient [[@Williams1966]](#bibliography), a number between 0 and 1 measuring the tendency for clusters to grow by the addition of clusters much smaller rather than by fusion with other clusters of comparable size."
           ,
            "Tree balance [[@Fernandez2020]](#bibliography), a number between 0 and 1 measuring the equality in the number of leaves in the branches concerned at each fusion in the hierarchical tree."
            )
)
knitr::kable(components, col.names = c("", ""))

## ----echo = FALSE-------------------------------------------------------------
methods <- data.frame(
  linkage   = c("single", "complete", "arithmetic, U", "arithmetic, W",
                "geometric, U/W", "harmonic, U/W", "versatile, U/W, $p$",
                "---", "ward", "centroid, U", "centroid, W",
                "flexible, U, $\\beta$", "---",
                "flexible, W, $\\beta$", "---"),
  hclust    = c("single", "complete", "average", "mcquitty",
                "---", "---", "---", "ward", "ward.D2",
                "centroid", "median", "---", "---", "---", "---"),
  agnes     = c("single", "complete", "average", "weighted",
                "---", "---", "---", "---", "ward",
                "---", "---",
                "gaverage, $\\beta$",
                "gaverage, $\\alpha_1$, $\\alpha_2$, $\\beta$, $\\gamma$",
                "flexible, $(1-\\beta)/2$",
                "flexible, $\\alpha_1$, $\\alpha_2$, $\\beta$, $\\gamma$")
)
knitr::kable(methods, col.names = c("`linkage`", "`hclust`", "`agnes`"))

