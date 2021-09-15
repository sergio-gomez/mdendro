TYPES.PROX <- c("distance", "similarity")
METHODS <- c("arithmetic", "single", "complete", "geometric", "harmonic",
    "versatile", "ward", "centroid", "flexible")
GROUPS <- c("variable", "pair")
MEASURES <- c("cor", "sdr", "ac", "cc", "tb")

linkage <- function(prox, type.prox = "distance", digits = NULL,
    method = "arithmetic", par.method = 0, weighted = FALSE, group = "variable")
{
  # Check parameters
  if (class(prox) != "dist") {
    stop("'prox' must be an object of class \"dist\"")
  }
  if (attr(prox, "Size") < 2L) {
    stop("'prox' must have at least 2 objects to cluster")
  }
  if (anyNA(prox)) {
    stop("NA-values in 'prox' not allowed")
  }
  storage.mode(prox) <- "double"
  type.prox <- match.arg(type.prox, TYPES.PROX)
  if (min(prox) < 0) {
    stop("Data must be non-negative")
  } else if (max(prox) == Inf) {
    stop("Data must be less than Inf")
  }
  if (is.null(digits)) {
    digits <- -1L
  }
  method <- match.arg(method, METHODS)
  if (((method == "ward") || (method == "centroid"))
      && (type.prox == "similarity")) {
    stop("'ward' and 'centroid' cannot be used with similarity data")
  }
  if ((method == "flexible") && ((par.method < -1) || (+1 < par.method))) {
    stop("'par.method' for 'flexible' must be between -1 and +1")
  }
  group <- match.arg(group, GROUPS)
  # Build agglomerative hierarchical clustering
  lnk <- rcppLinkage(prox=as.numeric(prox),
      isDistance=as.logical(type.prox == "distance"), digits=as.integer(digits),
      method=as.character(method), methodPar=as.numeric(par.method),
      isWeighted=as.logical(weighted),
      isVariable=as.logical(group == "variable"))
  # Permute the original observations
  ord <- integer(attr(prox, "Size"))
  jOrd <- 0L
  stack <- length(lnk$height)
  while (length(stack) > 0L) {
    top <- stack[1L]
    stack <- stack[-1L]
    if (top < 0L) {  # leaf
      jOrd <- jOrd + 1L
      ord[jOrd] <- -top
    } else {  # inner node
      children <- c()
      m <- lnk$merger[[top]]
      nClusters <- length(m)
      for (i in 1L : nClusters) {
        children <- c(children, m[i])
      }
      stack <- c(children, stack)
    }
  }
  # Set attributes for cophenetic matrix
  attr(lnk$coph, "Labels") <- attr(prox, "Labels")
  attr(lnk$coph, "Size") <- attr(prox, "Size")
  attr(lnk$coph, "call") <- match.call()
  attr(lnk$coph, "Diag") <- FALSE
  attr(lnk$coph, "Upper") <- FALSE
  class(lnk$coph) <- "dist"
  # Return object of class "linkage"
  structure(list(
      call = match.call(),
      digits = lnk$digits,
      merger = lnk$merger,
      height = lnk$height,
      range = lnk$range,
	  order = ord,
      coph = lnk$coph,
	  binary = as.logical(length(lnk$height) == (attr(prox, "Size") - 1L)),
      cor = lnk$cor,
      sdr = lnk$sdr,
      ac = lnk$ac,
      cc = lnk$cc,
      tb = lnk$tb),
    class = "linkage")
}

descplot <- function(prox, ..., type.prox = "distance", digits = NULL,
    method = "versatile", par.method = c(-1,0,+1), weighted = FALSE,
    group = "variable", measure = "cor", slope = 10)
{
  method <- match.arg(method, c("versatile", "flexible"))
  measure <- match.arg(measure, MEASURES)
  if (slope <= 0) {
    stop("'slope' must be a positive number")
  }
  y <- numeric(length(par.method))
  for (i in seq_along(par.method)) {
    lnk <- linkage(prox, type.prox=type.prox, digits=digits, method=method,
        par.method=par.method[i], weighted=weighted, group=group)
    y[i] <- lnk[[measure]]
  }
  if (method == "flexible") {
    plot(x=par.method, y=y, ylab=measure, xlab="par.method", ...)
  } else {  # (method == "versatile")
    x <- tanh(par.method / slope)
    plot(x=x, y=y, ylab=measure, xlab="par.method", xaxt="n", ...)
    axis(side=1, at=x, labels=par.method)
  }
}

summary.linkage <- function(object, ...) {
  # Print call
  cat("Call:\n", sep="")
  cl <- object$call
  cat(deparse(cl[[1L]]), "(prox = ", deparse(cl$prox), ",\n", sep="")
  type.prox <- match.arg(cl$type.prox, TYPES.PROX)
  cat("        type.prox = \"", type.prox, "\",\n", sep="")
  cat("        digits = ", object$digits, ",\n", sep="")
  method <- match.arg(cl$method, METHODS)
  cat("        method = \"", method, "\",\n", sep="")
  if ((method == "versatile") || (method == "flexible")) {
    par.method <- if (is.null(cl$par.method)) "0" else deparse(cl$par.method)
    cat("        par.method = ", par.method, ",\n", sep="")
  }
  if ((method != "single") && (method != "complete")) {
    weighted <- if (is.null(cl$weighted)) "FALSE" else deparse(cl$weighted)
    cat("        weighted = ", as.logical(weighted), ",\n", sep="")
  }
  group <- match.arg(cl$group, GROUPS)
  cat("        group = \"", group, "\")\n\n", sep="")
  # Print binary
  cat("Binary dendrogram: ", object$binary, "\n\n", sep="")
  # Print measures
  cat("Descriptive measures:\n", sep="")
  print(c(cor = object$cor, sdr = object$sdr, ac = object$ac, cc = object$cc,
      tb = object$tb), ...)
  invisible(object)
}

cophenetic.linkage <- function(x) {
  x$coph
}

as.hclust.linkage <- function(x, ...) {
  size <- attr(x$coph, "Size")
  mrg <- matrix(nrow=size-1L, ncol=2L)
  storage.mode(mrg) <- "integer"
  hgt <- numeric(size-1L)
  nMergers <- length(x$height)
  ties <- rep.int(0L, times=nMergers)  # number of ties before each merger
  for (k in 1L : nMergers) {
    m <- x$merger[[k]]
    kties <- k + ties[k]
    mrg[kties,1L] <- if (m[1L] < 0L) m[1L] else m[1L] + ties[m[1L]]
    mrg[kties,2L] <- if (m[2L] < 0L) m[2L] else m[2L] + ties[m[2L]]
    hgt[kties] <- x$height[k]
    nClusters <- length(m)
    for (i in seq_len(nClusters - 2L)) {
      mrg[kties+i,1L] <- kties + i - 1L
      mrg[kties+i,2L] <- if (m[i+2L] < 0L) m[i+2L] else m[i+2L] + ties[m[i+2L]]
      hgt[kties+i] <- x$height[k]
    }
    if (nClusters > 2L) {
      for (kk in k : nMergers) {
        ties[kk] <- ties[kk] + nClusters - 2L
      }
    }
  }
  structure(list(
      merge = mrg,
      height = hgt,
      order = x$order,
      labels = attr(x$coph, "Labels"),
      method = match.arg(x$call$method, METHODS),
      call = x$call,
      dist.method = NULL),
    class = "hclust")
}

as.dendrogram.linkage <- function(object, ...) {
  labels <- attr(object$coph, "Labels")
  if (is.null(labels)) {
    size <- attr(object$coph, "Size")
    labels <- as.character(seq_len(size))
  }
  hBottom <- .hgtBottom(object)
  z <- list()
  nMergers <- length(object$height)
  ties <- rep.int(0L, times=nMergers)  # number of ties before each merger
  for (k in 1L : nMergers) {
    x <- object$merger[[k]]
    kties <- k + ties[k]
    zk <- list()
    if (x[1L] < 0L) {  # leaf
      fstPoint <- 0
      zk[[1L]] <- -x[1L]
      attr(zk[[1L]], "label") <- labels[-x[1L]]
      attr(zk[[1L]], "members") <- 1L
      attr(zk[[1L]], "height") <- hBottom
      attr(zk[[1L]], "leaf") <- TRUE
    } else {  # inner node
      fstPoint <- attr(z[[as.character(x[1L])]], "midpoint")
      zk[[1L]] <- z[[as.character(x[1L])]]
      z[[as.character(x[1L])]] <- NULL
    }
    nMembers <- attr(zk[[1L]], "members")
    if (x[2L] < 0L) {  # leaf
      lstPoint <- 0
      zk[[2L]] <- -x[2L]
      attr(zk[[2L]], "label") <- labels[-x[2L]]
      attr(zk[[2L]], "members") <- 1L
      attr(zk[[2L]], "height") <- hBottom
      attr(zk[[2L]], "leaf") <- TRUE
    } else {  # inner node
      lstPoint <- attr(z[[as.character(x[2L])]], "midpoint")
      zk[[2L]] <- z[[as.character(x[2L])]]
      z[[as.character(x[2L])]] <- NULL
    }
    prevMembers <- nMembers
    nMembers <- nMembers + attr(zk[[2L]], "members")
    attr(zk, "members") <- nMembers
    attr(zk, "midpoint") <- (fstPoint + prevMembers + lstPoint) / 2
    attr(zk, "height") <- object$height[k]
    z[[as.character(k)]] <- zk
    nClusters <- length(x)
    for (i in seq_len(nClusters - 2L)) {
      zk <- list()
      zk[[1L]] <- z[[as.character(k)]]
      z[[as.character(k)]] <- NULL
      if (x[i+2L] < 0L) {  # leaf
        lstPoint <- 0
        zk[[2L]] <- -x[i+2L]
        attr(zk[[2L]], "label") <- labels[-x[i+2L]]
        attr(zk[[2L]], "members") <- 1L
        attr(zk[[2L]], "height") <- hBottom
        attr(zk[[2L]], "leaf") <- TRUE
      } else {  # inner node
        lstPoint <- attr(z[[as.character(x[i+2L])]], "midpoint")
        zk[[2L]] <- z[[as.character(x[i+2L])]]
        z[[as.character(x[i+2L])]] <- NULL
      }
      prevMembers <- nMembers
      nMembers <- nMembers + attr(zk[[2L]], "members")
      attr(zk, "members") <- nMembers
      attr(zk, "midpoint") <- (fstPoint + prevMembers + lstPoint) / 2
      attr(zk, "height") <- object$height[k]
      z[[as.character(k)]] <- zk
    }
    if (nClusters > 2L) {
      for (kk in k : nMergers) {
        ties[kk] <- ties[kk] + nClusters - 2L
      }
    }
  }
  structure(z[[as.character(nMergers)]], class = "dendrogram")
}

as.multidendrogram <- function(object) {
  labels <- attr(object$coph, "Labels")
  if (is.null(labels)) {
    size <- attr(object$coph, "Size")
    labels <- as.character(seq_len(size))
  }
  hBottom <- .hgtBottom(object)
  z <- list()
  nMergers <- length(object$height)
  for (k in 1L : nMergers) {
    x <- object$merger[[k]]
    zk <- list()
    nMembers <- 0L
    nClusters <- length(x)
    if (x[1L] < 0L) {  # leaf
      fstPoint <- 0
    } else {  # inner node
      fstPoint <- attr(z[[as.character(x[1L])]], "midpoint")
    }
    if (x[nClusters] < 0L) {  # leaf
      lstPoint <- 0
    } else {  # inner node
      lstPoint <- attr(z[[as.character(x[nClusters])]], "midpoint")
    }
    for (i in 1L : nClusters) {
      if (x[i] < 0L) {  # leaf
        zk[[i]] <- -x[i]
        attr(zk[[i]], "label") <- labels[-x[i]]
        attr(zk[[i]], "members") <- 1L
        attr(zk[[i]], "height") <- hBottom
        attr(zk[[i]], "leaf") <- TRUE
      } else {  # inner node
        zk[[i]] <- z[[as.character(x[i])]]
        z[[as.character(x[i])]] <- NULL
      }
      prevMembers <- nMembers
      nMembers <- nMembers + attr(zk[[i]], "members")
    }
    attr(zk, "members") <- nMembers
    attr(zk, "midpoint") <- (fstPoint + prevMembers + lstPoint) / 2
    attr(zk, "height") <- object$height[k]
    attr(zk, "range") <- object$range[k]
    z[[as.character(k)]] <- zk
  }
  structure(z[[as.character(nMergers)]], class = "dendrogram")
}

plot.linkage <- function (x, type = c("rectangle", "triangle"), center = FALSE,
    edge.root = FALSE,
    nodePar = NULL, edgePar = list(),
    leaflab = c("perpendicular", "textlike", "none"), dLeaf = NULL,
    xlab = "", ylab = "", xaxt = "n", yaxt = "s",
    horiz = FALSE, frame.plot = FALSE, xlim, ylim,
    col.rng = "lightgray", ...)
{
  xd <- as.multidendrogram(x)
  type <- match.arg(type)
  if (type == "triangle") {
    col.rng <- NULL
  }
  leaflab <- match.arg(leaflab)
  hgt <- attr(xd, "height")
  yBot <- .hgtBottom(x)
  if (edge.root && is.logical(edge.root)) {
    edge.root <- 0.0625 * (if (stats::is.leaf(xd)) 1 else abs(hgt - yBot))
  }
  mem.xd <- .memberDend(xd)
  type.prox <- match.arg(x$call$type.prox, TYPES.PROX)
  yTop <- if (type.prox == "distance") hgt + edge.root else hgt - edge.root
  hgtTop <- .hgtDend(xd, type.prox=type.prox, col.rng=col.rng)
  if (center) {
    x1 <- 0.5 ; x2 <- mem.xd + 0.5
  } else {
    x1 <- 1   ; x2 <- mem.xd
  }
  xl. <- c(x1 - 1/2, x2 + 1/2)
  yl. <- c(yBot,
      if (type.prox == "distance") max(yTop, hgtTop) else min(yTop, hgtTop))
  if (horiz) {## swap and reverse direction on `x':
    tmp <- xl.; xl. <- rev(yl.); yl. <- tmp
    tmp <- xaxt; xaxt <- yaxt; yaxt <- tmp
  }
  if (missing(xlim) || is.null(xlim)) xlim <- xl.
  if (missing(ylim) || is.null(ylim)) ylim <- yl.
  grDevices::dev.hold(); on.exit(grDevices::dev.flush())
  plot(0, xlim=xlim, ylim=ylim, type="n", xlab=xlab, ylab=ylab,
      xaxt=xaxt, yaxt=yaxt, frame.plot=frame.plot, ...)
  if (is.null(dLeaf)) {
    dLeaf <- 0.75 *
        (if (horiz) graphics::strwidth("w") else graphics::strheight("x"))
  }

  if (edge.root) {
### FIXME: the first edge + edgetext is drawn here, all others in plotNode()
### -----  maybe use trick with adding a single parent node to the top ?
    x0 <- plotNodeLimit(x1, x2, xd, center)$x
    if (horiz) {
      graphics::segments(hgt, x0, yTop, x0)
    } else {
      graphics::segments(x0, hgt, x0, yTop)
    }
    if (!is.null(et <- attr(xd, "edgetext"))) {
      my <- mean(hgt, yTop)
      if (horiz) graphics::text(my, x0, et) else graphics::text(x0, my, et)
    }
  }
  plotNode(x1, x2, xd, type.prox=type.prox, type=type, center=center,
      leaflab=leaflab, dLeaf=dLeaf, nodePar=nodePar, edgePar=edgePar,
      horiz=horiz, col.rng=col.rng)
}

.memberDend <- function(x) {
  r <- attr(x, "x.member")
  if (is.null(r)) {
    r <- attr(x, "members")
    if (is.null(r)) r <- 1L
  }
  r
}

.hgtDend <- function(node, type.prox, col.rng) {
  hgt <- attr(node, "height")
  if (!stats::is.leaf(node)) {
    rng <- if (is.null(col.rng)) 0 else attr(node, "range")
    hgt <- if (type.prox == "distance") hgt + rng else hgt - rng
    todo <- NULL # Non-leaf nodes to traverse after this one.
    repeat {
      ## For each child: count and add to todo list if it is not a leaf.
      while (length(node)) {
        child <- node[[1L]]
        node <- node[-1L]
        if (stats::is.leaf(child) || is.null(col.rng)) {
          rng <- 0
        } else {
          rng <- attr(child, "range")
        }
        chldHgt <- attr(child, "height")
        if (type.prox == "distance") {
          hgt <- max(hgt, chldHgt + rng)
        } else { # (type.prox == "similarity")
          hgt <- min(hgt, chldHgt - rng)
        }
        if (!stats::is.leaf(child)) {
          todo <- list(node=child, todo=todo)
        }
      }
      ## Advance to next node, terminating when no nodes left.
      if (is.null(todo)) {
        break
      } else {
        node <- todo$node
        todo <- todo$todo
      }
    }
  }
  hgt
}

.hgtBottom <- function(x) {
  type.prox <- match.arg(x$call$type.prox, TYPES.PROX)
  if (type.prox == "distance") 0 else 10^(ceiling(log10(max(x$coph))))
}

### the work horse: plot node (if pch) and lines to all children
plotNode <- function(x1, x2, subtree, type.prox, type, center, leaflab, dLeaf,
    nodePar, edgePar, horiz = FALSE, col.rng)
{
  wholetree <- subtree
  depth <- 0L
  llimit <- list()
  KK <- integer()
  kk <- integer()

  repeat {
    inner <- !stats::is.leaf(subtree) && (x1 != x2)
    yTop <- attr(subtree, "height")
    yRng <- if (is.null(col.rng)) 0 else attr(subtree, "range")
    bx <- plotNodeLimit(x1, x2, subtree, center)
    xTop <- bx$x
    depth <- depth + 1L
    llimit[[depth]] <- bx$limit

    ## handle node specific parameters in "nodePar":
    nPar <- attr(subtree, "nodePar")
    hasP <- !is.null(nPar)
    if (!hasP) nPar <- nodePar

    if (getOption("verbose")) {
      cat(if (inner) "inner node" else "leaf", ":")
      if (!is.null(nPar)) {
        cat(" with node pars\n"); utils::str(nPar)
      }
      cat(if (inner ) paste(" height", formatC(yTop), "; "), "(x1,x2)= (",
          formatC(x1, width=4), ",", formatC(x2, width=4), ")",
          "--> xTop=", formatC(xTop, width=8), "\n", sep="")
    }

    Xtract <- function(nam, L, default, indx) {
      rep(if (nam %in% names(L)) L[[nam]] else default, length.out=indx)[indx]
    }
    asTxt <- function(x) { # to allow 'plotmath' labels:
      if (is.character(x) || is.expression(x) || is.null(x)) {
        x
      } else {
        as.character(x)
      }
    }

    i <- if (inner || hasP) 1 else 2 # only 1 node specific par

    if (!is.null(nPar)) { ## draw this node
      pch <- Xtract("pch", nPar, default=1L:2, i)
      cex <- Xtract("cex", nPar, default=c(1,1), i)
      col <- Xtract("col", nPar, default=graphics::par("col"), i)
      bg <- Xtract("bg", nPar, default=graphics::par("bg"), i)
      graphics::points(if (horiz) cbind(yTop, xTop) else cbind(xTop, yTop),
          pch=pch, bg=bg, col=col, cex=cex)
    }

    if (leaflab == "textlike") {
      p.col <- Xtract("p.col", nPar, default="white", i)
    }
    lab.col <- Xtract("lab.col", nPar, default=graphics::par("col"), i)
    lab.cex <- Xtract("lab.cex", nPar, default=c(1,1), i)
    lab.font <- Xtract("lab.font", nPar, default=graphics::par("font"), i)
    lab.xpd <- Xtract("xpd", nPar, default=c(TRUE,TRUE), i)
    if (stats::is.leaf(subtree)) {
      ## label leaf
      if (leaflab == "perpendicular") { # somewhat like plot.hclust
        if (horiz) {
          X <- yTop + dLeaf * lab.cex
          Y <- xTop; srt <- 0; adj <- c(0, 0.5)
        } else {
          Y <- yTop - dLeaf * lab.cex
          X <- xTop; srt <- 90; adj <- 1
        }
        nodeText <- asTxt(attr(subtree, "label"))
        graphics::text(X, Y, nodeText, xpd=lab.xpd, srt=srt, adj=adj,
            cex=lab.cex, col=lab.col, font=lab.font)
      }
    } else if (inner) {
      segmentsHV <- function(x0, y0, x1, y1) {
        if (horiz) {
          graphics::segments(y0, x0, y1, x1, col=col, lty=lty, lwd=lwd)
        } else {
          graphics::segments(x0, y0, x1, y1, col=col, lty=lty, lwd=lwd)
        }
      }
      rectHV <- function(x0, y0, x1, y1, col.rng) {
        if (horiz) {
          graphics::rect(y0, x0, y1, x1, col=col.rng, lty=lty, lwd=lwd)
        } else {
          graphics::rect(x0, y0, x1, y1, col=col.rng, lty=lty, lwd=lwd)
        }
      }
      if (type == "rectangle") {
        k1 <- 1L
        k2 <- length(subtree)
        child1 <- subtree[[k1]]
        child2 <- subtree[[k2]]
        hasE <- !is.null(ePar <- attr(child1, "edgePar"))
        if (!hasE) {
          ePar <- edgePar
        }
        i <- if (!stats::is.leaf(child1) || hasE) 1 else 2
        ## define line attributes for rectHV():
        col <- Xtract("col", ePar, default=graphics::par("col"), i)
        lty <- Xtract("lty", ePar, default=graphics::par("lty"), i)
        lwd <- Xtract("lwd", ePar, default=graphics::par("lwd"), i)
        xk1 <- if (center) mean(bx$limit[k1:(k1 + 1L)])
               else bx$limit[k1] + .midDend(child1)
        xk2 <- if (center) mean(bx$limit[k2:(k2 + 1L)])
               else bx$limit[k2] + .midDend(child2)
        rectHV(xk1, yTop, xk2,
            if (type.prox == "distance") yTop+yRng else yTop-yRng, col.rng) # h
      }
      for (k in seq_along(subtree)) {
        child <- subtree[[k]]
        ## draw lines to the children and draw them recursively
        yBot <- attr(child, "height")
        if (getOption("verbose")) {
          cat("ch.", k, "@ h=", yBot, "; ")
        }
        if (is.null(yBot)) {
          yBot <- 0
        }
        xBot <- if (center) mean(bx$limit[k:(k + 1)])
                else bx$limit[k] + .midDend(child)

        ePar <- attr(child, "edgePar")
        hasE <- !is.null(ePar)
        if (!hasE) {
          ePar <- edgePar
        }
        i <- if (!stats::is.leaf(child) || hasE) 1 else 2
        ## define line attributes for segmentsHV():
        col <- Xtract("col", ePar, default=graphics::par("col"), i)
        lty <- Xtract("lty", ePar, default=graphics::par("lty"), i)
        lwd <- Xtract("lwd", ePar, default=graphics::par("lwd"), i)
        if (type == "triangle") {
          segmentsHV(xTop, yTop, xBot, yBot)
        } else { # rectangle
          segmentsHV(xBot, yTop, xBot, yBot) # v
        }
        vln <- NULL
        if (stats::is.leaf(child) && (leaflab == "textlike")) {
          nodeText <- asTxt(attr(child, "label"))
          if (getOption("verbose")) {
            cat('-- with "label"', format(nodeText))
          }
          hln <- 0.6 * graphics::strwidth(nodeText, cex=lab.cex) / 2
          vln <- 1.5 * graphics::strheight(nodeText, cex=lab.cex) / 2
          graphics::rect(xBot-hln, yBot, xBot+hln, yBot+2*vln, col=p.col)
          graphics::text(xBot, yBot+vln, nodeText, xpd=lab.xpd, cex=lab.cex,
              col=lab.col, font=lab.font)
        }
        if (!is.null(attr(child, "edgetext"))) {
          edgeText <- asTxt(attr(child, "edgetext"))
          if (getOption("verbose")) {
            cat('-- with "edgetext"', format(edgeText))
          }
          if (!is.null(vln)) {
            mx <-
                if (type == "triangle")
                  (xTop+xBot+((xTop-xBot)/(yTop-yBot))*vln)/2
                else xBot
            my <- (yTop+yBot+2*vln)/2
          } else {
            mx <- if (type == "triangle") (xTop+xBot)/2 else xBot
            my <- (yTop+yBot)/2
          }
          ## Both for "triangle" and "rectangle" : Diamond + Text

          p.col <- Xtract("p.col", ePar, default="white", i)
          p.border <- Xtract("p.border", ePar, default=graphics::par("fg"), i)
          ## edge label pars: defaults from the segments pars
          p.lwd <- Xtract("p.lwd", ePar, default=lwd, i)
          p.lty <- Xtract("p.lty", ePar, default=lty, i)
          t.col <- Xtract("t.col", ePar, default=col, i)
          t.cex <- Xtract("t.cex", ePar, default=1, i)
          t.font <- Xtract("t.font", ePar, default=graphics::par("font"), i)

          vlm <- graphics::strheight(c(edgeText, "h"), cex=t.cex) / 2
          hlm <- graphics::strwidth(c(edgeText, "m"), cex=t.cex) / 2
          hl3 <- c(hlm[1L], hlm[1L]+hlm[2L], hlm[1L])
          if (horiz) {
            graphics::polygon(my+c(-hl3, hl3), mx+sum(vlm)*c(-1L:1L, 1L:-1L),
                col=p.col, border=p.border, lty=p.lty, lwd=p.lwd)
            graphics::text(my, mx, edgeText, cex=t.cex, col=t.col, font=t.font)
          } else {
            graphics::polygon(mx+c(-hl3, hl3), my+sum(vlm)*c(-1L:1L, 1L:-1L),
                col=p.col, border=p.border, lty=p.lty, lwd=p.lwd)
            graphics::text(mx, my, edgeText, cex=t.cex, col=t.col, font=t.font)
          }
        }
      }
    }

    if (inner && length(subtree)) {
      KK[depth] <- length(subtree)
      if (storage.mode(kk) != storage.mode(KK)) {
        storage.mode(kk) <- storage.mode(KK)
      }

      ## go to first child
      kk[depth] <- 1L
      x1 <- bx$limit[1L]
      x2 <- bx$limit[2L]
      subtree <- subtree[[1L]]
    } else {
      repeat {
        depth <- depth - 1L
        if (!depth || (kk[depth] < KK[depth])) break
      }
      if (!depth) break
      length(kk) <- depth
      kk[depth] <- k <- kk[depth] + 1L
      x1 <- llimit[[depth]][k]
      x2 <- llimit[[depth]][k+1L]
      subtree <- wholetree[[kk]]
    }
  } ## repeat
  invisible()
}

plotNodeLimit <- function(x1, x2, subtree, center) {
  ## get the left borders limit[k] of all children k=1..K, and
  ## the handle point `x' for the edge connecting to the parent.
  inner <- !stats::is.leaf(subtree) && (x1 != x2)
  limit <- c(x1,
      if (inner) {
        K <- length(subtree)
        mTop <- .memberDend(subtree)
        limit <- integer(K)
        xx1 <- x1
        for (k in 1L : K) {
          m <- .memberDend(subtree[[k]])
          ## if(is.null(m)) m <- 1
          xx1 <- xx1 + (if (center) (x2-x1)*m/mTop else m)
          limit[k] <- xx1
        }
        limit
      } else { ## leaf
        x2
      })
  mid <- attr(subtree, "midpoint")
  center <- center || (inner && !is.numeric(mid))
  x <- if (center) mean(c(x1,x2)) else x1 + (if (inner) mid else 0)
  list(x = x, limit = limit)
}

.midDend <- function(x) {
  if (is.null(mp <- attr(x, "midpoint"))) 0 else mp
}
