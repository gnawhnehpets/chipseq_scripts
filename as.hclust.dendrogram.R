function (x, ...) 
{
     stopifnot(is.list(x), length(x) == 2)
     n <- length(ord <- as.integer(unlist(x)))
     unlist(x)
     iOrd <- sort.list(ord)
     iOrd
     if (!identical(ord[iOrd], seq_len(n))) 
          stop(gettextf("dendrogram entries must be 1,2,..,%d (in any order), to be coercable to \"hclust\"", 
                        n), domain = NA)
     stopifnot(n == attr(x, "members"))
     n.h <- n - 1L
     n.h
     labsu <- unlist(labels(x))
     labs <- labsu[iOrd]
     x <- .add.dendrInd(x)
     SIMP <- function(d) {
          if (is.leaf(d)) {
               -as.vector(d)
          }
          else {
               j <<- j + 1L
               height[j] <<- attr(d, "height")
               inds[[j]] <<- attr(d, ".indx.")
               attributes(d) <- NULL
               d[] <- lapply(d, SIMP)
               d
          }
     }
     height <- numeric(n.h)
     inds <- vector("list", n.h)
     j <- 0L
     xS <- SIMP(x)
     ii <- sort.list(height)
     merge <- matrix(NA_integer_, 2L, n.h)
     for (k in seq_len(n.h)) {
          if (k < n.h) {
               in.k <- inds[[ii[k]]]
               s <- xS[[in.k]]
          }
          else s <- xS
          if (getOption("as.hclust.dendr", FALSE)) {
               cat(sprintf("ii[k=%2d]=%2d -> s=xS[[in.k]]=", k, 
                           ii[k]))
               str(s)
          }
          stopifnot(length(s) == 2L, all(vapply(s, is.integer, 
                                                NA)))
          merge[, k] <- unlist(s)
          if (k < n.h) 
               xS[[in.k]] <- +k
     }
     structure(list(merge = t(merge), height = height[ii], order = ord, 
                    labels = labs, call = match.call(), method = NA_character_, 
                    dist.method = NA_character_), class = "hclust")
}