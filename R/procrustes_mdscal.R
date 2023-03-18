#' @title procrustes_mdscal
#' @description Fast multidimensional scaling with divide and conquer
#' strategy using Procrustes transformations.
#' @param X0 A matrix or data-frame containing \eqn{n} observations of the
#' explanatory variables, with dimension \eqn{p} of at least 2.
#' @param k Desired dimension in the reduced space, default is 2.
#' @param diss Type of metric to produce the dissimilarities between observations.
#' \texttt{minkowski} works with Minkowski metric of order \texttt{power},
#' and \texttt{gower} works with Gower metric, using specific weights for variables.
#' @param power Power of the Minkowski metric. Default is 2 (Euclidean).
#' @param wght A vector of \eqn{p} weights, for the original variables in the
#' Gower metric. Default is 1 for each variable.
#' @param l Maximum size of the sub-groups in the splitting procedure. Default is 400.
#' @param c Number of alignment points, common to all the sub-groups. Default is 4.
#' @param n_cores Number of CPU cores to be used in the parallelization. Default is 1.
#' This option is not available in Windows.
#' @param seed A seed for the random selection of points to split the data.
#' @details This functions implements DCO-MDS (\textit{Divide and conquer multidimensional scaling}))
#' with Procrustes transformations to align the different groups.
#' The scaling method is Kruskal distance scaling by minimizing Stress-1 with the SMACOF algorithm,
#' using splines as a mapping procedure of the dissimilarities.
#' The results are rotated to the principal axis. Two options
#' are available for the dissimilarities: Minkowski metric and Gower metric.
#' There is a parallelization of the main procedure of scaling and aligning the sub-groups.
#' @return A data matrix of \eqn{n} rows and \eqn{k} columns, with the coordinates
#' of points in the reduced space.
#' @examples
#' x <- MASS::mvrnorm(n=5000,mu=c(0,0,0,0),diag(c(1,2,3,4)))
#'
#' Z <- procrustes_mdscal(x, k=2,diss='minkowski', power=1.5, seed = 16497)
#' plot(Z, pch=20, main = 'Fast multidimensional scaling k=2')
#'
#' #' @references
#' Borg, I., & Groenen, P. J. F. (2005). Modern multidimensional scaling: theory and applications. Springer Science & Business Media.
#' Pedro Delicado y Cristian Pachon-Garcia. “Multidimensional Scaling for Big Data”. En: arXiv preprint arXiv:2007.11919 (2020).
#'
#' @rdname procrustes_mdscal
#' @export
procrustes_mdscal <- function(X0,k=2,diss,
                              power=NULL,wght=NULL,l=400,c=4,n_cores=1,seed=NULL){
  require(smacof)
  X0 <- as.matrix(X0) ; n <- nrow(X0)


  #definition of the dissimilarity function
  {
    diss <- match.arg(diss,choices=c('minkowski','gower'))
    switch(diss,
           minkowski = {
             if(is.null(power)){power <- 2}
             dissimilarity <- function(X0){
               dist(X0,method='minkowski',p=power)
             }
           },
           gower = {
             if(is.null(wght)){wght <- rep.int(1,ncol(X0))}
             ranges <- apply(X0,2,function(x) if(is.numeric(x)){max(x) - min(x)}
                             else {NA})
             dissimilarity <- function(X0){
               as.dist(StatMatch::gower.dist(X0,rngs=ranges,
                                             var.weights = wght))
             }
           })

  }

  #function to get partitions
  get_partitions <- function(n, l, c, k) {

    if (l-c <= 0) {
      stop("\"l\" must be greater than \"c\"")
    } else if (l-c <= c) {
      stop("\"l-c\" must be greater than \"c\"")
    } else if (l-c <= k) {
      stop ("\"l-c\" must be greater than \"r\"")
    }

    set.seed(seed); permutation <- sample(x = n, size = n, replace = FALSE)

    if (n<=l) {
      list_indexes <- list(permutation)
    } else {
      min_sample_size <- max(k+2, c)
      p <- 1 + ceiling((n-l)/(l-c))
      last_partition_sample_size <- n - (l + (p-2) * (l-c))

      if (last_partition_sample_size < min_sample_size & last_partition_sample_size > 0) {
        p <- p - 1
        last_partition_sample_size <- n - (l + (p-2) * (l-c))
      }

      first_partition <- permutation[1:l]
      last_partition <- permutation[(n-last_partition_sample_size+1):n]
      list_indexes <- split(x = permutation[(l+1):(n-last_partition_sample_size)],
                            f = 1:(p-2))
      names(list_indexes) <- NULL
      list_indexes[[p-1]] <- list_indexes[[1]]
      list_indexes[[p]] <- last_partition
      list_indexes[[1]] <- first_partition
    }

    return(list_indexes)
  }

  #function for alignment with procrustes transformation
  alignment <- function(X1, X2, Xa){

    n<-nrow(X1)

    #Procrustes transformation
    H <- diag(rep(1,n))-(1/n)*matrix(1,n,1) %*% t(matrix(1,n,1))
    C <- t(X1) %*% H %*% X2
    svd <- svd(C)

    R <- svd$v %*% t(svd$u) #rotation matrix
    s <- sum(diag(C %*% R))/sum(diag(t(X2) %*% H %*% X2)) #dilation factor
    t <- (1/n)*t(X1 - s * X2 %*% R) %*% matrix(1,n,1) #translation

    return(s * Xa %*% R + matrix(1,nrow(Xa),1) %*% t(t))
  }

  #create 'empty' configuration
  X <- matrix(0,nrow = nrow(X0),ncol = k)

  #list of indexes
  idx <- get_partitions(nrow(X0),l,c,k)

  #select alignment points
  idm <- idx[[1]][1:c]

  #first MDS
  X[idx[[1]],] <- smacof::mds(dissimilarity(X0[idx[[1]],]),
                              ndim=k,type='mspline',
                              principal=TRUE)$conf

  #scaling and alignment of the rest with parallel computation
  results_list <- parallel::mclapply(idx[-1],
                                     function(i) {

                                       Xj <- smacof::mds(dissimilarity(X0[c(idm,i),]),
                                                         ndim=k,type='mspline',
                                                         principal=TRUE)$conf

                                       alignment(X[idm,],Xj[1:c,],Xj[-(1:c),])
                                     },
                                     mc.cores = n_cores)

  #final configuration matrix
  X[Reduce(c,idx[-1]),] <- Reduce(rbind,results_list)
  X <- X %*% prcomp(X)$rotation
  colnames(X) <- paste('Z',1:k,sep="")

  return(X)
}
