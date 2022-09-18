#' @title ss.split
#' @description Returns the Within-groups sum of squares and
#' Between-group sum of squares of a data set, given a grouping variable.
#' @param x A matrix or data-frame containing \eqn{n} observations of the
#' explanatory variables,
#' with dimension \eqn{p} of at least 2.
#' @param grouping A factor specifying the class for each observation.
#' @details Given a data matrix \eqn{\mathbf{X} \in \mathbb{R}^{n \times p}} with
#' global sample mean \eqn{\overline{\mathbf{x}}}, containing
#' \eqn{i=\{1,2,...,k\}} different groups, with sample sizes
#' \eqn{n_1, n_2, ..., n_k} respectively,
#' and their \eqn{\overline{\mathbf{x}_i}} sample means,
#' the Total sum of squares:
#'
#'   \eqn{\mathbf{T} = \sum_{j=1}^n (\mathbf{x}_j -
#'   \overline{\mathbf{x}})(\mathbf{x}_j - \overline{\mathbf{x}})^T}
#'
#'   is decomposed into the Within sum of squares
#'
#'   \eqn{\mathbf{U} = \sum_{i=1}^k \sum_{j=1}^{n_i} (\mathbf{x} -
#'   \overline{\mathbf{x}}_i)(\mathbf{x} - \overline{\mathbf{x}}_i)^T}
#'
#'   and the between sum of squares
#'
#'   \eqn{\mathbf{H} = \sum_{i=1}^k n_i (\overline{\mathbf{x}}_i -
#'   \overline{\mathbf{x}})(\overline{\mathbf{x}}_i -
#'   \overline{\mathbf{x}})^T}.
#'
#'
#'  \eqn{\mathbf{T} = \mathbf{U} + \mathbf{H}}.
#'
#' @return A list containing three \eqn{p \times p} numeric arrays.
#' \code{total.SS} contains \eqn{\mathbf{T}}, \code{within.SS}
#' contains \eqn{\mathbf{U}}and \code{between.SS}
#' contains \eqn{\mathbf{H}}.
#' @examples
#' set.seed(16497)
#' x1 <- MASS::mvrnorm(n=10,mu=c(0,0,0),diag(c(1,1,1)))
#' x2 <- MASS::mvrnorm(n=12,mu=c(2,2,2),diag(c(1,1,1)))
#' x3 <- MASS::mvrnorm(n=10,mu=c(4,4,4),diag(c(1,1,1)))
#' x <- rbind.data.frame(x1,x2,x3)
#' grouping <- factor(c(rep(1,10), rep(2,12), rep(3,10)))
#' names(x) <- c('Var1','Var2','TheThird')
#'
#' split <- ss.split(x,grouping)
#' split$total.SS
#' split$within.SS
#' split$between.SS
#' @rdname ss.split
#' @export
ss.split<-function(x, grouping){

x <- split(data.frame(x),grouping)

k <- length(x)
p <- ncol(x[[1]])
if(p==1) {stop('At least two-dimensional data is required')}
n <- sapply(x,nrow)
N <- sum(n)

#Means by group
XR <- Reduce('+',
             lapply(x,function(i) nrow(i)*colMeans(i)))/sum(sapply(x,nrow))

#Within groups sum of squares
U <- Reduce('+',lapply(x,function(i) (nrow(i)-1)*cov(i)))

#Between groups sum of squares
H <- Reduce('+',lapply(x,function(i) nrow(i)*(colMeans(i)-XR) %*%
                         t((colMeans(i)-XR))))

rownames(H) <- colnames(H)

return(list('total.SS' = U+H,
            'within.SS' = U,
            'between.SS' = H))
}

