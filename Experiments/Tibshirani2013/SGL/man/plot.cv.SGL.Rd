\name{plot.cv.SGL}
\alias{plot.cv.SGL}
\title{plots the cross-validated error curve produced by cv.SGL}
\description{
Plots the cross-validated error curve, and confidence bounds for each \code{lambda} in our regularization path.}
\usage{
\method{plot}{cv.SGL}(x, ...)
}
\arguments{
  \item{x}{fitted \code{"cv.SGL"} object}
  \item{...}{additional arguments to be passed to plot}
}
\details{
A plot is produced}
\references{Simon, N. and Tibshirani, R. (2011)
  \emph{A Sparse-Group Lasso}, \cr
  \url{http://www-stat.stanford.edu/~nsimon/SGL.pdf}
  }
\author{Noah Simon and Rob Tibshirani\cr
Maintainer: Noah Simon \email{nsimon@stanford.edu}
}
\seealso{\code{SGL} and \code{cv.SGL}.}
\examples{
n = 100; p = 300; size.groups = 5
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 2*rnorm(n)
data = list(x = X, y = y)
cvFit = cvSGL(data, index, type = "linear")
plot(cvFit)
}
\keyword{models}
\keyword{regression}
