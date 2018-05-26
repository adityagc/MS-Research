\name{predict.SGL}
\alias{predict.SGL}
\title{Outputs predicted responses from an SGL model for new observations}
\description{
Outputs predicted response values for new user input observations at a specified \code{lambda} value}
\usage{
\method{predict}{SGL}(x, newX, lam)
}
\arguments{
  \item{x}{fitted \code{"SGL"} object}
  \item{newX}{covariate matrix for new observations whose responses we wish to predict}
  \item{lam}{the index of the lambda value for the model with which we desire to predict}
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
Fit = SGL(data, index, type = "linear")
X.new = matrix(rnorm(n * p), ncol = p, nrow = n)
predict(Fit, X.new, 5)
}
\keyword{models}
\keyword{regression}
