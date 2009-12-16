\name{summary.sample.data}
\alias{summary.sample.data}
\title{Summary of a CONIFERS sample.data object}

\description{
   These functions create and print summary results for a CONIFERS
   sample.data object.
}

\usage{
   \method{summary}{sample.data}( object, ... )
}

\arguments{
   \item{object}{an object of class \code{\link{sample.data}}.}
   \item{...}{not used by user.}
}

\value{Applying \code{summary} on an object of class
   \code{\link{sample.data}} prints out a simple summary of the
   \code{\link{sample.data}}. A \code{\link{sample.data}} object is
   the basic data type used for most CONIFERS functions.
}



\references{

Ritchie, M.W. 2008. User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
4.10. See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Martin W. Ritchie \email{mritchie@fs.fed.us} }

\seealso{\code{\link{sample.data}},
	     \code{\link{sp.sums}}
}

\examples{
library( rconifers )

## set the variant to the SWO variant
set.species.map(set.variant( 0 ))

## grow the data that was originally swo in the smc variant
# load and display CONIFERS example plots
data( plots.swo )
print( plots.swo )

# load and display CONIFERS example plants
data( plants.swo )
print( plants.swo )

# create the sample.data list object
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0 )
class(sample.swo.3)  <- "sample.data"

summary(sample.swo.3)

}

\keyword{models}



