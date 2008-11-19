

\keyword{datasets}
\name{plots}
\alias{plots}
\non_function{}
\title{Example plots data set for CONIFERS growth model}
\usage{data(plots)}

\description{
This a simple plots \code{\link{data.frame}} that contains basic plot attributes used in the CONIFERS growth model.

The \code{plots} data frame has 4 rows and 5 columns.
}

\format{
This data frame contains the following columns:
\describe{
  \item{plot}{A numeric plot identifier.}
  \item{elevation}{Plot elevation, in feet.}
  \item{slope}{Plot slope, in percent.}
  \item{aspect}{Plot aspect, in degrees from North.}
  \item{whc}{Water holding capacity of the plot, in inches.}
  \item{map}{Mean annual precipitation, in inches.}
}
}

\details{The example plot data represents plot level observations for
a field sample that can be used to generate future forest conditions
using the CONIFERS forest growth model. The data are in imperial
units. The required values are the \code{plot}, \code{elevation}, and
\code{slope}, and \code{aspect} measurements. Values that have not
been observed, or are missing (e.g. whc and map), are denoted with an
\code{NA} and can be imputed using the \code{\link{impute}} function}


\references{

Ritchie, M.W. 2008. User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
4.10. See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Martin W. Ritchie \email{mritchie@fs.fed.us} }

\seealso{    \code{\link{calc.max.sdi}},
	     \code{\link{impute}},
	     \code{\link{plants}},	
	     \code{\link{plots}} 
	     \code{\link{project}},
	     \code{\link{rand.seed}},
	     \code{\link{rconifers}},
	     \code{\link{sample.data}},	     		  
	     \code{\link{set.species.map}},		  
	     \code{\link{set.variant}},		  
 	     \code{\link{smc}}, 
 	     \code{\link{summary.sample.data}},
	     \code{\link{swo}}, 
	     \code{\link{thin}}
}


\examples{
library( rconifers )
data(plots)
print( plots )
}
