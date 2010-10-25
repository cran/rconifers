

\keyword{datasets}
\name{smc}
\alias{smc}
\non_function{}
\title{Species table for the SMC variant of the CONIFERS growth model}
\usage{data(smc)}

\description{
These are species codes used for the Stand Management Cooperative (SMC) variant of the CONIFERS growth model.

The \code{smc} data frame has 3 rows and 13 columns.
}

\format{
This data frame contains the following columns:
\describe{

\item{idx}{The index of the species code that is used within the
  CONIFERS source code.}
\item{code}{The alpha-numeric code that is used during data input,
  normally in the field data recorder. These are used in the
  summaries and outputs rather than the species names.}
\item{fsp}{The functional species code (found in the coefficients)
  that is used to project the plant record. Every species code is
  mapped to a single functional species code. See \code{\link{set.species.map}}.}
\item{name}{The long name for the plant species.}
\item{organon}{The code that should be used to export the species to
  the ORGANON growth and yield model. See \url{http://www.cof.orst.edu/cof/fr/research/organon/index.htm} for more details}
\item{cactos}{The code that should be used to export the species to
  the CACTOS growth and yield model. See \url{http://www.cnr.berkeley.edu/~wensel/cactos/cactoss.htm} for
  more details}
\item{fvs}{The code that should be used to export the species to
  the FVS growth and yield model. See \url{http://www.fs.fed.us/fmsc/index.php} for more
  details}
\item{endemic.mort}{ is a percent for the endemic mortality, where
(0.002) means that two tenths of a percent of the trees of that
functional species type will be "killed" for each year of projection.}
\item{max.sdi}{ is the maximum stand density index the functional
  species can obtain before competition induced mortality is applied.}
\item{browse.damage}{ is a percent for the browse damage,
  where (0.002) means that two tenths of a percent of the trees of that
  species type will be "killed" as a result of browse damage.}
\item{mechanical.damage}{ is a percent for the mechanical damage,
  where (0.002) means that two tenths of a percent of the trees of that
  species type will be "killed" as a result of mechanical damage.}
\item{genetic.worth.h}{ is the average expected level of gain,
expressed as percentage, for diameter at a given age. Values are
bounded by 0 and 20.}
\item{genetic.worth.d}{ is the average expected level of gain,
expressed as percentage, for diameter at a given age. Values are
bounded by 0 and 20.}

}
}

\details{
  The dataset is similar to the smc.txt file that is distributed
  with the GUI version of CONIFERS. A species lookup table is nothing
  more than a data.frame with specific column names. 
}

\references{

Ritchie, M.W. 2008 . User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
4.10. See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

Vaughn, Nicholas. 2007. An individual-tree model to predict the annual
growth of young stands of Douglas-fir (Pseudotsuga menziesii (Mirbel)
Franco) in the Pacific northwest. M.S. Thesis, University of
Washington. 91p.

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Martin W. Ritchie \email{mritchie@fs.fed.us} }

\seealso{    \code{\link{plants.smc}},
	     \code{\link{sample.data}},
	     \code{\link{set.variant}},
	     \code{\link{swo}}
}

\examples{
library( rconifers )

## set the variant to the SMC variant
##set.species.map(set.variant( 1 ))
set.variant( 1 )

# load the SMC species coefficients into R as a data.frame object
data( smc )

## grow the smc example data
# load and display CONIFERS example smc plots
data( plots.smc )
print( plots.smc )

# load and display CONIFERS example smc plants
data( plants.smc )
print( plants.smc )

# create the sample.data list object
sample.smc <- list( plots=plots.smc, plants=plants.smc, age=3, x0=0.0 )
class(sample.smc)  <- "sample.data"

}

\keyword{datasets}


