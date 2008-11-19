

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
  The dataset is similar to the smc.txt file that is distirbuted
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

## set the variant to the SMC variant
set.variant( 1 )

# load the Stand Management Cooperative species coefficients into R as a data.frame object
data( smc )

# set the species map 
sp.map <- list(idx=smc$idx,
                 fsp=smc$fsp,
                 code=as.character(smc$code),
                 em=smc$endemic.mort,
		 msdi=smc$max.sdi,
                 b=smc$browse.damage,
                 m=smc$mechanical.damage,
                 gwh=smc$genetic.worth.h,
                 gwd=smc$genetic.worth.d)
set.species.map( sp.map )

## grow the data that was originally swo in the smc variant
# load and display CONIFERS example plots
data( plots )
print( plots )

# load and display CONIFERS example plants
data( plants )
print( plants )

## change all plants to Douglas-fir
plants$sp.code <- "DF"

# randomly remove 10 crown.width observations to test
# the impute function
blanks <- sample( 1:nrow( plants ), 10, replace=FALSE )
plants[blanks,]$crown.width <- NA

# create the sample.data list object
sample.3 <- list( plots=plots, plants=plants, age=3, x0=0.0 )
class(sample.3)  <- "sample.data"

# fill in missing values
sample.3.imp <- impute( sample.3 )

# print the maximum stand density index for the current settings
print( calc.max.sdi( sample.3.imp ) )

# print a summary of the sample
print( sample.3.imp )

# now, project the sample forward for 20 years
# with all of the options turned off
sample.23 <- project( sample.3.imp, 
	  20, 
	  control=list(rand.err=0,rand.seed=0,endemic.mort=0,sdi.mort=0))

## print the projected summaries
print( sample.23 )

## plot the diagnostics from the fit a linear dbh-tht model 
## before thinning
opar <- par( mfcol=c(2,2 ) )
plot( lm( sample.23$plants$tht ~ sample.23$plants$dbh ) )
par( opar )

## thin the stand to capture mortality, redistribute growth, 
## and possibly generate revenue
## Proportional thin for selected tree species, does not remove shrubs
sample.23.t1 <- thin( sample.23, 
	     control=list(type=1, target=50.0, target.sp="DF" ) )
print( sample.23.t1 )

## Proportional thin across all tree species
sample.23.t2 <- thin( sample.23, 
	     control=list(type=2, target=50.0 ) )
print( sample.23.t2 )

## Thin from below, by dbh, all species
sample.23.t3 <- thin( sample.23, 
	     control=list(type=3, target=50.0 ) )
print( sample.23.t3 )

## Thin from below, by dbh for species "PM"
sample.23.t4 <- thin( sample.23, 
	     control=list(type=4, target=50.0, target.sp="PM" ) )
print( sample.23.t4 )

## plot the diagnostics from the fit a linear dbh-tht model
## after proportional thinning
opar <- par( mfcol=c(2,2 ) )
plot( lm( sample.23.t2$plants$tht ~ sample.23.t2$plants$dbh ) )
par( opar )

## print the differences, by species
print( sp.sums( sample.23.t4 ) - sp.sums( sample.23 ) )

}

\keyword{datasets}


