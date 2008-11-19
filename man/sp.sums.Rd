
\name{sp.sums}
\alias{sp.sums}
\title{A summary of a CONIFERS sample.data object by species}

\description{ This function returns a \code{\link{data.frame}} object
   that contain common summaries by species for
   \code{\link{sample.data}} object.  
}

\usage{
   sp.sums( x )
}

\arguments{
   \item{x}{an object of class \code{\link{sample.data}}.}
}

\format{

The \code{\link{data.frame}} object returned from sp.sums contains the
following columns for each species:

\describe{
\item{expf}{The total plants per acre.}
\item{tht}{The mean total height.}
\item{ba}{The total basal area for all plants.}
\item{qmd}{The quadratic mean diameter for all trees. See Details.}
}
}

\details{ This function returns a \code{\link{data.frame}} object that
	contains species level summary information. It is intended
	demonstration only and users are encouraged to examine the
	source code and modify. All results include all stems and not
	only those above breast height.
}

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

## set the variant to the SWO variant
set.variant( 0 )

# load the Southwest-Oregon species coefficients into R as a data.frame object
data( swo )

# set the species map 
sp.map <- list(idx=swo$idx,
                 fsp=swo$fsp,
                 code=as.character(swo$code),
                 em=swo$endemic.mort,
		 msdi=swo$max.sdi,
                 b=swo$browse.damage,
                 m=swo$mechanical.damage,
                 gwh=swo$genetic.worth.h,
                 gwd=swo$genetic.worth.d)
set.species.map( sp.map )

## grow the data that was originally swo in the smc variant
# load and display CONIFERS example plots
data( plots )
print( plots )

# load and display CONIFERS example plants
data( plants )
print( plants )

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

\keyword{models}
