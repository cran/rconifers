###
###	$Id: rconifers.r 613 2008-11-25 23:14:14Z hamannj $	
###
###            R interface package for conifers growth model
###
### Copyright 2003-2008 Jeff D. Hamann <jeff.hamann@forestinformatics.com>
###
### This file is part of the rconifers library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


## To submit to CRAN:
## $ ch to conifers directory
## $ R CMD CHECK rconifers
## $ R CMD BUILD rconifers
## $ ftp cran.r-project.org
## Connected to cran.wu-wien.ac.at.
## 220 Welcome to the CRAN FTP service.
## Name (cran.r-project.org:hamannj): anonymous
## 331 Please specify the password.
## Password:
## 230-Welcome, CRAN useR!
## 230-
## 230-If you have any unusual problems,
## 230-please report them via e-mail to <cran-sysadmin@statmath.wu-wien.ac.at>.
## 230-
## 230 Login successful.
## Remote system type is UNIX.
## Using binary mode to transfer files.
## ftp> cd incoming
## 250 Directory successfully changed.
## ftp> put rconifers_0.0-9.tar.gz
## local: rconifers_0.0-9.tar.gz remote: rconifers_0.0-9.tar.gz
## 229 Entering Extended Passive Mode (|||58381|)
## 150 Ok to send data.
## 100% |**************************************************************************************************************| 99741     146.79 MB/s    00:00 ETA
## 226 File receive OK.
## 99741 bytes sent in 00:04 (24.20 KB/s)
## ftp> quit
## 221 Goodbye.
## $



.First.lib <-
  function(lib, pkg)
{
  ## this package doesn't need any additional package...
  ## if it did, you'd modify the code below...
  ## but for now, simply display this little message
  ## or here?

  library.dynam("rconifers", pkg, lib)

  ## put package checks here...
###     if(!require(ts, quietly = TRUE))
###         stop("Package ts is needed.  Stopping")
    
  mylib <- dirname(.path.package("rconifers"))
  ver <- packageDescription("rconifers", lib = mylib)["Version"]

  vertxt <- paste("\n`rconifers' version:", ver, "\n")

  introtxt <-
    paste("\n`rconfiers' is a package that provides an alternative interface\n",
          "to the CONIFERS forest growth model.\n\n",
          "The user is able to define arbitrary silvicultural prescriptions,\n",
          "using standard R scripts, to predict near term (<20 years) future\n",
          "forest conditions under a variety of model configurations and\n",
          "management decisions.\n\n",
          "The output can be used in standard graphics, file input and output,\n",
          "and data analysis functions using R.\n\n",
          "Jeff D. Hamann, Forest Informatics, Inc. and \n",
          "Martin W. Ritchie, USFS PSW Research Station, Redding  Silviculture Lab\n\n",
          "See `library (help=rconifers)' for details.\n\n",
          sep = "")

  if(interactive() || getOption("verbose")) {
    cat(paste(vertxt, introtxt))
  }
  
  if( .Call( "r_set_variant", 0, PACKAGE="rconifers" ) ) {
    stop("Package could not be loaded. Could not load default coefficients. Stopping...")
  }
  
  ## set the species mappings
  data( swo )
  sp.map <- list(idx=swo$idx,
                 fsp=swo$fsp,
                 code=as.character(swo$code),
                 em=swo$endemic.mort,
                 msdi=swo$max.sdi,
                 b=swo$browse.damage,
                 m=swo$mechanical.damage,
                 gwh=swo$genetic.worth.h,
                 gwd=swo$genetic.worth.d)

  if( !.Call( "r_set_species_map", sp.map, verbose=FALSE,PACKAGE="rconifers" ) ) {
    stop("Package could not be loaded. Could not load swo species map. Stopping...")
  } else {
    introtxt <- paste("\nInitialized species map using data(swo)\nType help.start() for documentation\n", sep = "")
    if(interactive() || getOption("verbose")) {
      cat( introtxt)
    }
  }  
  
  ## do I register my routines here?

}

.Last.lib <- function( lib ) {
  .Call( "exit_conifers", PACKAGE="rconifers" )
  library.unload( "rconifers" )  
}

rand.seed <- function( control ) {
  val <- .Call( "r_reseed", control, PACKAGE="rconifers" )
  val
}


## this function will set the species codes (internal)
set.species.map <- function( sp.map=list(idx,fsp,code,em,msdi,b,m,gwh,gwd) ) {
  
  ## verify the lengths match
  if( length( sp.map$idx ) != length( sp.map$fsp ) ) {
    stop( "The lengths of the index and functional species vectors do not match." )
    return
  }

  val <- .Call( "r_set_species_map", sp.map, verbose=FALSE, PACKAGE="rconifers" )
  ## no need to return the number of species assigned
}

set.variant <- function( var=0 ).Call( "r_set_variant", var, PACKAGE="rconifers" )

## this function, which is used internally, to convert data from the c code
build.sample.data <- function( x ) {

  ## build the two data.frame objects and add them to the list
  ret.val <- vector( "list" )
  ret.val$x0 <- x[[1]]
  ret.val$age <- x[[2]]

  plots <-   as.data.frame( cbind( plot=x[[3]][[1]] ,
                                  elevation=x[[3]][[2]],
                                  slope=x[[3]][[3]],
                                  aspect=x[[3]][[4]],
                                  whc=x[[3]][[5]],
                                  map=x[[3]][[6]]
                                  )
                           )
  plants <-   as.data.frame( cbind( plot=x[[4]][[1]] ,
                                   sp.code=x[[4]][[2]],
                                   d6=x[[4]][[3]],
                                   dbh=x[[4]][[4]],
                                   tht=x[[4]][[5]],
                                   cr=x[[4]][[6]],
                                   n.stems=x[[4]][[7]],
                                   expf=x[[4]][[8]],
                                   crown.width=x[[4]][[9]],
                                   errors=x[[4]][[10]]
                                   )
                            )

  plants$d6 <- as.numeric( levels( plants$d6 ) )[plants$d6]
  plants$dbh <- as.numeric( levels( plants$dbh ) )[plants$dbh]
  plants$tht <- as.numeric( levels( plants$tht ) )[plants$tht]
  plants$cr <- as.numeric( levels( plants$cr ) )[plants$cr]
  plants$n.stems <- as.numeric( levels( plants$n.stems ) )[plants$n.stems]
  plants$expf <- as.numeric( levels( plants$expf ) )[plants$expf]
  plants$crown.width <- as.numeric( levels( plants$crown.width))[plants$crown.width]
  plants$errors <- as.numeric( levels( plants$errors))[plants$errors]

  ret.val$plots <- plots
  ret.val$plants <- plants

  class(ret.val)  <- "sample.data"
  ret.val

}

################################################################################
## file i/o interfaces for C code
################################################################################
project <- function( x,
                    years=1,
                    control=list(rand.err=0,
                      rand.seed=0,
                      endemic.mort=0,
                      sdi.mort=0) )
{

  if( class( x ) != "sample.data" ) {
    stop( "x is not a sample.data object." )
    return
  }
  
  x$plants$sp.code <- as.character( x$plants$sp.code )
  val <- build.sample.data( .Call( "r_project_sample",  x, years, control, PACKAGE="rconifers" ) )
  val
}


################################################################################
## file i/o interfaces for C code
################################################################################

thin <- function( x,
                 control=list(type=2,
                   target=50.0,
                   target.sp=NULL ) )
{

  if( class( x ) != "sample.data" ) {
    stop( "x is not a sample.data object." )
    return
  }

  if( control$type < 1 || control$type > 4 ) {
    stop( "control$type must be 1, 2, 3, or 4. See help." )
    return    
  }
  
  if( control$type %in% c(2,3) & !is.null( control$target.sp ) ) {
    stop( "If you want to thin all species, sp=NULL." )
    return    
  }

  if( control$type %in% c(1,4) & is.null( control$target.sp ) ) {
    stop( "If you want to thin a particular species, sp != NULL." )
    return    
  }

  ## make sure the sp.code values are not factors.
  x$plants$sp.code <- as.character( x$plants$sp.code )
  val <- build.sample.data( .Call( "r_thin_sample", x, control, PACKAGE="rconifers" ) )
  val
}

## fill in missing values wrapper function
impute <- function( x,
                   control=list(fpr=11.78,
                     min.dbh=5.6,
                     baf=40.0) )
{
  
  if( class( x ) != "sample.data" ) {
    stop( "x is not a sample.data object." )
    return
  }
  
  x$plants$sp.code <- as.character( x$plants$sp.code )
  val <- build.sample.data( .Call( "r_fill_in_missing_values",
                                  x,
                                  control,
                                  PACKAGE="rconifers" ) )
  val
  
}


calc.max.sdi <- function( x )
{
  if( class( x ) != "sample.data" ) {
    stop( "x is not a sample.data object." )
    return
  }
  
  x$plants$sp.code <- as.character( plants$sp.code )
  val <- .Call( "r_calc_max_sdi", x, PACKAGE="rconifers" )
  
  val
}


################################################################################
## print a few results of the whole system
print.sample.data <- function( x, digits = max( 3, getOption("digits") - 1 ),... ) {

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))
  
  cat("\n")
  cat( "sample contains", nrow(x$plots), "plots records\n" )
  cat( "sample contains", nrow(x$plants), "plant records\n" )

  ##cat( "years projected = ", x$prj.yrs, "\n")
  cat( "age = ", x$age, "\n")
  cat( "x0 = ", x$x0, "\n")
  cat( "max sdi = ", calc.max.sdi( x ), "\n")

  ## print out some species summaries...
  print( sp.sums( x ) )

  if( sum( x$plants$errors ) ) {

    ## iterate over the plant records and report how many of what kinds
    ## of errors there are in the plant list...
    ## the #defines in conifers.h contain the masks for the type
    ## of errors
    ## this also means that you have to have a dependent package (bitops)
    
    print( "plant records contain errors..." )
  }
  
  ##   print( stand.table( x ) )  
  invisible( x )
}


summary.sample.data <- function(object,...) {
  summary.sample.data <- object
  ##summary.sample.data$stand.table <- stand.table( object )
  summary.sample.data
}





## split and sum the plants by species
sp.sums <- function( x ) {

  if( class( x ) != "sample.data" ) {
    stop( "x is not a sample.data object." )
    return
  }

  x.by.sp <- split( x$plants, f=list(x$plants$sp.code ) )
  
  ## this function computes the summaries each species.
  sp.sums.f <- function( x, npts ) {
    expf <- sum( x$expf * x$n.stems ) / npts
    tht <- sum( x$expf * x$n.stems * x$tht ) / sum( x$expf * x$n.stems )
    ba <- sum( x$expf * x$n.stems * x$dbh^2*0.0054541539 ) / npts
    qmd <- sqrt( ba / expf / 0.0054541359 )
    sp.sums <- c(qmd,tht,ba,expf)
  }
  
  df <- as.data.frame( t(sapply( x.by.sp, sp.sums.f, nrow(x$plots) )) )
  names( df ) <- c("qmd","tht","ba","expf")

  df
}

  


