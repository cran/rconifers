
/* 	$Id: rconifers.c 674 2010-10-25 22:14:00Z hamannj $	 */

/* To submit to CRAN: */
/* $ ch to conifers directory */
/* $ R CMD CHECK rconifers */
/* $ R CMD BUILD rconifers */
/* $ ftp cran.r-project.org */
/* Connected to cran.wu-wien.ac.at. */
/* 220 Welcome to the CRAN FTP service. */
/* Name (cran.r-project.org:hamannj): anonymous */
/* 331 Please specify the password. */
/* Password:  */
/* 230-Welcome, CRAN useR! */
/* 230- */
/* 230-If you have any unusual problems, */
/* 230-please report them via e-mail to <cran-sysadmin@statmath.wu-wien.ac.at>. */
/* 230- */
/* 230 Login successful. */
/* Remote system type is UNIX. */
/* Using binary mode to transfer files. */
/* ftp> cd incoming */
/* 250 Directory successfully changed. */
/* ftp> put rconifers_0.0-9.tar.gz  */
/* local: rconifers_0.0-9.tar.gz remote: rconifers_0.0-9.tar.gz */
/* 229 Entering Extended Passive Mode (|||58381|) */
/* 150 Ok to send data. */
/* 100% |**************************************************************************************************************| 99741     146.79 MB/s    00:00 ETA */
/* 226 File receive OK. */
/* 99741 bytes sent in 00:04 (24.20 KB/s) */
/* ftp> quit */
/* 221 Goodbye. */
/* $  */



/* and send email to cran@r-project.org.  */
/* Please indicate the copyright situation (GPL, ...) in your submission.  */




#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "conifers.h"

/* R header files */
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>


/************************************************************************/
/* global declarations       						*/
/************************************************************************/
unsigned long current_variant;
//unsigned long use_genetic_gains;

double COEFFS_VERSION;
unsigned long N_COEFFS;
struct COEFFS_RECORD *COEFFS_PTR;

unsigned long N_SPECIES;
struct SPECIES_RECORD *SPECIES_PTR;

/************************************************************************/
/* function declarations       						*/
/************************************************************************/

/* helper functions */
SEXP get_list_element(SEXP list, char *str);
SEXP getvar(SEXP name, SEXP rho);

/* new functions for next version of the library */
SEXP r_init_variant( SEXP variant_sexp );
SEXP r_set_variant( SEXP variant_sexp );
/* SEXP r_get_mort_settings(); */
/* SEXP r_set_mort_settings(  SEXP map_sexp, SEXP verbose_sexp ); */
SEXP r_set_species_map(  SEXP map_sexp, SEXP verbose_sexp  );
SEXP r_reseed( SEXP ctl );
SEXP r_project_sample( SEXP data_sexp, SEXP n_years_sexp, SEXP ctl_sexp );
SEXP r_thin_sample( SEXP data_sexp,   SEXP ctl_sexp );
SEXP r_fill_in_missing_values( SEXP data_sexp, SEXP ctl_sexp ) ;
SEXP r_calc_max_sdi( SEXP data_sexp );

/* these functions are used to convert the plots between the two interfaces */
struct PLOT_RECORD *build_plot_array_from_sexp( SEXP plot_sexp, 
						unsigned long *n_plots );
SEXP build_sexp_from_plot_array( unsigned long n_plots, 
				 struct PLOT_RECORD *plots_ptr );

struct PLANT_RECORD *build_plant_array_from_sexp( SEXP plant_sexp, 
						  unsigned long *n_plants );
SEXP build_sexp_from_plant_array( unsigned long n_plants, 
				  struct PLANT_RECORD *plants_ptr );
SEXP build_return_data_sexp( double x0,
			     unsigned long age,
			     unsigned long n_plots,
			     struct PLOT_RECORD *plots_ptr,
			     unsigned long n_plants,
			     struct PLANT_RECORD *plants_ptr );


/* these are called when the library is loaded and unloaded */
SEXP init_conifers();
SEXP exit_conifers();


/************************************************************************/
/* function definitions							*/
/************************************************************************/

void R_init_rconifers(DllInfo *info)
{
   /* Register routines, allocate resources. */
/*    Rprintf( "Register routines, allocate resources\n" ); */

}


void R_unload_rconifers(DllInfo *info)
{
   /* Release resources. */
/*   Rprintf( "Release resources.\n" ); */
}


/************************************************************************/
/* data interface functions						*/
/************************************************************************/
SEXP r_reseed( SEXP ctl )
{
   unsigned long use_random_error = 0;
   unsigned long random_seed = 0;

   SEXP ans;
   PROTECT(ans = allocVector(INTSXP, 1));

   /* which one is the proper way to perform this operation */
   /* intitialize the config variables */
   use_random_error  = asInteger( get_list_element( ctl, "use.random.error" ) );
   random_seed  = asInteger( get_list_element( ctl, "random.seed" ) );

   if( random_seed == 0 )
   {
      Rprintf( "Using the clock to seed the number generator\n" );
      srand( (unsigned)time( NULL ) );
   }
   else
   {
      Rprintf( "Seeding the random number generator with %d\n", random_seed );
      srand( random_seed );
   }

   INTEGER(ans)[0] = 1;
   UNPROTECT( 1 );

   return ans;
}



/* this function initializes the coeffs */
SEXP r_set_variant( SEXP variant_sexp )
{
   SEXP ans;

   long variant = asInteger(variant_sexp);

   PROTECT(ans = allocVector(INTSXP, 1));

   if( variant == CONIFERS_SWO || variant == CONIFERS_SMC )
   {
      /* the coefficients section */
      if( COEFFS_PTR )
      {
	    free( COEFFS_PTR );
	    COEFFS_PTR = NULL;
	    N_COEFFS = 0;
      }   
      
      /* load the coeffs for the variant */
      COEFFS_PTR  = con_init_coeffs( variant, &N_COEFFS, &COEFFS_VERSION );
   }
   else
   {
      error( "That variant is invalid\n" );
      error( "Only variant zero (0) or one (1) is valid\n" );
      INTEGER(ans)[0] = 1;
      UNPROTECT(1);
      return ans;
   }


   if( !COEFFS_PTR )
   {
      error( "The coefficients were not initialized\n" );
      INTEGER(ans)[0] = 1;
      UNPROTECT(1);
      return ans;
   }

   if( COEFFS_PTR )
   {
      Rprintf( "Initialized %ld functional species coefficients for variant # %ld\n", N_COEFFS, variant );
      Rprintf( "coefficients version is %lf\n", COEFFS_VERSION );
   }

   /* set the global variable for the variant */
   current_variant = variant;

   /* you should perform a check to make sure that all the	*/
   /* species have a functional species in the coeffs array	*/

   /* return the current_variant */
   INTEGER(ans)[0] = current_variant;
   UNPROTECT(1);
   return ans;

}


/* i think you should convert this read from a data.frame object */
SEXP  r_set_species_map(  SEXP map_sexp, 
			  SEXP verbose_sexp )
{
   unsigned long i;

   int verbose;

   SEXP ans;

   SEXP idx_sexp;
   SEXP fsp_sexp;
   SEXP code_sexp;
   SEXP endemic_mort_sexp;
   SEXP max_sdi_sexp;
   SEXP browse_damage_sexp;
   SEXP mechanical_damage_sexp;
   SEXP genetic_worth_h_sexp;
   SEXP genetic_worth_d_sexp;

   /* protect in incomming list */
   PROTECT( map_sexp = AS_LIST( map_sexp ) );
   PROTECT(ans = allocVector(INTSXP, 1));

   /* added verbose argument to function */
   PROTECT( verbose_sexp = AS_INTEGER( verbose_sexp ) );
   verbose = asInteger(verbose_sexp);

   /* extract the list elements */
   idx_sexp = get_list_element( map_sexp, "idx" );
   fsp_sexp = get_list_element( map_sexp, "fsp" );
   code_sexp = get_list_element( map_sexp, "code" );
   endemic_mort_sexp = get_list_element( map_sexp, "em" );
   max_sdi_sexp = get_list_element( map_sexp, "msdi" );
   browse_damage_sexp = get_list_element( map_sexp, "b" );
   mechanical_damage_sexp = get_list_element( map_sexp, "m" );
   genetic_worth_h_sexp = get_list_element( map_sexp, "gwh" );
   genetic_worth_d_sexp = get_list_element( map_sexp, "gwd" );

   /* coerce the vectors */
   PROTECT( idx_sexp = coerceVector( idx_sexp, INTSXP ) );
   PROTECT( fsp_sexp = coerceVector( fsp_sexp, INTSXP ) );
   PROTECT( code_sexp = coerceVector( code_sexp, STRSXP ) );
   PROTECT( endemic_mort_sexp = coerceVector( endemic_mort_sexp, REALSXP ) );
   PROTECT( max_sdi_sexp = coerceVector( max_sdi_sexp, REALSXP ) );
   PROTECT( browse_damage_sexp = coerceVector( browse_damage_sexp, REALSXP ) );
   PROTECT( mechanical_damage_sexp = coerceVector( mechanical_damage_sexp, REALSXP ) );
   PROTECT( genetic_worth_h_sexp = coerceVector( genetic_worth_h_sexp, REALSXP ) );
   PROTECT( genetic_worth_d_sexp = coerceVector( genetic_worth_d_sexp, REALSXP ) );

   if( SPECIES_PTR )
   {
      free( SPECIES_PTR );
      N_SPECIES = 0;
   }
   
   N_SPECIES = length(idx_sexp);
   SPECIES_PTR = (struct SPECIES_RECORD *)calloc( 
      N_SPECIES, sizeof( struct SPECIES_RECORD ) );
   if( !SPECIES_PTR )
   {
      N_SPECIES = 0;
      INTEGER(ans)[0] = 0;
      UNPROTECT(12);
      return ans;
   }

   /* assign species mappings */
   for( i = 0; i < N_SPECIES; i++ )
   {
      SPECIES_PTR[i].idx = INTEGER( idx_sexp )[i];
      SPECIES_PTR[i].fsp_idx = INTEGER( fsp_sexp )[i];
      strcpy( SPECIES_PTR[i].sp_code, CHAR( STRING_ELT( code_sexp, i ) ) );

      SPECIES_PTR[i].endemic_mortality = REAL( endemic_mort_sexp )[i];
      SPECIES_PTR[i].max_sdi = REAL( max_sdi_sexp )[i];
      SPECIES_PTR[i].browse_damage = REAL( browse_damage_sexp )[i];
      SPECIES_PTR[i].mechanical_damage = REAL( mechanical_damage_sexp )[i];
      SPECIES_PTR[i].genetic_worth_h = REAL( genetic_worth_h_sexp )[i];
      SPECIES_PTR[i].genetic_worth_d = REAL( genetic_worth_d_sexp )[i];

      if( verbose )
      {
	 Rprintf( "%i,%i,%s,%f\n",
	 SPECIES_PTR[i].idx,
	 SPECIES_PTR[i].fsp_idx,
	 SPECIES_PTR[i].sp_code,
	 SPECIES_PTR[i].max_sdi
	    );
      }

   }

   /* now sort the species back to the "native" order (by index) */
   qsort(  (void*)SPECIES_PTR, 
	   (size_t)(N_SPECIES), 
	   sizeof( struct SPECIES_RECORD ),
	   compare_species_by_idx );

   INTEGER(ans)[0] = N_SPECIES;
   UNPROTECT(12);
   return ans;

}


SEXP exit_conifers()
{

   if( SPECIES_PTR )
   {
      Rprintf( "freeing the species array\n" );
      free( SPECIES_PTR );
      N_SPECIES = 0;
   }

   if( COEFFS_PTR )
   {
      Rprintf( "freeing the coefficients array\n" );
      free( COEFFS_PTR );
      N_COEFFS = 0;
      COEFFS_VERSION = 0.0;
   }

   return R_NilValue;

}

/* get the list element named str, or return NULL */
SEXP get_list_element(SEXP list, char *str)
{

   int i;
   SEXP elmt = R_NilValue;
   SEXP names = getAttrib(list, R_NamesSymbol);
   for (i = 0; i < length(list); i++) 
   {
      if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) 
      {
	 elmt = VECTOR_ELT(list, i);
	 break;
      }
   }

   return elmt;
}


SEXP getvar(SEXP name, SEXP rho)
{
   SEXP ans;
   if(!isString(name) || length(name) != 1) 
   {
      error("name is not a single string");
   }
   
   if(!isEnvironment(rho)) 
   {
      error("rho should be an environment");
   }

   ans = findVar(install(CHAR(STRING_ELT(name, 0))), rho);

   return(R_NilValue);
}


SEXP r_project_sample( 
   SEXP data_sexp,
   SEXP n_years_sexp,
   SEXP ctl_sexp ) 
{

   unsigned long return_code;
   unsigned long i;

   unsigned long endemic_mort;
   unsigned long sdi_mort;
   unsigned long rand_error;
   unsigned long rand_seed;
   unsigned long use_precip_in_hg = 0;
   unsigned long use_genetic_gains = 0;

   unsigned long n_plots;
   struct PLOT_RECORD *plots_ptr = NULL;

   unsigned long n_plants;
   struct PLANT_RECORD *plants_ptr = NULL;

   long nyrs = asInteger(n_years_sexp);
   unsigned long hcb_growth_on = 1;
   double x0;
   unsigned long age = 0;

   SEXP ret_val;

   /* now sort the species back to the "native" order (by index) */
   qsort(  (void*)SPECIES_PTR, 
	   (size_t)(N_SPECIES), 
	   sizeof( struct SPECIES_RECORD ),
	   compare_species_by_idx );


   /* intitialize the config variables, grab them from the ctl list */
   rand_error  = asInteger( get_list_element( ctl_sexp, "rand.err" ) );
   rand_seed  = asInteger( get_list_element( ctl_sexp, "rand.seed" ) );
   endemic_mort = asInteger( get_list_element( ctl_sexp, "endemic.mort" ) );
   sdi_mort  = asInteger( get_list_element( ctl_sexp, "sdi.mort" ) ); 
   use_genetic_gains  = asInteger( get_list_element( ctl_sexp, "genetic.gains" ) );

   
   /* all bet this damn variable is causing the troubles!!!! */
   /* this is now going into the data list object */
   x0 = asReal( get_list_element( data_sexp, "x0" ) );
   age = asInteger( get_list_element( data_sexp, "age" ) );

   if( x0 < 0.0 )
   {
      x0 = 0.0;
   }

/*    Rprintf( "value of x0 = %lf\n", x0 ); */
/*    Rprintf( "value of age = %ld\n", age ); */

/*    Rprintf( "build_plot_array_from_sexp..." ); */
      plots_ptr = build_plot_array_from_sexp( 
	 get_list_element( data_sexp, "plots" ), &n_plots );
/*       Rprintf( "n_plots = %ld\n", n_plots ); */
/*    Rprintf( "done\n" ); */
   
/*    Rprintf( "build_plant_array_from_sexp..." ); */
      plants_ptr = build_plant_array_from_sexp( 
      get_list_element( data_sexp, "plants" ), &n_plants );
/*    Rprintf( "n_plants = %ld\n", n_plants ); */
/*    Rprintf( "done\n" ); */
   
      /* a check to ensure the site index values for the plots are non-zero */
      if( current_variant == CONIFERS_SMC )
      {
	 for( i = 0; i < n_plots; i++ )
	 {
	    if( plots_ptr[i].site_30 <= 0.0 )
	    {

	       Rprintf( "Invalid plots for variant %ld. Check site index values\n", current_variant );
	       ret_val = build_return_data_sexp( x0,
						 age,
						 n_plots, plots_ptr, 
						 n_plants, plants_ptr  );  
	       
	       free( plots_ptr );
	       free( plants_ptr );
	       
	       /* unprotect the return value (list) */
	       UNPROTECT( 1 );
	       return ret_val;	   
	       
	    }
	 }
      }
      

   /* project the sample.data for 1 year, nyrs times */
   for( i = 0; i < nyrs; i++ )
   {
        if(hcb_growth_on)
        {
	      hcb_growth_on = FALSE;  /*  turn off hcb growth  */
        }
        else
        {
	      hcb_growth_on = TRUE;   /*  turn on hcb growth  */
        }

	return_code = CONIFERS_SUCCESS;
	

/*  	Rprintf( "age = %d, value of x0 = %lf, before\n", age, x0 ); */
	project_plant_list( &return_code,

			    N_SPECIES,
			    SPECIES_PTR,
			    N_COEFFS,
			    COEFFS_PTR,

			    n_plants,
			    plants_ptr,
			    n_plots,
			    plots_ptr,

			    &x0,

			    endemic_mort,
			    sdi_mort,
			    hcb_growth_on,
			    use_precip_in_hg,
			    rand_error,
			    current_variant,
			    use_genetic_gains);

/* 	Rprintf( "age = %d, value of x0 = %lf, after\n", age, x0 ); */

	/* if the project didn't work, you must print an error message and return the unprojected data */
	if( return_code != CONIFERS_SUCCESS )
	{
	   Rprintf( "unable to project, return_code = %ld, %lf, check conifers.h for list of return codes\n", return_code, x0 );
	   ret_val = build_return_data_sexp( x0,
					     age,
					     n_plots, plots_ptr, 
					     n_plants, plants_ptr  );  
	  
	   free( plots_ptr );
	   free( plants_ptr );
	   
	   /* unprotect the return value (list) */
	   UNPROTECT( 1 );
	   return ret_val;	   
	}
	

	age++;
   }


/*   Rprintf( "done\n" ); */

/*   Rprintf( "building return data sexp..." ); */
  ret_val = build_return_data_sexp( x0,
				    age,
				    n_plots, plots_ptr, 
				    n_plants, plants_ptr  );  
/*   Rprintf( "done\n" ); */
  
  free( plots_ptr );
  free( plants_ptr );

  /* unprotect the return value (list) */
  UNPROTECT( 1 );
  return ret_val;

}



/* this function seems to be working as intended */
SEXP build_return_data_sexp( double x0,
			     unsigned long age,
			     unsigned long n_plots,
			     struct PLOT_RECORD *plots_ptr,
			     unsigned long n_plants,
			     struct PLANT_RECORD *plants_ptr )
{

   /* r interface variables */
   SEXP ret_val;
   SEXP sexp_plots;
   SEXP sexp_plants;

   PROTECT( ret_val = allocVector( VECSXP, 4 ) );

/*    Rprintf( "value of x0 = %lf\n", x0 ); */
   SET_VECTOR_ELT( ret_val, 0, ScalarReal( x0 ) );
   SET_VECTOR_ELT( ret_val, 1, ScalarInteger( age ) );

   PROTECT( sexp_plots =  build_sexp_from_plot_array( n_plots, plots_ptr ) );
   SET_VECTOR_ELT( ret_val, 2, sexp_plots );

   PROTECT( sexp_plants =  build_sexp_from_plant_array( n_plants, plants_ptr ) );
   SET_VECTOR_ELT( ret_val, 3, sexp_plants );

   UNPROTECT( 4 );
   return ret_val;
   
}

struct PLOT_RECORD *build_plot_array_from_sexp( SEXP plot_sexp, 
						unsigned long *n_plots )
{

   unsigned long i;
   struct PLOT_RECORD* plots_ptr;
   
   /* plots s expression variables */
   SEXP plot_plot_sexp;
   SEXP plot_elev_sexp;
   SEXP plot_slp_sexp; 
   SEXP plot_asp_sexp; 
   SEXP plot_h20_sexp; 
   SEXP plot_map_sexp;
   SEXP plot_si30_sexp; //site index

   PROTECT( plot_sexp = AS_LIST( plot_sexp ) );

   plot_plot_sexp = get_list_element( plot_sexp, "plot" );
   plot_elev_sexp = get_list_element( plot_sexp, "elevation" );
   plot_slp_sexp  = get_list_element( plot_sexp, "slope" );
   plot_asp_sexp  = get_list_element( plot_sexp, "aspect" );
   plot_h20_sexp  = get_list_element( plot_sexp, "whc" );
   plot_map_sexp  = get_list_element( plot_sexp, "map" );
   plot_si30_sexp = get_list_element( plot_sexp, "si30"); //site index

   PROTECT( plot_plot_sexp = coerceVector( plot_plot_sexp, INTSXP ) );
   PROTECT( plot_elev_sexp = coerceVector( plot_elev_sexp, REALSXP ) );
   PROTECT( plot_slp_sexp  = coerceVector( plot_slp_sexp, REALSXP ) );
   PROTECT( plot_asp_sexp  = coerceVector( plot_asp_sexp, REALSXP ) );
   PROTECT( plot_h20_sexp  = coerceVector( plot_h20_sexp, REALSXP ) );
   PROTECT( plot_map_sexp  = coerceVector( plot_map_sexp, REALSXP ) );
   PROTECT( plot_si30_sexp = coerceVector( plot_si30_sexp, REALSXP ) ); //site index

   /* build the plots vector */
   *n_plots = length( plot_plot_sexp );
   plots_ptr = (struct PLOT_RECORD*)calloc( 
      (*n_plots), sizeof( struct PLOT_RECORD ) );

/*    Rprintf( "n_plots %d\n", (*n_plots) ); */

   /* assign the plot array */
   for( i = 0; i < (*n_plots); i++ )
   {
      plots_ptr[i].plot = INTEGER( plot_plot_sexp )[i];
      plots_ptr[i].elevation = REAL( plot_elev_sexp )[i];
      plots_ptr[i].slope = REAL( plot_slp_sexp )[i];
      plots_ptr[i].aspect = REAL( plot_asp_sexp )[i];
      plots_ptr[i].water_capacity = REAL( plot_h20_sexp )[i];
      plots_ptr[i].mean_annual_precip = REAL( plot_map_sexp )[i];
      plots_ptr[i].site_30 = REAL( plot_si30_sexp )[i];  //site index

      /* should this be a check on the minimim site index value */
      /* only applies for CONIFERS_SMC */
      if( current_variant == CONIFERS_SMC )
      {
	 if( ISNA( REAL( plot_si30_sexp )[i] ) ||
	     ISNAN( REAL( plot_si30_sexp )[i] )  ||
	     plots_ptr[i].site_30 <= 0.0 )
	 {
	    Rprintf(
	       "Invalid si30 for plot %ld, setting value to %lf\n",
	       plots_ptr[i].plot,
	       plots_ptr[i].site_30 );
	    
	    plots_ptr[i].site_30 = 0.0;
	 }	 
      }

/*       Rprintf( */
/* 	 "%ld %lf %lf %lf %lf %lf\n", */
/* 	 plots_ptr[i].plot, */
/* 	 plots_ptr[i].elevation, */
/* 	 plots_ptr[i].slope, */
/* 	 plots_ptr[i].aspect, */
/* 	 plots_ptr[i].water_capacity, */
/* 	 plots_ptr[i].mean_annual_precip ); */

   }

   /* this needs to match the number of PROTECT statements in the function */
   UNPROTECT( 8 );
   
   return plots_ptr;
}

SEXP build_sexp_from_plot_array( unsigned long n_plots, 
				 struct PLOT_RECORD *plots_ptr )
{

   unsigned long i;
   SEXP ret_val;

   /* plot variables */
   SEXP ret_plots_id;
   SEXP ret_plots_elev;
   SEXP ret_plots_slope;
   SEXP ret_plots_asp;
   SEXP ret_plots_h20;
   SEXP ret_plots_map;
   SEXP ret_plots_si30;  //site index

   PROTECT( ret_val = allocVector( VECSXP, 7 ) );

   /* plots */
   PROTECT( ret_plots_id    = allocVector( INTSXP, n_plots ) );
   PROTECT( ret_plots_elev  = allocVector( REALSXP, n_plots ) );
   PROTECT( ret_plots_slope = allocVector( REALSXP, n_plots ) );
   PROTECT( ret_plots_asp   = allocVector( REALSXP, n_plots ) );
   PROTECT( ret_plots_h20   = allocVector( REALSXP, n_plots ) );
   PROTECT( ret_plots_map   = allocVector( REALSXP, n_plots ) );
   PROTECT( ret_plots_si30  = allocVector( REALSXP, n_plots ) );  //site index

   for( i = 0; i < n_plots; i++ )
   {
      INTEGER(ret_plots_id)[i] = plots_ptr[i].plot;
      REAL(ret_plots_elev)[i] = plots_ptr[i].elevation;
      REAL(ret_plots_slope)[i] = plots_ptr[i].slope;
      REAL(ret_plots_asp)[i] = plots_ptr[i].aspect;
      REAL(ret_plots_h20)[i] = plots_ptr[i].water_capacity;
      REAL(ret_plots_map)[i] = plots_ptr[i].mean_annual_precip;
      REAL(ret_plots_si30)[i] = plots_ptr[i].site_30; //site index
   }

   SET_VECTOR_ELT( ret_val, 0, ret_plots_id );
   SET_VECTOR_ELT( ret_val, 1, ret_plots_elev );
   SET_VECTOR_ELT( ret_val, 2, ret_plots_slope );
   SET_VECTOR_ELT( ret_val, 3, ret_plots_asp );
   SET_VECTOR_ELT( ret_val, 4, ret_plots_h20 );
   SET_VECTOR_ELT( ret_val, 5, ret_plots_map );
   SET_VECTOR_ELT( ret_val, 6, ret_plots_si30  );

   /* unprotect the vectros from the plots */
   UNPROTECT( 7 );

   return ret_val;
}


/********************************************************************************/
/* new function for converting plants into array of plants for the model	*/
/********************************************************************************/


struct PLANT_RECORD *build_plant_array_from_sexp( SEXP plant_sexp, 
						  unsigned long *n_plants )
{

   unsigned long i;
   struct PLANT_RECORD* plants_ptr;
   
   char                    temp_sp_code[16];
   struct SPECIES_RECORD *sp_ptr;

   /* plants s expression variables */
   SEXP plant_plot_sexp;
   SEXP plant_sp_code_sexp;
   SEXP plant_d6_sexp;
   SEXP plant_dbh_sexp;
   SEXP plant_tht_sexp;
   SEXP plant_cr_sexp;
   SEXP plant_n_stems_sexp;
   SEXP plant_expf_sexp;
   SEXP plant_crown_width_sexp;

   PROTECT( plant_sexp = AS_LIST( plant_sexp ) );

   plant_plot_sexp = get_list_element( plant_sexp, "plot" );
   plant_sp_code_sexp = get_list_element( plant_sexp, "sp.code" );
   plant_d6_sexp = get_list_element( plant_sexp, "d6" );
   plant_dbh_sexp = get_list_element( plant_sexp, "dbh" );
   plant_tht_sexp = get_list_element( plant_sexp, "tht" );
   plant_cr_sexp = get_list_element( plant_sexp, "cr" );
   plant_n_stems_sexp = get_list_element( plant_sexp, "n.stems" );
   plant_expf_sexp = get_list_element( plant_sexp, "expf" );
   plant_crown_width_sexp = get_list_element( plant_sexp, "crown.width" );

   /* read the plants */
   PROTECT( plant_plot_sexp = coerceVector( plant_plot_sexp, INTSXP ) );
   PROTECT( plant_sp_code_sexp = coerceVector( plant_sp_code_sexp, STRSXP ) );
   PROTECT( plant_d6_sexp = coerceVector( plant_d6_sexp, REALSXP ) );
   PROTECT( plant_dbh_sexp = coerceVector( plant_dbh_sexp, REALSXP ) );
   PROTECT( plant_tht_sexp = coerceVector( plant_tht_sexp, REALSXP ) );
   PROTECT( plant_cr_sexp = coerceVector( plant_cr_sexp, REALSXP ) );
   PROTECT( plant_n_stems_sexp = coerceVector( plant_n_stems_sexp, INTSXP ) );
   PROTECT( plant_expf_sexp = coerceVector( plant_expf_sexp, REALSXP ) );
   PROTECT( plant_crown_width_sexp = coerceVector( plant_crown_width_sexp, REALSXP ) );

   /* build the plots vector */
   *n_plants = length( plant_plot_sexp );
   plants_ptr = (struct PLANT_RECORD*)calloc( 
      (*n_plants), sizeof( struct PLANT_RECORD ) );

/*    Rprintf( "n_plants %d\n", (*n_plants) ); */

    /* sort the species codes based on sp_code */
    qsort(  (void*)SPECIES_PTR,
            (size_t)(N_SPECIES),
            sizeof( struct SPECIES_RECORD ),
	    compare_species_by_sp_code );
    
   /* assign the plot array */
   for( i = 0; i < (*n_plants); i++ )
   {
      plants_ptr[i].plot = INTEGER( plant_plot_sexp )[i];

/*       plants_ptr[i].plant = INTEGER( plant_plant_sexp )[i]; */
      plants_ptr[i].plant = i+1;

      strcpy( temp_sp_code, CHAR( STRING_ELT( plant_sp_code_sexp, i ) ) );

      /* get the species code and look up the correct index */
      sp_ptr = get_species_entry_from_code(    N_SPECIES,
					       SPECIES_PTR,
					       temp_sp_code );
      if( !sp_ptr )
      {
	 Rprintf( "Couldn't find the species code for %s, %s in species map\n",
		  temp_sp_code, CHAR( STRING_ELT( plant_sp_code_sexp, i ) ) );
	 Rprintf( "Make sure you have the entry in your species map. See help\n" );
	 continue;
      }

/*       Rprintf( "initialize species code %ld %s %s %ld %s\n", */
/* 	       i, */
/* 	       temp_sp_code, */
/* 	       CHAR( STRING_ELT( plant_sp_code_sexp, i ) ), */
/* 	       sp_ptr->idx, */
/* 	       sp_ptr->sp_code ); */
      
      /* this is the index of the "unsorted" array */
      plants_ptr[i].sp_idx = sp_ptr->idx;
      plants_ptr[i].d6 = REAL( plant_d6_sexp )[i];
      plants_ptr[i].dbh = REAL( plant_dbh_sexp )[i];
      plants_ptr[i].tht = REAL( plant_tht_sexp )[i];
      plants_ptr[i].cr = REAL( plant_cr_sexp )[i];
      plants_ptr[i].n_stems = INTEGER( plant_n_stems_sexp )[i];
      plants_ptr[i].expf = REAL( plant_expf_sexp )[i];
      plants_ptr[i].crown_width = REAL( plant_crown_width_sexp )[i];

      /* these are calculated values */
      plants_ptr[i].d6_area = plants_ptr[i].d6*plants_ptr[i].d6*FC_I;
      plants_ptr[i].basal_area = plants_ptr[i].dbh*plants_ptr[i].dbh*FC_I;
      plants_ptr[i].crown_area = plants_ptr[i].crown_width * 
	  plants_ptr[i].crown_width * MY_PI / 4.0;

      /* perform some basic error checking here */
      /* see if you can use the ISNAN macro here */

      if( ISNA( REAL( plant_d6_sexp )[i] ) ||
	  ISNAN( REAL( plant_d6_sexp )[i] )  ||
	  plants_ptr[i].d6 < 0.0 )
      {
	    plants_ptr[i].d6 = 0.0;
      }

      if( ISNA( REAL( plant_dbh_sexp )[i] ) ||
	  ISNAN( REAL( plant_dbh_sexp )[i] )  ||
	  plants_ptr[i].dbh < 0.0 )
      {
	    plants_ptr[i].dbh = 0.0;
      }

      if( ISNAN( REAL( plant_tht_sexp )[i] )  || plants_ptr[i].expf < 0.0 )
      {
	    plants_ptr[i].tht = 0.0;
      }

      if( ISNAN( REAL( plant_cr_sexp )[i] )  || plants_ptr[i].cr < 0.0 )
      {
	    plants_ptr[i].cr = 0.0;
      }

      if( ISNAN( REAL( plant_expf_sexp )[i] )  || plants_ptr[i].expf < 0.0 )
      {
	    plants_ptr[i].expf = 0.0;
      }

      if( ISNAN( REAL( plant_crown_width_sexp )[i] )  || plants_ptr[i].crown_width < 0.0 )
      {
	    plants_ptr[i].crown_width = 0.0;
	    plants_ptr[i].crown_area = 0.0;
      }

/*       if( ISNAN( REAL( plant_crown_area_sexp )[i] )  || plants_ptr[i].crown_area < 0.0 ) */
/*       { */
/* 	 plants_ptr[i].crown_area = 0.0; */
/*       } */

/*       plants_ptr[i].crown_area = REAL( plant_crown_area_sexp )[i]; */

   }

   /* now sort the species back to the "native" order (by index) */
   qsort(  (void*)SPECIES_PTR,
	   (size_t)(N_SPECIES),
	   sizeof( struct SPECIES_RECORD ),
	   compare_species_by_idx );

   UNPROTECT( 10 );   /* plot lists */
   
   return plants_ptr;
}


/* this function includes the optional "errors" associated with each tree */
SEXP build_sexp_from_plant_array( 
   unsigned long n_plants, 
   struct PLANT_RECORD *plants_ptr )
{

   unsigned long i;
   SEXP ret_val;
   
   /* plot variables */
   SEXP ret_plants_plot;
   SEXP ret_plants_sp_code;
   SEXP ret_plants_d6;
   SEXP ret_plants_dbh;
   SEXP ret_plants_tht;
   SEXP ret_plants_cr;
   SEXP ret_plants_n_stems;
   SEXP ret_plants_expf;
   SEXP ret_plants_crown_width;
   SEXP ret_plants_errors;

   PROTECT( ret_val = allocVector( VECSXP, 10 ) );

   /* plants */
   PROTECT( ret_plants_plot = allocVector( INTSXP, n_plants ) );
   PROTECT( ret_plants_sp_code = allocVector(STRSXP,n_plants));
   PROTECT( ret_plants_d6 = allocVector( REALSXP, n_plants ) );
   PROTECT( ret_plants_dbh = allocVector( REALSXP, n_plants ) );
   PROTECT( ret_plants_tht = allocVector( REALSXP, n_plants ) );
   PROTECT( ret_plants_cr = allocVector( REALSXP, n_plants ) );
   PROTECT( ret_plants_n_stems = allocVector( INTSXP, n_plants ) );
   PROTECT( ret_plants_expf = allocVector( REALSXP, n_plants ) );
   PROTECT( ret_plants_crown_width = allocVector( REALSXP, n_plants ) );
   PROTECT( ret_plants_errors = allocVector( INTSXP, n_plants ) );

   /* now sort the species back to the "native" order (by index) */
   qsort(  (void*)SPECIES_PTR,
	   (size_t)(N_SPECIES),
	   sizeof( struct SPECIES_RECORD ),
	   compare_species_by_idx );

   for( i = 0; i < n_plants; i++ )
   {
      INTEGER(ret_plants_plot)[i] = plants_ptr[i].plot;
      SET_STRING_ELT(ret_plants_sp_code, i, mkChar( SPECIES_PTR[plants_ptr[i].sp_idx].sp_code ) );

      REAL(ret_plants_d6)[i] = plants_ptr[i].d6;
      REAL(ret_plants_dbh)[i] = plants_ptr[i].dbh;
      REAL(ret_plants_tht)[i] = plants_ptr[i].tht;
      REAL(ret_plants_cr)[i] = plants_ptr[i].cr;
      INTEGER(ret_plants_n_stems)[i] = plants_ptr[i].n_stems;
      REAL(ret_plants_expf)[i] = plants_ptr[i].expf;
      REAL(ret_plants_crown_width)[i] = plants_ptr[i].crown_width;

      INTEGER(ret_plants_errors)[i] = plants_ptr[i].errors;

   }

   SET_VECTOR_ELT( ret_val, 0, ret_plants_plot );
   SET_VECTOR_ELT( ret_val, 1, ret_plants_sp_code );
   SET_VECTOR_ELT( ret_val, 2, ret_plants_d6 );
   SET_VECTOR_ELT( ret_val, 3, ret_plants_dbh );
   SET_VECTOR_ELT( ret_val, 4, ret_plants_tht );
   SET_VECTOR_ELT( ret_val, 5, ret_plants_cr );
   SET_VECTOR_ELT( ret_val, 6, ret_plants_n_stems );
   SET_VECTOR_ELT( ret_val, 7, ret_plants_expf );
   SET_VECTOR_ELT( ret_val, 8, ret_plants_crown_width );
   SET_VECTOR_ELT( ret_val, 9, ret_plants_errors );

   /* unprotect the vectros from the plots */
   UNPROTECT( 10 );

   return ret_val;
}


SEXP r_thin_sample( 
   SEXP data_sexp,
   SEXP ctl_sexp ) 
{

   unsigned long p;
   unsigned long return_code;

   unsigned long n_plots;
   struct PLOT_RECORD *plots_ptr;
   struct PLOT_RECORD *plot_ptr;

   double x0;
   unsigned long age;
   unsigned long n_plants;
   struct PLANT_RECORD *plants_ptr;

   char                    temp_sp_code[SP_LENGTH];
   struct SPECIES_RECORD *sp_ptr;
   unsigned long		sp_idx;
   
   /* the thinning functions set these variables */
   double	n_plants_removed;
   double	ba_removed;
   long		thin_type;
   double	target;
/*    unsigned long n_ctl_args = length( ctl_sexp ); */

   SEXP ret_val;

   /* do you need to protect the data and control list going into the function? */


   /* intitialize the config variables, grab them from the ctl list */
   thin_type = asInteger( get_list_element( ctl_sexp, "type" ) );
   target  = asReal( get_list_element( ctl_sexp, "target" ) ); 

//   target  = asReal( get_list_element( ctl_sexp, "sp" ) ); 

   x0 = asReal( get_list_element( data_sexp, "x0" ) );
   age = asInteger( get_list_element( data_sexp, "age" ) );
   plots_ptr = build_plot_array_from_sexp( 
      get_list_element( data_sexp, "plots" ), &n_plots );

   plants_ptr = build_plant_array_from_sexp( 
      get_list_element( data_sexp, "plants" ), &n_plants );


   // if the species is not null, then the user wants to 
   // thin, specifically, a single species. otherwise, 
   // if the sp == NULL, the user wants to thin all species
   // so the first task is to find the species index...
   sp_idx = 0;
   if( thin_type == DO_EXPF_SP_THIN || thin_type == DO_EXPF_SP_THIN_FROM_BELOW )
   {
      if( CHAR(STRING_ELT(get_list_element( ctl_sexp, "target.sp" ), 0)) != NULL )
	 //if( temp_sp_code != NULL )
      {
	 strcpy( temp_sp_code, CHAR(STRING_ELT(get_list_element( ctl_sexp, "target.sp" ), 0)) );
	 
	 /* sort the species codes based on sp_code */
	 qsort(  (void*)SPECIES_PTR,
		 (size_t)(N_SPECIES),
		 sizeof( struct SPECIES_RECORD ),
		 compare_species_by_sp_code );
	 
	 /* get the species code and look up the correct index */
	 sp_ptr = get_species_entry_from_code(    N_SPECIES,
						  SPECIES_PTR,
						  temp_sp_code );
	 if( sp_ptr == NULL )
	 {
	    Rprintf( "Couldn't find the species code for target.sp = %s in the current species map\n",
		     temp_sp_code,
		     CHAR(STRING_ELT(get_list_element( ctl_sexp, "target.sp" ), 0)) );
	    Rprintf( "Make sure you have the entry in your species map. See help\n" );
	    ret_val = build_return_data_sexp( x0,
					      age,
					      n_plots, plots_ptr, 
					      n_plants, plants_ptr  );  
	    
	    free( plots_ptr );
	    free( plants_ptr );
	    UNPROTECT( 1 );
	   
	 /* sort the species codes based on sp_code */
	 qsort(  (void*)SPECIES_PTR,
		 (size_t)(N_SPECIES),
		 sizeof( struct SPECIES_RECORD ),
		 compare_species_by_idx );

 	    return ret_val;
	 } 
	 
	 sp_idx = sp_ptr->idx;

      }
   }


   /* sort the species codes based on sp_code */
   qsort(  (void*)SPECIES_PTR,
	   (size_t)(N_SPECIES),
	   sizeof( struct SPECIES_RECORD ),
	   compare_species_by_idx );
   
   /* if the code is NULL, then send NULL to the function */
   /* and simply pass in the numeric values for the */
   /* thinning function */
   plot_ptr = &plots_ptr[0];
   for( p = 0; p < n_plots; p++, plot_ptr++ )
   {
//      Rprintf( "thinning plot %d\n", plot_ptr->plot );
      
      thin_plot(  &return_code,
		  n_plants,
		  plants_ptr,
		  plot_ptr,
		  
		  /* globals (at least for this interface) */
		  N_SPECIES,
		  SPECIES_PTR,
		  N_COEFFS,
		  COEFFS_PTR,
		  
		  sp_idx,
		  
		  thin_type,
		  target,

		  &n_plants_removed,
		  &ba_removed );

/*       Rprintf( "thinning plot %ld for %ld, by type %ld, to a target of %lf, n_plants_removed = %lf, ba_removed = %lf, return_code = %ld\n", */
/* 	       plot_ptr->plot, */
/* 	       sp_idx, */
/* 	       thin_type, */
/* 	       target, */
/* 	       n_plants_removed, ba_removed, */
/* 	 return_code ); */

   }

   /* do we ever update x0, say before-after a thinning?	*/
   /* if so you 1) need to expose/export the ability to		*/
   /* update the value and 2) update x0 here			*/

   //x0 = asReal( get_list_element( data_sexp, "x0" ) );
/*    Rprintf( "value of x0 = %lf\n", x0 ); */

   
   //ret_val =  build_sexp_from_plant_array( n_plants, plants_ptr );
   /* build the output sample */
   /* we don't really need to modify the plot attribs (yet) and so	*/
   /* we don't need to update the values. we can simple copy the	*/
   /* original into the return structure in here (or in the R code)	*/
/*   Rprintf( "building return data sexp..." ); */
  ret_val = build_return_data_sexp( x0, 
				    age,
				    n_plots, plots_ptr, 
				    n_plants, plants_ptr  );  
/*   Rprintf( "done\n" ); */

   free( plots_ptr );
   free( plants_ptr );

   UNPROTECT( 1 );
   return ret_val;

}

SEXP r_fill_in_missing_values( 
   SEXP data_sexp,
   SEXP ctl_sexp ) 
{

/*    unsigned long p; */
   unsigned long return_code;

/* data variables */
   double x0;
   unsigned long age;
   unsigned long n_plots;
   struct PLOT_RECORD *plots_ptr;
   unsigned long n_plants;
   struct PLANT_RECORD *plants_ptr;

/* control variables */
   double fixed_plot_radius;
   double min_dbh;
   double baf;
//   unsigned long model_variant;

   SEXP ret_val;   

/* set the control variables */
   fixed_plot_radius  = asReal( get_list_element( ctl_sexp, "fpr" ) ); 
   min_dbh  = asReal( get_list_element( ctl_sexp, "min.dbh" ) ); 
   baf  = asReal( get_list_element( ctl_sexp, "baf" ) ); 

//   model_variant = asInteger( get_list_element( ctl_sexp, "model.variant" ) ); 

/* set the data variables */
   x0 = asReal( get_list_element( data_sexp, "x0" ) );
   age = asInteger( get_list_element( data_sexp, "age" ) );

   plots_ptr = build_plot_array_from_sexp( 
      get_list_element( data_sexp, "plots" ), &n_plots );

   plants_ptr = build_plant_array_from_sexp( 
      get_list_element( data_sexp, "plants" ), &n_plants );
   

/*    /\* perform a general check for valid pointers *\/ */
/*    if( plots_ptr == NULL || plants_ptr == NULL ) */
/*    { */

/*       if( plots_ptr == NULL ) */
/*       { */
/* 	 Rprintf( "unable to impute values. plots may have missing values for variant\n" ); */
/*       } */
      
/*       //UNPROTECT( 1 ); */
/*       //return data_sexp; */
/*       return R_NilValue; */

/*    } */


   /* if you have everything, then fill in the missing values */
   fill_in_missing_values( &return_code,
			   
			   N_SPECIES,
			   SPECIES_PTR,
			   N_COEFFS,
			   COEFFS_PTR,
			   
			   //model_variant, 
			   current_variant,
			   
			   n_plants,
			   plants_ptr,
			   n_plots,
			   plots_ptr,
			   
			   fixed_plot_radius,
			   min_dbh,
			   baf );

   if( return_code != CONIFERS_SUCCESS )
   {
      Rprintf( "unable to impute values, return_code = %ld\n", return_code );
   }

  ret_val = build_return_data_sexp( x0, 
				    age,
				    n_plots, plots_ptr, 
				    n_plants, plants_ptr  );  
/*   Rprintf( "done\n" ); */

   free( plots_ptr );
   free( plants_ptr );

   UNPROTECT( 1 );
   return ret_val;


}


SEXP r_calc_max_sdi( 
   SEXP data_sexp )
{

   unsigned long return_code;
   unsigned long n_plots;
   unsigned long n_plants;
   struct PLANT_RECORD *plants_ptr;

   double	max_sdi = 0.0;

   SEXP ans;

   PROTECT(ans = allocVector(REALSXP, 1));

   /* build the plots vector */
   n_plots = length( get_list_element( data_sexp, "plots" ) );
   plants_ptr = build_plant_array_from_sexp( 
      get_list_element( data_sexp, "plants" ), &n_plants );

   calc_max_sdi( &return_code,
		 
		 N_SPECIES,
		 SPECIES_PTR,

		 N_COEFFS,
		 COEFFS_PTR,
		 
		 n_plants,
		 plants_ptr,

		 n_plots,
		 &max_sdi );
   
   if( return_code != CONIFERS_SUCCESS )
   {
      Rprintf( "unable to compute max_sdi, return_code = %ld\n", return_code );
   }
   
   free( plants_ptr );

   REAL(ans)[0] = max_sdi;
   UNPROTECT( 1 );
   return ans;


}

