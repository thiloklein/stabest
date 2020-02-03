//#define ARMA_USE_SUPERLU
//#define DEBUG
//#define DEBUG2
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <limits> // for infinity
#include <random>
#include <cassert> // for checks
#include "aux_functions.h"
#include "truncnorm_approx.h"

//#define INF arma::datum::inf
#define INF std::numeric_limits<double>::infinity()
//#define FMIN(a, b) ((a) < (b) ? (a) : (b)) // I use std::min instead!
//#define FMAX(a, b) ((a) > (b) ? (a) : (b))
#define TRUNCNORM(mu, sigma, lower, upper) truncn2(mu, sigma, lower, upper)
//#define TRUNCNORM(mu, sigma, lower, upper) truncnorm_approx(runiform(gen), mu, sigma, lower, upper)

// lookup utility to get the cid of an arbitrary row 'idx', based on a lookup table that contains the start indices of the next college block for each cid
int get_cid(int idx, Rcpp::IntegerVector cid_lookup) {
  int cid = 0;
  while (idx >= cid_lookup[cid]) cid++;
  return cid;
}

using namespace Rcpp;
//using namespace arma;

#ifndef DEBUG
	#define ARMA_NO_DEBUG
#endif



// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List stabest_internal(
    Rcpp::NumericVector Vc, // students' valuation over colleges
    Rcpp::NumericVector Vs, // colleges' valuations over students
    arma::sp_mat Xc, // variables that drive students' preferences over colleges 
    arma::sp_mat XcXcInv, // pre-computed in R because sparse matrix inversion requires extra library in Rcpp
    arma::sp_mat Xs, // variables that drive colleges' preferences over students
    arma::sp_mat XsXsInv, // pre-computed in R because sparse matrix inversion requires extra library in Rcpp
    Rcpp::NumericMatrix betaR, // vector of coefficients driving students' prefs over colleges, 1 col and student.nCoef rows
    Rcpp::NumericMatrix gammaR, // vector of coefs driving colleges' prefs over students, 1 col and college.nCoef rows
    Rcpp::LogicalVector VacantSeats, // VacantSeats[i]=1 if college i has vacant seats, vector of length nColleges
    Rcpp::IntegerVector nObs_college, // nObs_college[c] gives the number of observations that college c has
    Rcpp::IntegerVector ID_eqCollege, // ID_ of equilibrium college for each student (R index, cycles)
    Rcpp::LogicalVector match, // =1 if observation is an observed match
    Rcpp::IntegerVector sid, // student id of an observation (R indexing)
    Rcpp::IntegerVector ID_cBetter, // ID_ of a better college (R indexing), 0 if NA
    Rcpp::IntegerVector ID_cWorse, // ID_ of a worse college (R indexing), 0 if NA, -1 max of all unacceptable ones, -2 = used to indicate unacceptable
    Rcpp::IntegerVector ID_sBetter, // ID_ of a better student (R indexing), 0 if NA
    Rcpp::IntegerVector ID_sWorse, // ID_ of a worse student (R indexing), 0 if NA, -1 max of all unacceptable ones, -2 = used to indicate unacceptable
    Rcpp::IntegerVector ID_nextCollege, // ID_ of the next obs with the current student and another college (R indexing), it cycles
    int niter,
    int thin,
    int demean
) {
  
  
  //-------------------------------------------------------------------------------------------------//
  //--- set up --------------------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------------------------//
  
  
  
  
  // preliminaries and dimension checks
#ifdef DEBUG
	std::cout << "Preliminary checks ...\n" << std::flush;
#endif
  int nColleges = VacantSeats.size();
  assert( nObs_college.size() == nColleges );
  assert( VacantSeats.size() == nColleges );
  unsigned int nData = Vc.size();
  assert( Vs.size() == nData );
  assert( Xc.n_rows == nData );
  assert( Xs.n_rows == nData );
  assert( ID_sBetter.size() == nData );
  assert( ID_cBetter.size() == nData );
  assert( ID_sWorse.size() == nData );
  assert( ID_cWorse.size() == nData );
  assert( ID_nextCollege.size() == nData );
  assert( sid.size() == nData );
  assert( match.size() == nData );
  int nStudents = Rcpp::max(sid);
  assert( ID_eqCollege.size() == nStudents );
  assert( Rcpp::max(ID_eqCollege) <= (int) nData ); // ID_ is R-style index, so <= rather than <
  assert( XcXcInv.n_cols == XcXcInv.n_rows );
  assert( (int) XcXcInv.n_cols == betaR.nrow() );
  assert( betaR.ncol() == 1);
  assert( XsXsInv.n_cols == XsXsInv.n_rows );
  assert( (int) XsXsInv.n_cols == gammaR.nrow() );
  assert( gammaR.ncol() == 1);
  
  clock_t begin = 0, begin2 = 0, draw_time = 0, loop_time=0;
  
  
  // count number of obs per student (needed for demeaning)

  Rcpp::IntegerVector nObs_student(nStudents);
  for (int s=0; s<nStudents; s++) nObs_student[s] = 0;
  for (int i=0; i<sid.size(); i++) nObs_student[sid[i]-1] = std::max(1, nObs_student[sid[i]-1]+1); // since we divide by it, we don't want a zero here.
#ifdef DEBUG
  printf("We have %d students and %d colleges.\n", nStudents, nColleges);
#endif
  
  
  
  
  // create lookup table for college start indices
  
  Rcpp::IntegerVector idx_college_start(nColleges);
  int cumsum = 0;
  for (int c=0; c<nColleges; c++) {
    cumsum += nObs_college[c];
    idx_college_start[c] = cumsum;
  }
  assert( idx_college_start(nColleges-1) == nData ); // important for the lookup utility so that it returns cid < nColleges
  
  
  
  
  // set up projection matrices
 #ifdef DEBUG
	std::cout << "Set up projection matrices ...\n" << std::flush;
#endif
  arma::mat Proj_Xc( XcXcInv * trans(Xc) ); // NB: although design matrix may be sparse, projection matrix is probably dense
  arma::mat Proj_Xs( XsXsInv * trans(Xs) );
 #ifdef DEBUG
	std::cout << "Set up dense versions of inv(X'X) ...\n" << std::flush;
#endif
  arma::mat XcXcInv_dense(XcXcInv); // quick & dirty, need dense matrix for cholesky decomp
  arma::mat XsXsInv_dense(XsXsInv); // quick & dirty, need dense matrix for cholesky decomp
  
  
  
  
  // set up coefficient matrices
#ifdef DEBUG
	std::cout << "Set up coefficient matrices ...\n" << std::flush;
#endif
  arma::colvec beta = Rcpp::as<arma::colvec>( betaR );
  arma::colvec beta_hat = Rcpp::as<arma::colvec>( betaR );
  arma::colvec gamma = Rcpp::as<arma::colvec>( gammaR );
  arma::colvec gamma_hat = Rcpp::as<arma::colvec>( gammaR );
  
  
  
  
  // set up sample paths
#ifdef DEBUG
	std::cout << "Set up sample paths ...\n" << std::flush;
#endif
  int Nsamples = floor(niter / thin);
  arma::mat betavalues = arma::zeros(Nsamples, Xc.n_cols);
  arma::mat gammavalues = arma::zeros(Nsamples, Xs.n_cols);
  
  
  
  
  // set up vectors for latent and mean utility
#ifdef DEBUG
	std::cout << "Set up vectors for latent and mean utility ...\n" << std::flush;
#endif
  arma::colvec Xc_beta = arma::zeros( Vs.size() );
  arma::colvec Xs_gamma = arma::zeros( Vc.size() );
  arma::colvec Vc_GroupMean = arma::zeros( nStudents );
  arma::colvec Vs_GroupMean = arma::zeros( nColleges );
  Rcpp::NumericVector Vs_min(nColleges);
  Rcpp::NumericVector Vs_maxUnacceptable(nColleges); // max of unacceptable valuations
  Rcpp::NumericVector Vc_maxUnacceptable(nStudents);
  
  // initialize utility vectors
#ifdef DEBUG
	std::cout << "Initialize utility vectors ...\n" << std::flush;
#endif
  Xc_beta = Rcpp::as<arma::colvec>( Vc );
  Xs_gamma = Rcpp::as<arma::colvec>( Vs );
  for (int c=0; c<nColleges; c++) Vs_min[c] = INF; // start with +INF to get the minimum across assigned students
  for (int c=0; c<nColleges; c++) Vs_maxUnacceptable[c] = - INF;
  for (int s=0; s<nStudents; s++) Vc_maxUnacceptable[s] = - INF;
  int idx_cs=0;
  int s = 0;
  for(int c=0; c <nColleges; c++){
    for(int j=0; j < nObs_college[c]; j++){
      s = sid[idx_cs] - 1; // data could be irregular so that we need a lookup for the student here, converted to C indexing
      if (match[idx_cs]) Vs_min[c] = std::min(Vs_min[c], Vs[idx_cs]); // update the minimum of the college's valuation over all matched students
      if (ID_cWorse[idx_cs]==-2) Vc_maxUnacceptable[s] = std::max(Vc_maxUnacceptable[s], Vc[idx_cs]); // if c is unacceptable to s, update max of unacceptable valuations
//#ifdef DEBUG
      //printf("c=%d, idx=%d, ID_cWorse=%d, Vs_maxUnacceptable=%f, Xs_gamma=%f, ", c, idx_cs, ID_sWorse[idx_cs], Vs_maxUnacceptable[c], Vs[idx_cs]);
      //printf("Vs_maxUnacceptable[c] > Vc[idx_cs] yields %d\n", Vs_maxUnacceptable[c] > Vs[idx_cs] );
      //printf("FMAX() = %f, std::max() = %f\n", FMAX(Vs_maxUnacceptable[c], Xs_gamma[idx_cs]), std::max(Vs_maxUnacceptable[c], Xs_gamma[idx_cs]));
//#endif
      if (ID_sWorse[idx_cs]==-2) Vs_maxUnacceptable[c] = std::max(Vs_maxUnacceptable[c], Vs[idx_cs]); // if s is unacceptable to c, update max of unacceptable valuations over students
      idx_cs += 1;
    }
  }
  for (int c=0; c<nColleges; c++) if (std::isinf(  Vs_min[c] ) ) Vs_min[c] = -INF; // now set to -INF if the college has no students (since Vs_min serves as a lower bound)
  
// #ifdef DEBUG
  // printf("\nInitial state of utility vectors:\n");
  // std::cout << "Xc_beta: " << Xc_beta << std::endl;
  // std::cout << "Xs_gamma: " << Xs_gamma << std::endl;
  // for (int c=0; c<nColleges; c++) printf("c=%d, Vs_min=%f, Vs_max_unacceptable=%f\n", c, Vs_min[c], Vs_maxUnacceptable[c]);
  // for (int s=0; s<nStudents; s++) printf("s=%d, Vc_max_unacceptable=%f\n", s, Vc_maxUnacceptable[s]);
// #endif
  
  
  
  
  // make custom progress report, every 10% of iterations
  
  Rcpp::Rcout << "Process " << getpid() << ": Drawing " << niter << " MCMC samples..." << std::endl;
  int chunk10percent = floor( (double) niter / (double) 10);
  
#ifdef DEBUG
	std::cout << "Start main iteration loop ...\n" << std::flush;
#endif
  
  
  //-------------------------------------------------------------------------------------------------//
  //--- start iteration loop ------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------------------------//
  
  for (int iter=0; iter<niter; iter++) {
    
    
    // progress report
    begin = clock();
    if (!fmod(iter, chunk10percent) ) {
      Rcpp::Rcout << "Iteration " << iter << " ..." ;
    }
    
    
    // reset / init temporary variables
    
    double Vcupperbar=INF, Vsupperbar=INF, Vclowerbar=-INF, Vslowerbar=-INF, 
      Vs_tmp=0, Vc_tmp=0, Vs_cs=0, Vc_cs=0, Vs_min_c=0;
    bool match_cs = false, HasVacantSeats_c=false;
    idx_cs = 0; // current position of college c and student s
    int idx_c_start = 0; // start position of college c
    int idx_eqCollege_s = 0, idx_eqCollege_sprime=0; // position of the equilibrium college of student s
    int idx_cprimes = 0; // position of college cprime and studen s
    int idx_tmp = 0; // temporary index
    int nObs_c = 0; // number of observations for this college
	  int cprime = 0; // cid of alternative college
    s = 0; // student ID
    arma::uword idx_cs_uword = 0; // index, converted to unsigned (for arma vectors)
    if (demean) {
      Vc_GroupMean *= 0;
      Vs_GroupMean *= 0;
    }
    
    
    for(int c=0; c <nColleges; c++){
      
      idx_c_start = idx_cs;
      HasVacantSeats_c = VacantSeats[c];
      nObs_c = nObs_college[c];
      
      for(int j=0; j < nObs_c; j++){
        
        
        // store variables to reduce memory access
        
        s = sid[idx_cs] - 1; // data could be irregular so that we need a lookup for the student here, converted to C indexing
        match_cs = match[idx_cs];
        idx_eqCollege_s = ID_eqCollege[s] - 1; // row index of student s's eq. college (ID_eqCollege[s]==0 indicates that student is not matched)
        Vs_cs = Vs[idx_cs];
        Vs_min_c = Vs_min[c];
        
        
        // reset bounds 
        
        Vcupperbar=INF, Vsupperbar=INF, Vclowerbar=-INF, Vslowerbar=-INF;

        std::cout.flush();
#ifdef DEBUG
        printf("\n**************************************************************************\n");
        printf("iter=%d, idx=%d, c=%d, s=%d, match_cs=%d, idx_eqCollege_s=%d\n", iter, idx_cs, c, s, match_cs, idx_eqCollege_s);
#endif
        
        
        
        //----------------------------------------------------------------------------------//
        // --- student s's valuation over college c ----------------------------------------//
        //----------------------------------------------------------------------------------//
        
        
        
        //--- equilibrium bounds -----------------------------------------------------------//
        
#ifdef DEBUG
          std::cout << "Find student's eq. bounds ..\n" << std::flush;
#endif
        if (idx_eqCollege_s>=0) { // if student s is matched to _some_ school
		
    		  if(match_cs == 0){  // case 1: non-equilibrium matches
    			  
				    // if s is not unacceptable to school c AND
    			  //   - student's valuation is higher than minimum of college's students (so college prefers student) [Eqn (A4)] OR
				    //   - c has a vacant seat ...
    			  if( ID_sWorse[idx_cs]!=-2 && (( Vs_cs > Vs_min_c ) || HasVacantSeats_c ) ){
    				  // ... then college's valuation has to be lower than that of student's equilibrium college
    				  Vcupperbar = std::min( Vcupperbar, Vc[idx_eqCollege_s] ); 
    			  }
    			  
    			} else {  // case 2: equilibrium matches
				    
    			  // check all other colleges: if student s could get in, set Vc so that he does not want to go there
    			  idx_cprimes = idx_cs;
					  cprime = c;
					  
					  // get record number of next college until we reach idx_cs again, converting R style ID_ index to C style 0 based index along the way
					  // we don't need to check college c again because student can only be matched to one college, so only nColleges - 1 remain
    			  while( (idx_cprimes = ID_nextCollege[idx_cprimes] - 1) != idx_cs ){ 
    			    
    				  cprime = get_cid(idx_cprimes, idx_college_start); // find college ID of next college, using the lookup utility
					    
						  // if s is not unacceptable to cprime AND
    					//   - student's valuation is higher than minimum of cprime's students (so college prefers student) [Eqn (A6)] OR
						  //   - cprime has a vacant seat ...
    					if( ID_sWorse[idx_cprimes]!=-2 && ( Vs[idx_cprimes] > Vs_min[cprime] || VacantSeats[cprime] ) ){
    					  
    					  // ... then student must value college i higher than iprime [Eqn (A11) = Eqn (A15_b), term 1]
    						Vclowerbar = std::max( Vclowerbar, Vc[idx_cprimes] ); 
    					  
    					} // end if
				    } // end loop over schools
    			} // non-equilibrium vs. equilibrium 
        } // end if
        
        
        
        //--- rank-order-list based bounds ------------------------------------------------------//
        
        
#ifdef DEBUG
          std::cout << "Determine student's rank order bounds ..\n" << std::flush;
#endif
        // recall that ID_cWorse and ID_cBetter uses R indexing, 0 indicates NA (no rank order bound)
        
        if( (idx_tmp = ID_cBetter[idx_cs])>0 ) // if a better college exists
          Vcupperbar = std::min( Vcupperbar, Vc[ idx_tmp - 1 ] ); // student j's valuation over the college that j ranks just above college i.
        if( (idx_tmp = ID_cWorse[idx_cs]) > 0 ) // if a worse ranked alternative college exists
          Vclowerbar = std::max( Vclowerbar, Vc[ idx_tmp - 1  ]); // student s's valuation over the college that s ranks just below college c.
        else if (idx_tmp==-1) // -1 indicates that the lower bound is the max of all unacceptable colleges (could by -INF if no unacceptable ones exist)
          Vclowerbar = std::max( Vclowerbar, Vc_maxUnacceptable[ s ]);
        
        
        
        //--- draw new valuations that respect the bounds --------------------------------------//
        
#ifdef DEBUG
          std::cout << "Drawing student's new latent valuation over college ..\n" << std::flush;
#endif
        if(Vclowerbar < Vcupperbar)
          Vc_cs = TRUNCNORM(Xc_beta[ idx_cs_uword ], 1.0, Vclowerbar, Vcupperbar); // mu, sigma, lower, upper
        else // Vcupperbar<Vclowerbar => Vcupper < +Infinity so we can use it as a new valuation
          Vc_cs = Vcupperbar;
        Vc[idx_cs] = Vc_cs;
        if (demean) Vc_GroupMean[s] += Vc_cs;
        
        
        // update the max of unacceptable valuations over colleges if c is unacceptable to s
        
        if (idx_tmp==-2) { // if college c is unacceptable to student s, need to update
          if ( Vc_cs > Vc_maxUnacceptable[ s ]) Vc_maxUnacceptable[ s ] = Vc_cs; // here it is easy to update
          else { // else we search through all unacceptable ones and find the new max
            Vc_tmp = -INF;
            idx_cprimes = idx_cs;
            for( int cprime=0; cprime < nColleges; cprime++ ){
              if ( ID_cWorse[idx_cprimes] == -2 ) Vc_tmp = std::max(Vc_tmp, Vc[idx_cprimes]);
              // find position of next occurence of college cprime and student s
              idx_cprimes = ID_nextCollege[idx_cprimes] - 1; // convert to C indexing
            }
            Vc_maxUnacceptable[ s ] = Vc_tmp;
          }
        }
        
        
#ifdef DEBUG
        printf("Student %d's valuation over college %d:\n", s, c);
        printf("Vcupperbar = %f, Vclowerbar = %f\n", Vcupperbar, Vclowerbar);
        printf("Vc_hat = %f, Vc = %f\n", Xc_beta[ idx_cs_uword ], Vc_cs);
        printf("Vc_maxUnacceptable = %f, Vc_GroupMean = %f\n", Vc_maxUnacceptable[ s ], Vc_GroupMean[s]);
        if (idx_tmp==-2) printf("College %d is unacceptable to student %d\n", c, s);
#endif
        
        
        //----------------------------------------------------------------------------------//
        // --- college c's valuation over student s ----------------------------------------//
        //----------------------------------------------------------------------------------//
        
#ifdef DEBUG
        printf("\nCollege %d's valuation over student %d:\n", c, s);
#endif
        //--- equilibrium bounds -----------------------------------------------------------//
        
#ifdef DEBUG
          std::cout << "Determine college's eq. bounds ..\n" << std::flush;
#endif
        //  the following if-else applies only to colleges that are at full capacity:
        if ( !HasVacantSeats_c ) {
			
          if(match_cs == 0){  // case 1: non-equilibrium matches
		  
            // if student prefers college to current match (Eqn (A3))
      		  if (idx_eqCollege_s>=0) { // case 1.1: student s is matched to _some_ school
      			  if( Vc_cs > Vc[idx_eqCollege_s] ){ // if student would prefer going to s instead of his eq. college ...
      				  // ... then students' valuation has to be lower than that of worst student attending the college
      				  Vsupperbar = std::min( Vsupperbar, Vs_min_c );
      				}
      			} else if (idx_tmp!=-2) { // case 1.2: student s is not matched to any school, then if idx_tmp!=-2 indicates that school c is not unacceptable to student s
					    Vsupperbar = std::min( Vsupperbar, Vs_min_c ); // .. we can also bound this valuation
				    }
				
          } else {  // case 2: equilibrium matches
		  
      		  for(  int idx_csprime=idx_c_start; idx_csprime < idx_c_start+nObs_c; idx_csprime++ ){ // go through all students
						
						// case 2.1: student sprime is matched to _some_ school
    				  if ( (idx_eqCollege_sprime = ID_eqCollege[sid[idx_csprime]-1] -1) >= 0 ) { 
					    
    					  if( match[idx_csprime]==0 ){ // if college not matched to student sprime ...
      						// ... but student jprime values i over his equilibrium school (Eqn (A9))
      						if( Vc[idx_csprime] > Vc[ idx_eqCollege_sprime ] ){ 
      						  // ... then school must value student s higher than sprime (Eqn (A10), term 1).
      						  Vslowerbar = std::max( Vslowerbar, Vs[idx_csprime] );
      						}
    					  } // end if: student sprime is not matched to college c
						  
    				  } else if (ID_cWorse[idx_csprime]!=-2 ) { // case 2.2: student sprime is not matched to _any_ school, then if school c is not unacceptable to student sprime ...
					      
						    Vslowerbar = std::max( Vslowerbar, Vs[idx_csprime] ); // ... we can also bound the valuation
						    
					    } // end if
				    } // end for: loop over students
    	    } // non-equilibrium vs. equilibrium 
        } // end of vacant seats
		
#ifdef DEBUG
        printf("Equilibrium bounds: Vsupperbar = %f, Vslowerbar = %f\n", Vsupperbar, Vslowerbar);
#endif
#ifdef DEBUG2
        printf("0 ");  
#endif
        
        //--- rank order list bounds ---------------------------------------------------------//
        
#ifdef DEBUG
          std::cout << "Determine college's rank order bounds ..\n" << std::flush;
#endif
        if( (idx_tmp = ID_sBetter[idx_cs]) ){
          // student j's valuation over the college that j ranks just above college i.
          Vsupperbar = std::min( Vsupperbar, Vs[ idx_tmp - 1 ] );
        }  
        if( (idx_tmp = ID_sWorse[idx_cs]) > 0 ){ // recall that ID_sWorse and ID_sBetter uses R indexing, 0 indicates NA (no rank order bound)
          Vslowerbar = std::max( Vslowerbar, Vs[ idx_tmp - 1  ]);  // colleges valuation over the student that c ranks just below student s.
        } else if (idx_tmp==-1) { // -1 indicates that the lower bound is the max of all unacceptable students (could by -INF if no unacceptable ones exist)
            Vslowerbar = std::max( Vslowerbar, Vs_maxUnacceptable[ c ]);
        }
#ifdef DEBUG2
        printf("1 ");  
#endif
        
        
        //--- draw new valuations that respect the bounds -----------------------------------//
        
#ifdef DEBUG
        std::cout << "Drawing college's new latent valuation of student ..\n" << std::flush;
#endif
        Vs_cs = (Vslowerbar < Vsupperbar) ? TRUNCNORM(Xs_gamma[ idx_cs_uword ], 1.0, Vslowerbar, Vsupperbar) : Vsupperbar; // mu, sigma, lower, upper
	      Vs[idx_cs] = Vs_cs;
        if (demean) Vs_GroupMean[c] += Vs_cs;
#ifdef DEBUG2
        printf("2 ");  
#endif
        // update Vs_min if student j and college i are matched
        
          
        if (match_cs) {
          if ( Vs_cs < Vs_min_c  )  Vs_min[c] = Vs_cs ; // here it is easy to update
          else { // now we need to find a new minimum
            // Vs_min(t)(i) = arma::min( Vs(t)( M(t)(iprimeuvec,studentIds(L(t)(iprime))))) ;
            Vs_tmp = INF;
            for(  int idx_csprime=idx_c_start; idx_csprime < idx_c_start+nObs_c; idx_csprime++ ) {
              Vs_tmp = match[idx_csprime] ? std::min( Vs_tmp, Vs[idx_csprime]  ) : Vs_tmp;
            }
            Vs_min[c] = Vs_tmp;
          }
        }

#ifdef DEBUG2
        printf("3 ");  
#endif
        // update the max of unranked Vs if s is unacceptable to c
        
        
        if (idx_tmp==-2) { // if student s is unacceptable to college c, need to update
          if ( Vs_cs > Vs_maxUnacceptable[ c ]) Vs_maxUnacceptable[ c ] = Vs_cs; // here it is easy to update
          else { // else we search through all unacceptable ones and find the new max
            Vs_tmp = -INF;
            for(  int idx_csprime=idx_c_start; idx_csprime < idx_c_start+nObs_c; idx_csprime++ )
              if ( ID_sWorse[idx_csprime] == -2 ) Vs_tmp = std::max(Vs_tmp, Vs[idx_csprime]);
            Vs_maxUnacceptable[ c ] = Vs_tmp;
          }
        }
        
#ifdef DEBUG
        printf("Rank order bounds:  Vsupperbar = %f, Vslowerbar = %f\n", Vsupperbar, Vslowerbar);
        printf("Vs_hat = %f, Vs = %f\n", Xs_gamma[ idx_cs_uword ], Vs_cs);
        printf("Vs_maxUnacceptable = %f, Vs_min = %f, Vs_GroupMean = %f\n", Vs_maxUnacceptable[ c ], Vs_min[ c ], Vs_GroupMean[ c ]);
        if (idx_tmp==-2) printf("Student %d is unacceptable to college %d\n", s, c);
#endif
        
        
        // increment indices
        idx_cs += 1;
        idx_cs_uword += 1;
        
      } // end obs belonging to college i (students, but not necessarily all of them)
    } // end college i
    
		
    // compute Vs and Vc group means, and subtract means
    if (demean) {
#ifdef DEBUG
          std::cout << "Compute Vs and Vc group means, and subtract means ..\n" << std::flush;
#endif
  	  for (int s=0; s<nStudents; s++) {
  		  Vc_GroupMean[s] /= nObs_student[s];
  		  Vc_maxUnacceptable[s] -= Vc_GroupMean[s];
  	  }
  	  idx_cs = 0;
      for (int c=0; c <nColleges; c++){
          
  		  Vs_GroupMean[c] /= nObs_college[c];
        Vs_min[c] -= Vs_GroupMean[c];
  		  Vs_maxUnacceptable[c] -= Vs_GroupMean[c];
          
        for(int j=0; j < nObs_college[c]; j++){
          Vc[idx_cs] -= Vc_GroupMean[sid[idx_cs]-1];
          Vs[idx_cs] -= Vs_GroupMean[c];
          idx_cs++;
        }
      }
    }
    
    
    //----------------------------------------------------------------------------------//
    //--- compute new parameter vectors ------------------------------------------------//
    //----------------------------------------------------------------------------------//
    begin2 = clock();
    loop_time += (begin2 - begin);
    
    // update parameter vector (problem: Vc and Vs are NumericVector, so we need to cast to an arma::colvec
    beta_hat = Proj_Xc * Rcpp::as<arma::colvec>( Vc );
    gamma_hat = Proj_Xs * Rcpp::as<arma::colvec>( Vs );
    
    // draw new bet and gamma
    beta = mvrnormArma(beta_hat, XcXcInv_dense); //beta = Proj*Y;
    gamma = mvrnormArma(gamma_hat, XsXsInv_dense); //beta = Proj*Y;
    
#ifdef DEBUG
    printf("\nUpdated parameter vectors:\n");
    std::cout << "beta_hat: " << beta_hat << std::endl;
    std::cout << "beta: " << beta << std::endl;
    std::cout << "gamma_hat: " << gamma_hat << std::endl;
    std::cout << "gamma: " << gamma << std::endl;
#endif
    
    // save draws
    if ( (iter % thin == 0) ) {
      betavalues.row(iter/thin) = trans(beta);
      gammavalues.row(iter/thin) = trans(gamma);
    }
    
    // update mean valuations
    Xs_gamma = Xs * gamma;
    Xc_beta = Xc * beta;
    
    draw_time += (clock() - begin2);
		
		if (!fmod(iter, chunk10percent) ) {
      Rcpp::Rcout << " estimated remaining time: " << ((double) draw_time + (double) loop_time) / (CLOCKS_PER_SEC * (iter+1)) * (niter - iter) << "s" << std::endl;
    }
  } // end iter
  
  
	// print timing information
	std::cout << "Total time spent in main iteration loop: " << ((double) draw_time + (double) loop_time) / CLOCKS_PER_SEC << "s" <<std::endl;
  std::cout << "        ... of which drawing valuations: " << (double) loop_time / CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "       ... of which updating parameters: " << (double) draw_time / CLOCKS_PER_SEC << "s" << std::endl << std::flush;

  
  
  return List::create(  
    // parameter draws
    Named("betadraws") = betavalues,
    Named("gammadraws") = gammavalues,
    //Named("Vc") = Vc,
    //Named("Vs") = Vs,
    Named("Vc_hat") = Xc_beta,
    Named("Vs_hat") = Xs_gamma
  );
}



