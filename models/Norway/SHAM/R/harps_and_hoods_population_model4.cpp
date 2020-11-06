// Including capelin and cod biomass as covariates on fecundity
// Martin Biuw and Tor Arne Øigård
// Modified from:
// Fitting state-space models to seal populations with scarce data
// Tor Arne Øigård and Hans J. Skaug
// ICES Journal of Marine Science, 2014
// Contact:  Martin Biuw (martin.biuw@hi.no) or Tor Arne Øigård (tor.arne.oigard@nr.no)

#include <TMB.hpp>
#include <iostream>

//Bound between 0 and 1
template <class Type>
Type ilogit(Type x){
  return Type(1.0)/(Type(1.0)+exp(-x));
}

//Bound between -1 and 1
template <class Type>
Type bound(Type x){
  return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}

//Logit Transform
template <class Type>
Type logitt(Type x){
  return log(x/(Type(1)-x));
}

// t distribution
template <class Type>
Type dt1(Type x){
  return Type(1.0)/M_PI/(Type(1.0)+x*x);
}


// Mixture density (normal + t distribution)
template <class Type>
Type nldens(Type x, Type mu, Type sd, Type p){
  Type z=(x-mu)/sd;
  return -log(Type(1.0)/sd*((Type(1.0)-p)*dnorm(z,Type(0.0),Type(1.0),false)+p*dt1(z)));
}

//posfun-function
template <class Type>
Type posfun(Type x, Type eps){
  if (x>=eps)
  {
    return x;
  }
  else
  {
    Type y=1.0-x/eps;
    return eps*(1./(1.+y+y*y));
  }
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  //-------------------------------------------------------------------
  // DATA SECTION
  // Input data:
  
  DATA_INTEGER(Amax);                      //Maximum age group
  DATA_ARRAY(Cdata);                       //Catch data
  DATA_INTEGER(Nc);                        //Number of years with catch data
  DATA_ARRAY(pup_production);              //Pup production estimates
  DATA_INTEGER(Np);                        //Number of years with pup production data
  DATA_ARRAY(Fecundity);                   //Observed fecundity rates
  DATA_INTEGER(Nf);                        //Number of years with fecundity data
  DATA_VECTOR(Pmat);                       //Collapsed maturity curve
  DATA_INTEGER(Npred);                     //Number of years to run projections
  DATA_ARRAY(priors);                      //Priors for parameters
  DATA_INTEGER(Npriors);                   //Number of priors
  DATA_VECTOR(CQuota);                     //Catch level in future projections
  DATA_VECTOR(cap);                        //Capelin biomass
  DATA_VECTOR(cod);                        //Cod biomass

  // Objects to hold derived data:
  
  array<Type>  Cmatrix(Nc+Npred+2,3);           //Catch data
  array<Type>  Pm(Nc+Npred+1,Amax);             //Birth curve
  vector<Type> N0(Nc+Npred+2);                  //Pup abundance
  vector<Type> N1(Nc+Npred+2);                  //Abundance of one year and older seals
  vector<Type> b(Npriors);                 	    //Concatenation of parameters
  vector<Type> Ft(Nc+Npred+2);                  //Fecundity rates
  vector<Type> logitFt(Nc+Npred+2);             //logit Fecundity rates
  
  
  //-------------------------------------------------------------------
  // PARAMETER SECTION
  // Input parameters:
  
  PARAMETER(logK);                         //Initial population size
  PARAMETER(Mtilde);                       //Natural adult mortality
  PARAMETER(M0tilde);                      //Natural pup mortality
  PARAMETER(ftilde);                       //Mean fecundity rate
//  PARAMETER(logBeta);                      //Intercept, linear predictor 
  PARAMETER(logBetaCap);                   //Slope, linear predictor for capelin
  PARAMETER(logBetaCod);                   //Slope, linear predictor for cod
  PARAMETER(logSdProc);                    //Process error, fecundity
  PARAMETER(logSdObs);                     //Observation error, fecundity
//  PARAMETER(logitp);                       //Mixture proportion for observation error
  PARAMETER_VECTOR(F);                     //Random effect vector (unobserved pregnancy rate)

  
  // Transformations:
  Type K = exp(logK);
  Type M = ilogit(Mtilde);
  Type M0 = ilogit(M0tilde);
  Type fmean = ilogit(ftilde);
//  Type beta = exp(logBeta);
  Type betaCap = exp(logBetaCap);
  Type betaCod = exp(logBetaCod);
  Type sdProc = exp(logSdProc);
  Type sdObs = exp(logSdObs);
//  Type p = ilogit(logitp);
  
  // Concatenation of parameters into b
  b(0) = K;
  b(1) = M;
  b(2) = M0;
  b(3) = fmean;
//  b(4) = beta; Equivalent to fmean?
  b(4) = betaCap;
  b(5) = betaCod;
  b(6) = sdProc;
  b(7) = sdObs;
//  b(8) = p;
  
  //Convert mortality to survival
  Type em = exp(-M);
  Type em0 = exp(-M0);
  
  
  //-------------------------------------------------------------------
  // PRELIMINARY CALCULATIONS 
  // Fecundity estimation:  

  Type counter = 0;
  Ft(CppAD::Integer(counter)) = fmean; 
  // Or should fmean be treated as the beta (intercept) parameter? How? 
  // By again including it in the linear predictor as logit(fmean)? 
  for(int i=1;i<Nc;i++){
  //  logitFt(CppAD::Integer(counter)) = beta+betaCap*cap(i-1)+betaCod*cod(i-1);
    logitFt(CppAD::Integer(counter)) = logit(fmean)+betaCap*cap(i-1)+betaCod*cod(i-1);
    Ft(CppAD::Integer(counter)) = exp(logitFt(CppAD::Integer(counter)))/(Type(1)+exp(logitFt(CppAD::Integer(counter))));  // Logistic transformation to restrict to (0,1) - EQ 7
    counter = counter + 1;
  }
  
  // Complete the fecundity-vector for projections
  // Hold fecundity constant after last prey/competition data

  for(int i=0;i<=(Npred+1);i++){
    Ft(CppAD::Integer(counter)) = fmean;
    logitFt(CppAD::Integer(counter)) = log((fmean+Type(1e-10))/(Type(1)+Type(1e-10)-fmean));  //Could use predefined functions from above instead....
    counter = counter + 1;
  }
  
  
  // Catch data
  Cmatrix(0,0) = 1945;
  Cmatrix(0,1) = 0;
  Cmatrix(0,2) = 0;
  for(int i=1;i<(Nc+1);i++){
    Cmatrix(i,0) = Cdata(i-1,0);
    Cmatrix(i,1) = Cdata(i-1,1);
    Cmatrix(i,2) = Cdata(i-1,2);
  }
  
  for(int i=(Nc+1);i<(Nc+Npred+1);i++){
    Cmatrix(i,0) = Cmatrix(i-1,0)+1;
    Cmatrix(i,1) = CQuota(0);
    Cmatrix(i,2) = CQuota(1);
  }
  
  
  // Maturity curves (version with constant curve)
  for(int i=0;i<(Nc);i++){
    for(int j=0;j<Amax;j++){
      Pm(i,j) = Pmat(j);
    }
  }
  
  for(int i=(Nc);i<(Nc+Npred+1);i++){
    for(int j=0;j<Amax;j++){
      Pm(i,j) = Pmat(j);
    }
  }
  
  
  //-------------------------------------------------------------------
  // POPULATION TRAJECTORY
  
  // Initiate of N in year 0 (1945) - EQ 1 and 2
  matrix<Type> N(Nc+Npred+2,Amax);
  for(int i=0;i<Amax;i++){
    N(0,i) = exp(-Type(i)*M);		// Adults
  }
  
  N(0,Amax-1) /= 1-em;			// Correct A+ group
  
  Type Nsum = 0;
  for(int i = 0;i<Amax;i++)
    Nsum = Nsum + N(0,i);
  
  for(int i = 0;i<Amax;i++)
    N(0,i) = K*N(0,i)/Nsum;		// Normalize vector to K
  
  N1(0) = K;                // Abundance of 1 year and older seals
  N0(0) = (1-em)/em0*K;			// To balance natural mortality of 1+ group
  
  // Calculate population trajectory
  for(int i=1;i<Nc+Npred+2;i++)
  {
    N(i,0) = (N0(i-1)-Cmatrix(i-1,1))*em0;			// 0-group from last year - EQ 3
    for(int j=1;j<Amax;j++)
    {
      N(i,j) = N(i-1,j-1)*(Type(1)-Cmatrix(i-1,2)/N1(i-1))*em;	// Pro-rata distribution of catches EQ 3 and 4
    }
    
    N(i,Amax-1) += N(i-1,Amax-1)*(Type(1)-Cmatrix(i-1,2)/N1(i-1))*em;      // A+ group correction - EQ 3
    
    // Ensures that N > 0
    Nsum = 0;
    for(int j = 0;j<Amax;j++)
      Nsum = Nsum + N(i,j);
    
    N1(i) = posfun(Nsum-Cmatrix(i,2),Type(1)) + Cmatrix(i,2);
    for(int j = 0;j<Amax;j++)
      N(i,j) = N1(i)*N(i,j)/Nsum;			// Scales up the full age structure
    
    // Recruitement equation
    for(int j=0;j<Amax;j++)
       N0(i) += 0.5*Ft(i)*Pm(i-1,j)*N(i,j);

    // Ensures that pup production is larger than pup catch
    N0(i) = posfun(N0(i)-Cmatrix(i,1),Type(1)) + Cmatrix(i,1);
  }
  
  //-------------------------------------------------------------------
  // LIKELIHOOD CONTRIBUTIONS:
  
  Type nll=Type(0);
  
  // Likelihood contribution from pup production estimates
  for(int i=0;i< Np;i++)
  {
    nll += -dnorm(N0(CppAD::Integer(pup_production(i,0)-Cdata(0,0)+1)),pup_production(i,1),pup_production(i,1)*pup_production(i,2),true);
  }
  
  // Likelihood contributions from process error in fecundity
  for(int i=0;i<Nc+1;i++)
  {
    nll += -dnorm(Ft(i), fmean, sdProc); // Is fmean correct here???
  }  
  // Likelihood contribution from fecundity estimates 
  // CURRENTLY DOES NOT account for SE VALUES FROM FECUNDITY OBSERVATIONS
  for(int i=0;i<Nf;i++)
  {
  // With gaussian errors only:
    nll += -dnorm(Fecundity(i,1), F(CppAD::Integer(Fecundity(i,0)-Cdata(0,0)+1)), sdObs); 
  
  // With mixture of gaussian and t distributed errors:
  //  nll += -nldens(Fecundity(i,1), F(CppAD::Integer(Fecundity(i,0)-Cdata(0,0)+1)), sdObs, p); 
  }
  
  // Likelihood contribution from prior distibutions
  for(int i=0;i< Npriors;i++)
  {
    nll += -dnorm(b(i),priors(i,0),priors(i,1),true);
  }
  
  ADREPORT(N0)
  ADREPORT(N1)
  ADREPORT(Ft)

  return nll;
}
