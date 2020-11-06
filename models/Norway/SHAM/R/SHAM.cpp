// Fitting state-space models to seal populations with scarce data
// Tor Arne Øigård and Hans J. Skaug
// ICES Journal of Marine Science, 2014
// Contact:  Tor Arne Øigård  (toro@imr.no) and Hans Julius Skaug (skaug@imr.no)

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
  //Input data
  DATA_INTEGER(Amax);                      //Maximum age group
  DATA_INTEGER(SSstart);                   //Start of state-space process
  DATA_ARRAY(Cdata);                       //Catch data
  DATA_INTEGER(Nc);                        //Number of years with catch data
  DATA_ARRAY(pup_production);              //Pup production estimates
  DATA_INTEGER(Np);                        //Number of years with pup production data
  DATA_ARRAY(Fecundity);                   //Observed fecundity rates
  DATA_INTEGER(Nf);                        //Number of years with fecundity data
  DATA_VECTOR(Pmat);                        //Birth ogives
  DATA_INTEGER(Npred);                     //Number of years to run projections
  DATA_ARRAY(priors);                      //Priors for parameters
  DATA_INTEGER(Npriors);                   //Number of priors
  DATA_INTEGER(xmax);                      //Length of state-space process
  DATA_VECTOR(CQuota);                     //Catch level in future projections
  DATA_VECTOR(cap);                        //Capelin biomass
  DATA_VECTOR(cod);                        //Cod biomass
  
 //Estimated parameters
  PARAMETER(logK);                         //Initial population size
  PARAMETER(Mtilde);                       //Natural adult mortality
  PARAMETER(M0tilde);                      //Natural pup mortality
  PARAMETER(ftilde);                       //Mean fecundity rate
  PARAMETER(atilde);                       //AR(1) parameter
  PARAMETER(b1tilde);
  PARAMETER(b2tilde);
  PARAMETER(logSigma);                     //sigma
  PARAMETER_VECTOR(u);                     //random effect vector
  
  //Transform estimated parameters
  Type K = exp(logK);
  Type M = ilogit(Mtilde);
  Type M0 = ilogit(M0tilde);
  Type fmean = ilogit(ftilde);
  Type a = bound(atilde);
  Type b1 = ilogit(b1tilde);
  Type b2 = ilogit(b2tilde);
  Type Sigma = exp(logSigma);
  
  //Adult and pup survival
  Type em = exp(-M);
  Type em0 = exp(-M0);

  array<Type>  Cmatrix(Nc+Npred+2,3);           //Catch data
  array<Type>  Pm(Nc+Npred+1,Amax);             //Birth ogive
  vector<Type> N0(Nc+Npred+2);                  //Pup abundance
  vector<Type> N1(Nc+Npred+2);                  //Abundance of one year and older seals
  vector<Type> Ft(Nc+Npred+2);                  //Fecundity rates
  vector<Type> logitFt(Nc+Npred+2);             //logit Fecundity rates
  vector<Type> b(Npriors);                 	//Concatination of parameters

  // Concatination of parameters into b
  b(0) = K;                                      // Initial population size
  b(1) = M;                                      // Natural adult mortality
  b(2) = M0;                                     // Natural pup mortality
  b(3) = fmean;                                  // Mean fecundity rate before state-space process        
  b(4) = a;                                      // AR 1 parameter
  b(5) = b1;
  b(6) = b2;
  b(7) = Sigma;                                  // Noise variance
 

  // Preliminary calculations - Preparations of Catch data
  Cmatrix(0,0) = 1945;
  Cmatrix(0,1) = 0;
  Cmatrix(0,2) = 0;
  for(int i=1;i<(Nc+1);i++){
    Cmatrix(i,0) = Cdata(i-1,0);
    Cmatrix(i,1) = Cdata(i-1,1);
    Cmatrix(i,2) = Cdata(i-1,2);
    }

  for(int i=(Nc+1);i<(Nc+Npred+2);i++){
    Cmatrix(i,0) = Cmatrix(i-1,0)+1;
    Cmatrix(i,1) = CQuota(0);
    Cmatrix(i,2) = CQuota(1);
    }
    
  // Preliminary calculations - Preparations of birth ogive
  // (version with constant curve)
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
 
 
  // Prepare the fecundity-vector - hold fecundity constant up until year stmod
  // Logistic transformation to restrict to (0,1) - EQ 7
  Type counter = 0;
  for(int i=0;i<=(SSstart-Cdata(0,0));i++){
    Ft(CppAD::Integer(counter)) = fmean;
    logitFt(CppAD::Integer(counter)) = log((fmean+Type(1e-10))/(Type(1)+Type(1e-10)-fmean));  //Could use predefined functions from above instead....
    counter = counter + 1;
    }

  // Logistic transformation to restrict to (0,1) - EQ 7
  Ft(CppAD::Integer(counter)) = fmean;
  logitFt(CppAD::Integer(counter)) = log((fmean+Type(1e-10))/(Type(1)+Type(1e-10)-fmean)); 
  counter = counter + 1;

  Type logit_fmean = log((fmean+Type(1e-10))/(Type(1)+Type(1e-10)-fmean)); 
 
  // Initialize AR(1) process
  vector<Type> x(xmax);
  x(0) = 0;
 
  // AR(1)-process
  for(int i=1;i<xmax;i++)  
   {
     x(i) = a*x(i-1) + b1*cap(i) - b2*cod(i) + Sigma*u(i-1);		               // AR1 update
     logitFt(CppAD::Integer(counter)) = logit_fmean+x(i);      // Part of Eq 7
     Ft(CppAD::Integer(counter)) = exp(logitFt(CppAD::Integer(counter)))/(Type(1)+exp(logitFt(CppAD::Integer(counter))));  // Logistic transformation to restrict to (0,1) - EQ 7
     counter = counter + 1;
   }


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
  
  N1(0) = K;                            // Abundance of 1 year and older seals
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
  
    // Reqruitement equation
    for(int j=0;j<Amax;j++)
      N0(i) += 0.5*Ft(i)*Pm(i-1,j)*N(i,j);
   

    // Ensures that pup production is larger than pup catch
    N0(i) = posfun(N0(i)-Cmatrix(i,1),Type(1)) + Cmatrix(i,1);	
  }


  // Likelihood contribution
  Type nll=Type(0);

  //Data contribution to likelihood
  
  // Likelihood contribution from pup production estimates - EQ 9// Pup production estimates
  for(int i=0;i< Np;i++)
  {
    nll += -dnorm(N0(CppAD::Integer(pup_production(i,0)-Cdata(0,0)+1)),pup_production(i,1),pup_production(i,1)*pup_production(i,2),true);
  }

  // Likelihood contribution from fecundity estimates - EQ 10
  for(int i=0;i<Nf;i++)
  {
    nll += -dnorm(Ft(CppAD::Integer(Fecundity(i,0)-Cdata(0,0)+1)),Fecundity(i,1),Fecundity(i,2),true);
  }

   // Likelihood contribution from prior distibutions - EQ 11
  for(int i=0;i< Npriors;i++)
  {
    nll += -dnorm(b(i),priors(i,0),priors(i,1),true);
  }
  
  // Random effect contribution
  for(int i=0;i<xmax-1;i++)
  {
    nll += -dnorm(u(i),Type(0),Type(1),true);
  }

  ADREPORT(N0)
  ADREPORT(N1)
  ADREPORT(Ft)

  return nll;
}


