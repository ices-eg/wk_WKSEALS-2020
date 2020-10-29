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
  DATA_ARRAY(Cdata);                       //Catch data
  DATA_INTEGER(Nc);                        //Number of years with catch data
  DATA_ARRAY(pupProductionData);          //Pup production estimates
  DATA_INTEGER(Np);                        //Number of years with pup production data
  DATA_VECTOR(Ftmp);                        //Observed fecundity rates
  DATA_ARRAY(Pmat);                        //Birth ogives
  DATA_INTEGER(Npred);                     //Number of years to run projections
  DATA_ARRAY(priors);                      //Priors for parameters
  DATA_INTEGER(Npriors);                   //Number of priors
  DATA_VECTOR(CQuota);                     //Catch level in future projections
  // Test
  //Estimated parameters
  PARAMETER(logK);                         //Initial population size
  PARAMETER(Mtilde);                       //Natural adult mortality
  PARAMETER(M0tilde);                      //Natural pup mortality


  //Transform estimated parameters
  Type K = exp(logK);
  Type M = ilogit(Mtilde);
  Type M0 = ilogit(M0tilde);

  //Adult and pup survival
  Type em = exp(-M);
  Type em0 = exp(-M0);

  array<Type>  Catch(Nc+Npred+2,3);             //Catch data
  array<Type>  P(Nc+Npred+1,Amax);              //Birth ogive
  //array<Type>  pupProd(Np,3);                 //Pup production estimates
  //array<Type>  Nout(Nc+Npred+2,3);
  array<Type>  N(Nc+Npred+2,Amax);              // Population matrix
  //vector<Type> years(Nc);                     //Years with catch data
  vector<Type> mub(Npriors);                    //Mean value of priors
  vector<Type> sdb(Npriors);                    //SD of priors
  vector<Type> N0(Nc+Npred+2);                  //Pup abundance
  vector<Type> N1(Nc+Npred+2);                  //Abundance of one year and older seals
  vector<Type> NTot(Nc+Npred+2);                //Total population
  //vector<Type> Ntot(Nc+1);                    //Total population
  vector<Type> Ft(Nc+Npred+2);
  vector<Type> b(Npriors);                      //Concatenation of parameters
  Type D1 = Type(0.0);
  Type DNmax = Type(0.0);
  Type N0CurrentYear = Type(0.0);

  Type N1CurrentYear = Type(0.0);
  Type NTotCurrentYear = Type(0.0);
  Type NTotmax = Type(0.0);
  Type NTotPred = Type(0.0);

  // Preliminary calculations - Preparations of Catch data and fecundity rates
  Catch(0,0) = 1945;
  Catch(0,1) = 0;
  Catch(0,2) = 0;
  Ft(0) = Type(0.0);
  for(int i=1;i<(Nc+1);i++){
    Catch(i,0) = Cdata(i-1,0);      //Year
    Catch(i,1) = Cdata(i-1,1);      //Pup catch
    Catch(i,2) = Cdata(i-1,2);      //1+ catch
    Ft(i) = Ftmp(i-1);
  }

  //CHECK END LIMIT LATER, MAYBE NC+NPRED+1 IS OK
  for(int i=(Nc+1);i<(Nc+Npred+2);i++){
    Catch(i,0) = Catch(i-1,0)+1;    //Year
    Catch(i,1) = CQuota(0);
    Catch(i,2) = CQuota(1);
    //Catch(i,1) = Type(0.0);
    //Catch(i,2) = Type(0.0);
    Ft(i) = Ftmp(Nc-1);
  }

  // Preliminary calculations - Preparations of birth ogive - does not look like this is needed, just call P = Pmat
  for(int i=0;i<(Nc);i++){
    for(int j=0;j<Amax;j++){
      P(i,j) = Pmat(i,j);
    }
  }


  for(int i=(Nc);i<(Nc+Npred+1);i++){
    for(int j=0;j<Amax;j++){
      P(i,j) = Pmat(Nc-1,j);
    }
  }

  /*for(int i=0;i<Nc;i++){
    years(i)=Catch(i+1);
  }
  */

  /* REMOVE THIS
  //Transform pup production estimates
  for(int i=0;i<Np;i++){
    pupProd(i,0) = pupProductionData(i,0)- years(0) + 1; //Year index
    pupProd(i,1) = pupProductionData(i,1);               //Mean
    pupProd(i,2) = pupProd(i,1)*pupProductionData(i,2);  //SD
  }
   */

  // Extract priors
  mub(0) = priors(0,0);
  mub(1) = priors(1,0);
  mub(2) = priors(2,0);
  sdb(0) = priors(0,1);
  sdb(1) = priors(1,1);
  sdb(2) = priors(2,1);

  //Rcout << "Priors:\n" << "mean = " << mub << endl << "sd   = " << sdb << endl;

  b(0) = K;
  b(1) = M;
  b(2) = M0;

 //----------------------------------------------
 // Initiate of N in year 0 (1945) - EQ 1 and 2
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


 //---------------------------------
 // Calculate population trajectory
 for(int i=1;i<Nc+Npred+2;i++)
 {
   N(i,0) = (N0(i-1)-Catch(i-1,1))*em0;			// 0-group from last year

   for(int j=1;j<Amax;j++)
   {
     N(i,j) = N(i-1,j-1)*(Type(1)-Catch(i-1,2)/N1(i-1))*em;	// Pro-rata distribution of catches
   }

   N(i,Amax-1) += N(i-1,Amax-1)*(Type(1)-Catch(i-1,2)/N1(i-1))*em;      // A+ group correction

   // Ensures that N > 0
   Nsum = Type(0.0);
   for(int j = 0;j<Amax;j++)
     Nsum = Nsum + N(i,j);

   N1(i) = posfun(Nsum-Catch(i,2),Type(1.0)) + Catch(i,2);

   for(int j = 0;j<Amax;j++)
     N(i,j) = N1(i)*N(i,j)/Nsum;			// Scales up the full age structure

   // Recruitment equation
   for(int j=0;j<Amax;j++)
     N0(i) += 0.5*Ft(i)*P(i-1,j)*N(i,j);

   // Ensures that pup production is larger than pup catch
   N0(i) = posfun(N0(i)-Catch(i,1),Type(1.0)) + Catch(i,1);

   // Calculate the D(1+) statistic
   NTotPred = Type(0.0);
   if(Catch(i,0)==Catch(Nc+1,0))
   {
     D1 = 1/(N1(i));
     N0CurrentYear = N0(i);
     N1CurrentYear = N1(i);
     NTotCurrentYear = (N0(i)+N1(i));
   }
   
   if(Catch(i,0)==Catch(Nc+Npred+1,0)){
     D1 *= (N1(i));
     DNmax = (N1(i)+N0(i));
     NTotPred = DNmax;
   }

 }
 
 // //array<Type> tmpNtot(Nc+Npred+2);
 // vector<double> tmpNtot(Nc+Npred+2);
 // 
 for(int i=0;i<Nc+Npred+2;i++)
 {
    NTot(i) = (N0(i) + N1(i));
}
 // 
 // for(int i =0; i<Nc+Npred+2; ++i){
 //   tmpNTot(i) = ((double)CppAD::Integer(NTot(i)));
 // }
 // 
 // NTotmax = Type(-999999);
 // for(int i=0;i<Nc+1;i++)
 // {
 //   if(NTotmax < tmpNTot(i)){
 //     NTotmax = tmpNTot(i);
 //   }
 // }
 // 
 
 NTotmax = max(NTot);
 
 DNmax = NTotPred/NTotmax;
 
 //NTotmax = NTot(CppAD::Integer(maxind));
 //DNmax /= NTot(CppAD::Integer(maxind));
 
 //--------------------------------
 // Likelihood contribution
 Type nll = Type(0.0);

 //Data contribution to likelihood

 // Likelihood contribution from pup production estimates
 for(int i=0;i< Np;i++)
 {
   nll += -dnorm(N0(CppAD::Integer(pupProductionData(i,0)-Cdata(0,0)+1)),pupProductionData(i,1),pupProductionData(i,1)*pupProductionData(i,2),true);
 }

 // Likelihood contribution from prior distibutions - EQ 11
 for(int i=0;i< Npriors;i++)
 {
   nll += -dnorm(b(i),priors(i,0),priors(i,1),true);
 }

 // Likelihood penalty to avoid negative N
 for(int i=0;i< Nc;i++)
   nll += .001*(fabs(N1(i)-10000)-(N1(i)-10000))*(fabs(N1(i)-10000)-(N1(i)-10000));


 /*
 // Likelihood contribution from fecundity estimates - EQ 10
 for(int i=0;i<Nf;i++)
 {
   nll += -dnorm(Ft(CppAD::Integer(Fecundity(i,0)-Cdata(0,0)+1)),Fecundity(i,1),Fecundity(i,2),true);
 }
 */

 /*
 // Random effect contribution
 for(int i=0;i<xmax-1;i++)
 {
    nll += -dnorm(u(i),Type(0),Type(1),true);
 }
 */

  //----------------------------------------------------------------------------
  // Report section, define which quantities should have uncertainty estimates
  ADREPORT(N0)
  ADREPORT(N1)
  ADREPORT(NTot)  
  ADREPORT(N0CurrentYear)
  ADREPORT(N1CurrentYear)
  ADREPORT(NTotCurrentYear)
  ADREPORT(D1)
  ADREPORT(DNmax)
  ADREPORT(NTotPred)
  ADREPORT(NTotmax)
  REPORT(Catch)
  return nll;

  }
