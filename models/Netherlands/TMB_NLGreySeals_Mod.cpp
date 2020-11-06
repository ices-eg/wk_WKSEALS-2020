#include <TMB.hpp>
#include <iostream>

//Inv.logit Transform
template <class Type>
Type ilogit(Type x){
  return exp(x)/(Type(1.0)+exp(x));
}

//Logit Transform
template <class Type>
Type logitt(Type x){
  return log(x/(Type(1)-x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  //Input data
  DATA_VECTOR(pup);       // Observed pup counts NL
  DATA_VECTOR(daynr);     // Day number of NL pup count
  DATA_IVECTOR(yr_pup);    // Year of pup count as integer
  DATA_VECTOR(yr_pup_v);    // Year of pup count as vector
  DATA_VECTOR(pup_uk);    // Model-estimates pup-counts UK
  DATA_VECTOR(tot_uk);    // Model-estimates total NS-UK population size
  DATA_VECTOR(molt);      // Observed moult counts UK
  DATA_IVECTOR(yr_molt);   // Year of moult count
  DATA_VECTOR(smr);       // Observed summer counts
  DATA_IVECTOR(yr_smr);    // Year of summer count
  DATA_INTEGER(tmax);     // Number of years, equal to length pup counts UK (include in code below?)
  DATA_INTEGER(K_smr);    // Length of summer counts (include in code below?)
  DATA_INTEGER(K_molt);   // Length of molt counts (include in code below?)
  DATA_INTEGER(K_pup);    // Length of pup counts (include in code below?) 

  //Estimated parameters
  PARAMETER(phi_a);        //Adult survival
  PARAMETER(phi_p);        //Juvenile survival
  PARAMETER(ad_ini);       //Initial number of adults at t0
  PARAMETER(alpha_pup);    //proportion of UK pups moving into NL
  PARAMETER(fec_a);        //fecundity
  PARAMETER(beta1);        //??
  PARAMETER(beta2);        //Change (i.e. annual forward shift in mean birth day)
  PARAMETER(mean_day);     //Mean birth date
  PARAMETER(pup_dur);      //duration of pups on land
  PARAMETER(corf2);        //Haulout probability (fixed, should not be estimated, how to?)
  PARAMETER(alpha_molt);   //Proportion of UK adults moving into NL during moult
  PARAMETER(alpha_smr);    //Proportion of UK adults moving into NL during summer
  PARAMETER(beta_fem);     //Fraction of females (fixed)

  //Transform estimated parameters
  Type frac_fem_sto = (1+beta_fem)/(2+beta_fem); // is this working?
  Type corf_sto = ilogit(corf2);
  Type surv_ad = 0.8+0.2*phi_a;
  Type fec_ad =  0.6+0.4*fec_a;
  Type surv_pup = phi_p * surv_ad; // is this working? 
 
  // Create empty matrixes and vectors
  array<Type>   Nat(7,tmax);         // Numbers at age (7,tmax) or (tmax,7)  // To have multiple populations you need 3 dimensions.
  vector<Type>  N_est1(K_pup);       // 0 year olds in line with observed pup counts
  vector<Type>  bp(K_pup);           // birth arrival probability 
  vector<Type>  lp(K_pup);           // birth departure probability curve
  vector<Type>  pup_prop(K_pup);     // pup presence probability:  bp[j]-lp[j]
  vector<Type>  N27(K_molt);         // Total number of 1-7 year olds, same length as moult counts 
  vector<Type>  N17(K_smr);          // Total population size, same length as summer counts
  
  // Preliminary calculations for first year
  Nat(1,0) = alpha_pup * pup_uk(0);
  Nat(2,0) = 0;
  Nat(3,0) = 0;
  Nat(4,0) = 0;
  Nat(5,0) = 0;
  Nat(6,0) = ad_ini;
  Nat(0,0) = fec_ad * Nat(6,0) * frac_fem_sto;
    
  //---------------------------------
  // Calculate population trajectory
  for(int i=1;i<tmax;i++)
  {
    Nat(0,i) = fec_ad * frac_fem_sto* Nat(6,i-1);  // I assume 0=age1
    Nat(1,i) = alpha_pup * pup_uk(i-1) + surv_pup * Nat(0,i-1);
    
    for(int j=2;j<6;j++) //j=2=age3
    {
      Nat(j,i) = surv_ad * Nat(j-1,i-1);
    }
  
    Nat(6,i) = surv_ad * (Nat(5,i-1) + Nat(6,i-1));
   }
    
  // Calculate pup probabilities
     for(int i=0;i<K_pup;i++){
             bp(i)= ilogit(-mean_day*beta1+beta1*daynr(i) -beta1*beta2*yr_pup_v(i));
             lp(i)= ilogit(-mean_day*beta1+beta1*daynr(i) -beta1*beta2*yr_pup_v(i)-beta1*pup_dur);
             pup_prop(i) = bp(i)-lp(i);
           }
   
  // Calculate expected pup counts
      for(int i=0;i<K_pup;i++)
        N_est1(i) = 0.001 + Nat(0,yr_pup(i))*pup_prop(i);
             
             
  
  // Calculate expected molt counts
        
        for(int i=0;i<K_molt;i++){
          Type tmp= Type(0.0);
          for(int ii=1;ii<=6;ii++){ //age 2 to 7
            tmp = tmp + Nat(ii,yr_molt(i));
          }
          Nat.block(1,yr_molt(i),6,1).sum();  // The more compacte R-style code just to check that it compiles
        N27(i) = corf_sto * surv_pup * Nat(0,yr_molt(i)) + pow(surv_ad,106.0/365.0) * tmp + alpha_molt * tot_uk(yr_molt(i)+1);
      }
          
  // Calculate expected summer counts
      for(int i=0;i<K_smr;i++){
        Type tmp= Type(0.0);
        for(int ii=0;ii<=6;ii++){
          tmp = tmp + Nat(ii,yr_smr(i));
        }
        
          N17(i) = corf_sto * surv_pup * Nat(0,yr_smr(i)) + pow(surv_ad,207.0/365.0) * tmp + alpha_smr * tot_uk(yr_smr(i)+1);

      }
  
  
   
  //--------------------------------
  // Likelihood contribution
      Type nll = Type(0.0);

  //Data contribution to likelihood

  // Likelihood contribution from pup counts
     for(int i=0;i< K_pup;i++)
     {
     nll += -dpois(pup(i),N_est1(i),true);
     }

  // Likelihood contribution from molt counts
     for(int i=0;i< K_molt;i++)
     {
     nll += -dpois(molt(i),N27(i),true);
     }

  // Likelihood contribution from summer counts
     for(int i=0;i< K_smr;i++)
     {
     nll += -dpois(smr(i),N17(i),true);
     }

  // Likelihood contribution from informative prior distibutions - EQ 11
     nll += -dbeta(phi_a,Type(1.6),Type(1.2),true);
     nll += -dbeta(phi_p,Type(2.87),Type(1.78),true);
     nll += -dbeta(fec_a,Type(2),Type(1.5),true);
     
  // The following priors are fixed; i.e. parameters should not be estimated by the model
     nll += -dnorm(corf2,Type(-1.0935),Type(21.51307),true);
     nll += -dbeta(beta_fem,Type(2),Type(8),true);
         
  //----------------------------------------------------------------------------
  // Report section, define which quantities should have uncertainty estimates
     ADREPORT(N_est1)
     ADREPORT(N17)
     ADREPORT(N27)  
     ADREPORT(Nat)
     return nll;

}
