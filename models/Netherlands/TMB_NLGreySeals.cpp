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
  DATA_VECTOR(yr.pup);    // Year of pup count
  DATA_VECTOR(pup.uk);    // Model-estimates pup-counts UK
  DATA_VECTOR(tot.uk);    // Model-estimates total NS-UK population size
  DATA_VECTOR(molt);      // Observed moult counts UK
  DATA_VECTOR(yr.molt);   // Year of moult count
  DATA_VECTOR(smr);       // Observed summer counts
  DATA_VECTOR(yr.smr);    // Year of summer count
  DATA_INTEGER(tmax);     // Number of years, equal to length pup counts UK (include in code below?)
  DATA_INTEGER(K.smr);    // Length of summer counts (include in code below?)
  DATA_INTEGER(K.molt);   // Length of molt counts (include in code below?)
  DATA_INTEGER(K.pup);    // Length of pup counts (include in code below?) 

  //Estimated parameters
  PARAMETER(phi_a);        //Adult survival
  PARAMETER(phi_p);        //Juvenile survival
  PARAMETER(ad_ini);       //Initial number of adults at t0
  PARAMETER(alpha.pup);    //proportion of UK pups moving into NL
  PARAMETER(fec_a);        //fecundity
  PARAMETER(beta1);        //??
  PARAMETER(beta2);        //Change (i.e. annual forward shift in mean birth day)
  PARAMETER(mean_day);     //Mean birth date
  PARAMETER(pup_dur);      //duration of pups on land
  //PARAMETER(corf2);        //Haulout probability (fixed)
  PARAMETER(alpha.molt);   //Proportion of UK adults moving into NL during moult
  PARAMETER(alpha.smr);    //Proportion of UK adults moving into NL during summer
  //PARAMETER(beta_fem);     //Fraction of females (fixed)

  //Transform estimated parameters
  Type frac.fem.sto = (1+beta_fem)/(2+beta_fem) // is this working?
  Type frac.fem.cut = cut(frac.fem.sto)
  Type corf.sto = ilogit(corf2)
  Type corf.cut = cut(corf.sto) // cut makes sure that prior does not contribute to likelihood function
  Type surv.ad = 0.8+0.2*phi_a
  Type fec.ad =  0.6+0.4*fec_a
  Type surv.pup = phi_p * surv.ad // is this working? 
 
  // Create empty matrixes and vectors
  array<Type>   Nat(7,tmax);       // Numbers at age (7,tmax) or (tmax,7)
  vector<Type>  N_est1(Kpup);      // 0 year olds in line with observed pup counts
  vector<Type>  bp(Kpup)           // birth arrival probability 
  vector<Type>  lp(Kpup)           // birth departure probability curve
  vector<Type>  pup_prop[j](Kpup)  // pup presence probability:  bp[j]-lp[j]
  vector<Type>  N27(Kmolt)         // Total number of 1-7 year olds, same length as moult counts 
  vector<Type>  N17(Ksmr)          // Total population size, same length as summer counts
  
  // Preliminary calculations for first year
  Nat(2,0) = alpha.pup * pup.uk[1];
  Nat(3,0) = 0;
  Nat(4,0) = 0;
  Nat(5,0) = 0;
  Nat(6,0) = 0;
  Nat(7,0) = ad_ini
  Nat(1,0) = fec.ad * Nat(7,0) * frac.fem.cut
    
  //---------------------------------
  // Calculate population trajectory
  for(int i=1;i<tmax;i++)
  {
    Nat(0,i) = fec.ad * frac.fem.cut * (Nat(7,i-1)  // I assume 0=age1
    Nat(1,i) = alpha.pup * pup.uk(t) + surv.pup * Nat(0,i-1)
    
    for(int j=2;j<5;j++) //j=2=age3
    {
      N(j,i) = surv.ad * Nat(j,i-1)
    }
  
    Nat(6,i) = surv.ad * (Nat(5,i-1) + Nat(6,i-1))
   }
    
  // Calculate pup probabilities
     for(int i=1;i<Kpup;i++)
             bp(i)= ilogit(-mean_day*beta1+beta1*daynr(i) -beta1*beta2*yr.pup(i))
             lp(i)= ilogit(-mean_day*beta1+beta1*daynr(i) -beta1*beta2*yr.pup(i)-beta1*pup_dur)
             pup_prop(i) = bp(i)-lp(i)
           }
   
  // Calculate expected pup counts
      for(int i=1;i<Kpup;i++)
        N_est1(i) <- 0.001+ Nat(0,yr.pup(i))*pup_prop(i)
             
  // Calculate expected molt counts
      for(int i=1;i<K.molt;i++)
        N27(i) <- corf.cut * surv.pup.cut * Nat.cut(1,yr.molt(ik)) + pow(surv.ad.cut,106/365) * sum(Nat.cut(2:7,yr.molt(i))) + alpha.molt * tot.uk(yr.molt(i)+1)
          
  // Calculate expected molt counts
      for(int i=1;i<K.molt;i++)
        N17(i) <- corf.cut * surv.pup.cut * Nat.cut(1,yr.smr(ik)) + pow(surv.ad.cut,207/365) * sum(Nat.cut(2:7,yr.smr(i))) + alpha.smr * tot.uk(yr.molt(i)+1)
  
  
   
  //--------------------------------
  // Likelihood contribution
      Type nll = Type(0.0);

  //Data contribution to likelihood

  // Likelihood contribution from pup counts
     for(int i=0;i< Kpup;i++)
     {
     nll += -dpoisson(N_est1(i),pup(i),true);
     }

  // Likelihood contribution from molt counts
     for(int i=0;i< K.molt;i++)
     {
     nll += -dpoisson(N27(i),molt(i),true);
     }

  // Likelihood contribution from summer counts
     for(int i=0;i< K.smr;i++)
     {
     nll += -dpoisson(N17(i),smr(i),true);
     }

  // Likelihood contribution from informative prior distibutions - EQ 11
     nll += -dbeta(phi_a,1.6.1.2,true);
     nll += -dbeta(phi_p,2.87,1.78,true);
     nll += -dbeta(fec_a,2,1.5,true);
     nll += -dbeta(phi_a,1.6.1.2,true);
   
  // The following priors are fixed; i.e. parameters should not be estimated by the model
     nll += -dnorm(corf2,-1.0935,21.51307,true);
     nll += -dbeta(beta_fem,2,8,true);
         
  //----------------------------------------------------------------------------
  // Report section, define which quantities should have uncertainty estimates
     ADREPORT(N_est1)
     ADREPORT(N17)
     ADREPORT(N27)  
     ADREPORT(Nat)
     return nll;

}
