

################################################################################
#                                                                              #
#                              R-CODE Brasseur et al.                          #
#     RAPID RECOVERY OF DUTCH GRAY SEAL COLONIES FUELLED BY IMMIGRATION        #
#                 Code written by Geert Aarts and Tim Gerrodette               #
#                                                                              #
################################################################################



#!  SET WORKING DIRECTORY, LOAD DATA AND READ LIBRARIES

  # Rx64 3.6.0
  
  # Winbugs installation (Windows 10)
  # Go to https://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/ 
  # Select Download WinBUGS
  # Unzip file and copy folder WinBUGS14 into e.g. c:/program_files/
  # No installation needed
  
  # Load library
    library(boot)
    library(R2WinBUGS)
    library(coda)
    library(fields)
    library(RODBC)
   
  # Set github directory
    setwd("C:/github/wk_WKSEALS-2020")
    
  # Load data
    load("data/Netherlands/NLGreySeals.rdata")
    
                
################################################################################

#!  BUGS MODEL FITTED TO ALL DATA

################################################################################

    sink("models/Netherlands/Bugs_model_puponly.txt")
    cat("model {
    
      # SUMMER model
           Nat[2,1] <- alpha.pup * pup.uk[1]
           Nat[3,1] <- 0
           Nat[4,1] <- 0
           Nat[5,1] <- 0
           Nat[6,1] <- 0
           Nat[7,1] <- ad_ini
           Nat[1,1] <- fec.ad * Nat[7,1] * frac.fem.cut
           
           for (t in 2:tmax){
             Nat[1,t] <- fec.ad * frac.fem.cut * (Nat[7,t])
             Nat[2,t] <- alpha.pup * pup.uk[t] + surv.pup * Nat[1,t-1]
             Nat[3,t] <- surv.ad * Nat[2,t-1]
             Nat[4,t] <- surv.ad * Nat[3,t-1]
             Nat[5,t] <- surv.ad * Nat[4,t-1]
             Nat[6,t] <- surv.ad * Nat[5,t-1]
             Nat[7,t] <- surv.ad * (Nat[6,t-1] + Nat[7,t-1])
             }
        
           for (j in 1:K.pup){
             pup[j] ~ dpois(N_est1[j])
             N_est1[j] <- 0.001+ Nat[1,yr.pup[j]]*pup_prop[j]
             logit(bp[j])<- -mean_day*beta1+beta1*daynr[j] -beta1*beta2*yr.pup[j]
             logit(lp[j])<- -mean_day*beta1+beta1*daynr[j] -beta1*beta2*yr.pup[j]-beta1*pup_dur
             pup_prop[j] <- bp[j]-lp[j]
           }
        
           Nat.cut[2,1] <- alpha.pup.cut * pup.uk[1]
           Nat.cut[3,1] <- 0
           Nat.cut[4,1] <- 0
           Nat.cut[5,1] <- 0
           Nat.cut[6,1] <- 0
           Nat.cut[7,1] <- ad_ini.cut
           Nat.cut[1,1] <- fec.ad.cut * Nat.cut[7,1] * frac.fem.cut
           
           for (t in 2:tmax){
             Nat.cut[1,t] <- fec.ad.cut * frac.fem.cut * (Nat.cut[7,t])
             Nat.cut[2,t] <- alpha.pup.cut * pup.uk[t] + surv.pup.cut * Nat.cut[1,t-1]
             Nat.cut[3,t] <- surv.ad.cut * Nat.cut[2,t-1]
             Nat.cut[4,t] <- surv.ad.cut * Nat.cut[3,t-1]
             Nat.cut[5,t] <- surv.ad.cut * Nat.cut[4,t-1]
             Nat.cut[6,t] <- surv.ad.cut * Nat.cut[5,t-1]
             Nat.cut[7,t] <- surv.ad.cut * (Nat.cut[6,t-1] + Nat.cut[7,t-1])
             }       
        
        # Molt-summer model
           for (k in 1:K.molt){
             molt[k] ~ dpois(N27[k])
             N27[k] <- corf.cut * surv.pup.cut * Nat.cut[1,yr.molt[k]] + pow(surv.ad.cut,106/365) * sum(Nat.cut[2:7,yr.molt[k]]) + alpha.molt * tot.uk[yr.molt[k]+1]
                 }
                 
           for (k in 1:K.smr){
             smr[k] ~ dpois(N17[k])
             N17[k] <- corf.cut * (surv.pup.cut * Nat.cut[1,yr.smr[k]] + pow(surv.ad.cut,207/365) * sum(Nat.cut[2:7,yr.smr[k]]) + alpha.smr * tot.uk[yr.smr[k]+1])
                }
                 
        # fixed priors (estimated from pup model)
          phi_a ~ dbeta(1.6,1.2)
          phi_p ~ dbeta(2.87, 1.78)
          ad_ini ~ dunif(1,100)
          alpha.pup ~  dunif(0.0001,0.02)
          fec_a ~ dbeta(2,1.5)
          beta2 ~ dunif(-2, 2)
          beta1 ~ dunif(0.0001,5)
          mean_day ~ dunif(10,80)
          pup_dur ~ dunif(26,57)
           
        # Fixed priors (correction factor and female:male ratio
          corf2 ~ dnorm(-1.0935,21.51307) # 1/(0.2156^2)     # informative prior
          beta_fem ~ dbeta(2,8)         # similar to female:male ratio = 1+gamma(2,0.1) as suggested in SCOS 2012, but with upper limit of 2:1
                
        # Priors for estimation (molt and summer model) 
          alpha.molt ~ dunif(0.00001,1)
          alpha.smr ~  dunif(0.00001,1)
                       
        # derived parameters
          frac.fem.sto <- (1+beta_fem)/(2+beta_fem)
          frac.fem.cut <- cut(frac.fem.sto)
          corf.sto <- exp(corf2)/(exp(corf2)+1)
          corf.cut <- cut(corf.sto) # cut makes sure that prior does not contribute to likelihood function
          surv.ad <- 0.8+0.2*phi_a
          fec.ad <-  0.6+0.4*fec_a
          surv.pup <- phi_p * surv.ad
          
        # Exclude parameters from estimation
          ad_ini.cut<-cut(ad_ini)
          alpha.pup.cut<-cut(alpha.pup)
          beta2.cut<-cut(beta2)
          beta1.cut<-cut(beta1)
          mean_day.cut<-cut(mean_day)
          pup_dur.cut<-cut(pup_dur)
          surv.ad.cut <- cut(surv.ad)
          fec.ad.cut <-  cut(fec.ad)
          surv.pup.cut <- cut(surv.pup)
                         
          }",
      fill=TRUE)
      sink()




#!  FIT PUP-ONLY MODEL

  # MCMC controls
  # Testing settings below should take ~2-3mins
    n.chains <- 3
    n.burn <- 10
    n.updates <- 100
    n.thin <- 2
    
  # Settings below should probably last about 10 x 5 = 2-3 hours
    #n.burn <- 100
    #n.updates <- 1000
    #n.thin <- 10
    
    
  # initialization for pup-only model
    inits.list.puponly <- function(){
         list(phi_a=0.6, phi_p=0.8, ad_ini=6, alpha.pup=0.010, fec_a=0.6, beta2=-1.3, beta1=0.3, mean_day=70, pup_dur=35, corf2=-1.1,alpha.molt=0.002, alpha.smr=0.002, beta_fem=0.20)  # blue
           }
     
  # parameters to monitor
    monitor.puponly         <- c("ad_ini","alpha.pup","alpha.molt","alpha.smr","beta1","beta2","mean_day","pup_dur","surv.ad","fec.ad","surv.pup","corf.cut","frac.fem.cut","beta_fem")
    
   
  # R2WinBUGS call; calling OpenBUGS requires BRugs package
  # WinBUGS chains are much better mixed than OpenBUGS chains
    old.wd<-getwd()
    setwd("models/Netherlands/")
    out <- bugs(data=data.list,
                 parameters.to.save=monitor.puponly,
                 model.file="Bugs_model_puponly.txt",
                 #inits=NULL,
                 inits=inits.list.puponly,       # reading of initial values not working
                 n.chains=n.chains,
                 n.thin=n.thin,
                 n.iter=(n.updates+n.burn)*n.thin,  # total number of MCMC iterations per chain before thinning and including burn-in
                 n.burnin=n.burn*n.thin,
                 debug=TRUE,
                 #over.relax=TRUE,          # default is FALSE; for grey seal data, this doesn't seem to make any difference
                 DIC=TRUE,
                 bugs.seed=sample.int(14,1),
                 bugs.directory="C:/Program_files/WinBUGS14",
                 program=c("WinBUGS"))
                 #bugs.directory="C:/Program Files (x86)/OpenBUGS/OpenBUGS322",
                 #program=c("OpenBUGS"))
 
  # Change working directory back
    setwd(old.wd)
    
  # Put results into data.frame          
    post <- data.frame(matrix(out$sims.array,ncol=dim(out$sims.array)[3],dimnames=list(NULL,names(out$sims.list))))
  
  # Save results
    save.image("models/Netherlands/results/Winbugs_model_results_with_data_2020_11_03.rdata")          
    #save.image("models/Netherlands/results/Winbugs_full_model_results_with_data_2020_11_03.rdata")          
    

################################################################################

#!  PLOTTING AND DIAGNOSTICS

################################################################################

#! LOAD DATA AND LIBRARIES
   
   # Set working directory
     setwd("C:/github/wk_WKSEALS-2020/")
    
   # Load data
     load("models/Netherlands/results/Winbugs_full_model_results_with_data_2020_11_03.rdata")        
   
   # Load libraries
     library(coda)
     library(fields)
     library(boot)

#! DIAGNOSTICS AND SUMMARY

   # diagnostics and summaries
   # summary table
    probs <- c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)
    cbind(par.num=1:ncol(post),mean=apply(post,2,mean),sd=apply(post,2,sd),t(apply(post,2,quantile,probs)))
     # or using coda for means, SD and SE; one advantage of this is the sophisticated SE estimation
    summary(mcmc(post))

  # trace
  # if n.chains>1, each chain is plotted in a different color
    layout <- c(ceiling(sqrt(ncol(post))),ceiling(ncol(post)/ceiling(sqrt(ncol(post)))))
    windows(11,8); par(mfrow=layout,mar=c(3,3,2,1))
    for (j in 1:ncol(post)) {
       plot(post[1:n.updates,j],type="n",ylim=c(min(post[,j]),max(post[,j])),main=names(post)[j],xlab="",ylab="")
       sapply(1:n.chains,function(i)lines(post[((i-1)*n.updates+1):(i*n.updates),j],col=i+1))}
  
  # the effective size is the number of independent samples which contain the same amount of information
    sapply(1:n.chains,function(i) effectiveSize(mcmc(post[((i-1)*n.updates+1):(i*n.updates),])))    

  # posterior histograms and smoothed densities
  # if n.chains>1, the density plot of each chain is plotted in a different color, with total in black
    windows(11,8); par(mfrow=layout,mar=c(3,5,2,1))
    for (j in 1:ncol(post)) {
       ymax <- max(sapply(1:n.chains,function(i)max(density(post[((i-1)*n.updates+1):(i*n.updates),j])$y)))
       hist(post[,j],prob=T,ylim=c(0,ymax),main=names(post)[j])
       sapply(1:n.chains,function(i)lines(density(post[((i-1)*n.updates+1):(i*n.updates),j]),col=i+1))
       lines(density(post[,j]),lwd=3)
    }

  # plot image of correlation matrix
  #crosscorr.plot.tg(mcmc(post))                # tg version has improved colors and axis labeling
    par(las=2)
    crosscorr.plot(mcmc(post))                # tg version has improved colors and axis labeling


#!  HISTOGRAM POSTERIOR AND PRIORS

  # Define titles
    titles<-list(ad_ini=expression(paste("Initial adults ", N["6+,t=0"])),
                 alpha.pup=expression(paste("Relative import ", alpha[pup])),
                 alpha.molt=expression(paste("Visitors molt ", alpha[molt])),
                 alpha.smr=expression(paste("Visitors summer ", alpha[summer])),
                 beta2=expression(paste("Annual shift birth day ", beta[2])),
                 beta1=expression(paste("variability birth day ", beta[1])),
                 mean_day=expression(paste("Mean birth day ", mu[t[birth]])),
                 pup_dur=expression(paste("Pup duration ", mu[t[duration]])),
                 surv.ad=expression(paste("Adult survival  ", phi[adult])),
                 fec.ad=expression(paste("Fecundity ", f["6+"])),
                 surv.pup=expression(paste("Pup survival  ", phi[pup])),
                 corf.cut=expression(paste("Haul-out fraction ", pi)),
                 frac.fem.cut=expression(paste("Fraction females ", rho[F])))
             
  # Define uniform priors           
     unif.prior <- c("ad_ini","alpha.pup","alpha.molt","alpha.smr","beta1","beta2","mean_day","pup_dur","alpha.molt","alpha.smr")

  # Define priors
    prior.functions<-list(ad_ini=function(x){ifelse(x<1|x>100,0,1)},
                          alpha.pup=function(x){ifelse((x<0.0001 | x>0.02),0,1)},
                          alpha.molt=function(x){ifelse((x<0.00001 | x>1),0,1)},
                          alpha.smr=function(x){ifelse((x<0.00001 | x>1),0,1)},
                          beta2=function(x){ifelse((x<(-2) | x>2),0,1)},
                          beta1=function(x){ifelse((x<(0.0001) | x>5),0,1)},
                          mean_day=function(x){ifelse((x<10 | x>80),0,1)},
                          pup_dur=function(x){ifelse((x<26 | x>57),0,1)},
                          surv.ad=function(x){ifelse(x<0.8,0,dbeta((x-0.8)/0.2,1.6,1.2))},
                          surv.pup=function(x,x2){dbeta(x2,1.6,1.2)*dbeta(x,1.6,1.2)},
                          fec.ad=function(x){ifelse(x<0.6,0,dbeta((x-0.6)/0.4,2,1.5))},
                          corf.cut=function(x){dnorm(logit(x),-1.0935,0.2156)},
                          beta_fem=function(x){dbeta(x,2,8)},
                          frac.fem.cut=function(x){x2=(2*x-1)/(1-x); dbeta(x2,2,8)})
                          
  # Define layout and plotting region
    layout <- c(4,4)
    #layout <- c(ceiling(sqrt(ncol(post))),ceiling(ncol(post)/ceiling(sqrt(ncol(post)))))
    windows(11,8); par(mfrow=layout,mar=c(3,4,3,0),family="serif",ps=16)

  # Define parameters to plot
    pars<-names(post)[is.element(names(post),c("deviance","beta_fem"))==FALSE]
    
  # Make plot
    for (j in pars) {          
       
     # Define xlimits  
       xlim <- range(density(post[,j])$x)
       if (names(post[j])=="surv.ad") xlim <- c(0.8,1.0)
       if (names(post[j])=="fec.ad") xlim <- c(0.6,1.0)
       if (names(post[j])=="surv.pup") xlim <- c(0,1)
       if (names(post[j])=="surv.pup.prior") xlim <- c(0,1)  
   
     # Construct and plot histogram
       the.hist<-hist(post[,j],prob=T,seq(xlim[1],xlim[2],length=20),plot=FALSE)
       if (is.element(j,c("ad_ini","beta1","surv.ad","frac.fem.cut"))) ylabel="Density"
       else (ylabel="")
       the.hist<-hist(post[,j],prob=T,seq(xlim[1],xlim[2],length=20),ylab=ylabel,xlab="",main=titles[[j]],xlim=xlim,ylim=c(0,1.15*max(the.hist$density)))
       sum.hist<-sum(the.hist$density)
   
     # Define x values used for plotting prior
       x <- seq(xlim[1],xlim[2],length=100)
 
     # Plot priors
       if (names(post[j])!="surv.pup") {
         standardize<-sum(prior.functions[[j]](x)/100)/(sum.hist/19) 
         if (is.element(j,c("frac.fem.cut","corf.cut"))) the.col<-"red"
         else (the.col="grey")
         lines(x,y=prior.functions[[j]](x)/standardize,lwd=3,col=the.col)} 
      
       else { a<-rbeta(100000,1.6,1.2)
              a<-0.8+0.2*a
              b<-rbeta(100000,2.87, 1.78)
              the.hist2<-hist(a*b,seq(xlim[1],xlim[2],length=20),plot=FALSE)
              lines(the.hist2$mids,the.hist2$density,lwd=3,col="grey")
              #x2 <- x*(0.8+0.2*x); 
              #standardize<-sum(prior.functions[[j]](x,x2)/100)/(sum.hist/19) 
              #lines(x2,y=prior.functions[[j]](x,x2)/standardize,lwd=3,col="gray")
              }
      
      # Plot posterior smoothed density
        lines(density(post[,j]),lwd=2)
  
  }
    
     
#!  CORRELATION DENSITY PLOT SURVIVAL PUP, ADULTS AND IMPORT

  # bivariate density plots of parameters
    dev.off()
    ncut <- 30                            # number of bins for variables
    ndx.par <- 1:ncol(post)               # all parameters
    #ndx.par<-which(is.element(names(post),c("alpha.pup","alpha.smr","alpha.molt","surv.ad","surv.pup")))
    ndx.par<-which(is.element(names(post),c("alpha.pup","surv.ad","surv.pup")))
    ndx <- combn(ndx.par,2)               # matrix of all combinations of indices to be plotted
    par(mfrow=c(floor(sqrt(ncol(ndx))),ceiling(ncol(ndx)/floor(sqrt(ncol(ndx))))),mar=c(5,5,2,1),family="serif",ps=16)
    for (k in 1:ncol(ndx)) {
      x <- seq(min(post[,ndx[1,k]]),max(post[,ndx[1,k]]),length=ncut)
      y <- seq(min(post[,ndx[2,k]]),max(post[,ndx[2,k]]),length=ncut)
      z <- as.matrix(table(cut(post[,ndx[1,k]],ncut),cut(post[,ndx[2,k]],ncut)))
      Pcor <- round(cor(post[,ndx[1,k]],post[,ndx[2,k]]),3)
      #contour(x,y,z,drawlabels=F,xlab=names(post)[ndx[1,k]],ylab=names(post)[ndx[2,k]],cex.lab=1.5,main=paste("Pearson cor =",Pcor)) # contour lines
      #image(x,y,z,col=tim.colors(64),xlab=names(post)[ndx[1,k]],ylab=names(post)[ndx[2,k]],main=paste("Pearson cor =",Pcor))  # color plot
      image(x,y,z=image.smooth(z)$z,col=tim.colors(64),xlab=titles[names(post)[ndx[1,k]]][[1]],ylab=titles[names(post)[ndx[2,k]]][[1]],main=paste("Pearson cor =",Pcor))  # smoothed color plot
      box()
    }


################################################################################

#!  DEFINE FUNCTION TO CALCULATE POPULATION MATRIX AND PLOT TREND OVER YEARS

################################################################################
   
# Define function to calculate population matrix     

  # Define values
    pup.uk<-data.list$pup.uk
    tot.uk<-data.list$tot.uk
    pup.nl<-data.list$pup
      
  # Define maximum number of years
    tmax=27      # was 26 + 1 year for plotting
   
  # Define function 
    Nat.function<-function(){
     
     Nat<-matrix(0,nrow=7,ncol=tmax)
      
       Nat[2,1] <- alpha.pup * pup.uk[1]
       Nat[3,1] <- 0
       Nat[4,1] <- 0
       Nat[5,1] <- 0
       Nat[6,1] <- 0
       Nat[7,1] <- ad_ini
       Nat[1,1] <- fec.ad * Nat[7,1] * frac.fem
     
       for (t in 2:tmax){
         Nat[2,t] <- alpha.pup * pup.uk[t] + surv.pup * Nat[1,t-1]
         Nat[3,t] <- surv.ad * Nat[2,t-1]
         Nat[4,t] <- surv.ad * Nat[3,t-1]
         Nat[5,t] <- surv.ad * Nat[4,t-1]
         Nat[6,t] <- surv.ad * Nat[5,t-1]
         Nat[7,t] <- surv.ad * (Nat[6,t-1] + Nat[7,t-1])
         Nat[1,t] <- fec.ad * frac.fem * (Nat[7,t])}
         
     return(Nat)}    
    
# Calculate population sizes
  
  tmax=length(c(1985:2013))
  
  # Create empty vector 
    all.total<-all.pups<-all.moult<-all.summer<- all.year1<-all.UKimport<-all.smr.visitors<-all.molt.visitors<-all.pup.ratio<-rep(NA,tmax)
   
  # Loop 
    for (i in 1:nrow(post)){
      
      # Get parameter values
         ad_ini     = post[i,"ad_ini"]
         alpha.pup  = post[i,"alpha.pup"]
         alpha.molt = post[i,"alpha.molt"]
         alpha.smr  = post[i,"alpha.smr"]
         beta2      = post[i,"beta2"]
         beta1      = post[i,"beta1"]
         mean_day   = post[i,"mean_day"]
         pup_dur    = post[i,"pup_dur"]
         surv.ad    = post[i,"surv.ad"]
         fec.ad     = post[i,"fec.ad"]
         surv.pup   = post[i,"surv.pup"]
         corf       = post[i,"corf.cut"]
         frac.fem   = post[i,"frac.fem.cut"]
 
      # Caclualte population matrix  
        Nat<-Nat.function()
      
      # Extract total number of individuals
        total<-colSums(Nat)
        
      # Extract pups  
        pups<-Nat[1,]
        
      # Extract moult  
        #moult <- corf * Nat[1,] + colSums(Nat[2:7,]) + alpha.molt * pup.uk/0.16
        #moult <- corf * Nat[1,] + colSums(Nat[2:7,])
        moult <- surv.pup*corf * Nat[1,] + (surv.ad^(106/365))*colSums(Nat[2:7,])
       
      # Extract summer       
        #summer <- corf * (colSums(Nat[1:7,]) + alpha.smr * pup.uk/0.16)
        #summer <- corf * (colSums(Nat[1:7,]))
        summer <- surv.pup*corf * Nat[1,] + (surv.ad^(207/365))*corf*colSums(Nat[2:7,])
        
      # Extract UK relative import of yearlings
        UKimport<-    alpha.pup * pup.uk[1:29]
        
      # Extract Summer temporary visitors
        smr.visitors<-alpha.smr * tot.uk[1:29]
        
      # Extract Molt temporary visitors  
        molt.visitors<-alpha.molt * tot.uk[1:29]
        
      # Calculate pup ratio
        pup.ratio<-pups/total 
        
      # Extract 1 year olds 
        year1 <- Nat[2,]  
        
      # Put total and pups in total table  
        all.total<-rbind(all.total,total)
        all.pups<-rbind(all.pups,pups)
        all.moult<-rbind(all.moult,moult) 
        all.summer<-rbind(all.summer,summer)
        all.UKimport<-rbind(all.UKimport,UKimport)
        all.smr.visitors<-rbind(all.smr.visitors,smr.visitors)
        all.molt.visitors<-rbind(all.molt.visitors,molt.visitors)
        all.pup.ratio<-rbind(all.pup.ratio,pup.ratio)
        all.year1<-rbind(all.year1,year1)
    
        }
        
      # Remove first row
        all.total<-all.total[-1,]
        all.pups<-all.pups[-1,]
        all.summer<-all.summer[-1,]
        all.moult<-all.moult[-1,]       
        all.UKimport<-all.UKimport[-1,]
        all.smr.visitors<-all.smr.visitors[-1,]
        all.molt.visitors<-all.molt.visitors[-1,]
        all.year1<-all.year1[-1,]
    
      # Put into list
        var.list<-list(total=all.total, pups=all.pups, summer=all.summer, moult=all.moult,UKimport=all.UKimport, smr.visitors=all.smr.visitors, molt.visitors=all.molt.visitors, year1=all.year1)
        
      # Put everything in a table
          for (i in names(var.list))
          {
          the.mean<-apply(var.list[i][[1]],2,mean)
          the.quant<-t(apply(var.list[i][[1]],2,function(x){quantile(x,c(0.025,0.25,0.75,0.975))}))       
          pop.table<-cbind(data.frame(name=i,years=1985:2013,mean=the.mean),the.quant)   
          if (i=="total") all.pop.table<-pop.table
          else all.pop.table<-rbind(all.pop.table,pop.table)
          }

       # Show relative contribution of import
         all.pop.table$mean[all.pop.table$name=="UKimport"]/all.pop.table$mean[all.pop.table$name=="year1"]

       # Write table
         #write.csv(all.pop.table,"models/Netherlands/results/population_variables_2020_11_03.csv")


################################################################################                    

#!  MAKE PLOT OF POPULATION TRENDS

################################################################################
  
    dev.off()
  
  # Define add 
    mm<-356/365
    years<-1985:2013
    molt.smr.years<-1986:2014 # molt and summer is assumed to take place after pups are born
    
    # Set par
      windows(width=11/2.5,height=16/2.5,rescale="fixed")
       par(family="serif",ps=14,mfrow=c(3,1),mai=c(0.4,0.55,0.05,0.05),omi=c(0.2,0,0.1,0))
        
    # Make plot of total population size and pups
      #all.mean<-apply(all.total,2,function(x){quantile(x,0.5)})
      all.mean<-all.pop.table$mean[all.pop.table$name=="total"]
      all.lower<-all.pop.table[all.pop.table$name=="total","2.5%"]
      all.upper<-all.pop.table[all.pop.table$name=="total","97.5%"]
      plot(years+mm,all.mean,xlim=c(1985,2013.8),ylim=c(0,4000),type="l",yaxs="i", xaxs="i",ylab="Number of animals", xlab="Year")
      polygon(c(years,rev(years),1985)+mm,c(all.lower,rev(all.upper),all.lower[1]),col="lightgrey",border="darkgrey")
      lines(years+mm,all.mean,lwd=2,lty=2)
      polygon(x=c(1986.25,1986.25,1987.75,1987.75,1986.25),y=c(3475,3675,3675,3475,3475),col="lightgrey",border="lightgrey")
      polygon(x=c(1986.25,1986.25,1987.75,1987.75,1986.25),y=c(3150,3350,3350,3150,3150),col="Tomato",border="Tomato")
      legend(x=1986,y=3900,legend=c("Estimated Dutch population size (all ages)","Estimated pup production","Pup counts"), lwd=c(2,2,1),lty=c(2,1,-1),col=c("black","darkred","black"),pch=c(NA,NA,17),bty="n") 
      #text(2013,3900,"a",cex=1.3)
      #pup.mean<-apply(all.pups,2,function(x){quantile(x,0.5)})
      pup.mean<-all.pop.table$mean[all.pop.table$name=="pups"]
      pup.lower<-all.pop.table[all.pop.table$name=="pups","2.5%"]
      pup.upper<-all.pop.table[all.pop.table$name=="pups","97.5%"]
      polygon(c(years,rev(years),1985)+mm,c(pup.lower,rev(pup.upper),pup.lower[1]),col="Tomato",border="Tomato")  # col=rgb(0,0,1,0.5),border=rgb(0,0,1,0.9)
      #polygon(c(molt.smr.years,rev(molt.smr.years),1986)+mm,c(summer.lower,rev(summer.upper),summer.lower[1]),col=rgb(1,0,0,1),border=rgb(1,0,0,0.9)) # col=rgb(1,0,0,0.5)
      lines(years+mm,pup.mean,lwd=2,col="darkred")
      points(pup.nl$seasonyear+mm,pup.nl$SumOfHG_pup,pch=17)
      text(2013,3900,"a",cex=1.3)
      box(lwd=2)
      
  # Make plot of moult
    plot(years+mm,all.mean,xlim=c(1985,2013.8),ylim=c(0,4000),type="l",yaxs="i",xaxs="i",ylab="Number of animals", xlab="Year",col=0)
    mm<-106/365
    #moult.mean<-apply(all.moult,2,function(x){quantile(x,0.5)})
    moult.mean<-all.pop.table$mean[all.pop.table$name=="moult"]
    moult.lower<-all.pop.table[all.pop.table$name=="moult","2.5%"]
    moult.upper<-all.pop.table[all.pop.table$name=="moult","97.5%"]
    #polygon(c(molt.smr.years,rev(molt.smr.years),1986)+mm,c(moult.lower,rev(moult.upper),moult.lower[1]),col=rgb(0,1,0,0.5),border=rgb(0,1,0,0.9))
    polygon(c(molt.smr.years,rev(molt.smr.years),1986)+mm,c(moult.lower,rev(moult.upper),moult.lower[1]),col="LightGreen",border="LightGreen")
    points(molt.smr.years+mm,moult.mean,type="l",lwd=2,col="darkgreen")
    points(Data$year[Data$Type=="Moult"]+mm,Data$TOT[Data$Type=="Moult"],pch=15)   
    points(counts$seasonyear[counts$Type=="Moult"]+mm,counts$SumOfHG_num[counts$Type=="Moult"],pch=1)
    mm=356/365
    lines(years+mm,all.mean,type="l",lty=2,lwd=2)
    #lines(years+mm,pup.mean,lty=2)
    #legend(x=1986,y=c(3900),legend=c("Expected molt count excluding visitors","Molt counts used for model fitting","Remaining molt counts","Estimated Dutch population size (all ages)","Estimated pup production"), lwd=c(2,1,1,1,1),col=c("darkgreen","black","black","black","black"),lty=c(1,-1,-1,3,2),pch=c(NA,15,1,NA,NA),bty="n") 
    polygon(x=c(1986.25,1986.25,1987.75,1987.75,1986.25),y=c(3475,3675,3675,3475,3475),col="LightGreen",border="LightGreen")
    legend(x=1986,y=c(3900),legend=c("Expected molt count excluding visitors","Molt counts used for model fitting","Remaining molt counts"), lwd=c(2,1,1),col=c("darkgreen","black","black"),lty=c(1,-1,-1),pch=c(NA,15,1),bty="n") 
    text(2013,3900,"b",cex=1.3)
    box(lwd=2)
  
    
  # Make plot for summer
    plot(years+mm,all.mean,xlim=c(1985,2013.8),ylim=c(0,4000),type="l",yaxs="i",xaxs="i",ylab="Number of animals", xlab="Year",col=0)
    mm<-207/365  
    #summer.mean<-apply(all.summer,2,function(x){quantile(x,0.5)})
    summer.mean<-all.pop.table$mean[all.pop.table$name=="summer"]
    summer.lower<-all.pop.table[all.pop.table$name=="summer","2.5%"]
    summer.upper<-all.pop.table[all.pop.table$name=="summer","97.5%"]
    all.summer.mean<- all.pop.table$mean[all.pop.table$name=="summer"] + all.pop.table$mean[all.pop.table$name=="smr.visitors"]
    polygon(c(molt.smr.years,rev(molt.smr.years),1986)+mm,c(summer.lower,rev(summer.upper),summer.lower[1]),col="LightBlue",border="LightBlue") # col=rgb(1,0,0,0.5)
    points(molt.smr.years+mm,summer.mean,type="l",lwd=2,col="darkblue")
    points(counts$seasonyear[counts$Type=="Other"]+mm,counts$SumOfHG_num[counts$Type=="Other"],pch=16)
    #lines(years+mm,all.summer.mean,lwd=4)
    #polygon(c(molt.smr.years,rev(molt.smr.years),1986)+mm,c(summer.mean,rev(all.summer.mean),summer.mean[1]),density=20,angle=-45) # col=rgb(1,0,0,0.5)
    mm=356/365
    lines(years+mm,all.mean,type="l",lty=2,lwd=2)
    mm<-356/365
    #lines(years+mm,pup.mean,lty=2)
    #legend(x=1986,y=c(3900),legend=c("Expected summer count excluding visitors","Summer counts","Estimated Dutch population size (all ages)","Estimated pup production"), lwd=c(2,1,1,1),col=c("darkred","black","black","black"),pch=c(NA,16,NA,NA),lty=c(1,-1,3,2),bty="n") 
    #legend(x=1986,y=c(3900),legend=c("Expected summer count excluding visitors","Summer counts"), lwd=c(2,1),col=c("darkred","black"),pch=c(NA,16),lty=c(1,-1),bty="n") 
    polygon(x=c(1986.25,1986.25,1987.75,1987.75,1986.25),y=c(3475,3675,3675,3475,3475),col="LightBlue",border="LightBlue")
    legend(x=1986,y=c(3900),legend=c("Expected summer count excluding visitors","Summer counts"), lwd=c(2,1),col=c("darkblue","black"),pch=c(NA,16),lty=c(1,-1),bty="n") 
    text(2013,3900,"c",cex=1.3)
    mtext("Year",side=1,line=2.5,cex=0.7)
    box(lwd=2)
  
################################################################################                    

#!  MAKE PLOT OF PUPS AND 1 YEAR OLDS

################################################################################

   # Set par
      windows(width=11/2.5,height=16/2.5,rescale="fixed")
       par(family="serif",ps=14,mfrow=c(3,1),mai=c(0.4,0.55,0.05,0.05),omi=c(0.2,0,0.1,0))
        
    # Make plot of total population size and pups
      #all.mean<-apply(all.total,2,function(x){quantile(x,0.5)})
      year1.mean<-all.pop.table$mean[all.pop.table$name=="year1"]
      year1.lower<-all.pop.table[all.pop.table$name=="year1","2.5%"]
      year1.upper<-all.pop.table[all.pop.table$name=="year1","97.5%"]
      plot(years+mm,year1.mean,xlim=c(1985,2013.8),ylim=c(0,600),type="l",yaxs="i", xaxs="i",ylab="Number of animals", xlab="Year")
      polygon(c(years,rev(years),1985)+mm,c(year1.lower,rev(year1.upper),year1.lower[1]),col="Gold",border="Gold")
      lines(years+mm,year1.mean,lwd=2,lty=1,col="black")  # Chocolate
      ukimport.mean<-all.pop.table$mean[all.pop.table$name=="UKimport"]
      ukimport.lower<-all.pop.table[all.pop.table$name=="UKimport","2.5%"]
      ukimport.upper<-all.pop.table[all.pop.table$name=="UKimport","97.5%"]
      polygon(c(years,rev(years),1985)+mm,c(ukimport.lower,rev(ukimport.upper),ukimport.lower[1]),col="Orange",border="Orange")
      lines(years+mm,ukimport.mean,lwd=2,lty=1,col="black")  # "OrangeRed"
      polygon(x=c(1986.25,1986.25,1987.75,1987.75,1986.25),y=c(510,540,540,510,510),col="Gold",border="Gold")
      polygon(x=c(1986.25,1986.25,1987.75,1987.75,1986.25),y=c(460,490,490,460,460),col="Orange",border="Orange")
      legend(x=1986,y=575,legend=c("Estimated total 1-year olds","Estimated UK import 1-year olds"), lwd=c(2,2),lty=c(1,1),col=c("black","black"),pch=c(NA,NA),bty="n") 
      text(2013,575,"a",cex=1.3)
      box(lwd=2)
            
    # Make plot of total population size and pups
      opvang<-read.csv("data/opvang.csv")
      pup.mean<-all.pop.table$mean[all.pop.table$name=="pups"]
      pup.lower<-all.pop.table[all.pop.table$name=="pups","2.5%"]
      pup.upper<-all.pop.table[all.pop.table$name=="pups","97.5%"]
      plot(years+mm,pup.mean,xlim=c(1985,2013.8),ylim=c(0,600),type="l",yaxs="i", xaxs="i",ylab="Number of animals", xlab="Year")
      polygon(c(years,rev(years),1985)+mm,c(pup.lower,rev(pup.upper),pup.lower[1]),col="Tomato",border="Tomato")  # col=rgb(0,0,1,0.5),border=rgb(0,0,1,0.9)
      #polygon(c(molt.smr.years,rev(molt.smr.years),1986)+mm,c(summer.lower,rev(summer.upper),summer.lower[1]),col=rgb(1,0,0,1),border=rgb(1,0,0,0.9)) # col=rgb(1,0,0,0.5)
      lines(years+mm,pup.mean,lwd=2,col="darkred")
      points(opvang$Year+mm,opvang$Total,pch=1,lwd=2,cex=2)
      points(pup.nl$seasonyear+mm,pup.nl$SumOfHG_pup,pch=17)     
      polygon(x=c(1986.25,1986.25,1987.75,1987.75,1986.25),y=c(510,540,540,510,510),col="Tomato",border="Tomato")     
      legend(x=1986,y=575,legend=c("Estimated pup production","Pup counts","Grey seals in rehabilitation centres"), lwd=c(2,1,2),lty=c(1,-1,-1),col=c("darkred","black","black"),pch=c(NA,17,1),pt.cex=c(1,1,2),bty="n") 
      text(2013,575,"b",cex=1.3)
      box(lwd=2)
    

################################################################################                    

#!  MAKE PANEL PLOT OF PUP NUMBERS

################################################################################
    par(family="Times")
    
  # Define parameters
    # Get parameter values
           ad_ini     = median(post[,"ad_ini"])
           alpha.pup  = median(post[,"alpha.pup"])
           alpha.molt = median(post[,"alpha.molt"])
           alpha.smr  = median(post[,"alpha.smr"])
           beta2      = median(post[,"beta2"])
           beta1      = median(post[,"beta1"])
           mean_day   = median(post[,"mean_day"])
           pup_dur    = median(post[,"pup_dur"])
           surv.ad    = median(post[,"surv.ad"])
           fec.ad     = median(post[,"fec.ad"])
           surv.pup   = median(post[,"surv.pup"])
           corf       = median(post[,"corf"])
           frac.fem   = median(post[,"frac.fem"])
  
  
  # Create prediction table
    prgr<-expand.grid(daynr=1:120,years=c(1985:2013))
    prgr$pred_obs<-rep(NA,nrow(prgr))
    prgr$dat="model.pred"  
    
  # Create similar table for model data
    pups<-pup.nl[,c("daynr","seasonyear","SumOfHG_pup")]
    names(pups)<-c("daynr","years","pred_obs")
    pups$dat<-rep("model.data",nrow(pups))
    
  # Table for omited data
    allp<-counts[counts$Type=="Pup",]
    allp<-allp[,c("daynr","seasonyear","SumOfHG_pup")]
    names(allp)<-c("daynr","years","pred_obs")
    include<-is.element(paste(allp$daynr,allp$years),paste(pups$daynr,pups$year))==FALSE
    allp<-allp[include,]
    allp$dat<-rep("omit",nrow(allp))
  
  # Combine
    prgr<-rbind(prgr,pups,allp)
    rm(pups,allp)
      
  # Some processing
   
     # Change year id    
       prgr$yearid<-prgr$years-1985+1 
       
     # Define ordering of factor
       prgr$dat<-factor(prgr$dat,levels=c("model.data","omit","model.pred"))
      
     # Only select data from 1986 and beyond
       prgr<-prgr[prgr$years>=1989,]
     
     # Create date time variable
       prgr$year_fact<-factor(prgr$years,levels=c(2009:2013,2004:2008,1999:2003,1994:1998, 1989:1993))
    
  # Calculate number of Pups for each day and year
    bp<- inv.logit(-mean_day*beta1+beta1*prgr$daynr -beta1*beta2*prgr$yearid)
    lp<- inv.logit(-mean_day*beta1+beta1*prgr$daynr -beta1*beta2*prgr$yearid-beta1*pup_dur)
    pup_prop<-bp-lp
    eta<- pup.mean[prgr$yearid]*pup_prop
    prgr$pred_obs[prgr$dat=="model.pred"]<-eta[prgr$dat=="model.pred"]
   
  # Make plot    
    library(lattice)
    xyplot(pred_obs~daynr|factor(year_fact),scales=list(x=list(at=c(31,62,93),labels=c("1-Dec","1-Jan","1-Feb"),rot=60),fontfamily="serif",cex=1),fontfamily="serif",ps=14, data=prgr,type=c("p","p","l"), pch=c(17,1,0), col="black", groups=dat,distribute.type=TRUE,xlab=list("Date",fontfamily="serif",cex=1.15), ylab=list("Number of pups",fontfamily="serif",cex=1.15),
              strip = strip.custom(bg = grey(0.9)),par.strip.text=list(fontfamily="serif",ps=18),main="")
  

################################################################################                    

#!  SUMMARY STATISTICS

################################################################################

#!  CONSTRUCT SUMMARY TABLE OF PARAMETERS 
    
  # Construct summary table    
    sum.table<-cbind(par.num=1:ncol(post),mean=apply(post,2,mean),sd=apply(post,2,sd),t(apply(post,2,quantile,c(0.025,0.5,0.975))))
    sum.table<-sum.table[c("surv.ad","surv.pup","fec.ad","ad_ini","alpha.pup","beta1","beta2","mean_day","pup_dur","alpha.smr","alpha.molt","corf.cut","frac.fem.cut"),]   
    
  # Show table  
    sum.table<-round(sum.table,4)
    sum.table
    
  # Write table  
    write.csv(sum.table,"results/summary_table_all_parameters_2014_06_23.csv")
    
    
#!  CALCULATE MEANS AND SDs OF PRIORS

  # Mean and sd for adult survival
    alpha.=1.6
    beta.=1.2
    #rand<-0.8+0.2*rbeta(1e+8,1.6,1.2)
    mean(rand)
    sd(rand)
    0.8+0.2*(alpha./(alpha.+beta.))
    0.2*sqrt(alpha.*beta./((alpha.+beta.)^2 * (alpha.+beta.+1)))
 
  # For Pup survival
    alpha.=2.87
    beta.=1.78
    (alpha./(alpha.+beta.))
    sqrt(alpha.*beta./((alpha.+beta.)^2 * (alpha.+beta.+1)))
 
  # For Fecundity
    alpha.=2
    beta.=1.5
    0.6+0.4*(alpha./(alpha.+beta.))
    0.4*sqrt(alpha.*beta./((alpha.+beta.)^2 * (alpha.+beta.+1)))
 


#!  SUMMARY STATISTICS OF FORWARD SHIFT
    
  # Calculate date of birth in 1985
    as.POSIXct(strptime("1985-11-01 00:00",format="%Y-%m-%d %H:%M"))+mean(post$mean_day)*60*60*24
   
  # Calculate peak in pup numbers n 1985
    as.POSIXct(strptime("1985-11-01 00:00",format="%Y-%m-%d %H:%M"))+mean(post$mean_day)*60*60*24 + mean(post$pup_dur)*60*60*24/2
    
  # Calculate date of birth in 2013
    as.POSIXct(strptime("2013-11-01 00:00",format="%Y-%m-%d %H:%M"))+mean(post$mean_day)*60*60*24 + mean(post$beta2)*(2013-1985)*60*60*24
   
  # Calculate peak in pup numbers in 2013
    as.POSIXct(strptime("2013-11-01 00:00",format="%Y-%m-%d %H:%M"))+mean(post$mean_day)*60*60*24 + + mean(post$pup_dur)*60*60*24/2 + mean(post$beta2)*(2013-1985)*60*60*24
   
  # Forward shift from 1985 to 2013
    2013-1985
    mean(post$beta2)*(2013-1985)
  
  

