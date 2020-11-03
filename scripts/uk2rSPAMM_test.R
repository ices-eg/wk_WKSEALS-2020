###############################################################
## Testing conversion of UK pup production data 
## for input into rSPAMM. 
## It runs, but is naive and with lots of 'dummy' data.
## M. Biuw 2020-11-04 
###############################################################


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Read UK pup production data:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

uk <- read.table('../../../data/UK/PupData_1984_2018.txt', header=T)
uk$North.Sea[which(uk$North.Sea==-1)] <- NA
uk$Inner.Hebrides[which(uk$Inner.Hebrides==-1)] <- NA
uk$Outer.Hebrides[which(uk$Outer.Hebrides==-1)] <- NA
uk$Orkney[which(uk$Orkney==-1)] <- NA



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create rSPAMM data list 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

uk.dat <- list()
uk.dat$Amax <- 20
uk.dat$Cdata <- min(uk$Year, na.rm=T):max(uk$Year, na.rm=T)
uk.dat$Cdata <- cbind(V1=uk.dat$Cdata, V2=rep(0, length(uk.dat$Cdata)), V3=rep(0, length(uk.dat$Cdata)))
uk.dat$Nc <- nrow(uk.dat$Cdata)

## Use mean CV from harpeast for pup counts
## Either sum of all regions:
uk.dat$pupProductionData <- cbind(uk$Year,
                                  apply(uk[,-1], 1, sum),
                                  rep(mean(data$pupProductionData[,3]), nrow(uk)))

## Or only single region (e.g. Outer Hebrides)
# uk.dat$pupProductionData <- cbind(uk$Year, 
#                                   uk$Outer.Hebrides, 
#                                   rep(mean(data$pupProductionData[,3]), nrow(uk)))

uk.dat$pupProductionData <- uk.dat$pupProductionData[which(!is.na(uk.dat$pupProductionData[,2])),]
uk.dat$Np <- nrow(uk.dat$pupProductionData)

## Use posterior "mean" fecundity from Thomas et al. (2018)
uk.dat$Ftmp <- rep(0.9, nrow(uk.dat$Cdata))

## Use mean ogives from harpeast:

uk.dat$Pmat <- matrix(apply(data$Pmat, 2, mean), nrow=1)

for(i in 2:nrow(uk.dat$Cdata)) uk.dat$Pmat <- rbind(uk.dat$Pmat, uk.dat$Pmat[1,])

uk.dat$Npred=0

## Assume starting population size of 4 * first estimated pup production
uk.dat$priors <- cbind(V1=c(uk.dat$pupProductionData[1,2]*4,
                                    0.05, 0.52),
                                  V2=c(0.2*uk.dat$pupProductionData[1,2]*4,
                                    0.05, 0.05))

## Or wild priors:
# uk.dat$priors <- cbind(V1=c(300000,
#                             0.2, 0.3),
#                        V2=c(0.2*300000,
#                             0.5, 0.5))


uk.dat$Npriors <- as.integer(3)
uk.dat$CQuota <- rep(0, 2)
uk.dat$fecundity <- data.frame(V1=uk.dat$Cdata[,1],
                               V2=uk.dat$Ftmp,
                               V3=rep(0.06, nrow(uk.dat$Cdata)))
uk.dat$Pper <- data.frame(Pstart=min(uk.dat$Cdata[,1]), Pstop=max(uk.dat$Cdata[,1]))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Set initial parameter values: 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Initial values taken from mean priors as a first try:

uk.par <- list(logK=log(as.numeric(uk.dat$priors[1,1])),
               Mtilde=logit(as.numeric(uk.dat$priors[2,1])), 
               M0tilde=logit(as.numeric(uk.dat$priors[3,1])))

ukobj <- run.model(data = uk.dat, par = uk.par)

res <- model.results(data = uk.dat,optobject = ukobj)

partab <- par.table(results=res, dat=uk.dat) 

plotRes(res, uk.dat, plotNlims=F, plotProjections = F)

