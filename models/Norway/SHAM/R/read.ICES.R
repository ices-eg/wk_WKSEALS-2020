read.ICES <- function(filename="Data/_10193_10193.xml") {
  require(XML)
  require(plyr)
  tst <- xmlParse(filename)
  tst <- xmlToList(tst)
  tst.met <- unlist(tst[-grep('FishData', names(tst))])
  tst <- tst[grep('FishData', names(tst))]
  
  tst <- lapply(tst, function(x) as.data.frame(x, stringsAsFactors=F))
  
  tst.df <- tst[[1]]
  for(i in 2:length(tst)) {
    tst.df <- rbind.fill(tst.df, tst[[i]])
  }
  
  for(i in 1:length(tst.df)) tst.df[,i] <- as.numeric(tst.df[,i])
  
  list(Meta=tst.met, Data=tst.df)
}


