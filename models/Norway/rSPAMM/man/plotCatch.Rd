% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outputFunctions.R
\name{plotCatch}
\alias{plotCatch}
\title{Plot the reported catch data}
\usage{
plotCatch(
  catch = cdata,
  width = 9,
  height = 7,
  position = "dodge",
  grDev = FALSE
)
}
\arguments{
\item{catch}{Reported catcj data}

\item{width}{Figure width}

\item{height}{Figure height}

\item{position}{Position of bars: If Pup catch and 1+ catch next to each other use position = "dodge" (default). On top of each other use position = "stack"}

\item{grDev}{Logical parameter to decide wether to open a OS independent graphical window}
}
\value{
plot Returns a plot of predicted population size for different population components
}
\description{
Plot the reported catch data.
}
\examples{
plotCatch(data$Cdata)
}
\keyword{model}
\keyword{population}
