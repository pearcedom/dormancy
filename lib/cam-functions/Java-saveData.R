# Save the results (Aest, Sest) as CSV formatted files in the directory "results".
# The files will be named by the current time.
# Author:  Fan Meng (mengfan@vt.edu)
# Version: 4.26.2012
# This software is under license BSD 3.

saveData <- function(Aest, Sest) {
  now      <- format(Sys.time(),"%Y-%m-%d-%H-%M")
  AestName <- paste("../results/", now, "-Aest.csv", sep="")
  SestName <- paste("../results/", now, "-Sest.csv", sep="")
  write.csv(Aest, file=AestName, row.names=FALSE)
  write.csv(Sest, file=SestName, row.names=FALSE)
}