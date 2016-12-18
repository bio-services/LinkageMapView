# edit file
# make sure all positions in the file are numeric

editlgdf  <- function (df) {
    # make sure all positions are numeric
    if (!is.numeric(df[,2]))
    { for (i in 1:nrow(df)) {
      message(c("Some positions are not numeric - in linkage group: ",df[i,1]))
      if (!is.numeric(df[i,2])) {
        message(c("I see ",  df[i,2] , " in row " , i))

      }
    }

      stop("Some positions are not numeric") }
}
