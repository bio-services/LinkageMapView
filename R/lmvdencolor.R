#' LinkageMapView density color function
#'
#' lmvdencolor is a helper function which you can use to create a data frame
#' of colors to be used as the sectcoldf input parameter on the lmv.linkage.plot command.
#' The colors will be used to color the linkage group based on the density
#' of position/marker.This function is called with default values when the
#' denmap = TRUE parameter is specified for lmv.linkage.plot and no sectcoldf parameter
#' is found.
#'
#' @param df Required, a data frame with the first two columns in this order:
#'      \enumerate{
#'        \item Linkage group name.
#'        \item Position - must be in numerical order ascending within
#'              linkage group name.  If the maximum position in any linkage
#'              group is < 1000, the density will be calculated
#'              for each position.  Otherwise the number of positions included
#'              for each density calculation will be:
#'              ceiling(maximum position of an linkage group/1000)
#'      }
#'
#' @param wsize Optional, default = 30.  # of postions in the sliding window
#'              for calculating postitions/marker. If the maximum position in
#'              any linkage group is >= 1000, the default sliding window size
#'              will be adjusted by the same ratio as the number of positions
#'              included for each density calculate.
#'
#' @param bias Optional, default = 5. a positive number. Higher values give
#'        more widely spaced colors at the high end.
#' @seealso \code{\link[grDevices]{colorRamp}}
#'
#' @param colorin Optional, a vector of colors to use where the first value is
#'        the color for the lowest density and the last value is the color
#'        for the highest density.  Default is:
#'        rev(colorRampPalette(RColorBrewer::brewer.pal(8, "Spectral"))(25))
#'
#' @return a data frame that can be used as sectcoldf input on the lmv.linkage.plot
#'        function to color the chromosome for a density map.
#'
#' @seealso \code{\link{lmv.linkage.plot}}
#'
#' @export
#'
#' @examples
#'
#' # add a column to a linkage group data frame to specify colors for
#' # line segments in lmv.linkage.plot using default colors from RColorBrewer
#' # Spectral palette.  Then just plot the returned colors out to see how
#' # they look.
#'
#' data(oat)
#'
#' sectcoldf <- lmvdencolor(oat)
#'
#' # see colors produced
#'
#' image(seq_along(oat[,2]), 1, as.matrix(seq_along(oat[,2])),
#'  col=sectcoldf$col, axes=FALSE, xlab="", ylab="")


lmvdencolor <- function(df,
                        wsize = 30,
                        bias = 5,
                        colorin = colorRampPalette(RColorBrewer::brewer.pal(8, "Spectral"))(25)) {
  editlgdf(df)
  mapthese <- unique(df[, 1])

  chr <- vector(mode = "character")
  s <- vector(mode = "double")
  e <- vector(mode = "double")
  col <- vector(mode = "character")
  dens <- vector(mode = "double")

  # assume if wsize=30 it defaulted.  If we need to adjust sectsize then we will
  # adjust wsize also since it defaulted.  If user passed it, use user value.
  if (wsize == 30) {
    adjwsize <- TRUE
  }
  else {
    adjwsize <- FALSE
  }

  # determine size of section to color.  Increase size until < 1000 sections for
  # largest linkage group (in case input is in Mbp instead of cM)
  sectsize <- 1
  if (max(df[, 2]) > 1000) {
    nbrsects <- max(df[, 2])
    sectsize <- 1
    while (nbrsects > 1000) {
      sectsize <- sectsize * 10
      # adjust wsize in the same way
      if (adjwsize) {
        wsize <- wsize * 10
      }
      nbrsects <- max(df[, 2]) / sectsize
    }
  }
  # calculate density within linkage group in sectsize position units
  dfix <- 0
  for (i in 1:length(mapthese)) {
    thislg <- df[, 2][df[, 1] == mapthese[i]]
    k <- min(thislg) + sectsize / 2
    while (k <= max(thislg) - sectsize / 2) {
      dfix <- dfix + 1
      chr[dfix] <- mapthese[i]
      s[dfix] <- k - sectsize / 2
      # if < 1 sectionleft to end go to end
      if ((k + (sectsize/2 * 3)) > max(thislg)) {
        e[dfix] <- max(thislg)
      }
      else {
        e[dfix] <- k + sectsize / 2
      }
      if (sum(thislg < (k + wsize / 2) &
              thislg > (k - wsize / 2)) > 0) {
        # calculate actual window size
        awsize <-
          min(k + wsize / 2, max(thislg)) - max(k - wsize / 2, min(thislg))
        dens[dfix] <- awsize / (sum(thislg <= (k + wsize / 2) &
                                      thislg >= (k - wsize / 2)))
      }
      else {
        # means we have a window size region with no markers
        # so cM/locus is infinity
        message(c("There are sections of linkage group ", mapthese[i],
                  " with no markers found between ",
                  k - wsize / 2, " and ",  k + wsize / 2,". Your window size is ", wsize, "."))
        message("These sections will be colored at the highest density.")
        message("Or use function lmvdencolor to set an alternative window size. ")
        dens[dfix] <- NA  # will set to 1 after scaling
      }
      col[dfix] = "default"
      k <- k + sectsize
    }
  }
  sectcoldf <-
    data.frame(chr, s, e, col, dens, stringsAsFactors = FALSE)

  # scale all densities to range of 0 to 1
  densrange <- range(sectcoldf$dens, na.rm = TRUE) #ignore NAs
  if (diff(densrange) > 0) {
    densvals  <- (sectcoldf$dens - densrange[1]) / diff(densrange)
  }
  else
  {
    densvals <- rep(0, length(sectcoldf$dens))
  }
  # set NA to 1 as highest position / marker
  densvals[is.na(densvals)] <- 1.0

  f <- colorRamp(colorin, bias = bias)  # create color ramp function

  sectcoldf$col <- rgb(f(densvals) / 255)
  sectcoldf

}
