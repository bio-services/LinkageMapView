# spread out labels with different fontsizes

# variation of spread.labs from package TeachingDemos
# but handles a passed in cex vector for adjusting font sizes

spreadcexlabs <-
  function(x,
           mindiff,
           strh,
           maxiter = 7500,
           stepsize = 1 / 10,
           min = -Inf,
           max = Inf,
           cex = par("cex")) {
    # assumes ordered input due to variable sizes of labels based on cex

    cex = rep(cex, length.out = length(x))
    mingap <- mindiff - strh

    reqdiff <- vector(length = length(x) - 1)

    reqdiff[1:(length(x) - 1)] <-
      0.5 * strh * cex[1:(length(x) - 1)] + 0.5 * strh * cex[2:(length(x))] + mingap

    df <-
      x[-1] - x[-length(x)]  # get differences between points

    stp <- mindiff * stepsize

    i <- 1
    while (any(df < reqdiff)) {
      tmp <- c(df < reqdiff, FALSE)
      if (tmp[1] &&
          (x[1] - stp) < min) {
        # don't move bottom set
        tmp2 <- as.logical(cumprod(tmp))
        tmp <- tmp & !tmp2
      }
      x[tmp] <- x[tmp] - stp  # move all non-fitters down


      tmp <- c(FALSE, df < reqdiff)
      if (tmp[length(tmp)] &&
          (x[length(x)] + stp) > max) {
        # don't move top
        tmp2 <- rev(as.logical(cumprod(rev(tmp))))
        tmp <- tmp & !tmp2
      }
      x[tmp] <- x[tmp] + stp   # move new all non-fitters up

      df <- x[-1] - x[-length(x)]
      i <- i + 1
      if (i > maxiter) {
        #                warning("Maximum iterations reached")
        break
      }
    }

    return(x)
  }
