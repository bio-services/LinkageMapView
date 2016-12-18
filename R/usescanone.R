# extract qtl from scanone output of r/qtl package

usescanone <- function(qtlscanone, qtldf=NULL, mapthese, fg, maxdec) {

  if (is.null(qtldf)) {
    # make a dataframe to hold qtl info
    qtldf <- data.frame(
      chr = character(),
      qtl = character(),
      so = numeric(),
      si = numeric(),
      ei = numeric(),
      eo = numeric(),
      col = character(),
      stringsAsFactors = FALSE
    )
  }
  # extract qtl info from scanone df
  qtlchr <- unique(qtlscanone[, 1])
  for (i in 1:length(qtlchr)) {
    if (qtlchr[i] %in% mapthese) {
      qtlone <- qtl::bayesint(qtlscanone, chr = qtlchr[i])
      qtldf <- rbind(
        qtldf,
        data.frame(
          chr = qtlchr[i],
          qtl = paste(qtlchr[i], "@" , round(qtlone[2,2], digits = maxdec)),
          so = qtlone[1, 2],
          si = qtlone[2, 2],
          ei = qtlone[2, 2],
          eo = qtlone[3, 2],
          col=fg
        )
      )
    }
  }
  # remove factors if any
  fas <- sapply(qtldf, is.factor)
  qtldf[fas] <- lapply(qtldf[fas], as.character)
  return (qtldf)
}
