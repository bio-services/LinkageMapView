# extract one linkage group

getlg <- function(df, lg, dupnbr,roundpos) {
    # if dupnbr is true user wants only first marker name
    # at a position to show with (### more) after the marker name
    # to indicate how many more markers at that position
    # We'll use (### more) as the second label so that the first
    # label isn't corrupted for connecting homologous markers
    # based on matching marker names.  That means we need to
    # set maxnbrcolsfordups <- 2 in lmv.linkage.plot.

    templg <- subset(df, df[, 1] == lg)
    # make sure linkage group is in order by position
    thislg <- templg[order(templg$position),]
    thislg$position <- round(thislg$position, roundpos)
    if (dupnbr) {
        keeprows <- vector()
        dupcount <- vector(length = nrow(thislg))
        keeprows <- append(keeprows, 1)  # always keep first row
        dupcount[1] <- 0
        i <- 0
        for (i in 2:nrow(thislg)) {
            if (thislg$position[i] == thislg$position[i - 1]) {
                dupcount[i] <- dupcount[i - 1] + 1
                if (dupcount[i] == 1) {
                    keeprows <- append(keeprows, i) # for (## more label)
                }
            }
            else {
                if (dupcount[i - 1] > 0) {
                    thislg$locus[i - dupcount[i - 1]] <-
                        paste(" (",
                              dupcount[i - 1],
                              " more)",
                              sep = "")
                }
                keeprows <- append(keeprows, i)
                dupcount[i] <- 0
            }
        }
        # finish last row
        if (dupcount[i] > 0) {
            thislg$locus[i - dupcount[i - 1]] <-
                paste(" (",
                      dupcount[i],
                      " more)",
                      sep = "")
        }
        retdf <- thislg[unique(keeprows), ]
    }
    else {
        retdf <- thislg
    }
    return(retdf)
}
