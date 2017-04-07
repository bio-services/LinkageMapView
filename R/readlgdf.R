# read passed dataframe

readlgdf <-
  function(df,
           mapthese)
  {
    if (ncol(df) < 3) {
      for (i in 1:ncol(df)) {
        stop("Less than 3 columns found on input data frame - see help for mapthis parameter")
      }

    }

    # assign my own column names so I can reference by name
    colnames(df)[1:3] <- c("group", "position", "locus")

    # make group a character field

    df$group <- as.character(df$group)

    if (!missing(mapthese)) {
      if (!all(mapthese %in% df$group)) {
        stop ("chrnames to map not found in input data frame")
      }
    }

    # make sure data frame passed has no factors
    if (!is.null(df)) {
      fas <- sapply(df, is.factor)
      df[fas] <- lapply(df[fas], as.character)
    }

    return (df)
  }
