# read text or csv input

readlgtext <-
  function(fn,
           mapthese,
           header = TRUE,
           stringsAsFactors = FALSE)
  {
    # get file extension
    fnparts <- strsplit(fn, ".", fixed = TRUE)[[1]]
    if (length(fnparts) > 1) {
      ext = fnparts[length(fnparts)]
    } else {
      ext <- ""
    }

    # figure out how many rows and columns on input file
    # and issue message

    if (ext == 'csv') {
      nbrcols <- count.fields(fn, sep = ",")
    }
    else {
      nbrcols <- count.fields(fn)
    }

    if (any(nbrcols < 2)) {
      for (i in 1:length(nbrcols)) {
        if (nbrcols[i] < 2) {
          message(c("line ", i, " has only ", nbrcols[i], " columns."))
        }
      }
      stop("Less than 2 columns found on input file")
    }
    if (any(nbrcols < 3)) {
      for (i in 1:length(nbrcols)) {
        if (nbrcols[i] < 3) {
          message(c("row ", i, " is missing locus name."))
        }
      }
      warning(" Locus name missing on some lines")
    }

    if (ext == 'csv') {
      df <- read.csv(
        fn,
        header = header,
        stringsAsFactors = FALSE,
        strip.white = TRUE
      )
    }
    else {
      df <-
        read.table(
          fn,
          header = header,
          stringsAsFactors = FALSE,
          colClasses = c("character", "numeric", "character"),
          strip.white = TRUE
        )

    }
    # assign my own column names so I can reference by name
    colnames(df) <- c("group", "position", "locus")

    # make group a character field

    df$group <- as.character(df$group)

    if (!missing(mapthese)) {
      if (!all(mapthese %in% df$group)) {
        stop ("chrnames to map not found in input file")
      }
    }

    return (df)
  }
