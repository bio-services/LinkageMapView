#'  LinkageMapView plotting function
#'
#'  lmv is the main function to produce linkage group maps and has many
#'      parameters to customize the pdf output.
#'
#' @param mapthis Required, either a 'cross' object from r/qtl, a csv or txt
#'        file or a data frame with the following 3 columns in this order:
#'      \enumerate{
#'        \item Linkage group name. This will be the title for the linkage
#'              group unless overridden - see lgtitles.
#'        \item Position - must be in numerical order ascending within
#'              linkage group name.
#'        \item Locus - marker name at this position.
#'      }
#'
#' @param outfile Required, name for the output pdf file.
#'
#' @param mapthese Optional vector of linkage group names to print.
#'        The default, NULL, will print all linkage groups in mapthis.
#'
#' @param autoconnadj If TRUE (the default), locus with the same name
#'        (homologs) on adjacent linkage groups will be connected with a line.
#'
#' @param bg Background coloer for the pdf.  Default is "transparent".
#'
#' @param cex.main The magnification to be used for linkage group titles.
#'        The default is par("cex.main").
#'
#' @param col.main The color to be used for linkage group titles.
#'        Defaults to par("col.main").
#'
#' @param conndf An optional dataframe containing markers to be connected
#'        with lines (homologs).  If autoconnadj = TRUE, these lines will
#'        appear as well as those with the same name in adjacent linkage
#'        groups.  Required columns:
#'        \itemize{
#'        \item fromchr Linkage group for the line start.
#'        \item fromlocus Locus name for the line start.
#'        \item tochr Linkage group for the line end.
#'        \item tolocus Locus name for the line end.
#'        }
#'
#' @param dupnbr If TRUE, only the first marker name at a position will print
#'        with (## more) afterwards indicating the number of duplicate markers
#'        at that position.  dupnbr should be left to the default, FALSE,
#'        if showonly provided.
#'
#' @param family Font family for all text.  Default is "Helvitica".
#'
#' @param fg Foreground color for the pdf.  Default is black.
#'
#' @param font.main An integer which specifies which font to use for title text.
#'        The default is par("font.main"). 1 is plain text.  2 is bold.
#'        3 is italic. 4 is bold italic.
#'
#' @param header A boolean indicating if the input file has a header row.
#'        Default is TRUE.
#'
#' @param labdist Distance in inches from the chromosome to the position
#'        and locus labels. The default is 0.3 inches.
#'
#' @param lcex A numerical value giving the amount by which position labels
#'        should be magnified.  The default is par("cex").
#'        See also rcex for locus labels.
#'
#' @param lcol The color for the position labels.  The default is par("col").
#'        See also rcol for locus labels.
#'
#' @param lfont An integer which specifies which font to use for the position
#'        labels.  The default is par("font").  See also rfont for locus labels.
#'
#' @param lgperrow An integer specifying how many linkage groups to plot in
#'        one row.  As many rows as needed to plot all requested linkage
#'        groups will be plotted.
#'
#' @param lgtitles Optional vector of titles for the linkage groups.  These
#'        will override the default, which is that the linkage group names
#'        in the input print as titles.  This may be useful if in mapthese you
#'        have indicated to print the same linkage group more than once for the
#'        purpose of showing homologous markers without having lines cross.
#'
#' @param lgw Width of chromosome in inches.  Default is 0.25 inches.
#'
#' @param lg.col Linkage group color.  The color of the chromosomes.
#'        The default is the background color (bg).
#'
#' @param lg.lwd Linkage group linewidth. The width of the line around
#'        the chromosome.  Defaults to par("lwd").
#'
#' @param markerformatlist An optional list containing the following vectors:
#'        \itemize{
#'          \item locus Required.  A vector of loci for which the following
#'          should be applied.
#'          \item col Optional.  The color for these locus labels.  This color
#'          will override rcol.  See also rsegcol.
#'          \item cex Optional.  A numerical vlue giving the amount by which
#'          these locus labels should be magnified.  This value will override
#'          rcex.
#'          \item font Optional. An integer which specifies which font to use
#'          for these locus labels. This value will override rfont.
#'        }
#'
#' @param maxnbrcolsfordups Indicates the number of columns across the page for
#'        locus labels appearing at duplicate positions.  The default is 3.
#'
#' @param pdfheight Height of the output file in inches.  Defaults to 7.
#'
#' @param pdfwidth Width of the output file in inches.  Defaults to 7.
#'
#' @param pointsize The default point size to be used.  Defaults to 12.
#'
#' @param posonleft A vector of boolean (TRUE or FALSE) the length of the
#'        number of linkage groups to be plotted. if FALSE, print positions on
#'        right hand side of linkage group and locus names on left hand side
#'        of linkage group.  Default is TRUE.
#'
#' @param prtlgtitles If FALSE do not print linkage group titles.
#'        Default is TRUE.
#'
#' @param qtldf An optional data frame containing QTL information for plotting.
#'        The data frame, if provided, must contain:
#'        \itemize{
#'          \item chr Linkage group name for QTL.
#'          \item qtl Name (label) for QTL.
#'          \item so Start of outer interval. Numeric.
#'          \item si Start of inner interval. Numeric.
#'          \item ei End of inner interval. Numeric.
#'          \item eo End of outer interval. Numeric.
#'          \item col Color for QTL.
#'        }
#'
#' @param qtlscanone Optional scanone dataframe from package r/qtl.  If provided,
#'        all QTLs in the dataframe will be drawn by calculating their
#'        start and end with the r/qtl function bayesint with defaults.
#'
#' @param revthese Optional vector of linkage group names to reverse.  The end
#'        position becomes position 0 and position 0 becomes the end position.
#'
#' @param rcex A numerical value giving the amount by which locus labels
#'        should be magnified.  The default is par("cex").
#'        See also lcex for position labels.
#'
#' @param rcol The color for the locus labels.  The default is par("col").
#'        See also lcol for position labels.
#'
#' @param rfont An integer which specifies which font to use for the locus
#'        labels.  The default is par("font").
#'        See also lfont for position labels.
#'
#' @param roundpos Number of positions after the decimal point to print for
#'        positions.  Default is 1
#'
#' @param rsegcol Color of the segments across the chromosome and to the label.
#'        TRUE, the default, indicates the color should be the
#'        same as the label.
#'
#' @param ruler A single boolean (TRUE OR FALSE).  If TRUE, a cM ruler is
#'        drawn on the left hand side of the page and the position labels
#'        are not printed on any linkage group.  The default is FALSE.
#'
#' @param sectcoldf Optional data frame containing the following named columns
#'       indicating sections of the chromosome to be colored:
#'         \itemize{
#'           \item chr - matches from input file or cross object
#'            linkage group name
#'           \item s - start position in cM
#'           \item e - end position in cM
#'           \item col - color for section
#'          }
#'
#' @param showonly Optional vector of marker names.  If provided, only these
#'        marker names will be printed.
#'
#' @param title Title to be passed to pdf as metadata.  This title does not
#'        appear except in the pdf metadata.  Defaults to
#'        "LinkageMapView R output".
#'
#' @export
#'
#' @examples
#'
#' ## take a cross object from r/qtl and produce linkage map
#' ## on chr 1,4,6,15
#'
#' library(qtl)
#' data(hyper)
#' lmv(hyper,"hyper.pdf",mapthese=c(1,4,6,15))
#'
#' ## color some of the markers for emphasis
#'
#' library(qtl)
#' data(hyper)
#'
#' # make a list to pass label options
#' flist <- list()
#' locus <- c("D1Mit123","D1Mit105","D6Mit273","D15Mit56","D15Mit156")
#' col   <- c("red")
#' flist[[1]] <- list(locus=locus,col=col)
#' lmv(hyper,"hyperred.pdf",mapthese=c(1,4,6,15),markerformatlist=flist)
#'
#' ## change some of the pdf options and chromosome color
#' ## changing title color (col.main) to same as foreground pdf color
#'
#' library(qtl)
#' data(hyper)
#'
#' lmv(hyper,"hyperlg.pdf",
#' mapthese=c(1,4,6,15),
#' bg="black",fg="white",col.main="white",
#' pdfheight=8,title="myhyper",lg.col="tan")
#'
#' ## change all label colors and fonts
#'
#' library(qtl)
#' data(hyper)
#'
#' lmv(hyper,"hypercol.pdf",mapthese=c(1,4,6,15),
#' lcol="blue",lfont=2,lcex=1.2,rcol="red",rfont=3,rcex=2)
#'
#' ## make a dataframe to pass sections of chr to col
#' ## use a ruler instead of printing positions as labels
#' ## only allow one column for duplicate markers at same position
#' ## (default is 3)
#'
#' library(qtl)
#' data(hyper)
#'
#' chr = c(1, 4, 6, 15)
#' s = c(82,35,9.8,7.7)
#' e = c(94,47,21.9,13.1)
#' col = c("pink","blue","blue","green")
#' sectcoldf <-  data.frame(chr, s, e, col,stringsAsFactors = FALSE)
#'
#' lmv(hyper,"hyperruler.pdf",mapthese=c(1,4,6,15),
#' ruler=TRUE,maxnbrcolsfordups = 1, sectcoldf=sectcoldf)
#'
#' ## plot qtls also out of a r/qtl scanone object
#' ## plot marker names on left (instead of right) of chr 4 and 7
#'
#' library(qtl)
#' data(hyper)
#'
#' # create scanone df for testing
#' hyper <-
#'  calc.genoprob(hyper,
#'                step = 2.0,
#'                map.function = "haldane",
#'                stepwidth = "fixed")
#' hyper.scanone <- scanone(hyper)
#'
#' lmv(hyper,
#'    "testrqtlhyper2.pdf",mapthese=c(1,4,6,7,15),
#'    qtlscanone = hyper.scanone,
#'    posonleft = c(TRUE,FALSE,TRUE,FALSE,TRUE))
#'
#' ## plot a carrot comparative linkage map
#' ## kindly provided by Massimo Iorizzo:
#' ## Cavagnaro et al. BMC Genomics 2014, 15:1118
#'
#' # make a df to pass qtl info
#' qtldf <- data.frame(
#'   chr = character(),
#'   qtl = character(),
#'   so = numeric(),
#'   si = numeric(),
#'   ei = numeric(),
#'   eo = numeric(),
#'   col = character(),
#'   stringsAsFactors = FALSE
#' )
#' qtldf <- rbind(qtldf,
#'                data.frame(
#'                  chr = "70349LG3",
#'                  qtl = "RTPE-Q1",
#'                  so = 36.6,
#'                  si = 37,
#'                  ei = 37,
#'                  eo = 38,
#'                  col="red"
#'                ))

#' # make a list to pass label options
#' flist <- list()

#' locus <- c("BSSR-094", "K0149", "K0627", "K2161", "ESSR-087", "ESSR-057")
#' font  <- c(2)   #bold
#' flist[[1]] <- list(locus = locus, font = font)

#' locus <- c("F3H", "FLS1")
#' font  <- c(4)   #bold italic
#' flist[[2]] <- list(locus = locus, font = font)

#' locus <- c("P3", "P1", "Raa1")
#' font  <- c(3)   #italic
#' col <- c("red")
#' flist[[3]] <- list(locus = locus, font = font, col = col)
#' filename <- system.file("extdata", "carrot.csv", package="LinkageMapView")

#' lmv(
#'   mapthis = filename,
#'   outfile = "carrot.pdf",
#'   ruler = TRUE,
#'   lgtitle = c("2170", "70349", "10117"),
#'   maxnbrcolsfordups = 1,
#'   markerformatlist = flist,
#'   lg.col = "lightblue1",
#'   pdfwidth =10,
#'   revthese = c("70349LG3"),
#'   qtldf=qtldf
#' )


lmv <- function(mapthis,
                outfile,
                mapthese = NULL,
                autoconnadj = TRUE,
                bg = "transparent",
                cex.main = par("cex.main"),
                col.main = par("col.main"),
                conndf = NULL,
                dupnbr = FALSE,
                family = "Helvetica",
                fg = "black",
                font.main = par("font.main"),
                header=TRUE,
                labdist = .3,
                lcex = par("cex"),
                lcol = par("col"),
                lfont = par("font"),
                lgperrow = NULL,
                lgtitles = NULL,
                lgw = 0.25,
                lg.col = NULL,
                lg.lwd = par("lwd"),
                markerformatlist = NULL,
                maxnbrcolsfordups = 3,
                pdfwidth = NULL,
                pdfheight = NULL,
                pointsize = 12,
                posonleft = NULL,
                prtlgtitles = TRUE,
                qtldf = NULL,
                revthese = NULL,
                rcex = par("cex"),
                rcol = par("col"),
                rfont = par("font"),
                roundpos = 1,
                rsegcol = TRUE,
                ruler = FALSE,
                sectcoldf = NULL,
                qtlscanone = NULL,
                showonly = NULL,
                title = "LinkageMapView R output"
                )


{
  pgx <- 0.5   # where on page x-axis to draw
  pgy <- 0.5   # where on page y-axis to draw

  oldpar <-
    par(no.readonly = TRUE)   # save parameters for restoring
  on.exit(par(oldpar))

  # ----------------- Begin edit arguments -----------------------------------
  if (!is.null(roundpos)) {
    if (!is.numeric(roundpos)) {
      stop("roundpos - number of places after decimal point - must be numeric")
    }
  }

  # read input for further edits ---

  if ("cross" %in% class(mapthis)) {
    if (is.null(mapthese)) {
      mapthese <- qtl::chrnames(mapthis)
    }
    lgin <- readlgcross(mapthis, mapthese)
  }
  else if ("character" == class(mapthis)) {
    if (is.null(mapthese)) {
      lgin <- readlgtext(mapthis,header=header)
      mapthese <- unique(lgin$group)
    } else{
      lgin <- readlgtext(mapthis, mapthese, header=header)
    }
  }
  else if ("data.frame" == class(mapthis)) {
    if (is.null(mapthese)) {
      lgin <- readlgdf(mapthis)
      mapthese <- unique(lgin$group)
    } else{
      lgin <- readlgdf(mapthis, mapthese)
    }
  }
  else {
    stop("first parameter, mapthis, must be a filename or an r/qtl cross object")
  }

  if (!is.null(revthese)) {
    if (!all(revthese %in% mapthese)) {
      stop ("Requested chrnames to reverse not in those requested to plot")
    }
  }

  if (!is.null(lgperrow)) {
    if (!(is.numeric(lgperrow) && floor(lgperrow) == lgperrow)) {
      stop ("lgperrow must be NULL or an integer")
    }
  }
  else {
    lgperrow <- length(mapthese)
  }

  nbrrows <- ceiling(length(mapthese) / lgperrow)

  # if input has null or NA Locus names convert existing locus names
  # to showonly list so position labels won't show for the null/NA
  # ones

  if(any(lgin$locus == "") || any(is.na(lgin$locus))) {
    notnull <- lgin$locus[which(lgin$locus != "")]
    showonly <- notnull[which(!is.na(notnull))]
  }

  if (!is.null(showonly)) {
    if (!all(showonly %in% lgin$locus)) {
      stop (c("Requested showonly locus not found:", showonly[!(showonly %in% lgin$locus)]))
    }
    else {
      # dups really screw up the code, correct for the user
      showonly <- unique(showonly)
    }
  }
  # make sure qtl data frame passed has no factors
  if (!is.null(qtldf)) {
    fas <- sapply(qtldf, is.factor)
    qtldf[fas] <- lapply(qtldf[fas], as.character)
  }
  # make sure qtlscanone is df and add it to (or create) qtldf
  if (!is.null(qtlscanone)) {
    if (!("data.frame" %in% class(qtlscanone) & "scanone" %in% class(qtlscanone))) {
      stop (c("qtlscanone should be a data.frame for r/qtl")) }
    else {
      qtldf <- usescanone(qtlscanone,qtldf,mapthese,fg,maxdec=roundpos)
    }
  }

  # make sure sectcoldf data frame passed has no factors
  if (!is.null(sectcoldf)) {
    fas <- sapply(sectcoldf, is.factor)
    sectcoldf[fas] <- lapply(sectcoldf[fas], as.character)
  }

  # make sure vector of positioning requests is valid
  if (!is.null(posonleft)) {
    if (!(all(posonleft %in% c(T, F)) &&
          length(posonleft) == length(mapthese))) {
      stop(
        "Position on vector must be same length as # of linkage groups to map and only TRUE or FALSE"
      )
    }
  }
  else {
    posonleft <- rep(TRUE, length.out = length(mapthese))
  }

  # make sure ruler is valid
  if (!(is.null(ruler)) && !(length(ruler) == 1)) {
    stop("ruler must be a single TRUE or FALSE")
  }


  if (!(ruler == T) && !(ruler == F))              {
    stop("ruler may only be TRUE or FALSE")
  }

  # make sure dupnbr is valid
  if (!(is.null(dupnbr)) && !(length(dupnbr) == 1)) {
    stop("dupnbr must be a single TRUE or FALSE")
  }


  if (!(dupnbr == T) && !(dupnbr == F))              {
    stop("dupnbr may only be TRUE or FALSE")
  }

  if (dupnbr == T) {
    maxnbrcolsfordups <- 2
  }

  if (!(prtlgtitles == T) && !(prtlgtitles == F))              {
    stop("prtlgtitles may only be TRUE or FALSE")
  }
  # ----------------- End edit arguments -----------------------------------


  pdf.options(
    bg = bg,
    title = title,
    family = family,
    pointsize = pointsize,
    fg = fg
  )
  on.exit(pdf.options(reset = TRUE), add = TRUE)

  # pdf size doesn't matter here - just for reqdim plotting
  # pdf will be reallocated before actually drawing

  pdf(outfile, width = 30, height = 30)
  on.exit(dev.off(), add = TRUE)


  lg <- list()
  dim <- list()
  lrcol <- list()
  llcol <- list()
  lrfont <- list()
  llfont <- list()
  lrcex <- list()
  llcex <- list()
  qtl <- list()
  solist <- list()
  sectcol <- list()

  # make a list of linkage groups to process
  # reverse position numbers if requested
  for (i in 1:length(mapthese)) {
    lg[[i]] <- getlg(lgin, mapthese[i], dupnbr)
    editlgdf(lg[[i]])  # display message and stop if not in correct format
    if (lg[[i]]$group[1] %in% revthese) {
      lg[[i]]$locus <- rev(lg[[i]]$locus)
      lg[[i]]$position <- revpos(lg[[i]]$position, roundpos)
    }
    else {
      lg[[i]]$position <- round(lg[[i]]$position, roundpos)
    }

    # add to the list for this linkage group the qtls
    # and reverse postions if necessary

    if (!is.null(qtldf))
    {
      qtl[[i]] <- subset(qtldf, qtldf$chr == lg[[i]][1, 1])
      for (qtlix in 1:nrow(qtl[[i]])) {
        if (nrow(subset(qtldf, qtldf$chr == lg[[i]][1, 1])) != 0) {
          if (qtl[[i]]$chr[qtlix] %in% revthese) {
            # reverse positions
            so <-
              qtl[[i]]$so[qtlix]
            si <-
              qtl[[i]]$si[qtlix]
            ei <-
              qtl[[i]]$ei[qtlix]
            eo <-
              qtl[[i]]$eo[qtlix]
            qtl[[i]]$eo[qtlix] <-
              max(lg[[i]]$position) - so - min(lg[[i]]$position)
            qtl[[i]]$ei[qtlix] <-
              max(lg[[i]]$position) - si - min(lg[[i]]$position)
            qtl[[i]]$si[qtlix] <-
              max(lg[[i]]$position) - ei - min(lg[[i]]$position)
            qtl[[i]]$so[qtlix] <-
              max(lg[[i]]$position) - eo - min(lg[[i]]$position)
          }
        }
      }
    }

    # if user requested sections to be colored make list for this lg
    if (!is.null(sectcoldf))
    {
      sectcol[[i]] <- subset(sectcoldf, sectcoldf$chr == lg[[i]][1, 1])
    }
  }

  # -- Begin determine smallest and largest position for each row --------

  miny <- vector(length = nbrrows)
  maxy <- vector(length = nbrrows)

  for (nr in 1:nbrrows) {
    fromlg <- (nr - 1) * lgperrow + 1
    tolg <- min(c(length(mapthese),  nr * lgperrow))

    miny[nr] <- min(lg[[fromlg]]$position)
    maxy[nr] <- max(lg[[fromlg]]$position)

    for (i in fromlg:tolg)     {
      if (max(lg[[i]]$position) > maxy[nr]) {
        maxy[nr] <- max(lg[[i]]$position)
      }
      if (min(lg[[i]]$position) < miny[nr]) {
        miny[nr] <- min(lg[[i]]$position)
      }
    }
  }

  # -- End determine smallest and largest position for each row ----------

  # -- Begin loop for lg - apply formats & get required size -------------
  for (i in 1:length(mapthese)) {
    #if any of these options specified - apply first
    lrcol[[i]] <- rcol
    lrcol[[i]] = rep(lrcol[[i]], length.out = length(lg[[i]]$locus))

    lrfont[[i]] <- rfont
    lrfont[[i]] = rep(lrfont[[i]], length.out = length(lg[[i]]$locus))

    lrcex[[i]] <- rcex
    lrcex[[i]] <-
      rep(lrcex[[i]], length.out = length(lg[[i]]$locus))

    llcol[[i]] <- lcol
    llcol[[i]] = rep(llcol[[i]], length.out = length(lg[[i]]$position))

    llfont[[i]] <- lfont
    llfont[[i]] = rep(llfont[[i]], length.out = length(lg[[i]]$position))

    llcex[[i]] <- lcex
    llcex[[i]] <-
      rep(llcex[[i]], length.out = length(lg[[i]]$position))

    # next apply specific requst by locus
    if (!is.null(markerformatlist)) {
      for (ml in 1:length(markerformatlist)) {
        if (!is.null(markerformatlist[[ml]]$col) &&
            !is.na(markerformatlist[[ml]]$col)) {
          lrcol[[i]][which(lg[[i]]$locus %in% markerformatlist[[ml]]$locus)] <-
            markerformatlist[[ml]]$col
        }
        if (!is.null(markerformatlist[[ml]]$cex) &&
            !is.na(markerformatlist[[ml]]$cex)) {
          lrcex[[i]][which(lg[[i]]$locus %in% markerformatlist[[ml]]$locus)] <-
            markerformatlist[[ml]]$cex
        }
        if (!is.null(markerformatlist[[ml]]$font) &&
            !is.na(markerformatlist[[ml]]$font)) {
          lrfont[[i]][which(lg[[i]]$locus %in% markerformatlist[[ml]]$locus)] <-
            markerformatlist[[ml]]$font
        }
      }
    }
    if (!is.null(qtldf)) {
      qtldfone <- qtl[[i]]
      if (nrow(qtl[[i]]) == 0) {
        qtldfone <- NULL
      }
    }
    else {
      qtldfone <- NULL
    }


    dim[[i]] <- reqdim(
      lg[[i]],
      c(miny[nr], maxy[nr]),
      maxnbrcolsfordups = maxnbrcolsfordups,
      pdfwidth = 30,
      labdist = labdist,
      lcol = llcol[[i]],
      lfont = llfont[[i]],
      lcex = llcex[[i]],
      rcol = lrcol[[i]],
      rfont = lrfont[[i]],
      rcex = lrcex[[i]],
      cex.main = cex.main,
      qtldf = qtldfone,
      ruler = ruler,
      prtlgtitles = prtlgtitles,
      showonly = showonly
    )

  }

  # -- End loop for lg - apply formats & get required size -------------

  pgxlg <- vector(length = length(mapthese))
  width <- vector(length = length(mapthese))

  totwidth <- rep(0, nbrrows)
  totheight <- rep(0, nbrrows)

  for (nr in 1:nbrrows) {
    fromlg <- (nr - 1) * lgperrow + 1
    tolg <- min(c(length(mapthese),  nr * lgperrow))

    #    calculate relative widths of plot areas
    #    and y range for all to be mapped

    # give a little for left and right margin first

    dim[[fromlg]]$reqwidth <-
      dim[[fromlg]]$reqwidth + 0.3
    dim[[tolg]]$reqwidth <-
      dim[[tolg]]$reqwidth + 0.3

    for (i in fromlg:tolg)    {
      totwidth[nr] <- dim[[i]]$reqwidth + totwidth[nr]
      if (dim[[i]]$reqheight > totheight[nr]) {
        totheight[nr] <- dim[[i]]$reqheight
      }
    }
  }

  allrowwidth <- max(totwidth)
  allrowheight <- sum(totheight)

  message(c("Required pdfwidth = ", allrowwidth))
  message(c("Required pdfheight = ", allrowheight))

  # determine relative height of each row for layout
  relheight <- vector(length = nbrrows)

  for (nr in 1:nbrrows) {
    fromlg <- (nr - 1) * lgperrow + 1
    tolg <- min(c(length(mapthese),  nr * lgperrow))

    maxrowheight <- 0
    for (i in fromlg:tolg)    {
      if (dim[[i]]$reqheight > maxrowheight) {
        maxrowheight <- dim[[i]]$reqheight
      }
    }
    relheight[nr] <- maxrowheight / allrowheight
  }

  # if user did not specify size, set to required size
  if (is.null(pdfwidth)) {
    pdfwidth <- ceiling(allrowwidth)
  }
  if (is.null(pdfheight)) {
    pdfheight <- ceiling(allrowheight)
  }

  message(c("Using pdfwidth = ", pdfwidth))
  message(c("Using pdfheight = ", pdfheight))

  if (ruler) {
    leftmar <- 2
  }
  else {
    leftmar <- 0
  }
  if (prtlgtitles) {
    cextitle = cex.main
  } else {
    cextitle = 0
  }

  # turn off pdf used just for sizing and start the real one
  dev.off()
  pdf.options(
    bg = bg,
    title = title,
    family = family,
    pointsize = pointsize,
    fg = fg
  )
  pdf(outfile, width = pdfwidth, height = pdfheight)
  layout(c(seq(1, nbrrows)), heights = relheight)
  par(mar = c(1, leftmar, cextitle +2, 0))

  # --- Begin loop for nbr rows to draw linkage groups -----------------
  for (nr in 1:nbrrows) {
    plot(
      .5,
      .5,
      xlim = c(0, 1),
      ylim = c(maxy[nr], miny[nr]),
      type = "n",
      cex = 1,
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      xaxs = "i",
      # don't pad x axis on each side
      bty = "n"
    )

    pin <- par("pin")[1]

    if (ruler) {
      lastlab <- floor(maxy[nr] / 5) * 5
      tickat <- seq(0, floor(maxy[nr]))
      axlab <- vector()
      for (lab in 0:maxy[nr]) {
        if (!lab %% 5) {
          axlab <- c(axlab, lab)
        }
        else {
          axlab <- c(axlab, NA)
        }
      }

      axis(
        side = 2,
        at = tickat,
        labels = axlab,
        col.axis = col.main
      )

    }

    # if last row put footnote
    if (nr == nbrrows) {
      mtext(
        "Rendered by LinkageMapView",
        side = 1,
        cex = 0.5,
        col = fg,
        adj = 0
      )
    }
    # set up list for returned values from drawone
    # needed for placement of segments connecting markers
    # between linkage groups

    yrlabwidth <- list()
    adjyr <- list()
    adjyl <- list()
    dups <- list()


    widthused <- 0

    fromlg <- (nr - 1) * lgperrow + 1
    tolg <- min(c(length(mapthese),  nr * lgperrow))

    # ----------- Begin loop for each lg in this row -----------------------
    for (i in fromlg:tolg) {
      pctwidth <- dim[[i]]$reqwidth / totwidth[nr]
      width <- pctwidth * pin
      if (posonleft[i]) {
        if (!ruler) {
          llablen <- dim[[i]]$maxlenllab
          lspace <-
            strwidth("M", units = "inches") * max(lcex)
          llabdist <- labdist
        }
        else {
          if (i > 1 && !posonleft[i - 1]) {
            llablen <- 0
            # make sure titles do not overlap when mirror no labels come together
            lspace <-
              max(strwidth(lg[[i]][1, 1]) * cex.main ,
                  strwidth(lg[[i - 1]][1, 1]) * cex.main) / 2 + strwidth("M", units = "inches") * cex.main
            llabdist <- 0
          }
          else {
            llablen <- 0
            lspace <- 0
            llabdist <- 0
          }
        }
      }
      else {
        llablen <- dim[[i]]$maxlenrlab
        lspace <-
          strwidth("M", units = "inches") * max(rcex)
        llabdist <- labdist
      }
      pgxlg[i] <-
        (sum(llablen,
             lgw / 2,
             llabdist,
             lspace)) / pin  + widthused / pin

      # give a little for left margin
      if (i == 1) {
        pgxlg[i] <- pgxlg[i] + 0.3 / pdfwidth
      }
      widthused <- widthused + width

      if (!is.null(qtldf)) {
        qtldfone <- qtl[[i]]
        if (nrow(qtl[[i]]) == 0) {
          qtldfone <- NULL
        }
      }
      else {
        qtldfone <- NULL
      }

      if (!is.null(sectcoldf)) {
        sectcoldfone <- sectcol[[i]]
        if (nrow(sectcol[[i]]) == 0) {
          sectcoldfone <- NULL
        }
      }
      else {
        sectcoldfone <- NULL
      }


      if (!is.null(lgtitles)) {
        lgtitleone <- lgtitles[i]
      }
      else {
        lgtitleone <- NULL
      }

      dolist <-
        drawone(
          lg[[i]],
          dim[[i]],
          totwidth[nr],
          c(miny[nr], maxy[nr]),
          maxnbrcolsfordups = maxnbrcolsfordups,
          pdfwidth = pin,
          lgw = lgw,
          lg.col = lg.col,
          lg.lwd = lg.lwd,
          bg = bg,
          fg = fg,
          pgx = pgxlg[i] ,
          labdist = labdist ,
          lcol = llcol[[i]],
          lfont = llfont[[i]],
          lcex = llcex[[i]],
          rcol = lrcol[[i]],
          rfont = lrfont[[i]],
          rcex = lrcex[[i]],
          rsegcol = rsegcol,
          cex.main = cex.main,
          font.main = font.main,
          col.main = col.main,
          qtldf = qtldfone,
          posonleft = posonleft[i],
          ruler = ruler,
          prtlgtitles = prtlgtitles,
          lgtitles = lgtitleone,
          showonly = showonly,
          sectcoldf = sectcoldfone
        )
      yrlabwidth[[i]] <- dolist$yrlabwidth
      adjyr[[i]] <- dolist$adjyr
      adjyl[[i]] <- dolist$adjyl
      dups[[i]] <- dolist$dups
      solist[[i]] <- dolist$solist

    }

    # ----------- End loop for each lg in this row -----------------------

    # connect markers across linkage groups if requested

    # make sure connect marker data frame passed has no factors

    if (autoconnadj == TRUE)
    {
      #only pass the markers on this row

      autoconndf <- autoconn(lg[fromlg:tolg], fg, lgperrow)
      if (!is.null(conndf)) {
        fas <- sapply(conndf, is.factor)
        conndf[fas] <- lapply(conndf[fas], as.character)
        allconndf <- rbind(conndf, autoconndf)
        # get rid of duplicates if user specified and automatically
        # since auto is at the end, the user specified will be kept
        allconndf <- allconndf[!duplicated(allconndf[,1:4]),]
      }
      else {
        allconndf <- autoconndf
      }
    }
    else {
      if (!is.null(conndf)) {
        fas <- sapply(conndf, is.factor)
        conndf[fas] <- lapply(conndf[fas], as.character)
        allconndf <- conndf
      }
      else{
        allconndf <- data.frame() # empty
      }
    }

    if (nrow(allconndf) > 0) {
      for (i in 1:nrow(allconndf)) {
        # determine index for from and to chr

        fromi <- match(allconndf$fromchr[i], mapthese)
        if (is.na(fromi)) {
          stop(
            c(
              "Connect marker from position not found for chr = ",
              allconndf$fromchr[i],
              " and locus = ",
              allconndf$fromlocus[i]
            )
          )
        }
        toi <- match(allconndf$tochr[i], mapthese)
        if (is.na(toi)) {
          stop(
            c(
              "Connect marker to position not found for chr = ",
              allconndf$tochr[i],
              " and locus = ",
              allconndf$tolocus[i]
            )
          )
        }
        if (toi < fromi) {
          temp <- toi
          toi <- fromi
          fromi <- temp
        }


        # user requested connections must be in same row
        if (ceiling(fromi / lgperrow) != ceiling(toi / lgperrow)) {
          stop(
            c(
              "Connect marker from chr ",
              allconndf$fromchr[i],
              " must be on the same row as to chr ",
              allconndf$tochr[i],
              " and you have lgperrow = ",
              lgperrow
            )
          )
        }

        # determine x and y for from and to marker
        # look up marker in lgin to get position
        fpos <-
          lg[[fromi]]$position[lg[[fromi]]$locus == allconndf$fromlocus[i]]
        if (identical(fpos, numeric(0))) {
          stop(
            c(
              "Connect marker from position not found for chr = ",
              allconndf$fromchr[i],
              " and locus = ",
              allconndf$fromlocus[i]
            )
          )
        }
        tpos <-
          lg[[toi]]$position[lg[[toi]]$locus == allconndf$tolocus[i]]
        if (identical(tpos, numeric(0))) {
          stop(
            c(
              "Connect marker to position not found for chr = ",
              allconndf$tochr[i],
              " and locus = ",
              allconndf$tolocus[i]
            )
          )
        }

        if (!is.null(showonly)) {
          fy <- match(fpos, solist[[fromi]]$newllab)
          ty <- match(tpos, solist[[toi]]$newllab)
        }
        else
        {
          fy <- match(fpos, lg[[fromi]]$position)
          ty <- match(tpos, lg[[toi]]$position)
        }
        # determine from
        if (posonleft[fromi]) {
          connfxpos <-
            ((yrlabwidth[[fromi]][fy] + lgw / 2 + labdist / 2) / pin) + pgxlg[fromi]
          if (!is.null(showonly)) {
            connfypos <-
              adjyr[[fromi]][match(fpos, solist[[fromi]]$newllab[dups[[fromi]]$rkeep])]
          }
          else {
            connfypos <-
              adjyr[[fromi]][match(fpos, lg[[fromi]]$position[dups[[fromi]]$rkeep])]
          }
        }

        else {
          if (!ruler) {
            connfxpos <-
              ((
                dim[[fromi]]$maxlenllab + lgw / 2 + labdist + strwidth("M", units = "inches")
              ) / pin) + pgxlg[fromi]
            if (!is.null(showonly)) {
              connfypos <-
                adjyl[[fromi]][match(fpos, solist[[fromi]]$newllab[dups[[fromi]]$lkeep])]
            }
            else {
              connfypos <-
                adjyl[[fromi]][match(fpos, lg[[fromi]]$position[dups[[fromi]]$lkeep])]
            }

          } else
          {
            connfxpos <- ((lgw / 2) / pin) + pgxlg[fromi]
            connfypos <-
              fpos

          }
        }

        if (!posonleft[toi]) {
          conntxpos <-
            -((yrlabwidth[[toi]][ty]  + lgw / 2 + labdist / 2) / pin) + pgxlg[toi]
          if (!is.null(showonly)) {
            conntypos <-
              adjyr[[toi]][match(tpos, solist[[toi]]$newllab[dups[[toi]]$rkeep])]
          }
          else {
            conntypos <-
              adjyr[[toi]][match(tpos, lg[[toi]]$position[dups[[toi]]$rkeep])]
          }
        }

        else {
          if (!ruler) {
            conntxpos <-
              -((
                dim[[toi]]$maxlenllab  +  lgw / 2 + labdist + strwidth("M", units = "inches")
              ) / pin) + pgxlg[toi]
            if (!is.null(showonly)) {
              conntypos <-
                adjyl[[toi]][match(tpos, solist[[toi]]$newllab[dups[[toi]]$lkeep])]
            }
            else {
              conntypos <-
                adjyl[[toi]][match(tpos, lg[[toi]]$position[dups[[toi]]$lkeep])]
            }
          }
          else {
            conntxpos <- -((lgw / 2) / pin) + pgxlg[toi]
            conntypos <-
              tpos

          }
        }
        if (is.null(allconndf$col[i])) {
          allconndf$col[i] = fg
        }
        segments(connfxpos,
                 connfypos,
                 conntxpos,
                 conntypos,
                 col = allconndf$col[i])
      }
    }

  }
  # --- End loop for nbr rows to draw linkage groups -----------------

}
