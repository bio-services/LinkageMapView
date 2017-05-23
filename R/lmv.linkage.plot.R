#'  LinkageMapView plotting function
#'
#'  lmv.linkage.plot is the main function to produce linkage group maps and has many
#'      parameters to customize the pdf output.
#'
#' @param mapthis Required, either a 'cross' object from r/qtl, a csv or txt
#'        file or a data frame with the following 3 columns in this order:
#'      \enumerate{
#'        \item Required, linkage group name. This will be the title for the linkage
#'              group unless overridden - see lgtitles.
#'        \item Required, position - must be in numerical order ascending within
#'              linkage group name.
#'        \item Required, locus - marker name at this position.
#'        \item Optional, segcol - color for the line across the chromosome
#'              at this marker.  See also segcol parameter.
#'      }
#'
#' @param outfile Required, name for the output pdf file.
#'
#' @param mapthese Optional vector of linkage group names to print.
#'        The default, NULL, will print all linkage groups in mapthis.
#'
#' @param at.axis Optional. The points at which tick-marks are to be drawn on the ruler.
#'        Non-finite (infinite, NaN or NA) values are omitted.
#'        By default (when NULL) tickmark locations are computed.
#'  #' @seealso \code{\link[graphics]{axis}}
#'
#' @param autoconnadj If TRUE (the default), locus with the same name
#'        (homologs) on adjacent linkage groups will be connected with a line.
#'
#' @param cex.axis The magnification to be used for axis (ruler) text.
#'        The default is par("cex.axis").
#'
#' @param cex.lgtitle The magnification to be used for linkage group titles.
#'        The default is par("cex.main").
#'
#' @param cex.main The magnification to be used for main title.
#'        The default is par("cex.main").
#'
#' @param col.axis The color to be used for axis (ruler) text.
#'        Defaults to par("col.axis").
#'
#' @param col.lgtitle The color to be used for linkage group titles.
#'        Defaults to par("col.main").
#'
#' @param col.main The color to be used for the main title.
#'        Defaults to par("col.main").
#'
#' @param conndf An optional data frame containing markers to be connected
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
#' @param denmap If TRUE, you are requesting a density map which means no locus
#'        or position labels will be printed and the following parameters are
#'        set:
#'        ruler = TRUE
#'        autoconndf = FALSE
#'        conndf = NULL
#'        See also sectcoldf parameter
#'
#' @param dupnbr If TRUE, only the first marker name at a position will print
#'        with (## more) afterwards indicating the number of duplicate markers
#'        at that position.  dupnbr should be left to the default, FALSE,
#'        if showonly provided.
#'
#' @param font.axis An integer which specifies which font to use for the
#'        axis (ruler) text.
#'        The default is par("font.axis"). 1 is plain text.  2 is bold.
#'        3 is italic. 4 is bold italic.
#'
#' @param font.lgtitle An integer which specifies which font to use for the
#'        linkage group titles text.
#'        The default is par("font.main"). 1 is plain text.  2 is bold.
#'        3 is italic. 4 is bold italic.
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
#' @param labels.axis	Optional. This can either be a logical value specifying
#'        whether (numerical) annotations are to be made at the tickmarks on the ruler,
#'        or a character or expression vector of labels to be placed at the tickpoints.
#'        If this is not logical, at should also be supplied and of the same length.
#'        The default is TRUE.
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
#'        See also cex.lgtitle, col.lgtitle, font.lgtitle
#'
#' @param lgw Width of chromosome in inches.  Default is 0.25 inches.
#'
#' @param lg.col Linkage group color.  The color of the chromosomes.
#'        The default is the background color (pdf.bg).
#'
#' @param lg.lwd Linkage group linewidth. The width of the line around
#'        the chromosome.  Defaults to par("lwd").
#'
#' @param lty.axis Optional. Line type for both the axis line and the tick marks.
#'
#' @param lwd.axis Optional.  Line width for the axis line.  The default is 1.
#' @param lwd.ticks.axis Optional.  Line width for the axis tick marks.
#'        Default is lwd.axis
#'
#' @param main An optional title for the linkage group map.  See also
#'        cex.main, col.main, and font.main.
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
#' @param pdf.bg Background color for the pdf.  Default is "transparent".
#'
#' @param pdf.family Font family for all text.  Default is "Helvetica".
#'
#' @param pdf.fg Foreground color for the pdf.  Default is black.
#'
#' @param pdf.height Height of the output file in inches.  Defaults to the size
#' necessary to fit all linkage groups with other options specified.
#'
#' @param pdf.pointsize The default point size to be used.  Defaults to 12.
#'
#' @param pdf.title Title to be passed to pdf as metadata.  This title does not
#'        appear except in the pdf metadata.  Defaults to
#'        "LinkageMapView R output".
#'
#' @param pdf.width Width of the output file in inches.  Defaults to the size
#' necessary to fit all linkage groups with other options specified.
#'
#' @param posonleft A vector of boolean (TRUE or FALSE) the length of the
#'        number of linkage groups to be plotted. If FALSE, print positions on
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
#' @param qtlscanone Optional scanone data frame from package r/qtl.  If provided,
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
#' @param ruler A single boolean (TRUE OR FALSE).  If TRUE, an axis is
#'        drawn on the left hand side of the page and the position labels
#'        are not printed on any linkage group.  The default is FALSE.
#'
#' @param sectcoldf Optional data frame containing the following named columns
#'       indicating sections of the chromosome to be colored:
#'         \itemize{
#'           \item Required, chr - matches from input file or cross object
#'            linkage group name
#'           \item Required, s - start position in cM
#'           \item Required, e - end position in cM
#'           \item Required, col - color for section
#'           \item Optional, dens - a numeric cm / marker value used to print
#'            the density map legend.
#'          }
#'        For a density map, use the lmvdencolor function to populate sectcoldf.
#'        When denmap = TRUE and no sectcoldf parameter is supplied, lmvdencolor
#'        is called with defaults fully populating the sectcoldf data frame.
#'        See also the denmap parameter.
#'
#'        @seealso \code{\link{lmvdencolor}}
#'
#' @param segcol Optional. Name of the column in mapthis that contains colors for the line
#'        segments across the chromosome.  If specified, this overrides rsegcol.
#'
#' @param showonly Optional vector of marker names.  If provided, only these
#'        marker names will be printed.
#'
#' @param units Units of the position values supplied in mapthis.  The default
#'        value is cM (centimorgan) but any value can be provided.  The
#'        value provided is only used for a ruler (y axis) title and the density
#'        map legend text.
#'
#' @param ylab Optional.  Title for the y-axis (ruler).  The default value is
#'        "Density (<units>)"  See units parameter.
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
#' outfile = file.path(tempdir(), "hyper.pdf")
#' lmv.linkage.plot(hyper,outfile,mapthese=c(1,4,6,15))
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
#'
#' outfile = file.path(tempdir(), "hyperred.pdf")
#' lmv.linkage.plot(hyper,outfile,mapthese=c(1,4,6,15),markerformatlist=flist)
#'
#' ## change some of the pdf options and chromosome color
#' ## changing linkage group title color (col.lgtitle) to same as
#' ## foreground pdf color
#'
#' library(qtl)
#' data(hyper)
#'
#' outfile = file.path(tempdir(), "hyperlg.pdf")
#' lmv.linkage.plot(hyper,outfile,
#' mapthese=c(1,4,6,15),
#' pdf.bg="black",pdf.fg="white",col.lgtitle="white",
#' pdf.height=8,pdf.title="myhyper",lg.col="tan")
#'
#' ## change all label colors and fonts
#'
#' library(qtl)
#' data(hyper)
#'
#' outfile = file.path(tempdir(), "hypercol.pdf")
#' lmv.linkage.plot(hyper,outfile,mapthese=c(1,4,6,15),
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
#' outfile = file.path(tempdir(), "hyperruler.pdf")
#' lmv.linkage.plot(hyper,outfile,mapthese=c(1,4,6,15),
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
#' outfile = file.path(tempdir(), "testrqtlhyper2.pdf")
#' lmv.linkage.plot(hyper,
#'    outfile, mapthese=c(1,4,6,7,15),
#'    qtlscanone = hyper.scanone,
#'    posonleft = c(TRUE,FALSE,TRUE,FALSE,TRUE))
#'
#' \dontrun{
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
#' filename <- system.file("extdata", "Carrot.csv", package="LinkageMapView")

#' outfile = file.path(tempdir(), "carrot.pdf")
#' lmv.linkage.plot(
#'   mapthis = filename,
#'   outfile = outfile,
#'   ruler = TRUE,
#'   lgtitle = c("2170", "70349", "10117"),
#'   maxnbrcolsfordups = 1,
#'   markerformatlist = flist,
#'   lg.col = "lightblue1",
#'   pdf.width =10,
#'   revthese = c("70349LG3"),
#'   qtldf=qtldf
#' )
#' }
#'
#' ## do a density map with default colors
#' data(oat)
#'
#' outfile = file.path(tempdir(), "oat_Mrg01.pdf")
#' lmv.linkage.plot(oat,outfile,mapthese=c("Mrg01","Mrg02"),denmap=TRUE)
#'
#' \dontrun{
#' ## do a density map and provide your own colors with lmvdencolor helper
#' data(oat)
#' ##
#' outfile = file.path(tempdir(), "oat_Mrg01_YlGn.pdf")
#'
#' sectcoldf <- lmvdencolor(oat,colorin =
#' colorRampPalette(RColorBrewer::brewer.pal(8, "YlGn"))(5))
#'
#' lmv.linkage.plot(oat,outfile,denmap=TRUE,sectcoldf=sectcoldf)
#' }


lmv.linkage.plot <- function(mapthis,
                outfile,
                mapthese = NULL,
                at.axis = NULL,
                autoconnadj = TRUE,
                cex.axis = par("cex.axis"),
                cex.lgtitle = par("cex.main"),
                cex.main = par("cex.main"),
                col.axis = par("col.axis"),
                col.lgtitle = par("col.main"),
                col.main = par("col.main"),
                conndf = NULL,
                denmap = FALSE,
                dupnbr = FALSE,
                font.axis = par("font.axis"),
                font.lgtitle = par("font.main"),
                font.main = par("font.main"),
                header = TRUE,
                labdist = .3,
                labels.axis = TRUE,
                lcex = par("cex"),
                lcol = par("col"),
                lfont = par("font"),
                lgperrow = NULL,
                lgtitles = NULL,
                lgw = 0.25,
                lg.col = NULL,
                lg.lwd = par("lwd"),
                lty.axis = "solid",
                lwd.axis = 1,
                lwd.ticks.axis = lwd.axis,
                main = NULL,
                markerformatlist = NULL,
                maxnbrcolsfordups = 3,
                pdf.bg = "transparent",
                pdf.family = "Helvetica",
                pdf.fg = "black",
                pdf.width = NULL,
                pdf.height = NULL,
                pdf.pointsize = 12,
                pdf.title = "LinkageMapView R output",
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
                segcol = NULL,
                qtlscanone = NULL,
                showonly = NULL,
                units = "cM",
                ylab = paste("Density (",units,")",sep=""))


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

  # edit and convert font from text if necessary

  font.axis <- convertfont("font.axis", font.axis)
  font.lgtitle <- convertfont("font.lgtitle", font.lgtitle)
  font.main <- convertfont("font.main", font.main)
  lfont <- convertfont("lfont", lfont)
  rfont <- convertfont("rfont", rfont)
  if (!is.null(markerformatlist$font)) {
    markerformatlist$font <-
      convertfont("markerformatlist$font", markerformatlist$font)
  }

  # read input for further edits ---

  if ("cross" %in% class(mapthis)) {
    if (is.null(mapthese)) {
      mapthese <- qtl::chrnames(mapthis)
    }
    lgin <- readlgcross(mapthis, mapthese)
  }
  else if ("character" %in% class(mapthis)) {
    if (is.null(mapthese)) {
      lgin <- readlgtext(mapthis, header = header)
      mapthese <- unique(lgin$group)
    } else{
      lgin <- readlgtext(mapthis, mapthese, header = header)
    }
  }
  else if ("data.frame" %in% class(mapthis)) {
    if (is.null(mapthese)) {
      lgin <- readlgdf(as.data.frame(mapthis))
      mapthese <- unique(lgin$group)
    } else{
      lgin <- readlgdf(mapthis, mapthese)
    }
  }
  else {
    stop("first parameter, mapthis, must be a data frame, a filename, or an r/qtl cross object")
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

  if (denmap) {
    autoconnadj <- FALSE
    rsegcol <- FALSE
    ruler <- TRUE
    showonly <- NULL
    markerformatlist <- NULL
    conndf <- NULL
  }

  # user can specify colors for segments
  if (!is.null(segcol)) {
    # make sure segcol column exists in input
    if (!(segcol %in% colnames(lgin))) {
      stop (c("segcol column ", segcol, " not found in mapthis"))
    }
    else {
      rsegcol <- FALSE
    }
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

  # else if input has null or NA Locus names convert existing locus names
  # to showonly list so position labels won't show for the null/NA
  # ones

  else {
    if (any(lgin$locus == "") || any(is.na(lgin$locus))) {
      notnull <- lgin$locus[which(lgin$locus != "")]
      showonly <- notnull[which(!is.na(notnull))]
    }
  }

  # make sure qtl data frame passed has no factors
  if (!is.null(qtldf)) {
    fas <- sapply(qtldf, is.factor)
    qtldf[fas] <- lapply(qtldf[fas], as.character)
  }
  # make sure qtlscanone is df and add it to (or create) qtldf
  if (!is.null(qtlscanone)) {
    if (!("data.frame" %in% class(qtlscanone) &
          "scanone" %in% class(qtlscanone))) {
      stop (c("qtlscanone should be a data.frame for r/qtl"))
    }
    else {
      qtldf <-
        usescanone(qtlscanone, qtldf, mapthese, pdf.fg, maxdec = roundpos)
    }
  }


  if (is.null(sectcoldf) && denmap == TRUE) {
    # use default density map coloring
    sectcoldf <- lmvdencolor(lgin)
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
# ------------- Begin set par options --------------------

  # use top par("mar" for linkage group titles)
  # use top par("oma" for main title)
  # use left par("mar" for ruler)
  # use left par("oma" for ruler units)

  if (ruler) {
    leftmar <- 1 + ceiling(cex.axis)
    if (!is.null(units)) {
      leftmar <- leftmar + 1
    }
  }
  else {
    leftmar <- 1
  }

  if (prtlgtitles) {
    topmar <- ceiling(cex.lgtitle) + 1
  }
  else {
    topmar <- 1
  }

  if (!is.null(main)) {
    topoma <- 1+ ceiling(cex.main)
  }
  else {
    topoma <- 1
  }

  # ------------- End set par options --------------------

  pdf.options(
    bg = pdf.bg,
    title = pdf.title,
    family = pdf.family,
    pointsize = pdf.pointsize,
    fg = pdf.fg
  )
  on.exit(pdf.options(reset = TRUE), add = TRUE)

  # pdf size doesn't matter here - just for reqdim plotting
  # pdf will be reallocated before actually drawing

  pdf(outfile, width = 30, height = 30)
  on.exit(dev.off(), add = TRUE)

  par(mar = c(0, leftmar, topmar, 0),oma = c(1, 0, topoma, 1),new=FALSE)


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
    lg[[i]] <- getlg(lgin, mapthese[i], dupnbr,roundpos)
    editlgdf(lg[[i]])  # display message and stop if not in correct format
    if (lg[[i]]$group[1] %in% revthese) {
      lg[[i]]$locus <- rev(lg[[i]]$locus)
      lg[[i]]$position <- revpos(lg[[i]]$position, roundpos)
      #reverse segcol column if it exists
      if (!is.null(segcol)) {
        lg[[i]][[eval(segcol)]] <- rev(lg[[i]][[eval(segcol)]])
      }
      # if user requested sections to be colored make list for this lg
      if (!is.null(sectcoldf))
      {
        sectcol[[i]] <- rev(subset(sectcoldf, sectcoldf$chr == lg[[i]][1, 1]))
      }
    }
    else {
      lg[[i]]$position <- round(lg[[i]]$position, roundpos)
      # if user requested sections to be colored make list for this lg
      if (!is.null(sectcoldf))
      {
        sectcol[[i]] <- subset(sectcoldf, sectcoldf$chr == lg[[i]][1, 1])
      }
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
      denmap = denmap,
      maxnbrcolsfordups = maxnbrcolsfordups,
      pdf.width = 30,
      labdist = labdist,
      lcol = llcol[[i]],
      lfont = llfont[[i]],
      lcex = llcex[[i]],
      rcol = lrcol[[i]],
      rfont = lrfont[[i]],
      rcex = lrcex[[i]],
      cex.lgtitle = cex.lgtitle,
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
      dim[[fromlg]]$reqwidth + 0.5
 #   dim[[tolg]]$reqwidth <-
#    dim[[tolg]]$reqwidth + 0.3

    for (i in fromlg:tolg)    {
      totwidth[nr] <- dim[[i]]$reqwidth + totwidth[nr]
      if (dim[[i]]$reqheight > totheight[nr]) {
        totheight[nr] <- dim[[i]]$reqheight
      }
    }
  }

  allrowwidth <- max(totwidth)
  if (denmap & !is.null(sectcoldf$dens)) {
    # add a one 1/2 inch row for density map legend
    allrowheight <- sum(totheight) + 1.5
    relheight <- vector(length = (nbrrows + 1))
    relheight[(nbrrows + 1)] <- 1.5 / allrowheight
  }
  else {
    allrowheight <- sum(totheight)
    relheight <- vector(length = nbrrows)
  }

  message(c("Required pdf.width = ", allrowwidth))
  message(c("Required pdf.height = ", allrowheight))

  # determine relative height of each row for layout


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
  if (is.null(pdf.width)) {
    pdf.width <- ceiling(allrowwidth)
  }
  if (is.null(pdf.height)) {
    pdf.height <- ceiling(allrowheight)
  }

  message(c("Using pdf.width = ", pdf.width))
  message(c("Using pdf.height = ", pdf.height))

  # turn off pdf used just for sizing and start the real one
  dev.off()
  pdf.options(
    bg = pdf.bg,
    title = pdf.title,
    family = pdf.family,
    pointsize = pdf.pointsize,
    fg = pdf.fg
  )

  pdf(outfile, width = pdf.width, height = pdf.height)

  par(mar = c(0, leftmar, topmar, 0),oma = c(1, 0.2, topoma, 0.2),new=FALSE)

  if (denmap & !is.null(sectcoldf$dens) ) {
    layout(c(seq(1, nbrrows + 1)), heights = relheight)
  }
  else {
    layout(c(seq(1, nbrrows)), heights = relheight)
  }

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
      axis(
        side = 2,
        col.axis = col.axis,
        cex.axis = cex.axis,
        font.axis = font.axis,
        at = at.axis,
        labels = labels.axis,
        lty = lty.axis,
        lwd = lwd.axis,
        lwd.ticks = lwd.ticks.axis
      )
      mtext(ylab, side = 2, line=leftmar-1)
    }

    # if last row put footnote
    if (nr == nbrrows) {
      mtext(
        "Rendered by LinkageMapView",
        side = 1,
        cex = 0.5,
        col = pdf.fg,
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
              max(strwidth(lg[[i]][1, 1]) * cex.lgtitle ,
                  strwidth(lg[[i - 1]][1, 1]) * cex.lgtitle) / 2 + strwidth("M", units = "inches") * cex.lgtitle
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

      # give a little for left margin on first one on row
      if (i == (nr-1)*lgperrow +1) {
        pgxlg[i] <- pgxlg[i] + 0.3 / pdf.width
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

      if (i == 1) {
        #pass title to drawone
        main <- main
      }
      else {
        main = NULL
      }

      dolist <-
        drawone(
          lg[[i]],
          dim[[i]],
          totwidth[nr],
          c(miny[nr], maxy[nr]),
          denmap = denmap,
          maxnbrcolsfordups = maxnbrcolsfordups,
          pdf.width = pin,
          pdf.fg = pdf.fg,
          lgw = lgw,
          lg.col = lg.col,
          lg.lwd = lg.lwd,
          pgx = pgxlg[i] ,
          labdist = labdist ,
          lcol = llcol[[i]],
          lfont = llfont[[i]],
          lcex = llcex[[i]],
          rcol = lrcol[[i]],
          rfont = lrfont[[i]],
          rcex = lrcex[[i]],
          rsegcol = rsegcol,
          main = main,
          cex.main = cex.main,
          font.main = font.main,
          col.main = col.main,
          cex.lgtitle = cex.lgtitle,
          font.lgtitle = font.lgtitle,
          col.lgtitle = col.lgtitle,
          qtldf = qtldfone,
          posonleft = posonleft[i],
          ruler = ruler,
          prtlgtitles = prtlgtitles,
          lgtitles = lgtitleone,
          segcol = segcol,
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

      autoconndf <- autoconn(lg[fromlg:tolg], pdf.fg, lgperrow)
      if (!is.null(conndf)) {
        fas <- sapply(conndf, is.factor)
        conndf[fas] <- lapply(conndf[fas], as.character)
        allconndf <- rbind(conndf, autoconndf)
        # get rid of duplicates if user specified and automatically
        # since auto is at the end, the user specified will be kept
        allconndf <- allconndf[!duplicated(allconndf[, 1:4]), ]
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
          allconndf$col[i] = pdf.fg
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

  # if density map legend to be displayed
  if (denmap &
      !is.null(sectcoldf$dens)) {
    # plot legend is last row
    leg <- sectcoldf[order(sectcoldf$dens), ]
    if (max(leg$dens > 100)) {
      leg$dens <- round(leg$dens, digits = 0)
    }
    else {
      leg$dens <- round(leg$dens, digits = 1)
    }
    uleg <- leg[!duplicated(leg$dens), ]
    # reduce to reasonable number of buckets 1/4 inch wide bucket minimum
    if ((pdf.width / .25) < nrow(uleg)) {
      nbrbuckets <- round(pdf.width / .25, digits = 0)
      bplotdens <-
        uleg$dens[seq(1, length(uleg$dens), length(uleg$dens) / nbrbuckets)]
      bplotcol <-
        uleg$col[seq(1, length(uleg$col), length(uleg$col) / nbrbuckets)]
    } else {
      bplotdens <- uleg$dens
      bplotcol <- uleg$col
    }
    if (!bplotdens[length(bplotdens)] == uleg$dens[length(uleg$dens)]) {
      # include largest density
      bplotdens <- append(bplotdens, uleg$dens[length(uleg$dens)])
      bplotcol <- append(bplotcol, uleg$col[length(uleg$col)])
    }
    par(mar = c(5, leftmar, 1, 1))
    barplot(
      rep(1, length(bplotcol)),
      col = bplotcol,
      space = 0,
      axes = F,
      xlab = paste("Density (", units, "/Locus)", sep = ""),
      names = bplotdens,
      cex.names = .75,
      cex.lab = .75
    )

  }
}
