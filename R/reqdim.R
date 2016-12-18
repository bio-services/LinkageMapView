# calculate required dimension for one linkage group
# based on being passed label columns position and locus
# width of chromosome and length of distance to labels

reqdim <- function(df,
                   yrange,
                   pdfwidth = 12,
                   maxnbrcolsfordups = 3,
                   lgw = 0.25,
                   pgx = 0.5 ,
                   labdist = 0.1 ,
                   rcex = par("cex"),
                   lcex = par("cex"),
                   rfont = par("font"),
                   lfont = par("font"),
                   rcol = par("col"),
                   lcol = par("col"),
                   cex.main = par("cex.main"),
                   qtldf = NULL,
                   ruler = FALSE,
                   prtlgtitles = TRUE,
                   showonly = showonly)
{
  if (ruler) {
    leftmar <- 1
  } else {
    leftmar <- 0
  }
  if (prtlgtitles) {
    cextitle = cex.main
  } else {
    cextitle = 0
  }
  par(mar = c(1, leftmar, cextitle + 2, 0))      # only need margin at top for title

  # set up points to label - the points will be plotted invisibly
  # x1 is on the left side of the chromosome, x2 on the right
  # these will be recalculated after required dimensions are determined

  lgwpct <- lgw / pdfwidth
  x1 <-
    rep(pgx - lgwpct / 2, length(df$position))
  x2 <-
    rep(pgx + lgwpct / 2, length(df$position))

  y <- df$position
  rlab <- df$locus
  llab <- df$position

  # ylim reversed so 0 is at top of y-axis
  # don't print the points: type="n"
  # and don't print axis: xaxt and yaxt = "n"
  # and don't print axis labels: yaxt="n",xlab=""
  # and don't print the box around the plot: bty="n"
  # in other words don't print anything but establish the points
  # for the markers


  plot(
    x2,
    y,
    xlim = c(0, 1),
    ylim = rev(range(y)),
    type = "n",
    cex = 1,
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    xaxs = "i",
    bty = "n"

  )

  if (!is.null(showonly)) {
    solist <- show(
      showonly,
      llab,
      rlab,
      rcex = rcex,
      lcex = lcex,
      rfont = rfont,
      lfont = lfont,
      rcol = rcol,
      lcol = lcol
    )
    llab <- solist$newllab
    rlab <- solist$newrlab
    rcex <- solist$newrcex
    lcex <- solist$newlcex
    rfont <- solist$newrfont
    lfont <- solist$newlfont
    rcol <- solist$newrcol
    lcol <- solist$newlcol
  }

  # find dups to figure out how many columns to reserve space for
  dups <- fsdups(llab, maxnbrcolsfordups)

  # Determine width of each right label including dups

  yrlabwidth <- vector(length = length(llab))
  yrlabwidth[setdiff(dups$rkeep, dups$frkeep)] <-
    strwidth(rlab[setdiff(dups$rkeep, dups$frkeep)], units = "inches") *
    rcex[setdiff(dups$rkeep, dups$frkeep)] + labdist
  yrlabwidth[dups$frkeep] <-
    strwidth(rlab[dups$frkeep], units = "inches") *
    rcex[dups$frkeep] + labdist * 1.2
  # now add space for duplicates in cols
  if (maxnbrcolsfordups > 1) {
    for (i in 1:length(llab)) {
      for (m in 1:(maxnbrcolsfordups - 1)) {
        if (!is.na(dups$yd[m, i])) {
          yrlabwidth[i] <-
            yrlabwidth[i] + strwidth(rlab[(i + m)], units = "inches") * rcex[(i + m)] + (strwidth(" ", units =
                                                                                                    "inches")  * rcex[(i + m)]) / 2
        }
      }
    }
  }

  # If qtl provided, calculate width required
  # Since we don't know pdf dimensions yet, we
  # can't determine if qtl will overlap each other.
  # so we assume they will and add enough space to
  # the required width for each qtl

  # adding in strwidth("M",units=inches) takes care of
  # 0.5 of a character width default between points and labels
  # and leaves a little room between linkage groups and at edge of pdf

  if (!ruler) {
    reqwidth = sum(
      max(strwidth(llab, units = "inches") * lcex),
      max(yrlabwidth),
      lgw,
      labdist * 2,
      strwidth("M", units = "inches")  * max(lcex),
      strwidth("M", units  = "inches")  * max(rcex),
      nrow(qtldf) * (lgw / 3 + strheight("M", units = "inches") *
        3)
    )

    reqheight = max(sum(strheight(llab, units = "inches") * lcex * 1.4),
                    sum(strheight(rlab, units = "inches") * rcex * 1.4),
                    sum(strheight(llab, units = "inches") * rcex * 1.4)) # because positions spread like markers
  }
  else {
    reqwidth = sum(max(yrlabwidth),
                   lgw / 2,
                   labdist,
                   strwidth("M", units  = "inches")  * max(rcex),
                   nrow(qtldf) * (lgw / 3 + strheight("M", units = "inches") *
                                    3))

    reqheight = sum(strheight(rlab, units = "inches") * rcex * 1.4)

  }

  # give a margin at top and bottom for chromosome ends
  reqheight = reqheight + lgw
  # add in margins at top and bottom
  reqheight = reqheight + strheight("M", units = "inches") * (cex.main +
                                                                2)

  list(
    reqwidth = reqwidth,
    reqheight = reqheight,
    maxlenrlab = max(yrlabwidth),
    maxlenllab = max(strwidth(llab, units = "inches") * lcex)
  )

}
