# draw linkage group

drawone <-
  function(df,
           dim,
           totwidth,
           yrange,
           maxnbrcolsfordups = 3,
           pdfwidth = 12,
           lgw = 0.25,
           lg.col = NULL,
           lg.lwd = par("lwd"),
           bg = bg,
           fg = fg,
           pgx = 0.5 ,
           labdist = 0.1 ,
           rcex = par("cex"),
           lcex = par("cex"),
           rfont = par("font"),
           lfont = par("font"),
           rcol = par("col"),
           lcol = par("col"),
           rsegcol = TRUE,
           cex.main = cex.main,
           font.main = font.main,
           col.main = col.main,
           qtldf = NULL,
           posonleft = TRUE,
           ruler = FALSE,
           prtlgtitles = TRUE,
           lgtitles = NULL,
           showonly = NULL,
           sectcoldf = NULL) {
    y <- df$position
    rlab <- df$locus
    llab <- df$position

    pctwidth <- dim$reqwidth / totwidth
    width <- pctwidth * pdfwidth

    lgwpct <- lgw / pdfwidth

    # if user requested to have positions show on right
    # set up adjustments


    if (posonleft) {
      labdistpct <- labdist / pdfwidth
      x1 <- rep(pgx - lgwpct / 2, length(df$position))
      x2 <- rep(pgx + lgwpct / 2, length(df$position))
      rpos <- 4
      lpos <- 2
      posmult <- 1    # adjuster for left or right positioning
    }
    else {
      labdistpct <- -labdist / pdfwidth
      x1 <- rep(pgx + lgwpct / 2, length(df$position))
      x2 <- rep(pgx - lgwpct / 2, length(df$position))
      rpos <- 2
      lpos <- 4
      posmult <- -1
    }

    points(
      x2,
      y,
      type = "n",
      cex = 1,
      xlab = "",
      ylab = ""
    )

    # put title on top of chromosome
    if (prtlgtitles) {
      if (!is.null(lgtitles)) {
        lgtext = lgtitles
      }
      else {
        lgtext = df[1, 1]
      }
      mtext(
        lgtext,
        at = pgx,
        line = floor(cex.main),
        cex = cex.main,
        col = col.main,
        font = font.main
      )
    }

    # eliminate all but showonly labels if requested

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
    else {
      solist <- NULL
    }

    # find and save dup locations before calling spreadcexlabs
    dups <- fsdups(llab, maxnbrcolsfordups)

    # spread the labels as necessary
    # heights are negative because y axis is reversed

    if (length(dups$rkeep) > 1) {
      adjyr <- spreadcexlabs(
        llab[dups$rkeep],
        max(strheight(rlab)) * -1.4,
        strh = -max(strheight(rlab)),
        min = min(yrange),
        max = max(yrange),
        cex = rcex[dups$rkeep],
        maxiter = 99999
      )

    }
    else {
      adjyr <- llab[dups$rkeep]
    }

    pos = rpos

    # label the first columns except those that are dups(
    if (length(setdiff(dups$rkeep, dups$frkeep)) > 0) {
      text(
        x2[setdiff(dups$rkeep, dups$frkeep)] + labdistpct,
        adjyr[setdiff(1:length(adjyr), dups$fykeep)],
        labels = rlab[setdiff(dups$rkeep, dups$frkeep)],
        pos = pos,
        col = rcol[setdiff(dups$rkeep, dups$frkeep)],
        cex = rcex[setdiff(dups$rkeep, dups$frkeep)],
        font = rfont[setdiff(dups$rkeep, dups$frkeep)]
      )
    }

    # label the first columns of dups inset for visibility
    if (length(dups$frkeep > 0)) {
      text(
        x2[dups$frkeep] + labdistpct * 1.2,
        adjyr[dups$fykeep],
        labels = rlab[dups$frkeep],
        pos = pos,
        col = rcol[dups$frkeep],
        cex = rcex[dups$frkeep],
        font = rfont[dups$frkeep]
      )
    }

    # now fill out duplicates in columns
    if (maxnbrcolsfordups > 1) {
      for (i in 1:length(llab)) {
        for (m in 1:(maxnbrcolsfordups - 1)) {
          if (!is.na(dups$yd[m, i])) {
            if (m == 1) {
              rx <- x2[i] + labdistpct * 1.2
            }
            ry <- adjyr[dups$yd[m, i]]
            rx <-
              rx + posmult * (strwidth(rlab[(i + m - 1)]) * rcex[(i + m - 1)] + (strwidth(" ")  * rcex[(i + m -
                                                                                                          1)]) / 2)

            text(
              rx,
              ry,
              labels = rlab[(i + m)],
              pos = pos,
              col = rcol[(i + m)],
              cex = rcex[(i + m)],
              font =
                rfont[(i + m)]
            )
          }
        }
      }
    }

    # xpd = NA to turn off clipping and arcs can go into margins
    par(xpd = NA)
    if (!is.null(lg.col)) {
      rect(pgx - lgwpct / 2, min(y), pgx + lgwpct / 2, max(y), col = lg.col)
      symbols(
        x = pgx,
        y = min(y),
        circles = lgwpct / 2,
        bg = lg.col,
        add = TRUE,
        fg = lg.col,
        inches = FALSE
      )
      symbols(
        x = pgx,
        y = max(y),
        circles = lgwpct / 2,
        bg = lg.col,
        add = TRUE,
        fg = lg.col,
        inches = FALSE
      )
    }

    # color sections
    if (!is.null(sectcoldf)) {
      if (nrow(sectcoldf) > 0) {
        for (sc in 1:nrow(sectcoldf)) {
          rect(pgx - lgwpct / 2,
               sectcoldf$s,
               pgx + lgwpct / 2,
               sectcoldf$e,
               col = sectcoldf$col)
        }
      }
    }

    segments(pgx - lgwpct / 2, min(y), pgx - lgwpct / 2, max(y), lwd = lg.lwd)
    segments(pgx + lgwpct / 2, min(y), pgx + lgwpct / 2, max(y), lwd = lg.lwd)

    plotrix::draw.arc(
      x = pgx,
      y = min(y),
      radius = lgwpct / 2,
      deg1 = 0,
      deg2 = 180,
      lwd = lg.lwd
    )
    plotrix::draw.arc(
      x = pgx,
      y = max(y),
      radius = lgwpct / 2,
      deg1 = -180,
      deg2 = 0,
      lwd = lg.lwd
    )

    if (rsegcol) {
      segcol = rcol[setdiff(dups$rkeep, dups$frkeep)]
    }
    else {
      segcol = fg
    }
    # segments for nondups from chr to marker

    segments(x2[setdiff(dups$rkeep, dups$frkeep)] + labdistpct,
             adjyr[setdiff(1:length(adjyr), dups$fykeep)],
             x2[setdiff(dups$rkeep, dups$frkeep)],
             llab[setdiff(dups$rkeep, dups$frkeep)],
             col = segcol)
    #connect across chromosome
    segments(x1[setdiff(dups$rkeep, dups$frkeep)],
             llab[setdiff(dups$rkeep, dups$frkeep)],
             x2[setdiff(dups$rkeep, dups$frkeep)],
             llab[setdiff(dups$rkeep, dups$frkeep)], col = segcol)

    if (rsegcol) {
      segcol = rcol[dups$frkeep]
    }
    else {
      segcol = fg
    }
    #segments for dups
    if (length(dups$frkeep) > 0) {
      segments(x2[dups$fykeep] + labdistpct * 1.2, adjyr[dups$fykeep],
               x2[dups$fykeep],
               llab[dups$frkeep],
               col = segcol)

      #connect across chromosome
      segments(x1[dups$frkeep],
               llab[dups$frkeep],
               x2[dups$frkeep],
               llab[dups$frkeep], col = segcol)

      #segments connecting dups vertically
      if (length(y[dups$rkeep]) > 1) {
        for (i in 2:length(y[dups$rkeep])) {
          if (llab[dups$rkeep][i] == llab[dups$rkeep][i - 1]) {
            segments(x2[dups$rkeep][i] + labdistpct * 1.2,
                     adjyr[i],
                     x2[dups$rkeep][i] + labdistpct * 1.2,
                     adjyr[(i - 1)])
          }
        }
      }
    }

    # now draw lines across chromosome for non-shown markers if any

    if (!is.null(showonly)) {
      segments(
        rep(x1[1], length.out = length(setdiff(y, (
          llab[dups$rkeep]
        )))),
        setdiff(y, (llab[dups$rkeep])),
        rep(x2[1], length.out = length(setdiff(y, (
          llab[dups$rkeep]
        )))),
        setdiff(y, (llab[dups$rkeep])),
        col = fg
      )
    }

    pos = lpos
    if (!ruler) {
      text(
        x1[dups$lkeep] - labdistpct,
        adjyr[dups$flkeep],
        labels = llab[dups$lkeep],
        pos = pos,
        col = lcol[dups$lkeep],
        cex = lcex[dups$lkeep],
        font =
          lfont[dups$lkeep]
      )

      segments(x1[dups$lkeep] - labdistpct, adjyr[dups$flkeep], x1[dups$lkeep], llab[dups$lkeep])
    }


    # figure locus labwidth to pass back for connecting markers
    # and for drawing qtls

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

    # save yrlabwidth before adding in any qtls to pass back for connecting markers
    returnlabwidth <- yrlabwidth

    # draw qtls
    if (!is.null(qtldf)) {
      # save calculated start and end (biggest of actual positions or label length)
      # for determining if labels will overlap
      sdf <- vector(length = nrow(qtldf))
      edf <- vector(length = nrow(qtldf))
      # determine x position of qtl

      if (nrow(qtldf) > 0) {
        for (ql in 1:nrow(qtldf)) {
          # if there are labels on y axis where qtl needs to be drawn
          # first decide if QTL line (so to eo) or text label is longest
          oneinch <- grconvertY(0, from = "user", to = "inches") -
            grconvertY(1, from = "user", to = "inches")
          lablen <- strwidth(qtldf$qtl[ql], units = "inches")
          lablenY <- lablen / oneinch
          if (lablenY > qtldf$eo[ql] - qtldf$so[ql]) {
            # use lablenY start is middle - half of label length
            qtlstart <-
              ((qtldf$si[ql] + qtldf$ei[ql]) / 2) - (lablenY / 2)
            qtlend <-
              ((qtldf$si[ql] + qtldf$ei[ql]) / 2) + (lablenY / 2)
          }
          else {
            qtlstart <- qtldf$so[ql]
            qtlend <- qtldf$eo[ql]
          }
          for (f in 1:length(sdf)) {
            if (qtlstart > sdf[f] && qtlstart < edf[f]) {
              qtlstart <- sdf[f]
            }
            if (qtlend < edf[f] && qtlend > sdf[f]) {
              qtlend <- edf[f]
            }
          }
          sdf[ql] <- qtlstart
          edf[ql] <- qtlend
          # adjust in case labels overlap


          if (any(adjyr >= (qtlstart + strheight("M") * max(rcex)) &
                  adjyr <= (qtlend   - strheight("M") * max(rcex)))) {
            qtlxpos <-
              max(yrlabwidth[yrlabwidth > 0][which(adjyr >= (qtlstart + strheight("M") * max(rcex)) &
                                                     adjyr <= (qtlend - strheight("M") * max(rcex)))]) / pdfwidth
            yrlabwidth[yrlabwidth > 0][which(adjyr >= (qtlstart + strheight("M") * max(rcex)) &
                                               adjyr <= (qtlend  - strheight("M") * max(rcex)))] <-
              yrlabwidth[yrlabwidth > 0][which(adjyr >= (qtlstart + strheight("M") * max(rcex)) &
                                                 adjyr <= (qtlend  - strheight("M") * max(rcex)))] + lgw / 3 + strheight("M", units = "inches") *
              3
          }
          else {
            qtlxpos <- 0
          }

          if (is.null(qtldf$col[ql])) {
            qtldf$col[ql] <- par("col")
          }
          # draw qtl
          segments(
            x2[1] + posmult * (qtlxpos + strwidth("M") + lgwpct / 6),
            qtldf$so[ql],
            x2[1] + posmult * (qtlxpos + strwidth("M") + lgwpct / 6),
            qtldf$eo[ql],
            col = qtldf$col[ql]
          )
          # make inner region 1/3 size of linkage group width
          # 1/3 is arbitrary, just looks good
          rect(
            x2[1] + posmult * (qtlxpos + strwidth("M")),
            qtldf$si[ql],
            x2[1] + posmult * (qtlxpos + strwidth("M") + lgwpct / 3),
            qtldf$ei[ql],
            col = qtldf$col[ql],
            border = qtldf$col[ql]
          )
          text(
            x2[1] + posmult * (qtlxpos + strwidth("MM") + lgwpct / 3),
            (qtldf$si[ql] + qtldf$ei[ql]) / 2,
            label = qtldf$qtl[ql],
            # col = qtldf$col[ql],
            srt = 270
          )
        }
      }
    }
    return (
      list(
        adjyr = adjyr,
        adjyl = adjyr[dups$flkeep],
        yrlabwidth = returnlabwidth,
        dups = dups,
        solist = solist
      )
    )
  }
