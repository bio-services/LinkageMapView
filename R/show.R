# remove all the left and right labels except those in the show only list

show <- function(showonly,
                 llab,
                 rlab,
                 rcex = par("cex"),
                 lcex = par("cex"),
                 rfont = par("font"),
                 lfont = par("font"),
                 rcol = par("col"),
                 lcol = par("col"))
{
  # when invoked all vectors except showonly should be the
  # same length - see lmv.linkgae.plot.R where they are set

  newllab <- vector()
  newrlab <- vector()
  newrcex <- vector()
  newlcex <- vector()
  newrfont <- vector()
  newlfont <- vector()
  newrcol <- vector()
  newlcol <- vector()

  # match(showonly, rlab) returns the index where there is a match
  # since the result is in showonly order, the index in rlab might
  # not be in order, sort them.  that also removes the NA since the
  # showonly might not be in this linkage group

  newllab <- append(newllab, llab[sort(match(showonly, rlab))])
  newrlab <- append(newrlab, rlab[sort(match(showonly, rlab))])
  newrcex <- append(newrcex, rcex[sort(match(showonly, rlab))])
  newlcex <- append(newlcex, lcex[sort(match(showonly, rlab))])
  newrfont <- append(newrfont, rfont[sort(match(showonly, rlab))])
  newlfont <- append(newlfont, lfont[sort(match(showonly, rlab))])
  newrcol <- append(newrcol, rcol[sort(match(showonly, rlab))])
  newlcol <- append(newlcol, lcol[sort(match(showonly, rlab))])


  # always show the first and last position label

  if (!(rlab[1] %in% showonly)) {
    newrlab <- append("", newrlab)
    newllab <- append(llab[1], newllab)
    newrcex <- append(rcex[1], newrcex)
    newlcex <- append(lcex[1], newlcex)
    newrfont <- append(rfont[1], newrfont)
    newlfont <- append(lfont[1], newlfont)
    newrcol <- append(rcol[1], newrcol)
    newlcol <- append(lcol[1], newlcol)
  }

  if (!(rlab[length(rlab)] %in% showonly)) {
    newrlab <- append(newrlab, "")
    newllab <- append(newllab, llab[length(llab)])
    newrcex <- append(newrcex, rcex[length(llab)])
    newlcex <- append(newlcex, lcex[length(llab)])
    newrfont <- append(newrfont, rfont[length(llab)])
    newlfont <- append(newlfont, lfont[length(llab)])
    newrcol <- append(newrcol, rcol[length(llab)])
    newlcol <- append(newlcol, lcol[length(llab)])
  }

  return (
    list(
      newllab = newllab,
      newrlab = newrlab,
      newrcex = newrcex,
      newlcex = newlcex,
      newrfont = newrfont,
      newlfont = newlfont,
      newrcol = newrcol,
      newlcol = newlcol
    )
  )
}
