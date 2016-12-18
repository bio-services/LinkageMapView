# reverse the positions

revpos <- function(pos, maxdec)
{
  newpos <- vector(length = length(pos))
  newpos[1] <- pos[1]
  for (i in 2:length(pos)) {
    newpos[i] <-
      newpos[(i - 1)] + pos[(length(pos) - i + 2)] - pos[(length(pos) - i + 1)]
  }
  return (round(newpos, digits = maxdec))
}
