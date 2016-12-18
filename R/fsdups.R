# find duplicate locations
# input is:
#     vector of y plot points sorted
#     maxnbrcolsfordups

# returns
#  lkeep - vector of left label indexes - nondups
#  rkeep - vector of right label indexes - nondups but possible you
#          may need more than one y axis label location
#          for all dups if more than maxnbrcolsfordups
#  frkeep - vector of right label indexes for the first marker of dups
#  fykeep - vector of right label indexes for subsequent first column dups
#           that is if 5 dups, and maxnbrcolfordups is 3, then the 4th
#           dup is in the first column of the second row of dups
#  flkeep - adjyr index for left labels - left labels and right labels
#           for the same position must be across from each other.
#  yd    - provide matrix of index entries for where to place dup labels


fsdups <- function(y, maxnbrcolsfordups) {
  k <- length(y)   # original length of y

  # For storing dup y locations
  yd <- matrix(nrow = maxnbrcolsfordups - 1, ncol = length(y))
  j <- 0           # dup counter

  lkeep <-
    vector(mode = "integer")  # keep track of which to keep labels for
  rkeep <-
    vector(mode = "integer")  # the rest go in columns no need to spread

  # keep a vector of indexes for the first column of duplicates
  # where more than one row of duplicates

  frkeep <- vector(mode = "integer")
  fykeep <- vector(mode = "integer")
  flkeep <- vector(mode = "integer")

  #always keep first one
  lkeep <- append(lkeep, 1)
  rkeep <- append(rkeep, 1)
  flkeep <- append(flkeep, 1)
  if (k > 1) {
    for (dupi in 2:k) {
      if (y[dupi] == y[dupi - 1]) {
        if (j == 0) {
          frkeep <- append(frkeep, dupi - 1)
          fykeep <- append(fykeep, length(rkeep))
        }
        if (j < maxnbrcolsfordups - 1) {
          j <- j + 1
          yd[j, (dupi - j)] <-
            length(rkeep)  # save what will be adjyr index for dup
        }
        else {
          j <- 0
          rkeep <- append(rkeep, dupi)
          frkeep <- append(frkeep, dupi)
          fykeep <-
            append(fykeep, length(rkeep))  # first one in subsequent row for dup
        }
      }
      else {
        if (j == 0 && dupi > 2 && y[(dupi - 1)] == y[(dupi - 2)]) {
          frkeep <- append(frkeep, dupi - 1)
          fykeep <- append(fykeep, length(rkeep))
        }
        j <- 0
        lkeep <- append(lkeep, dupi)
        rkeep <- append(rkeep, dupi)
        flkeep <-
          append(flkeep, length(rkeep))  # adjyr index for left labels
      }
    }
  }
  list(
    lkeep = lkeep[!duplicated(lkeep)],
    rkeep = rkeep[!duplicated(rkeep)],
    frkeep = frkeep[!duplicated(frkeep)],
    fykeep = fykeep[!duplicated(fykeep)],
    flkeep = flkeep[!duplicated(flkeep)],
    yd = yd
  )
}
