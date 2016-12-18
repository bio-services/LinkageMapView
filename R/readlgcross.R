# read cross input from qtl package

readlgcross <- function(cross, mapthese)
{

  if (!requireNamespace("qtl", quietly = TRUE)) {
    stop("Please install package qtl")
  }
  # if list not provided create list of all using chrnames function of qtl
  if (missing(mapthese)) {
    mapthese <- qtl::chrnames(cross)
  }
  else {
    if (!all(mapthese %in% qtl::chrnames(cross))) {
      stop ("chrnames to map not found in cross object")
    }
  }
  df <- qtl::pull.map(cross, mapthese, as.table = TRUE)

  # remove factors if any
  fas <- sapply(df, is.factor)
  df[fas] <- lapply(df[fas], as.character)

  df$locus <- rownames(df)

  # assign my own column names so I can reference by name
  colnames(df) <- c("group", "position", "locus")
  return (df)
}
