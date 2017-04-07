# convert text font entries to the standard R form:

# "plain text" to 1
# "bold" to 2
# "italic" to 3
# "bold italic" to 4
#

convertfont  <- function (parm, fonttext) {
  # make sure valid
  validfont <-
    c("plain text", "bold", "italic", "bold italic", "1", "2", "3", "4")
  retfont <- vector(length = length(fonttext))

  for (i in 1:length(fonttext)) {
    if (!(fonttext[i] %in% validfont)) {
      message("Invalid font entered for paramenter ",
              parm, ":", fonttext[i])
      if (length(fonttext) > 1) {
        message ("at row number ", i)
      }
      message("Please use one of the following standard R fonts or see pdf.family parameter: ")
      message("1 or 'plain text'")
      message("2 or 'bold'")
      message("3 or 'italic'")
      message("4 or 'bold italic'")
      stop("")
    }
    else {
      if (fonttext[i] == "1" | fonttext[i] == "plain text") {
        retfont[i] <- 1
      }


      else {
        if (fonttext[i] == "2" | fonttext[i] == "bold") {
          retfont[i] <-  2
        }


        else {
          if (fonttext[i] == "3" | fonttext[i] == "italic") {
            retfont[i] <- 3
          }


          else {
            if (fonttext[i] == "4" | fonttext[i] == "bold italic") {
              retfont[i] <- 4
            }
          }
        }
      }
    }
  }

  retfont
}
