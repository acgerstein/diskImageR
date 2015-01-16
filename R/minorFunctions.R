#' Used to determine whether something is a letter

#' @param x a single element

#' @export

is.letter <-
function(x) grepl("[[:alpha:]]", x)

#' Used to determine whether something is a number

#' @param x a single element

#' @export

is.number <-
function(x) grepl("[[:digit:]]", x)
