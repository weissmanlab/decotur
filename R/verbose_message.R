#' .expand_pattern_scores
#'
#' This helper function takes compressed pattern scores
#' and restores the human-readable gene name table.

#'@param verbose Boolean from user input
#'@param ... Message to print


.verbose_message <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    message(...)
  }
}
