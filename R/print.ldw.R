#' Print method for ldw objects
#'
#' @param x An object of class \code{ldw}.
#' @param ... Further arguments (unused).
#' @export
#' @method print ldw
print.ldw = function(x, ...) {
  cat("Levenshtein distance in words:", x$levdist_in_words, "\n\n")
  cat("Alignment:", trimws(x$operations_list), "\n\n")
  if (length(x$substitutions_list) > 0)
    cat("Substitutions:", paste(x$substitutions_list, collapse = "; "), "\n")
  if (length(x$omissions_list) > 0)
    cat("Omissions:    ", paste(x$omissions_list, collapse = "; "), "\n")
  if (length(x$additions_list) > 0)
    cat("Additions:    ", paste(x$additions_list, collapse = "; "), "\n")
  cat("\n👍 Tip: assign the result to a variable, e.g. x <- ldw(...))\n\nThen you can access the components directly, e.g. x$substitutions_list, x$levdist_matrix, etc.\n", sep = "")
  invisible(x)
}

