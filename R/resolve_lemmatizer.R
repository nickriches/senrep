#' Resolve a lemmatizer argument into a usable form
#'
#' Internal helper for \code{ldw()}. Not exported.
#'
#' @param lemmatizer See \code{ldw()} documentation.
#' @return NULL, a function, or a list with source/target elements.
#' @keywords internal
resolve_lemmatizer = function(lemmatizer) {

  if (is.null(lemmatizer)) return(NULL)

  if (is.function(lemmatizer)) return(lemmatizer)

  if (is.character(lemmatizer) && length(lemmatizer) == 1) {
    if (lemmatizer == "textstem") {
      return(function(tokens) {
        if (!requireNamespace("textstem", quietly = TRUE)) {
          warning("textstem not installed — lemma matching disabled.", call. = FALSE)
          return(tokens)
        }
        textstem::lemmatize_words(tolower(tokens))
      })
    }
    stop("Unrecognised lemmatizer string: '", lemmatizer,
         "'. Use 'textstem', a function, NULL, or a precomputed lemma list.")
  }

  if (is.list(lemmatizer) && all(c("source", "target") %in% names(lemmatizer))) {
    return(lemmatizer)
  }

  stop("lemmatizer must be NULL, a function, the string 'textstem', ",
       "or a list with 'source' and 'target' lemma vectors.")
}


