#' Resolve a lemmatizer argument into a usable form
#'
#' Internal helper for \code{ldw()}. Not exported.
#'
#' @param model See \code{ldw()} documentation.
#' @return a list with source/target elements.
#' @keywords internal


make_lemmatizer_udpipe = function(model) {
  function(tokens) {
    if (!requireNamespace("udpipe", quietly = TRUE)) {
      warning("udpipe not installed — lemma matching disabled.", call. = FALSE)
      return(tokens)
    }
    ann = udpipe::udpipe_annotate(
      model,
      x = tokens,
      tokenizer = "vertical"   # tells udpipe input is pre-tokenized, one token per element
    )
    ann = as.data.frame(ann)
    if (nrow(ann) != length(tokens)) {
      warning("udpipe token count mismatch — lemma matching disabled for this input.", call. = FALSE)
      return(tokens)
    }
    ann$lemma
  }
}
