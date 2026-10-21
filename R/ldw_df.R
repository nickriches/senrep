#' Apply ldw() across paired vectors of sentences
#'
#' @param sources Character vector of reference sentences.
#' @param targets Character vector of comparison sentences, same length as \code{sources}.
#' @param ... Further arguments passed to \code{ldw()} (e.g. \code{lemmatizer}, \code{greedy_substitutions}).
#' @return A data frame with one row per pair: \code{source}, \code{target},
#'   \code{levdist_in_words}, \code{operations_list} (character), and list-columns
#'   \code{substitutions_list}, \code{omissions_list}, \code{additions_list}.
#' @examples
#' # Basic use with two sentence pairs
#' sources = c("the cat was chased by the dog", "she quickly ran home")
#' targets = c("the dog chased the cat", "she ran home")
#' results_df = ldw_df(sources, targets)
#' results_df$levdist_in_words
#'
#' # Passing extra arguments through to ldw(), e.g. a lemmatizer
#' \dontrun{
#' results_df = ldw_df(sources, targets, lemmatizer = "textstem")
#' }
#'
#' # Accessing list-columns for a single row
#' results_df$substitutions_list[[1]]
#' results_df$omissions_list[[2]]
#'
#' # Typical usage from a dataframe of transcripts
#' transcripts = data.frame(
#'   original = c("the cat was chased by the dog", "she quickly ran home"),
#'   spoken   = c("the dog chased the cat", "she ran home"),
#'   stringsAsFactors = FALSE
#' )
#' out = ldw_df(transcripts$original, transcripts$spoken)
#' out[, c("source", "target", "levdist_in_words")]
#' @export
ldw_df = function(sources, targets, ...) {
  if (length(sources) != length(targets)) {
    stop("sources and targets must be the same length.")
  }

  results = Map(function(s, t) ldw(s, t, ...), sources, targets)

  data.frame(
    source = sources,
    target = targets,
    levdist_in_words = vapply(results, function(r) r$levdist_in_words, numeric(1)),
    operations_list  = vapply(results, function(r) trimws(r$operations_list), character(1)),
    substitutions_list = I(lapply(results, function(r) r$substitutions_list)),
    omissions_list      = I(lapply(results, function(r) r$omissions_list)),
    additions_list      = I(lapply(results, function(r) r$additions_list)),
    stringsAsFactors = FALSE
  )
}
