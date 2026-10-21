#' Levenshtein Distance in Words with Human-like Alignment
#'
#' Computes word-level Levenshtein distance between two strings and
#' backtraces the optimal alignment into substitutions, omissions, and
#' additions.
#'
#' Levenshtein distance itself only defines the *cost* of the cheapest
#' alignment; when several alignments share that cost, which one is reported
#' depends on how ties are broken during the backtrace. \code{ldw()} breaks
#' ties between a substitution and an omission/addition in this order:
#' \enumerate{
#'   \item lemma match (if a \code{lemmatizer} is supplied): prefer substitution;
#'   \item character similarity at or above \code{sim_threshold}: prefer
#'     substitution;
#'   \item otherwise, follow \code{greedy_substitutions}.
#' }
#' The default (\code{greedy_substitutions = FALSE}) means two unrelated words
#' are only aligned as a substitution when no equally cheap alignment exists,
#' which tends to match how a human would score a sentence repetition.
#'
#' @param source Character string, the reference/original sentence.
#' @param target Character string, the comparison sentence (e.g. the
#'   child's repetition).
#' @param greedy_substitutions Logical. When a tie can't be resolved by lemma
#'   match or character similarity, prefer substitution (\code{TRUE}, the
#'   conventional diagonal-first backtrace) or omission/addition
#'   (\code{FALSE}, default).
#' @param sim_threshold Numeric between 0 and 1. Minimum normalised character
#'   similarity (1 - edit distance / length of the longer word) at which a
#'   tied pair of words is treated as a substitution. Default 0.4.
#' @param lemmatizer Controls lemma-based tie-breaking. One of:
#'   \itemize{
#'     \item \code{NULL} (default): lemma tier disabled entirely.
#'     \item \code{"textstem"}: uses \code{textstem::lemmatize_words()}.
#'       Requires the \pkg{textstem} package.
#'     \item A function taking a character vector of tokens and returning a
#'       character vector of lemmas of the same length and order (e.g. built
#'       via a udpipe-based wrapper for non-English languages).
#'     \item A list with named elements \code{source} and \code{target}, each
#'       a precomputed character vector of lemmas matching the token count
#'       and order of the split \code{source}/\code{target} strings. Useful
#'       for passing in output from \pkg{spacyr} or other external taggers.
#'   }
#' @return An object of class \code{"ldw"}: a list with elements
#'   \code{levdist_in_words} (integer distance), \code{levdist_matrix} (the
#'   DP matrix, numeric), \code{levdist_path} (the same matrix for display,
#'   with each cell on the chosen backtrace prefixed by \code{*}),
#'   \code{backtrace_path} (a logical matrix, \code{TRUE} for cells on the
#'   path), \code{substitutions_list},
#'   \code{omissions_list},
#'   \code{additions_list}, and \code{operations_list} (a single annotated
#'   string, with substitutions as \code{[a->b]}, omissions as \code{*a*} and
#'   additions as \code{+a+}).
#' @examples
#' ldw("the cat was chased by the dog", "the dog chased the cat")
#' ldw("the cat was chased by the dog", "the dog chased the cat",
#'     greedy_substitutions = TRUE)
#' \dontrun{
#' ldw("the cat was chased by the dog", "the dog chased the cat",
#'     lemmatizer = "textstem")
#' }
#' @export
ldw = function(source, target,
               greedy_substitutions = FALSE,
               sim_threshold = 0.4,
               lemmatizer = NULL) {

  # --- Validate inputs ---
  if (!is.character(source) || length(source) != 1 || is.na(source)) {
    stop("`source` must be a single, non-missing character string.")
  }
  if (!is.character(target) || length(target) != 1 || is.na(target)) {
    stop("`target` must be a single, non-missing character string.")
  }
  if (!is.logical(greedy_substitutions) || length(greedy_substitutions) != 1 ||
      is.na(greedy_substitutions)) {
    stop("`greedy_substitutions` must be TRUE or FALSE.")
  }
  if (!is.numeric(sim_threshold) || length(sim_threshold) != 1 ||
      is.na(sim_threshold) || sim_threshold < 0 || sim_threshold > 1) {
    stop("`sim_threshold` must be a single number between 0 and 1.")
  }

  # --- Tokenise (tolerant of leading/trailing/multiple spaces) ---
  tokenise = function(x) {
    x = trimws(x)
    if (x == "") character(0) else strsplit(x, "\\s+")[[1]]
  }
  source = tokenise(source)
  target = tokenise(target)
  n_s = length(source)
  n_t = length(target)

  # --- Resolve lemmatizer argument into usable form ---
  lem = resolve_lemmatizer(lemmatizer)
  source_lemmas = NULL
  target_lemmas = NULL

  if (is.function(lem)) {
    source_lemmas = tryCatch(lem(source), error = function(e) source)
    target_lemmas = tryCatch(lem(target), error = function(e) target)
  } else if (is.list(lem)) {
    if (length(lem$source) != n_s || length(lem$target) != n_t) {
      stop("Precomputed lemma vectors must match token counts of source/target.")
    }
    source_lemmas = lem$source
    target_lemmas = lem$target
  }

  # --- Build DP matrix ---
  m = matrix(0, nrow = n_s + 1, ncol = n_t + 1)
  rownames(m) = c("", source)
  colnames(m) = c("", target)
  m[, 1] = 0:n_s
  m[1, ] = 0:n_t

  for (sc in seq_len(n_s)) {
    for (tc in seq_len(n_t)) {
      cost = if (source[sc] == target[tc]) 0 else 1
      m[sc + 1, tc + 1] = min(m[sc, tc + 1] + 1,   # omission
                              m[sc + 1, tc] + 1,   # addition
                              m[sc, tc] + cost)    # match / substitution
    }
  }

  i = n_s + 1
  j = n_t + 1
  ldw_val = as.integer(m[i, j])

  # --- Backtrace ---
  substitutions = character(0)
  omissions = character(0)
  additions = character(0)
  op = ""

  # Cells visited by the backtrace (start at bottom-right, end at top-left)
  on_path = matrix(FALSE, nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m))
  on_path[i, j] = TRUE

  char_sim = function(w1, w2) {
    1 - utils::adist(w1, w2)[1, 1] / max(nchar(w1), nchar(w2))
  }

  while (i > 1 || j > 1) {
    can_diag = i > 1 && j > 1
    cost = if (can_diag) (if (source[i - 1] == target[j - 1]) 0 else 1) else Inf

    diag_optimal = can_diag && (m[i - 1, j - 1] + cost) == m[i, j]
    up_optimal   = i > 1 && (m[i - 1, j] + 1) == m[i, j]
    left_optimal = j > 1 && (m[i, j - 1] + 1) == m[i, j]

    if (diag_optimal && cost == 0) {
      move = "match"
    } else if (diag_optimal && !up_optimal && !left_optimal) {
      move = "sub"
    } else if (diag_optimal) {
      # Tie between substitution and omission/addition
      w1 = source[i - 1]
      w2 = target[j - 1]
      lemma_hit = !is.null(source_lemmas) &&
        isTRUE(source_lemmas[i - 1] == target_lemmas[j - 1])
      prefer_sub = lemma_hit ||
        char_sim(w1, w2) >= sim_threshold ||
        greedy_substitutions
      move = if (prefer_sub) "sub" else if (up_optimal) "omit" else "add"
    } else if (up_optimal) {
      move = "omit"
    } else {
      move = "add"
    }

    if (move == "match") {
      op = paste0(source[i - 1], " ", op)
      i = i - 1; j = j - 1
    } else if (move == "sub") {
      substitutions = c(substitutions, paste(source[i - 1], "->", target[j - 1]))
      op = paste0("[", source[i - 1], "->", target[j - 1], "] ", op)
      i = i - 1; j = j - 1
    } else if (move == "omit") {
      omissions = c(omissions, source[i - 1])
      op = paste0("*", source[i - 1], "* ", op)
      i = i - 1
    } else {
      additions = c(additions, target[j - 1])
      op = paste0("+", target[j - 1], "+ ", op)
      j = j - 1
    }
    on_path[i, j] = TRUE
  }

  # Display copy of m with the backtrace marked by "*" (numbers padded to a
  # common width so columns line up)
  m_path = matrix(paste0(ifelse(on_path, "*", " "),
                         formatC(m, width = max(nchar(m)), format = "d")),
                  nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m))
  m_path = noquote(m_path)

  result = list(
    levdist_in_words   = ldw_val,
    levdist_matrix     = m,
    levdist_path       = m_path,
    backtrace_path     = on_path,
    substitutions_list = rev(substitutions),
    omissions_list     = rev(omissions),
    additions_list     = rev(additions),
    operations_list    = trimws(op)
  )
  class(result) = "ldw"
  result
}
