#' Levenshtein Distance in Words with Human-like Alignment
#'
#' Computes word-level Levenshtein distance between two strings and
#' backtraces the optimal alignment into substitutions, omissions, and
#' additions. Ties between operations are resolved using an optional
#' lemmatizer (for linguistically-informed substitution preference),
#' falling back to character similarity, then to the \code{greedy_substitutions}
#' flag.
#'
#' @param source Character string, the reference/original sentence.
#' @param target Character string, the comparison sentence.
#' @param greedy_substitutions Logical. When ties can't be resolved by lemma
#'   match or character similarity, prefer substitution (TRUE) or
#'   defer to omission/addition (FALSE, default).
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
#' @return An (invisible) list with elements \code{levdist_in_words} (integer
#'   distance), \code{levdist_matrix} (the DP matrix), \code{substitutions_list},
#'   \code{omissions_list}, \code{additions_list}, and \code{operations_list}
#'   (a single annotated string).
#' @examples
#' ldw("the cat was chased by the dog", "the dog chased the cat")
#' ldw("the cat was chased by the dog", "the dog chased the cat",
#'     lemmatizer = "textstem")
#' @export
ldw = function(source, target,
               greedy_substitutions = TRUE,
               sim_threshold = 1,
               lemmatizer = NULL) {

  source = strsplit(source, " ")[[1]]
  target = strsplit(target, " ")[[1]]
  if (length(source) == 0) return(length(target))
  if (length(target) == 0) return(length(source))

  # --- Resolve lemmatizer argument into usable form ---
  lem = resolve_lemmatizer(lemmatizer)

  if (is.null(lem)) {
    source_lemmas = NULL
    target_lemmas = NULL
  } else if (is.function(lem)) {
    source_lemmas = tryCatch(lem(source), error = function(e) source)
    target_lemmas = tryCatch(lem(target), error = function(e) target)
  } else if (is.list(lem)) {
    if (length(lem$source) != length(source) || length(lem$target) != length(target)) {
      stop("Precomputed lemma vectors must match token counts of source/target.")
    }
    source_lemmas = lem$source
    target_lemmas = lem$target
  }

  # --- Build DP matrix ---
  m <- matrix(0, nrow = length(source) + 1, ncol = length(target) + 1)
  rownames(m) = c("", source)
  colnames(m) = c("", target)
  for (i in 1:(length(source) + 1)) m[i, 1] <- i - 1
  for (j in 1:(length(target) + 1)) m[1, j] <- j - 1

  for (sc in 1:(length(source))) {
    i = sc + 1
    for (tc in 1:(length(target))) {
      j = tc + 1
      cost = ifelse(source[sc] == target[tc], 0, 1)
      m[i, j] = min(m[i-1, j] + 1,
                    m[i, j-1] + 1,
                    m[i-1, j-1] + cost)
    }
  }

  ldw_val = m[i, j]
  substitutions = c(); additions = c(); omissions = c()
  op = ""

  sim_threshold = 0.4  # tune against real examples

  while (i > 1 | j > 1) {
    up       = if (i > 1) m[i-1, j] else Inf
    left     = if (j > 1) m[i, j-1] else Inf
    diagonal = if (i > 1 & j > 1) m[i-1, j-1] else Inf
    cost     = if (i > 1 & j > 1) ifelse(source[i-1] == target[j-1], 0, 1) else Inf

    is_match     = (i > 1 & j > 1) && cost == 0 && diagonal == m[i, j]
    diag_optimal = (i > 1 & j > 1) && (diagonal + cost) == m[i, j]
    up_optimal   = (i > 1) && (up + 1) == m[i, j]
    left_optimal = (j > 1) && (left + 1) == m[i, j]

    if (is_match) {
      op = paste(rownames(m)[i], op)
      i = i - 1; j = j - 1

    } else if (diag_optimal && !up_optimal && !left_optimal) {
      substitutions = c(substitutions, paste(rownames(m)[i], "->", colnames(m)[j]))
      op = paste0("[", rownames(m)[i], "->", colnames(m)[j], "] ", op)
      i = i - 1; j = j - 1

    } else if (diag_optimal) {
      w1 = rownames(m)[i]; w2 = colnames(m)[j]

      lemma_hit = !is.null(source_lemmas) && source_lemmas[i-1] == target_lemmas[j-1]
      sim = tryCatch({
        d = adist(w1, w2)[1, 1]
        1 - d / max(nchar(w1), nchar(w2))
      }, error = function(e) 0)

      prefer_sub = if (lemma_hit) TRUE else if (sim >= sim_threshold) TRUE else greedy_substitutions

      if (prefer_sub) {
        substitutions = c(substitutions, paste(w1, "->", w2))
        op = paste0("[", w1, "->", w2, "] ", op)
        i = i - 1; j = j - 1
      } else if (up_optimal) {
        omissions = c(omissions, rownames(m)[i])
        op = paste0("*", rownames(m)[i], "* ", op)
        i = i - 1
      } else if (left_optimal) {
        additions = c(additions, colnames(m)[j])
        op = paste0("+", colnames(m)[j], "+ ", op)
        j = j - 1
      }

    } else if (up_optimal) {
      omissions = c(omissions, rownames(m)[i])
      op = paste0("*", rownames(m)[i], "* ", op)
      i = i - 1

    } else if (left_optimal) {
      additions = c(additions, colnames(m)[j])
      op = paste0("+", colnames(m)[j], "+ ", op)
      j = j - 1
    }
  }

  substitutions = rev(substitutions)
  omissions = rev(omissions)
  additions = rev(additions)

  # return(invisible(list(
  #   levdist_in_words = ldw_val,
  #   levdist_matrix = m,
  #   substitutions_list = substitutions,
  #   omissions_list = omissions,
  #   additions_list = additions,
  #   operations_list = op
  # )))

  # at the end of ldw(), replace the final return with:
  result = list(
    levdist_in_words = ldw_val,
    levdist_matrix = m,
    substitutions_list = substitutions,
    omissions_list = omissions,
    additions_list = additions,
    operations_list = op
  )

  class(result) = "ldw"
  return(result)


}
