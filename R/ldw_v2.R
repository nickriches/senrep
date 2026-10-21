#' Perform a Levenshtein Distance-in-words analysis.
#'
#' @param source The sentence that is being repeated
#' @param target The repetition attempt
#' @return A list. In order to view the list, the result of the function needs to be passed to a variable, e.g. result = ldw("The cat was chased by the dog", "The cat chase the dog"). Then you can access the various different objects in the list.
#'
#' result$levdist_in_words shows the Levenshtein Distance in Words.
#' result$levdist_matrix shows the Levenshtein Distance matrix.
#' result$substitutions_list shows the substitutions, as a character vector
#' result$omissions shows the omissions, as a character vector
#' result$additions shows the additions, as a character vector
#' results$operations_list shows the operations within a sentence.
#'
#' @examples
#' @export




ldw = function(source, target, greedy_substitutions = TRUE){

  source = strsplit(source, " ")[[1]]
  target = strsplit(target, " ")[[1]]

  if(length(source) == 0) return(length(target))
  if(length(target) == 0) return(length(source))

  m <- matrix(0, nrow = length(source) + 1, ncol = length(target) + 1)

  rownames(m) = c("", source)
  colnames(m) = c("", target)

  # Initialize first row and column
  for (i in 1:(length(source) + 1)) {
    m[i, 1] <- i - 1
  }
  for (j in 1:(length(target) + 1)) {
    m[1, j] <- j - 1
  }

  for (sc in 1:(length(source))){
    i = sc + 1
    for (tc in 1:(length(target))){
      j = tc + 1
      cost = ifelse(source[sc] == target[tc], 0, 1)

      m[i, j] = min(m[i-1, j] + 1,      # deletion
                    m[i, j-1] + 1,      # insertion
                    m[i-1, j-1] + cost) # substitution
    }
  }

  ldw = m[i,j]

  substitutions = c()
  additions = c()
  omissions = c()
  op = ""

  # up = m[i-1, j]; left = m[i, j-1]; diagonal = m[i-1, j-1]

  # Extracting operations


  while(i > 1 | j > 1) {
    up       = if (i > 1) m[i-1, j] else Inf
    left     = if (j > 1) m[i, j-1] else Inf
    diagonal = if (i > 1 & j > 1) m[i-1, j-1] else Inf
    cost     = if (i > 1 & j > 1) ifelse(source[i-1] == target[j-1], 0, 1) else Inf

    is_match     = (i > 1 & j > 1) && cost == 0 && diagonal == m[i, j]
    diag_optimal = (i > 1 & j > 1) && (diagonal + cost) == m[i, j]
    up_optimal   = (i > 1) && (up + 1) == m[i, j]
    left_optimal = (j > 1) && (left + 1) == m[i, j]

    if (is_match) {
      # free diagonal move — no operation to record
      op = paste(rownames(m)[i], op)
      i = i - 1; j = j - 1

    } else if (diag_optimal && !up_optimal && !left_optimal) {
      # substitution is the ONLY optimal move — unambiguous
      substitutions = c(substitutions, paste(rownames(m)[i], "->", colnames(m)[j]))
      op = paste0("[", rownames(m)[i], "->", colnames(m)[j], "] ", op)
      i = i - 1; j = j - 1

    } else if (diag_optimal && greedy_substitutions) {
      # tied, but substitution is preferred
      substitutions = c(substitutions, paste(rownames(m)[i], "->", colnames(m)[j]))
      op = paste0("[", rownames(m)[i], "->", colnames(m)[j], "] ", op)
      i = i - 1; j = j - 1

    } else if (up_optimal) {
      # tied (or best), and substitution deprioritised — take omission
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



  return(invisible
         (list
           (levdist_in_words = ldw,
             levdist_matrix = m,
             substitutions_list = substitutions,
             omissions_list = omissions,
             additions_list = additions,
             operations_list = op)
         )
  )

}





reverse = function(x){
  x = strsplit(x, " ")[[1]]
  x = rev(x)
  x = paste(x, collapse = " ")
  return(x)
}
