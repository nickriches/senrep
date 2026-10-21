
library(textstem)

ldw = function(source, target, greedy_substitutions = TRUE){
  source = strsplit(source, " ")[[1]]
  target = strsplit(target, " ")[[1]]
  if(length(source) == 0) return(length(target))
  if(length(target) == 0) return(length(source))
  m <- matrix(0, nrow = length(source) + 1, ncol = length(target) + 1)
  rownames(m) = c("", source)
  colnames(m) = c("", target)
  for (i in 1:(length(source) + 1)) m[i, 1] <- i - 1
  for (j in 1:(length(target) + 1)) m[1, j] <- j - 1
  for (sc in 1:(length(source))){
    i = sc + 1
    for (tc in 1:(length(target))){
      j = tc + 1
      cost = ifelse(source[sc] == target[tc], 0, 1)
      m[i, j] = min(m[i-1, j] + 1,
                    m[i, j-1] + 1,
                    m[i-1, j-1] + cost)
    }
  }
  ldw = m[i,j]
  substitutions = c(); additions = c(); omissions = c()
  op = ""

  # Helper: word similarity via normalized edit distance
  word_sim = function(a, b) {
    d = adist(a, b)[1,1]
    1 - d / max(nchar(a), nchar(b))
  }

  # Helper: lemma match check
  lemma_match = function(a, b) {
    la = lemmatize_words(tolower(a))
    lb = lemmatize_words(tolower(b))
    la == lb
  }

  SIM_THRESHOLD = 0.4  # tune against real examples

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
      op = paste(rownames(m)[i], op)
      i = i-1; j = j-1

    } else if (diag_optimal && !up_optimal && !left_optimal) {
      # unambiguous substitution
      substitutions = c(substitutions, paste(rownames(m)[i], "->", colnames(m)[j]))
      op = paste0("[", rownames(m)[i], "->", colnames(m)[j], "] ", op)
      i = i-1; j = j-1

    } else if (diag_optimal) {
      # genuine tie — resolve via lemma, then similarity, then flag
      w1 = rownames(m)[i]; w2 = colnames(m)[j]

      lemma_hit = tryCatch(lemma_match(w1, w2), error = function(e) FALSE)
      sim = tryCatch(word_sim(w1, w2), error = function(e) 0)

      prefer_sub = if (lemma_hit) {
        TRUE
      } else if (sim >= SIM_THRESHOLD) {
        TRUE
      } else {
        greedy_substitutions
      }

      if (prefer_sub) {
        substitutions = c(substitutions, paste(w1, "->", w2))
        op = paste0("[", w1, "->", w2, "] ", op)
        i = i-1; j = j-1
      } else if (up_optimal) {
        omissions = c(omissions, rownames(m)[i])
        op = paste0("*", rownames(m)[i], "* ", op)
        i = i-1
      } else if (left_optimal) {
        additions = c(additions, colnames(m)[j])
        op = paste0("+", colnames(m)[j], "+ ", op)
        j = j-1
      }

    } else if (up_optimal) {
      omissions = c(omissions, rownames(m)[i])
      op = paste0("*", rownames(m)[i], "* ", op)
      i = i-1

    } else if (left_optimal) {
      additions = c(additions, colnames(m)[j])
      op = paste0("+", colnames(m)[j], "+ ", op)
      j = j-1
    }
  }

  substitutions = rev(substitutions)
  omissions = rev(omissions)
  additions = rev(additions)
  return(invisible(list(
    levdist_in_words = ldw,
    levdist_matrix = m,
    substitutions_list = substitutions,
    omissions_list = omissions,
    additions_list = additions,
    operations_list = op
  )))
}
