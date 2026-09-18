#' .get_scores_pa_closepairs
#'
#' Computes DeCoTUR scores for unique directional discordance patterns
#' across a specified set of close pairs.
#'
#' Traits with identical directional discordance patterns are collapsed
#' before scoring. The returned scores therefore correspond to pairs of
#' unique patterns rather than individual trait pairs.
#'
#' @param pa_matrix Presence-absence matrix; rows are traits and columns
#'   are samples.
#' @param close_pairs Matrix with one close pair per row and two columns
#'   giving sample indices.
#' @param classweights Numeric vector giving one weight per close pair.
#' @param verbose Logical; show progress information.
#'
#' @return A data frame containing scores for pairs of unique
#'   directional discordance patterns. The mapping from traits to patterns
#'   is stored in the `"pattern_map"` attribute.
#'
.get_scores_pa_closepairs <- function(pa_matrix,
    close_pairs, classweights, verbose = FALSE) {

  gene_names <- rownames(pa_matrix)
  if (is.null(gene_names)) {gene_names <- as.character(seq_len(nrow(pa_matrix)))  }

  if (length(classweights) != nrow(close_pairs)) {stop("classweights must have one entry per close pair")  }

  ## ------------------------------------------------------------
  ## Directional changes across close pairs
  ##
  ## D =  1 : 0 -> 1
  ## D = -1 : 1 -> 0
  ## D =  0 : no change
  ## ------------------------------------------------------------

  x1 <- pa_matrix[, close_pairs[, 1], drop = FALSE]
  x2 <- pa_matrix[, close_pairs[, 2], drop = FALSE]

  D <- (x1 == 0 & x2 == 1) - (x1 == 1 & x2 == 0)

  rm(x1, x2)

  ## ------------------------------------------------------------
  ## Collapse identical directional-change patterns
  ## ------------------------------------------------------------

  pattern_key <- apply(D, 1, paste, collapse = ",")

  unique_key <- unique(pattern_key)
  pattern_id <- match(pattern_key, unique_key)

  npatterns <- length(unique_key)

  .verbose_message(
    verbose,
    nrow(pa_matrix), " traits collapse to ",
    npatterns, " unique close-pair patterns ",
    sprintf("(%.1fx compression).", nrow(pa_matrix) / npatterns)
  )

  ## One representative row per unique pattern
  representative <- match(seq_len(npatterns), pattern_id)
  Du <- D[representative, , drop = FALSE]

  ## Mapping from original traits to patterns
  pattern_map <- data.frame(Trait = gene_names,
    Pattern = pattern_id, stringsAsFactors = FALSE)

  rm(D, pattern_key, unique_key)

  ## ------------------------------------------------------------
  ## Score unique patterns
  ## ------------------------------------------------------------

  up <- Du == 1
  down <- Du == -1

  ## Weight columns explicitly: each column represents one close pair
  up_w <- sweep(up, 2, classweights, `*`)
  down_w <- sweep(down, 2, classweights, `*`)

  positive <- up_w %*% t(up) + down_w %*% t(down)

  negative <- up_w %*% t(down) + down_w %*% t(up)

  unweighted_positive <- up %*% t(up) + down %*% t(down)

  unweighted_negative <- up %*% t(down) + down %*% t(up)

  ## Discordance count is a property of the pattern
  discordances <- rowSums(up) + rowSums(down)

  ## Include the diagonal because distinct traits can share a pattern.
  inds <- which(upper.tri(positive, diag = TRUE), arr.ind = TRUE)

  k <- cbind(inds[, 1], inds[, 2])

  pattern_scores <- data.frame(Pattern1 = inds[, 1],
    Pattern2 = inds[, 2], PositiveAssociation = as.numeric(positive[k]),
    NegativeAssociation = as.numeric(negative[k]),
    Score = as.numeric(pmax(positive[k], negative[k])),
    UnweightedPositive = as.numeric(unweighted_positive[k]),
    UnweightedNegative = as.numeric(unweighted_negative[k]),
    Discordance1 = discordances[inds[, 1]],
    Discordance2 = discordances[inds[, 2]],
    stringsAsFactors = FALSE
  )

  attr(pattern_scores, "pattern_map") <- pattern_map

  pattern_scores
}



