#' .expand_pattern_scores
#'
#' This helper function takes compressed pattern scores
#' and restores the human-readable gene name table.

#'@param pattern_scores In get_scores_pa_closepairs
#'@param pattern_map In get_scores_pa_closepairs

expand_pattern_scores <- function(pattern_scores, pattern_map) {

  ## Split trait names once rather than repeatedly searching the map
  traits_by_pattern <- split(pattern_map$Trait,pattern_map$Pattern)

  out <- vector("list", nrow(pattern_scores))

  for (i in seq_len(nrow(pattern_scores))) {

    p1 <- pattern_scores$Pattern1[i]
    p2 <- pattern_scores$Pattern2[i]

    g1 <- traits_by_pattern[[as.character(p1)]]
    g2 <- traits_by_pattern[[as.character(p2)]]

    if (p1 == p2) {

      ## Distinct genes sharing the same pattern.
      ## A single-gene pattern produces no gene pair.
      if (length(g1) < 2L) {next}

      pairs <- utils::combn(g1, 2)

      n <- ncol(pairs)

      d <- data.frame(
        Trait1 = pairs[1, ],
        Trait2 = pairs[2, ],
        stringsAsFactors = FALSE
      )

    } else {

      ## Cartesian product of genes belonging to the two patterns
      d <- expand.grid(
        Trait1 = g1,
        Trait2 = g2,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )

      n <- nrow(d)
    }

    ## Add the pattern-pair result to every represented gene pair
    d$PositiveAssociation <- rep(pattern_scores$PositiveAssociation[i], n)

    d$NegativeAssociation <- rep(pattern_scores$NegativeAssociation[i], n)

    d$Score <- rep(pattern_scores$Score[i], n)

    d$UnweightedPositive <- rep(pattern_scores$UnweightedPositive[i], n)

    d$UnweightedNegative <- rep(pattern_scores$UnweightedNegative[i], n)

    if ("sig" %in% names(pattern_scores)) {d$sig <- rep(pattern_scores$sig[i], n)}

    out[[i]] <- d
  }

  out <- out[!vapply(out, is.null, logical(1))]

  if (length(out) == 0L) {
    return(
      data.frame(
        Trait1 = character(),
        Trait2 = character(),
        PositiveAssociation = numeric(),
        NegativeAssociation = numeric(),
        Score = numeric(),
        UnweightedPositive = numeric(),
        UnweightedNegative = numeric(),
        stringsAsFactors = FALSE
      )
    )
  }

  do.call(rbind, out)
}
