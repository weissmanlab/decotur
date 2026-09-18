#' get_scores
#'
#' Computes DeCoTUR scores and identifies significantly associated trait pairs.
#'
#' This function identifies close pairs of samples, computes close-pair weights,
#' collapses traits with identical directional discordance patterns, calculates
#' DeCoTUR scores and significance at the pattern level, and returns the
#' significant associations expanded to individual trait pairs.
#'
#' @param pa_matrix A presence-absence matrix. Rows are traits and columns are
#'   samples. Row names should contain trait names.
#' @param distance_matrix A pairwise distance matrix. Row and column order
#'   should correspond to the columns of `pa_matrix`.
#' @param closepair_method Method used to select close pairs. One of
#'   `"distance"`, `"fraction"`, or `"fixednumber"`.
#' @param closepair_params A list of parameters for the selected close-pair
#'   method. For `"distance"`, supply the distance cutoff and `show_hist`.
#'   For `"fraction"`, supply the distance fraction and `show_hist`.
#'   For `"fixednumber"`, supply the number of close pairs, random seed,
#'   and `show_hist`.
#' @param downweight Logical; whether to downweight close-pair classes.
#' @param verbose Logical; whether to display progress messages.
#'
#' @return A data frame containing significantly associated trait pairs and
#'   their DeCoTUR scores.
#'
#' @export
get_scores <- function(
    pa_matrix,
    distance_matrix,
    closepair_method,
    closepair_params,
    downweight = TRUE,
    verbose = FALSE
) {

  .verbose_message(verbose, "Starting analysis.")

  ## ------------------------------------------------------------
  ## 1. Obtain close pairs
  ## ------------------------------------------------------------

  if (closepair_method == "distance") {

    if (length(closepair_params) != 2) {
      stop(
        "Incorrect number of closepair_params for method 'distance'. ",
        "Expected: distance cutoff and show_hist."
      )
    }

    close_pairs <- get_closepairs_distance(
      distance_matrix,
      closepair_params[[1]],
      closepair_params[[2]],
      verbose
    )

  } else if (closepair_method == "fraction") {

    if (length(closepair_params) != 2) {
      stop(
        "Incorrect number of closepair_params for method 'fraction'. ",
        "Expected: distance fraction and show_hist."
      )
    }

    close_pairs <- get_closepairs_fraction(
      distance_matrix,
      closepair_params[[1]],
      closepair_params[[2]],
      verbose
    )

  } else if (closepair_method == "fixednumber") {

    if (length(closepair_params) != 3) {
      stop(
        "Incorrect number of closepair_params for method 'fixednumber'. ",
        "Expected: number of close pairs, seed, and show_hist."
      )
    }

    close_pairs <- get_closepairs_fixednumber(
      distance_matrix,
      closepair_params[[1]],
      closepair_params[[2]],
      closepair_params[[3]],
      verbose
    )

  } else {

    stop(
      "Unidentified close-pair method. Use 'distance', ",
      "'fraction', or 'fixednumber'."
    )
  }

  .verbose_message(
    verbose,
    "Obtained ", nrow(close_pairs), " close pairs."
  )

  ## ------------------------------------------------------------
  ## 2. Save close-pair distances BEFORE reindexing
  ## ------------------------------------------------------------

  # close_pairs currently indexes the original distance matrix.
  # filter_traits_by_closepairs() subsequently reindexes these pairs
  # relative to the subsetted presence-absence matrix.
  dpds <- as.numeric(distance_matrix[close_pairs])

  ## ------------------------------------------------------------
  ## 3. Calculate close-pair weights
  ## ------------------------------------------------------------

  classes <- get_closepair_classes(close_pairs)

  if (downweight) {
    class_weights <- get_classweights(classes)
  } else {
    class_weights <- rep(1, nrow(close_pairs))
  }

  .verbose_message(verbose, "Obtained close-pair weights.")

  ## ------------------------------------------------------------
  ## 4. Filter traits and samples
  ## ------------------------------------------------------------

  pam <- filter_traits_by_closepairs(
    pa_matrix,
    close_pairs
  )

  pa_matrix <- pam[[1]]
  close_pairs <- pam[[2]]

  n_traits <- nrow(pa_matrix)
  n_tests <- choose(n_traits, 2)

  .verbose_message(
    verbose,
    "Retained ", n_traits, " traits representing ",
    format(n_tests, big.mark = ","),
    " trait-pair tests."
  )

  ## ------------------------------------------------------------
  ## 5. Compute compressed pattern scores
  ## ------------------------------------------------------------

  scores <- .get_scores_pa_closepairs(
    pa_matrix = pa_matrix,
    close_pairs = close_pairs,
    classweights = class_weights,
    verbose = verbose
  )

  pattern_map <- attr(scores, "pattern_map")

  ## ------------------------------------------------------------
  ## 6. Prepare close-pair distances for null model
  ## ------------------------------------------------------------

  dpds <- as.numeric(distance_matrix[close_pairs])

  if (any(!is.finite(dpds))) {stop("Close-pair distances contain non-finite values.")}

  if (any(dpds < 0)) {
    .verbose_message(verbose, sum(dpds < 0), " negative close-pair distances were set to zero.")
    dpds[dpds < 0] <- 0
  }
    
  if (sum(dpds) == 0) {stop("Close-pair distances must contain at least one positive value.")}

  ndpds <- dpds / sum(dpds)

  ## ------------------------------------------------------------
  ## 7. Compute significance at the pattern level
  ## ------------------------------------------------------------

  .verbose_message(verbose, "Computing significance.")

  # The null probability vector depends on a trait pair only through
  # the product of its two discordance counts. Therefore, calculate
  # the critical value only once for each unique product.
  dprod <-
    scores$Discordance1 *
    scores$Discordance2

  unique_dprod <- sort(unique(dprod))

  q <- ndpds^2

  # Bonferroni correction is based on the number of actual trait-pair
  # hypotheses, not the number of compressed pattern pairs.
  alpha <- 0.05 / n_tests

  critical_values <- numeric(length(unique_dprod))

  for (i in seq_along(unique_dprod)) {

    probs <- unique_dprod[i] * q

    probs[!is.finite(probs)] <- 0
    probs[probs < 0] <- 0
    probs[probs > 1] <- 1

    critical_values[i] <- qpbinom_modified(
      1 - alpha,
      probs,
      method = "RefinedNormal"
    )
  }

  scores$CriticalValue <-
    critical_values[match(dprod, unique_dprod)]

  scores$Observed <-
    scores$UnweightedPositive +
    scores$UnweightedNegative

  scores$sig <-
    scores$Observed > scores$CriticalValue

  .verbose_message(
    verbose,
    "Found ",
    sum(scores$sig),
    " significant close-pair pattern combinations."
  )

  ## ------------------------------------------------------------
  ## 8. Keep significant patterns and expand to trait pairs
  ## ------------------------------------------------------------

  significant_scores <- scores[
    scores$sig,
    ,
    drop = FALSE
  ]

  if (nrow(significant_scores) == 0) {

    .verbose_message(
      verbose,
      "No significant trait associations found."
    )

    return(
      data.frame(
        Trait1 = character(),
        Trait2 = character(),
        PositiveAssociation = numeric(),
        NegativeAssociation = numeric(),
        Score = numeric(),
        UnweightedPositive = numeric(),
        UnweightedNegative = numeric(),
        sig = logical(),
        stringsAsFactors = FALSE
      )
    )
  }

  scores <- .expand_pattern_scores(
    significant_scores,
    pattern_map
  )

  .verbose_message(
    verbose,
    "Returning ",
    format(nrow(scores), big.mark = ","),
    " significant trait pairs."
  )

  scores
}
