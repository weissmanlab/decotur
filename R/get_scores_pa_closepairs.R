#' get_scores_pa_closepairs
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
.get_scores_pa_closepairs <- function(pa_matrix, close_pairs, classweights, verbose = FALSE, chunksize = 5000) {
  gene_names <- rownames(pa_matrix)
  if (is.null(gene_names)) {
    gene_names <- as.character(seq_len(nrow(pa_matrix)))
  }
  if (length(classweights) != nrow(close_pairs)) {
    stop("classweights must have one entry per close pair")
  }
  ntraits <- nrow(pa_matrix)
  npairs  <- nrow(close_pairs)
  ## ------------------------------------------------------------
  ## Build directional-pattern keys in chunks
  ## ------------------------------------------------------------
  .verbose_message(verbose, "Building directional-change pattern keys in chunks...")
  pattern_key <- rep("", ntraits)
  starts <- seq(1, npairs, by = chunksize)
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "  Pattern chunks [:bar] :percent eta: :eta",
      total = length(starts),
      clear = FALSE
    )
    pb$tick(0)
  }
  for (start in starts) { 
    end <- min(start + chunksize - 1, npairs)
    inds <- start:end
    x1 <- pa_matrix[, close_pairs[inds, 1], drop = FALSE]
    x2 <- pa_matrix[, close_pairs[inds, 2], drop = FALSE]
    Dblock <- (x1 == 0 & x2 == 1) - (x1 == 1 & x2 == 0)
    ## Make one string per trait for this chunk
    block_key <- apply(Dblock, 1, paste0, collapse = "," )
    ## Append chunk key to cumulative key
    pattern_key <- paste0(pattern_key, "|", block_key)
    rm(x1, x2, Dblock, block_key)
    if (verbose) {
      pb$tick()
    }
  }
  ## ------------------------------------------------------------
  ## Collapse identical full patterns
  ## ------------------------------------------------------------
  unique_key <- unique(pattern_key)
  pattern_id <- match(pattern_key, unique_key)
  npatterns <- length(unique_key)
  .verbose_message(verbose, ntraits, " traits collapse to ", npatterns, " unique close-pair patterns ", sprintf("(%.1fx compression).", ntraits / npatterns))
  representative <- match(seq_len(npatterns), pattern_id)
  pattern_map <- data.frame(Trait = gene_names, Pattern = pattern_id, stringsAsFactors = FALSE)
  rm(unique_key)
  ## ------------------------------------------------------------
  ## Reconstruct D only for representative traits
  ## ------------------------------------------------------------
  .verbose_message(verbose, "Reconstructing directional-change matrix for ", npatterns, " representative traits..." )
  rep_pa <- pa_matrix[representative, ,drop = FALSE]
  x1 <- rep_pa[,close_pairs[, 1], drop = FALSE]
  x2 <- rep_pa[,close_pairs[, 2], drop = FALSE]
  Du <- (x1 == 0 & x2 == 1) - (x1 == 1 & x2 == 0)
  rm(x1, x2, rep_pa, pattern_key)
  ## ------------------------------------------------------------
  ## Score unique patterns
  ## ------------------------------------------------------------
  up <- Du == 1
  down <- Du == -1
  up_w <- sweep(up,2,classweights,`*`)
  down_w <- sweep(down,2,classweights,`*`)
  positive <- up_w %*% t(up) + down_w %*% t(down)
  negative <- up_w %*% t(down) + down_w %*% t(up)
  unweighted_positive <- up %*% t(up) + down %*% t(down)
  unweighted_negative <- up %*% t(down) + down %*% t(up)
  discordances <- rowSums(up) + rowSums(down)
  ## Include diagonal because distinct traits can share a pattern
  inds <- which(upper.tri(positive, diag = TRUE), arr.ind = TRUE)
  k <- cbind(inds[, 1], inds[, 2])
  pattern_scores <- data.frame(Pattern1 = inds[, 1], Pattern2 = inds[, 2], PositiveAssociation = as.numeric(positive[k]),
    NegativeAssociation = as.numeric(negative[k]), Score = as.numeric(pmax(positive[k], negative[k])),
    UnweightedPositive = as.numeric(unweighted_positive[k]), UnweightedNegative = as.numeric(unweighted_negative[k]),
    Discordance1 = discordances[inds[, 1]], Discordance2 = discordances[inds[, 2]], stringsAsFactors = FALSE)
  attr(pattern_scores, "pattern_map") <- pattern_map
  pattern_scores
}

