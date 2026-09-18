#' get_scores
#'
#' This function takes a presence-absence matrix, a distance matrix, a choice of
#' how to determine close pairs, and a list of relevant parameters for that choice,
#' and returns the close pair scores. There is an option for whether or not to
#' downweight bushes. There is also an option to return significance levels.
#' Finally there is an option to run (relatively) speed-optimized or memory-optimized versions of
#' this process.
#'
#'@param pa_matrix A presence-absence matrix for the trait in question (i.e. gene, allele, phenotype, etc.). Rows are trait, columns are sample. Rownames should be trait names, colnames should be sample names.
#'@param distance_matrix A distance matrix. Rownames and colnames should be sample names and should be in the same order as the columns in pa_matrix.
#'@param closepair_method Either 'distance', 'fraction', or 'fixednumber'.
#'@param closepair_params A list. For 'distance', the distance cutoff and show_hist. For 'fraction', the fraction and show_hist. For 'fixednumber', the number of close pairs, a random seed, and show_hist. For 'auto', the which_valley, nbins, maxvalleyheight, and show_hist. verbose is inherited from the main function call.
#'@param blocksize The number of traits to consider at once in the computation. Changing this number may help speed up the computation. Must be less than or equal to the number of traits.
#'@param downweight TRUE to downweight bushes. Default TRUE.
#'@param withsig TRUE to return significance. Default TRUE.
#'@param verbose TRUE for progress reports. Default TRUE.
#'@param version 'speed' for speed-optimized, 'memory' for memory optimized. Default 'speed'; if you run into memory problems, try decreasing blocksize first.
#'@export
#'@examples
#'get_scores(pa_matrix, distance_matrix, 'fraction', list(0.1, TRUE), 10, TRUE, TRUE, TRUE)

get_scores <- function(pa_matrix, distance_matrix, closepair_method, closepair_params, blocksize, downweight = TRUE, withsig = TRUE, verbose = TRUE, version = 'speed'){
  if(verbose){print('Starting function.')}
  if(closepair_method == 'distance'){
    if(length(closepair_params) != 2){
      stop('Incorrect number of closepair_params. (Should be 2).')
    }
    distance_cutoff <- closepair_params[[1]]
    show_hist <- closepair_params[[2]]
    close_pairs <- get_closepairs_distance(distance_matrix, distance_cutoff, show_hist, verbose)
  } else if(closepair_method == 'fraction'){
    if(length(closepair_params) != 2){
      stop('Incorrect number of closepair_params. (Should be 2).')
    }
    distance_fraction <- closepair_params[[1]]
    show_hist <- closepair_params[[2]]
    close_pairs <- get_closepairs_fraction(distance_matrix, distance_fraction, show_hist, verbose)
  } else if(closepair_method == 'fixednumber'){
    number_close_pairs <- closepair_params[[1]]
    seed <- closepair_params[[2]]
    show_hist <- closepair_params[[3]]
    close_pairs <- get_closepairs_fixednumber(distance_matrix, number_close_pairs, seed, show_hist, verbose)
  } else{stop('Unidentified close pair method (should be one of: distance, fraction, fixednumber).')}
  if(verbose){
    print(paste0('Obtained ',  dim(close_pairs)[1], ' close pairs.'))
  }
  classes <- get_closepair_classes(close_pairs)
  if(verbose){print('Obtained close pair classes.')}
  if(downweight){
    class_weights <- get_classweights(classes)
  } else{
    class_weights <- rep(1, length(classes))
  }
  if(verbose){print('Obtained close pair class weights.')}
  ## close_pairs still contains indices into the original distance matrix
  dpds <- as.numeric(distance_matrix[close_pairs])
  pam <- filter_traits_by_closepairs(pa_matrix, close_pairs)
  pa_matrix <- pam[[1]]
  close_pairs <- pam[[2]]
  scores <- get_scores_pa_closepairs(pa_matrix, close_pairs, class_weights, blocksize, verbose, version, withsig)
  if(verbose){print('Obtained scores.')}
  if(withsig){
    if(verbose){print('Computing discordance information.')}
    discdat <- get_discordances(pa_matrix, close_pairs, verbose) # need to fix this function to get the right row_names
    if(verbose){print('Obtained discordance information.')}
    if(verbose){print('Computing significance')}
    dpds[which(dpds < 0)] <- 1/(2*1000000)# hardcoded. apologies. shouldn't really matter.
    disc1 <- discdat$disc[match(scores$Trait1, discdat$trait)]
    disc2 <- discdat$disc[match(scores$Trait2, discdat$trait)]
    sum <- scores$UnweightedPositive + scores$UnweightedNegative
    m <- rep(1, dim(close_pairs)[1])
    ndpds <- dpds/sum(dpds)
    ndpds[which(ndpds == 0)] <- min(ndpds[which(ndpds > 0)])/100
    # The null model is a Poisson Binomial with rate:
    # gene1_discordance*gene2_discordance*normalized close pair core distances^2
    # So I need a matrix where each column is a gene pair and each row is a close pair
    # wait I want to compare scores with an rpbinom parameter, so I want the transverse
    # I will have to double-check this, but right now I am just debugging the speed of this function
    blocksize <- min(3000, dim(pa_matrix)[1])
    numblocks <- length(disc1) %/% blocksize
    leftover <- length(disc1) %% blocksize
    #t <- proc.time()
    res <- c()
    for(i in 1:numblocks){
      start <- (i-1)*blocksize+1
      end <- i*blocksize
      spbmat <- ( disc1[start:end] * disc2[start:end] ) %*% t(ndpds^2)
      spbmat[which(is.na(spbmat))] <- 0
      spbmat[which(is.nan(spbmat))] <- 0
      spbmat[which(spbmat < 0)] <- 0
      spbmat[which(spbmat > 1)] <- 1
      res <- c(res, apply(spbmat, 1, function(x){
        qpbinom_modified(1-0.05/length(disc1), x, method = 'RefinedNormal')
      }))
    }
    if(leftover > 0){
      start <- i*blocksize + 1
      end <- length(disc1)
      spbmat <- ( disc1[start:end] * disc2[start:end] ) %*% t(ndpds^2)
      spbmat[which(is.na(spbmat))] <- 0
      spbmat[which(is.nan(spbmat))] <- 0
      spbmat[which(spbmat < 0)] <- 0
      spbmat[which(spbmat > 1)] <- 1
      res <- c(res, apply(spbmat, 1, function(x){
        qpbinom_modified(1-0.05/length(disc1), x, method = 'RefinedNormal')
      }))
    }
    #print(proc.time() - t)
    # then we have
    scores$sig <- sum > res
    #pbmat <- (disc1 * disc2) %*% t(ndpds^2)
    # running into memory issues here, which is why I used the block method above
    if(verbose){print('Obtained significance')}
  }
  return(scores)
}

