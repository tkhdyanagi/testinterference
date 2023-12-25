#' Error handling
#'
#' @returns NULL if there are no errors
#'
#' @noRd
#'
errorhandling <- function(Y,
                          Z,
                          A,
                          hypothesis,
                          method,
                          design,
                          prob,
                          focal_unit,
                          num_focal_unit,
                          num_randomization,
                          strata,
                          zmatrix,
                          kappa,
                          cores) {

  # Y, Z, A --------------------------------------------------------------------

  # Binary Z
  if (any(Z != 0 & Z != 1)) {

    stop(paste("Z must not contain the values other than 0 and 1."))

  }

  if (hypothesis == "Fisher") {

    # Numeric
    if (!is.numeric(Y) | !is.numeric(Z)) {

      stop(paste("The elements of Y and Z must be numeric."))

    }

    # Dimensions
    if (length(Y) != length(Z)) {

      stop(paste("The lengths of Y and Z must be the same."))

    }

    # NA
    if (any(is.na(Y)) | any(is.na(Z))) {

      stop(paste("Y and Z must not contain NA."))

    }

  } else {

    # Numeric
    if (!is.numeric(Y) | !is.numeric(Z) | !is.numeric(A)) {

      stop(paste("The elements of Y, Z, and A must be numeric."))

    }

    # Dimensions
    if (length(Y) != length(Z) |
        length(Y) != ncol(A)   |
        length(Y) != nrow(A)) {

      stop(paste("The lengths of Y, Z, and of the row and column of A must be the same."))

    }

    # NA
    if (any(is.na(Y)) | any(is.na(Z)) | any(is.na(A))) {

      stop(paste("Y, Z, and A must not contain NA."))

    }

    # Binary A
    if (any(A != 0 & A != 1)) {

      stop(paste("A must not contain the values other than 0 and 1."))

    }

    # Diagonal elements of A
    if (any(diag(A) != 0)) {

      stop(paste("The diagonal elements of A must be 0."))

    }
  }

  # hypothesis -----------------------------------------------------------------

  if (hypothesis != "Fisher"    &
      hypothesis != "SUTVA"     &
      hypothesis != "exposure1" &
      hypothesis != "exposure2") {

    stop(paste("The argument 'hypothesis' must be one of 'Fisher', 'SUTVA', 'exposure1', and 'exposure2'."))

  }

  # method ---------------------------------------------------------------------

  if (method != "3-net"  &
      method != "2-net"  &
      method != "random" &
      method != "manual") {

    stop(paste("The argument 'method' must be one of '3-net', '2-net', 'random', and 'manual'."))

  }

  # design ---------------------------------------------------------------------

  if (design != "Bernoulli"  &
      design != "complete"   &
      design != "stratified" &
      design != "other") {

    stop(paste("The argument 'design' must be one of 'Bernoulli', 'complete', 'stratified', and 'other'."))

  }

  # prob -----------------------------------------------------------------------

  if (design == "Bernoulli") {

    if (!is.numeric(prob)) {

      stop(paste("The argument 'prob' must be a probability or an n-dimensional vector of probabilities when `design = 'Bernoulli'`."))

    }

    if (any(prob <= 0) | any(prob >= 1)) {

      stop(paste("The argument 'prob' must be a probability or an n-dimensional vector of probabilities when `design = 'Bernoulli'`."))

    }
  }


  # focal_unit -----------------------------------------------------------------

  if (method == "manual") {

    if (is.null(focal_unit) | !is.logical(focal_unit)) {

      stop(paste("The argument 'focal_unit' must be an n-dimensional logical vector when `method = 'manual'`."))

    }

    if (length(focal_unit) != length(Y)) {

      stop(paste("The argument 'focal_unit' must be an n-dimensional logical vector when `method = 'manual'`."))

    }
  }

  # num_focal_unit -------------------------------------------------------------

  # num_randomization ----------------------------------------------------------

  if (!is.numeric(num_randomization)) {

    stop(paste("The argument 'num_randomization' must be a large positive integer."))

  }

  if (is.numeric(num_randomization) & num_randomization <= 0) {

    stop(paste("The argument 'num_randomization' must be a large positive integer."))

  }

  # strata ---------------------------------------------------------------------

  if (design == "stratified") {

    if (!is.numeric(strata)) {

      stop(paste("The argument 'strata' must be an n-dimensional numerical vector when `design == 'stratified'`."))

    }

    if (length(strata) != length(Y)) {

      stop(paste("The argument 'strata' must be an n-dimensional numerical vector when `design == 'stratified'`."))

    }
  }

  # zmatrix --------------------------------------------------------------------

  if (design == "other") {

    if (!is.matrix(zmatrix)) {

      stop(paste("The argument 'zmatrix' must be a large numerical matrix when `design == 'other'`."))

    }

    if (nrow(zmatrix) != length(Y)) {

      stop(paste("The number of rows of 'zmatrix' must equal to the sample size n when `design == 'other'`."))

    }
  }

  # kappa ----------------------------------------------------------------------

  if (hypothesis == "exposure2") {
    if (!is.null(kappa)) {

      if (!is.numeric(kappa)) {

        stop(paste("The argument 'kappa' must be NULL or a positive integer no less than 2."))

      }

      if (kappa < 2) {

        stop(paste("The argument 'kappa' must be NULL or a positive integer no less than 2."))

      }
    }
  }

  # cores ----------------------------------------------------------------------

  max_cores <- parallel::detectCores()

  if (!is.numeric(cores)) {

    stop(paste("The argument 'cores' must be a positive integer not greater than the number of available CPU cores", max_cores))

  }

  if (cores > max_cores) {

    stop(paste("The argument 'cores' must be a positive integer not greater than the number of available CPU cores", max_cores))

  }

  # Return ---------------------------------------------------------------------

  return(NULL)

}
