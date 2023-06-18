#' Check whether focal assignments are actually focal
#'
#' @param Z The n-dimensional binary treatment assignment vector
#' @param A The n times n binary adjacency matrix
#' @param hypothesis
#' A character specifying the null hypothesis of interest.
#' Options are "Fisher", "SUTVA", "exposure1", and "exposure2".
#' @param focal_unit
#' The n-dimensional logical vector indicating whether each unit is focal.
#' @param focal_asgmt
#' The n times (1 + num_randomization) matrix of possibly focal assignments
#' @param E0_func The null exposure mapping
#'
#' @returns The matrix of treatment assignments that are actually focal.
#'
#' @noRd
#'
check_focal_asgmt <- function(Z,
                              A,
                              hypothesis,
                              focal_unit,
                              focal_asgmt,
                              E0_func) {

  # Variable definitions -------------------------------------------------------

  # Sample size
  n <- length(Z)

  # ID of focal units
  id_focal_unit <- which(focal_unit)

  # Number of focal units
  num_focal_unit <- sum(focal_unit)

  # Number of possibly focal assignments
  num_focal_asgmt <- ncol(focal_asgmt)

  # Dimension of null exposure mapping
  E0_dim <- 1 + (hypothesis == "exposure2") * 1

  # Realized null exposures for focal units
  E0_obs_focal <- E0_func(Z = Z, A = A, i = id_focal_unit)

  # Focal units' null exposures under possibly focal assignments
  E0_rand_focal <- array(NA, dim = c(num_focal_unit, num_focal_asgmt, E0_dim))

  for (id_asgmt in 1:num_focal_asgmt) {

    E0_rand_focal[, id_asgmt, ] <- E0_func(Z = focal_asgmt[, id_asgmt],
                                           A = A,
                                           i = id_focal_unit)

  }

  # Compute level set ----------------------------------------------------------

  level <- matrix(TRUE, nrow = num_focal_unit, ncol = num_focal_asgmt)

  if (E0_dim == 1) {

    level <- (E0_rand_focal[, , 1] == E0_obs_focal)

  } else {

    for (e0 in 1:E0_dim) {

      level <- (level) & (E0_rand_focal[, , e0] == E0_obs_focal[, e0])

    }
  }

  # Find actually focal assignments ----------------------------------------------

  id_actual_focal_asgmt <- which(apply(level, MARGIN = 2, FUN = all))

  actual_focal_asgmt <- focal_asgmt[, id_actual_focal_asgmt]

  actual_focal_asgmt <- unique(actual_focal_asgmt, MARGIN = 2)

  # Return ---------------------------------------------------------------------

  return(actual_focal_asgmt)

}
