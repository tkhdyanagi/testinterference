#' Find focal assignments
#'
#' @param Z The n-dimensional binary treatment assignment vector
#' @param A The n times n binary adjacency matrix
#' @param hypothesis
#' A character specifying the null hypothesis of interest.
#' Options are "Fisher", "SUTVA", "exposure1", and "exposure2".
#' @param design
#' A character specifying the randomization design of the original experiment.
#' Options are "Bernoulli", "complete", "stratified", and "other".
#' @param prob
#' NULL,
#' a scalar indicating the homogeneous probability of being treatment eligible,
#' or an n-dimensional vector indicating the heterogeneous probabilities
#' for treatment eligibility (i.e., the i-th element of `prob` indicates
#' the probability that unit i becomes treatment eligible).
#' This argument is used only when `design = "Bernoulli"`.
#' Default is NULL.
#' @param strata
#' The n-dimensional vector indicating the stratum to which each unit belongs.
#' This argument is used only when `design = "stratified"`.
#' @param focal_unit
#' The n-dimensional logical vector indicating whether each unit is focal.
#' @param num_randomization
#' A positive scalar specifying the number of randomization.
#' Default is 9999.
#' @param E0_func The null exposure mapping
#' @param E1_func The coarser exposure mapping
#'
#' @importFrom foreach %dopar%
#'
#' @returns The n times (1 + num_randomization) matrix of focal assignments
#'
#' @noRd
#'
find_focal_asgmt <- function(Z,
                             A,
                             hypothesis,
                             design,
                             prob,
                             strata,
                             focal_unit,
                             num_randomization,
                             E0_func,
                             E1_func) {

  # Variable definitions -------------------------------------------------------

  # Sample size
  n <- length(Z)

  # Matrix of focal assignments
  focal_asgmt <- matrix(NA, nrow = n, ncol = 1 + num_randomization)

  # Set the first column as the realized treatment assignment
  focal_asgmt[, 1] <- Z

  # Number of focal assignments
  num_focal_asgmt <- ncol(focal_asgmt)

  # ID of focal units
  id_focal_unit <- which(focal_unit)

  # Number of focal units
  num_focal_unit <- sum(focal_unit)

  # Support of strata
  supp_strata <- sort(unique(strata))

  # Hypothesis: "Fisher" -------------------------------------------------------

  if (hypothesis == "Fisher") {

    if (design == "Bernoulli") {

      # Find focal assignments
      if (length(prob) == 1) {

        focal_asgmt[, -1] <- matrix(data = sample(x = 0:1,
                                                  prob = c(1 - prob, prob),
                                                  size = n * num_randomization,
                                                  replace = TRUE),
                                    nrow = n,
                                    ncol = num_randomization,
                                    byrow = TRUE)

      } else if (length(prob) == n) {

        focal_asgmt[, -1] <- t(apply(X = cbind(1 - prob, prob),
                                     FUN = sample,
                                     MARGIN = 1,
                                     x = 0:1,
                                     size = num_randomization,
                                     replace = TRUE))

      }
    }

    if (design == "complete" | design == "stratified") {

      # Randomization within each stratum ("st" stands for stratum)
      for (st in supp_strata) {

        # Units in this stratum
        st_unit <- (strata == st)

        # Number of units in this stratum
        num_st_unit <- sum(st_unit)

        if (num_st_unit == 1) {

          # Set focal assignments for this unit
          focal_asgmt[st_unit, -1] <- Z[st_unit]

        } else if (num_st_unit > 1) {

          # Support of randomization in this stratum
          supp_Z_st <- Z[st_unit]

          # Randomization in this stratum
          focal_asgmt[st_unit, -1] <-
            matrix(data = rep(supp_Z_st, num_randomization),
                   nrow = num_st_unit,
                   ncol = num_randomization,
                   byrow = FALSE) %>%
            apply(FUN = sample, MARGIN = 2)

        }
      }
    }
  }

  # Hypothesis: "SUTVA" --------------------------------------------------------

  if (hypothesis == "SUTVA") {

    # Units subject to randomization
    rand_unit <- !focal_unit

    # Number of units subject to randomization
    num_rand_unit <- sum(rand_unit)

    # Set focal assignments for focal units
    focal_asgmt[focal_unit, ] <- Z[focal_unit]

    # Bernoulli randomization --------------------------------------------------

    if (design == "Bernoulli") {

      # Find focal assignments for units subject to randomization
      if (length(prob) == 1) {

        focal_asgmt[rand_unit, -1] <-
          matrix(data = sample(x = 0:1,
                               prob = c(1 - prob, prob),
                               size = num_rand_unit * num_randomization,
                               replace = TRUE),
                 nrow = num_rand_unit,
                 ncol = num_randomization,
                 byrow = TRUE)

      } else if (length(prob) == n) {

        focal_asgmt[rand_unit, -1] <-
          t(apply(X = cbind(1 - prob[rand_unit], prob[rand_unit]),
                  FUN = sample,
                  MARGIN = 1,
                  x = 0:1,
                  size = num_randomization,
                  replace = TRUE))

      }
    }

    # Complete / stratified randomization --------------------------------------

    if (design == "complete" | design == "stratified") {

      # Randomization within each stratum ("st" stands for stratum)
      for (st in supp_strata) {

        # Units subject to randomization in this stratum
        st_rand_unit <- (rand_unit) & (strata == st)

        # Number of units subject to randomization in this stratum
        num_st_rand_unit <- sum(st_rand_unit)

        if (num_st_rand_unit == 1) {

          # Set focal assignments for this unit
          focal_asgmt[st_rand_unit, -1] <- Z[st_rand_unit]

        } else if (num_st_rand_unit > 1) {

          # Support of randomization in this stratum
          supp_Z_st_rand <- Z[st_rand_unit]

          # Find focal assignments for units subject to randomization
          focal_asgmt[st_rand_unit, -1] <-
            matrix(data = rep(supp_Z_st_rand, num_randomization),
                   nrow = num_st_rand_unit,
                   ncol = num_randomization,
                   byrow = FALSE) %>%
            apply(FUN = sample, MARGIN = 2)

        }
      }
    }
  }

  # Hypothesis: "exposure1" ----------------------------------------------------

  if (hypothesis == "exposure1") {

    # List of treatment eligible units in each focal unit's neighborhood -------

    list_Z1 <- foreach::foreach(id_unit = id_focal_unit, .combine = "append", .inorder = TRUE) %dopar% {

      if (Z[id_unit] == 1) {

        list(c(id_unit, which(A[id_unit, ] == 1 & Z == 1)))

      } else {

        list(which(A[id_unit, ] == 1 & Z == 1))

      }
    }

    # Randomization ------------------------------------------------------------

    for (id_rand in 2:num_focal_asgmt) {

      # Random selection of treatment eligible units fixed at Z = 1
      id_fix_Z1 <- NULL
      for (i in 1:num_focal_unit) {
        if (length(list_Z1[[i]]) == 1) {

          id_fix_Z1 <- c(id_fix_Z1, list_Z1[[i]])

        } else {

          id_fix_Z1 <- c(id_fix_Z1, sample(x = list_Z1[[i]], size = 1))

        }
      }

      # n-dimensional logical vector for treatment eligible units fixed at Z = 1
      fix_Z1 <- 1:n %in% id_fix_Z1

      # Units subject to randomization
      rand_unit <- !fix_Z1

      # Number of units subject to randomization
      num_rand_unit <- sum(rand_unit)

      # Set focal assignments for treatment eligible units fixed at Z = 1
      focal_asgmt[fix_Z1, id_rand] <- 1

      # Bernoulli randomization ------------------------------------------------

      if (design == "Bernoulli") {

        # Find focal assignments for units subject to randomization
        if (length(prob) == 1) {

          focal_asgmt[rand_unit, id_rand] <-
            sample(x = 0:1,
                   prob = c(1 - prob, prob),
                   size = num_rand_unit,
                   replace = TRUE)

        } else if (length(prob) == n) {

          focal_asgmt[rand_unit, id_rand] <-
            apply(X = cbind(1 - prob[rand_unit], prob[rand_unit]),
                  FUN = sample,
                  MARGIN = 1,
                  x = 0:1,
                  size = 1,
                  replace = FALSE)

        }
      }

      # Complete / stratified randomization ------------------------------------

      if (design == "complete" | design == "stratified") {

        # Randomization within each stratum ("st" stands for stratum)
        for (st in supp_strata) {

          # Units subject to randomization in this stratum
          st_rand_unit <- (rand_unit) & (strata == st)

          # Number of units subject to randomization in this stratum
          num_st_rand_unit <- sum(st_rand_unit)

          if (num_st_rand_unit > 0) {

            # Support of randomization in this stratum
            supp_Z_st_rand <- Z[st_rand_unit]

            # Find focal assignments for units subject to randomization
            focal_asgmt[st_rand_unit, id_rand] <- sample(supp_Z_st_rand)

          }
        }
      }
    }
  }

  # Hypothesis: "exposure2" ----------------------------------------------------

  if (hypothesis == "exposure2") {

    # Set focal assignments for focal units
    focal_asgmt[focal_unit, ] <- Z[focal_unit]

    # List of treatment eligible peers of each focal unit
    list_Z1 <- foreach::foreach(id_unit = id_focal_unit, .combine = "append", .inorder = TRUE) %dopar% {

      list(which(A[id_unit, ] == 1 & Z == 1))

    }

    # Random selection of treatment eligible units fixed at Z = 1
    id_fix_Z1 <- NULL
    for (i in 1:num_focal_unit) {
      if (length(list_Z1[[i]]) == 1) {

        id_fix_Z1 <- c(id_fix_Z1, list_Z1[[i]])

      } else {

        id_fix_Z1 <- c(id_fix_Z1, sample(list_Z1[[i]], size = 1))

      }
    }

    # n-dimensional logical vector for treatment eligible units fixed at Z = 1
    fix_Z1 <- 1:n %in% id_fix_Z1

    # Set focal assignments for treatment eligible units fixed at Z = 1
    focal_asgmt[fix_Z1, ] <- 1

    # Units subject to randomization
    rand_unit <- (!focal_unit) & (!fix_Z1)

    # Number of units subject to randomization
    num_rand_unit <- sum(rand_unit)

    # Bernoulli randomization --------------------------------------------------

    if (design == "Bernoulli") {

      # Find focal assignments for units subject to randomization
      if (length(prob) == 1) {

        focal_asgmt[rand_unit, -1] <-
          matrix(data = sample(x = 0:1,
                               prob = c(1 - prob, prob),
                               size = num_rand_unit * num_randomization,
                               replace = TRUE),
                 nrow = num_rand_unit,
                 ncol = num_randomization,
                 byrow = TRUE)

      } else if (length(prob) == n) {

        focal_asgmt[rand_unit, -1] <-
          t(apply(X = cbind(1 - prob[rand_unit], prob[rand_unit]),
                  FUN = sample,
                  MARGIN = 1,
                  x = 0:1,
                  size = num_randomization,
                  replace = TRUE))

      }
    }

    # Complete / stratified randomization --------------------------------------

    if (design == "complete" | design == "stratified") {

      # Randomization within each stratum ("st" stands for stratum)
      for (st in supp_strata) {

        # Units subject to randomization in this stratum
        st_rand_unit <- (rand_unit) & (strata == st)

        # Number of units subject to randomization in this stratum
        num_st_rand_unit <- sum(st_rand_unit)

        if (num_st_rand_unit == 1) {

          # Set focal assignments for this unit
          focal_asgmt[st_rand_unit, -1] <- Z[st_rand_unit]

        } else if (num_st_rand_unit > 1) {

          # Support of randomization in this stratum
          supp_Z_st_rand <- Z[st_rand_unit]

          # Find focal assignments for units subject to randomization
          focal_asgmt[st_rand_unit, -1] <-
            matrix(data = rep(supp_Z_st_rand, num_randomization),
                   nrow = num_st_rand_unit,
                   ncol = num_randomization,
                   byrow = FALSE) %>%
            apply(FUN = sample, MARGIN = 2)

        }
      }
    }
  }

  # Return ---------------------------------------------------------------------

  return(focal_asgmt)

}
