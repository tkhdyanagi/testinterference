#' Perform randomization tests
#'
#' @param Y The n-dimensional outcome vector
#' @param Z The n-dimensional binary treatment assignment vector
#' @param A The n times n binary adjacency matrix
#' @param hypothesis
#' A character specifying the null hypothesis of interest.
#' Options are "Fisher", "SUTVA", "exposure1", and "exposure2".
#' @param focal_unit
#' The n-dimensional logical vector indicating whether each unit is focal.
#' @param focal_asgmt The matrix of focal assignments
#' @param kappa A positive integer no less than 2
#' @param E0_func The null exposure mapping
#' @param E1_func The coarser exposure mapping
#'
#' @returns A list containing the following elements:
#' \item{pval}{A vector of p-values}
#' \item{stat}{A matrix of test statistics computed with focal assignments}
#'
#' @importFrom foreach %dopar%
#'
#' @noRd
#'
randomization_test <- function(Y,
                               Z,
                               A,
                               hypothesis,
                               focal_unit,
                               focal_asgmt,
                               kappa,
                               E0_func,
                               E1_func) {

  # Variable definitions -------------------------------------------------------

  # Sample size
  n <- length(Z)

  # ID of focal units
  id_focal_unit <- which(focal_unit)

  # Number of focal units
  num_focal_unit <- sum(focal_unit)

  # Number of focal assignments
  num_focal_asgmt <- ncol(focal_asgmt)

  # Outcome vector for focal units
  Y_focal <- Y[focal_unit]

  # Treatment assignment vector for focal units
  Z_focal <- Z[focal_unit]

  # Dimensions of exposures
  E0_dim <- 1 + (hypothesis == "exposure2") * 1
  E1_dim <- 2 - (hypothesis == "Fisher") * 1

  # Dimension of overlapping exposures in E0 and E1
  overlap_dim <- 0 + (hypothesis == "exposure2") * 1

  # Function to compute test statistics ----------------------------------------

  stat_func <- function(asgmt) {

    # Compute exposures for focal units
    E0_focal <- E0_func(Z = asgmt, A = A, i = id_focal_unit)
    E1_focal <- E1_func(Z = asgmt, A = A, i = id_focal_unit)

    # Partition focal units into kappa groups based on E1 values
    # In other words, construct the set "mathcal{S}_j(d)"
    group <- rep(0, num_focal_unit)

    if (hypothesis == "Fisher") {

      group <-
        1 * (E1_focal == 0) +
        2 * (E1_focal == 1)

    } else if (hypothesis == "SUTVA") {

      group <-
        1 * (E1_focal[, 1] == Z_focal & E1_focal[, 2] == 0) +
        2 * (E1_focal[, 1] == Z_focal & E1_focal[, 2] == 1)

    } else if (hypothesis == "exposure1") {

      group <-
        1 * (E1_focal[, 1] == 0 & E1_focal[, 2] == 1) +
        2 * (E1_focal[, 1] == 1 & E1_focal[, 2] == 0) +
        3 * (E1_focal[, 1] == 1 & E1_focal[, 2] == 1)

    } else if (hypothesis == "exposure2") {

      for (k in 1:kappa) {

        group <- group +
          k * (E1_focal[, 1] == Z_focal & E1_focal[, 2] == k)

      }
    }

    # Support of groups
    supp_group <- sort(unique(group))

    # Number of groups
    num_group <- length(supp_group)

    # Compute test statistics --------------------------------------------------

    if (num_group == 1) {

      rep(NA, 3)

    } else {

      # List of outcomes within each group
      Y_group <- list()
      for (g in supp_group) {

        Y_group <- append(Y_group, list(Y_focal[group == g]))

      }

      # KW statistic
      KW <- as.numeric(stats::kruskal.test(x = Y_group)[1])

      # ACD statistic
      diff <- NULL

      for (g1 in 1:num_group) {
        for (g2 in g1:num_group) {

          diff <- c(diff, abs(mean(Y_group[[g1]]) - mean(Y_group[[g2]])))

        }
      }

      ACD <- mean(diff)

      # OLS statistic
      reg <- stats::lm(Y_focal ~ E1_focal)

      OLS <- aod::wald.test(
        Sigma = aod::vcov(reg),
        b     = aod::coef(reg),
        Terms = (2 + overlap_dim):(1 + E1_dim))$result$chi2[1]

      # Test statistics
      c(KW, ACD, OLS)

    }
  }

  # Compute test statistics under realized treatment assignment -----------------

  stat_obs <- stat_func(asgmt = Z)

  if (any(is.na(stat_obs))) {

    stop(paste("Statistics cannot be computed with the given variables."))

  }

  # Compute test statistics under each focal assignment ------------------------

  packages <- c("aod", "stats")

  stat <- foreach::foreach(id_asgmt = 1:num_focal_asgmt, .combine = "rbind", .inorder = TRUE, .errorhandling = "remove", .packages = packages) %dopar% {

    # Set a focal assignment
    if (id_asgmt == 1) {

      asgmt <- Z

    } else {

      asgmt <- focal_asgmt[, id_asgmt]

    }

    # Compute test statistics
    stat_func(asgmt = asgmt)

  }

  # Compute p-values ---------------------------------------------------------

  # Vectors of test statistics
  KW_vec  <- stat[, 1]
  ACD_vec <- stat[, 2]
  OLS_vec <- stat[, 3]

  # Realized test statistics
  KW_obs  <-  KW_vec[1]
  ACD_obs <- ACD_vec[1]
  OLS_obs <- OLS_vec[1]

  # P-values
  pval_KW  <- mean(KW_vec  >=  KW_obs, na.rm = TRUE)
  pval_ACD <- mean(ACD_vec >= ACD_obs, na.rm = TRUE)
  pval_OLS <- mean(OLS_vec >= OLS_obs, na.rm = TRUE)

  # Vector of p-values
  pval <- c(pval_KW, pval_ACD, pval_OLS)

  # Return ---------------------------------------------------------------------

  # Names
  names(pval) <- colnames(stat) <- c("KW", "ACD", "OLS")

  return(list(pval = pval,
              stat = stat))

}
