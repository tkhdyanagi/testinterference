#' Randomization Test for the Interference Structure
#'
#' Randomization test for the specification of spillover effects
#' utilizing hierarchical relationships between different exposures
#'
#' @details
#' The `testinterference()` function performs randomization tests for
#' the null hypotheses about the interference structure
#' utilizing hierarchical relationships between different exposures.
#' The null hypotheses of interest are the following.
#' (1) Fisher's sharp null hypothesis:
#' One's outcome does not depend on one's own treatment.
#' (2) Stable Unit Treatment Value Assumption (SUTVA):
#' One's outcome is determined by one's own treatment.
#' (3) Exposure 1:
#' One's outcome is determined by whether
#' there is at least one treated unit in one's neighborhood including oneself.
#' (4) Exposure 2:
#' One's outcome is determined by
#' one's own treatment and whether there is at least one treated peer.
#' The testing procedures are developed by Hoshino and Yanagi (2023).
#'
#' @param Y The n-dimensional outcome vector
#' @param Z The n-dimensional binary treatment assignment vector
#' @param A The n times n, binary, possibly asymmetric, adjacency matrix.
#' The diagonal elements must be zero.
#' This argument can be NULL if `hypothesis = "Fisher"`.
#' Default is NULL.
#' @param hypothesis
#' A character specifying the null hypothesis of interest.
#' Options are "Fisher", "SUTVA", "exposure1", and "exposure2".
#' "Fisher" stands for
#' Fisher's sharp null hypothesis in the canonical randomization test.
#' "SUTVA", "exposure1", and "exposure2" correspond to
#' the null hypotheses a, b, and c in Table 2 of Hoshino and Yanagi (2023).
#' Default is "SUTVA".
#' If `hypothesis = "exposure2"`, the user can specify the argument `kappa`.
#' @param method
#' A character specifying how to find focal units.
#' Options are "3-net", "2-net", "random", and "manual",
#' which stand for the 3-net (MIS-G) method, the 2-net (MIS-A) method,
#' the random selection, and the manual selection.
#' Default is "3-net".
#' If `method = "random"`, the user can specify the argument `num_focal_unit`.
#' If `method = "manual"`, the user must specify the argument `focal_unit`.
#' @param design
#' A character specifying the randomization design of the original experiment.
#' Options are "Bernoulli", "complete", "stratified", and "other",
#' which stand for Bernoulli randomization, complete randomization,
#' stratified randomization, and other experimental designs.
#' If `design = "Bernoulli"`, the user must specify the argument `prob`.
#' If `design = "stratified"`,the user must specify the argument `strata`.
#' If `design = "other"`, the user must specify the argument `zmatrix`.
#' @param prob
#' NULL,
#' a scalar indicating the homogeneous probability of being treatment eligible,
#' or an n-dimensional vector indicating the heterogeneous probabilities
#' for treatment eligibility (i.e., the i-th element of `prob` indicates
#' the probability that unit i becomes treatment eligible).
#' This argument is used only when `design = "Bernoulli"`.
#' Default is NULL.
#' @param focal_unit
#' NULL
#' or an n-dimensional logical vector specifying whether each unit is focal.
#' This argument is used only when `method = "manual"`.
#' Default is NULL.
#' @param num_focal_unit
#' NULL or a positive scalar specifying the number of focal units.
#' This argument is used only when `method = "random"`.
#' Default is NULL.
#' @param num_randomization
#' A large positive integer specifying the number of randomization.
#' Default is 9999.
#' This argument is not used when `design = "other"`.
#' @param strata
#' NULL or an n-dimensional numerical vector indicating
#' the stratum to which each unit belongs.
#' This argument is used only when `design = "stratified"`.
#' Default is NULL.
#' @param zmatrix
#' NULL or a large matrix of realizations of treatment assignments.
#' This argument must be specified by the user only when `design = "other"`.
#' The number of rows must equal to the sample size n.
#' The number of columns is the number of realizations given by the user.
#' Default is NULL.
#' @param kappa
#' NULL or a positive integer no less than 2.
#' This argument is used only when `hypothesis = "exposure2"`.
#' If `kappa = NULL`,
#' kappa is automatically chosen to maximize the number of focal units.
#' Default is NULL.
#' @param cores
#' A positive integer specifying the number of cores
#' to use for parallel computing.
#' Default is 1.
#'
#' @returns A list containing the following elements:
#' \item{pval}{The vector of p-values obtained from KW, ACD, and OLS statistics}
#' \item{Simes}{The vector of testing results by Simes’ correction under significance levels 10%, 5%, and 1%.
#' TRUE (resp. FALSE) indicates the rejection (resp. acceptance) of the null hypothesis.}
#' \item{stat}{A matrix of test statistics computed with focal assignments}
#' \item{focal_unit}{The logical vector indicating whether each unit is focal}
#' \item{focal_asgmt}{The matrix of focal assignments}
#' \item{num_focal_unit}{The number of focal units}
#' \item{num_focal_asgmt}{The number of focal assignments}
#'
#' @importFrom purrr %>%
#' @importFrom utils globalVariables
#'
#' @examples
#' data <- datageneration(n = 200)
#' test <- testinterference(Y = data$Y,
#'                          Z = data$Z,
#'                          A = data$A,
#'                          hypothesis = "SUTVA",
#'                          method = "3-net",
#'                          design = "complete",
#'                          num_randomization = 999,
#'                          cores = 1)
#'
#' @references Hoshino, T., & Yanagi, T. (2023).
#' Randomization Test for the Specification of Interference Structure.
#' arXiv preprint arXiv:2301.05580.
#'
#' @export
#'
testinterference <- function(Y                 = NULL,
                             Z                 = NULL,
                             A                 = NULL,
                             hypothesis        = "SUTVA",
                             method            = "3-net",
                             design            = "complete",
                             prob              = NULL,
                             focal_unit        = NULL,
                             num_focal_unit    = NULL,
                             num_randomization = 9999,
                             strata            = NULL,
                             zmatrix           = NULL,
                             kappa             = NULL,
                             cores             = 1) {

  # Error handling -------------------------------------------------------------

  error <- errorhandling(Y                 = Y,
                         Z                 = Z,
                         A                 = A,
                         hypothesis        = hypothesis,
                         method            = method,
                         design            = design,
                         prob              = prob,
                         focal_unit        = focal_unit,
                         num_focal_unit    = num_focal_unit,
                         num_randomization = num_randomization,
                         strata            = strata,
                         zmatrix           = zmatrix,
                         kappa             = kappa,
                         cores             = cores)

  # Variables used for parallel computing --------------------------------------

  cl <- parallel::makeCluster(cores)

  doParallel::registerDoParallel(cl, cores = cores)

  # Variable definitions -------------------------------------------------------

  # Sample size
  n <- length(Z)

  # Names
  names(Y) <- names(Z) <- rownames(A) <- colnames(A) <- 1:n

  # Set kappa
  if (hypothesis == "Fisher" | hypothesis == "SUTVA") {

    kappa <- 2

  }

  if (hypothesis == "exposure1") {

    kappa <- 3

  }

  if (hypothesis == "exposure2" & is.null(kappa)) {

    # Compute out-degrees no less than 2
    out_degree <- rowSums(A) %>%
      setdiff(0:1)

    # Set kappa as the mode of out-degrees no less than 2
    kappa <- which.max(table(out_degree)) %>%
      names() %>%
      as.integer()

  }

  # Set strata = 1 if design == "complete"
  if (design == "complete") {

    strata <- rep(1, n)

  }

  # Exposure mappings ----------------------------------------------------------

  if (hypothesis == "Fisher") {

    E0_func <- constant
    E1_func <- exposure_a

  } else if (hypothesis == "SUTVA") {

    E0_func <- exposure_a
    E1_func <- exposure_c

  } else if (hypothesis == "exposure1") {

    E0_func <- exposure_b
    E1_func <- exposure_c

  } else if (hypothesis == "exposure2") {

    E0_func <- exposure_c
    E1_func <- exposure_d

  }

  # Realized exposures
  E0_obs <- E0_func(Z = Z, A = A, i = 1:n)
  E1_obs <- E1_func(Z = Z, A = A, i = 1:n)

  # Find focal units --------------------------------------------------------

  if (hypothesis == "Fisher") {

    focal_unit <- rep(TRUE, n)

  } else {

    if (method %in% c("3-net", "2-net", "random")) {

      focal_unit <- find_focal_unit(Z              = Z,
                                    A              = A,
                                    hypothesis     = hypothesis,
                                    method         = method,
                                    num_focal_unit = num_focal_unit,
                                    kappa          = kappa,
                                    E0_obs         = E0_obs)

    }
  }

  if (sum(focal_unit) == 0) {

    stop(paste("No focal units found. Perhaps the null hypothesis is not compatible with the network A given by the user."))

  }

  # Find focal assignments --------------------------------------------------

  if (design != "other") {

    focal_asgmt <- find_focal_asgmt(Z                 = Z,
                                    A                 = A,
                                    hypothesis        = hypothesis,
                                    design            = design,
                                    prob              = prob,
                                    strata            = strata,
                                    focal_unit        = focal_unit,
                                    num_randomization = num_randomization,
                                    E0_func           = E0_func,
                                    E1_func           = E1_func)

  } else {

    focal_asgmt <- zmatrix

  }

  # Check for focal assignments ------------------------------------------------

  focal_asgmt <- check_focal_asgmt(Z           = Z,
                                   A           = A,
                                   hypothesis  = hypothesis,
                                   focal_unit  = focal_unit,
                                   focal_asgmt = focal_asgmt,
                                   E0_func     = E0_func)

  # Randomization test ---------------------------------------------------------

  RT <- randomization_test(Y           = Y,
                           Z           = Z,
                           A           = A,
                           hypothesis  = hypothesis,
                           focal_unit  = focal_unit,
                           focal_asgmt = focal_asgmt,
                           kappa       = kappa,
                           E0_func     = E0_func,
                           E1_func     = E1_func)

  # Simes' correction ----------------------------------------------------------

  # Function for Simes' correction
  Simes_func <- function(pval_vec, alpha) {

    # Sort p-values
    pval_sort <- sort(pval_vec)

    # Number of p-values
    num_pval <- length(pval_sort)

    # Reject or not under the significance level alpha
    reject <- FALSE
    for (i in 1:num_pval) {
      if (pval_sort[i] <= i * alpha / num_pval) {

        reject <- TRUE

      }
    }

    # Return
    return(reject)

  }

  # Results of Simes' correction when alpha = 0.1, 0.05, 0.01
  Simes010 <- Simes_func(pval_vec = RT$pval, alpha = 0.10)
  Simes005 <- Simes_func(pval_vec = RT$pval, alpha = 0.05)
  Simes001 <- Simes_func(pval_vec = RT$pval, alpha = 0.01)

  # Vector of Simes' correction results
  Simes <- c(Simes010, Simes005, Simes001)

  # Names
  names(Simes) <- c("10%", "5%", "1%")

  # End parallel computing -----------------------------------------------------

  parallel::stopCluster(cl)

  # Return ---------------------------------------------------------------------

  return(list(pval            = RT$pval,
              Simes           = Simes,
              stat            = RT$stat,
              focal_unit      = focal_unit,
              focal_asgmt     = focal_asgmt,
              num_focal_unit  = sum(focal_unit),
              num_focal_asgmt = ncol(focal_asgmt)))

}
