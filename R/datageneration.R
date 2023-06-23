#' Generate Artificial Data by Simulation
#'
#' @param n An even number specifying the sample size
#' @param design
#' A character specifying the randomization design.
#' Options are "Bernoulli", "complete", and "stratified",
#' which stand for Bernoulli randomization, complete randomization, and
#' stratified randomization.
#' Default is "complete".
#' @param dgp Specify data generating process 1 or 2 for the outcome equation
#' @param beta_d
#' A scalar specifying the value
#' of a coefficient for the direct effect in the outcome equation.
#' There is no direct effect if beta_d is set to zero.
#' @param beta_s
#' A scalar specifying the value
#' of a coefficient for spillovers in the outcome equation.
#' There are no spillovers if beta_s is set to zero.
#' @param Y0
#' NULL or an n-dimensional vector of the untreated potential outcomes.
#' Default is NULL.
#' @param A
#' An n times n binary adjacency matrix or
#' a character specifying how to generate the adjacency matrix.
#' For the latter case, options are "Erdos-Renyi" and "pairs".
#' Default is "Erdos-Renyi".
#'
#' @returns A list containing the following elements.
#' \item{Y}{An n-dimensional outcome vector}
#' \item{Z}{An n-dimensional treatment assignment vector}
#' \item{A}{An n times n binary adjacency matrix}
#' \item{strata}{An n-dimensional vector indicating strata}
#'
#' @importFrom purrr %>%
#'
#' @examples
#' data <- datageneration(n = 200)
#'
#' @export
#'
datageneration <- function(n,
                           design = "complete",
                           dgp    = 1,
                           beta_d = 1,
                           beta_s = 0,
                           Y0     = NULL,
                           A      = "Erdos-Renyi") {

  # The untreated potential outcome --------------------------------------------

  if (is.null(Y0)) {

    Y0 <- stats::rnorm(n)

  }

  # Network --------------------------------------------------------------------

  if (A == "Erdos-Renyi") {

    A <- igraph::erdos.renyi.game(n = n, p.or.m = 3/n, type = "gnp") %>%
      igraph::get.adjacency() %>%
      as.matrix()

  } else if (A == "pairs") {

    A <- matrix(0, nrow = n, ncol = n)

    edges <- cbind(1:n,
                   c(rbind(seq(2,     n, by = 2),
                           seq(1, n - 1, by = 2))))

    A[edges] <- 1

  }

  # Treatment assignment -------------------------------------------------------

  if (design == "Bernoulli") {

    Z <- sample(x = 0:1, prob = c(0.5, 0.5), size = n, replace = TRUE)

    strata <- NULL

  }

  if (design == "complete") {

    n1 <- n0 <- 0.5 * n

    supp_Z <- c(rep(1, n1), rep(0, n0))

    Z <- sample(supp_Z)

    strata <- NULL

  }

  if (design == "stratified") {

    strata <- sample(x = 1:2, prob = c(0.5, 0.5), size = n, replace = TRUE)

    Z <- rep(NA, n)

    for (s in 1:2) {

      s_size <- sum(strata == s)

      n1 <- ceiling(0.5 * s_size)
      n0 <-   floor(0.5 * s_size)

      supp_Z <- c(rep(1, n1), rep(0, n0))

      Z[strata == s] <- sample(supp_Z)

    }
  }

  # Outcome --------------------------------------------------------------------

  if (dgp == 1) {

    # True exposure mapping
    exposure <- (A %*% Z)

  } else if (dgp == 2) {

    # Number of treatment eligible peers
    num_t_peer <- (A %*% Z)

    # True exposure mapping
    exposure <-
      1 * (num_t_peer == 1) +
      2 * (num_t_peer == 2) +
      ifelse(num_t_peer >= 3, 1 / num_t_peer, 0)

  }

  # Outcome
  Y <- Y0 + beta_d * Z + beta_s * exposure

  # Return ---------------------------------------------------------------------

  return(list(Y      = Y,
              Z      = Z,
              A      = A,
              strata = strata))

}
