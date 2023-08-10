#' Find focal units
#'
#' @param Z The n-dimensional binary treatment assignment vector
#' @param A The n times n binary adjacency matrix
#' @param hypothesis
#' A character specifying the null hypothesis of interest.
#' Options are "Fisher", "SUTVA", "exposure1", and "exposure2".
#' @param method
#' A character specifying how to specify focal units.
#' Options are "MIS" and "random".
#' @param num_focal_unit
#' A positive scalar specifying the number of focal units.
#' This argument is used only when `method = "random"`.
#' @param kappa
#' NULL or a positive integer no less than 2.
#' This argument is used only when `hypothesis = "exposure2"`.
#' @param E0_obs The n-dimensional vector of null exposures
#'
#' @returns The n-dimensional logical vector indicating whether each unit is focal.
#'
#' @noRd
#'
find_focal_unit <- function(Z,
                            A,
                            hypothesis,
                            method,
                            num_focal_unit,
                            kappa,
                            E0_obs) {

  # Variable definitions -------------------------------------------------------

  # Sample size
  n <- length(Z)

  # Specify N(kappa) -----------------------------------------------------------

  if (hypothesis == "SUTVA") {

    N_kappa <- (rowSums(A) > 0)

  } else if (hypothesis == "exposure1") {

    N_kappa <- (rowSums(A) > 0) & (E0_obs == 1)

  } else if (hypothesis == "exposure2") {

    N_kappa <- (rowSums(A) == kappa) & (E0_obs[, 2] == 1)

  }

  # Variable definitions -------------------------------------------------------

  # Units who belong to N(kappa)
  id_N_kappa <- which(N_kappa)

  # Number of units in N(kappa)
  num_N_kappa <- sum(N_kappa)

  # Method of MIS --------------------------------------------------------------

  if (method == "MIS") {

    # Common-friend graph ------------------------------------------------------

    # Adjacency matrix restricted on N(kappa)
    A_N_kappa <- A[N_kappa, N_kappa]

    # Zero matrix
    A_common <- matrix(0, nrow = num_N_kappa, ncol = num_N_kappa)

    # Matrix B with (i, j)-th element = "if unit j belongs to i's neighborhood"
    B_N_kappa <- A_N_kappa
    diag(B_N_kappa) <- 1

    # Common-friend adjacency matrix
    for (i in 1:num_N_kappa) {
      A_common[i, ] <- (B_N_kappa %*% B_N_kappa[i, ] > 0) * 1
    }

    # Names
    rownames(A_common) <- colnames(A_common) <- id_N_kappa

    # Common-friend graph
    G_common <- igraph::graph_from_adjacency_matrix(adjmatrix = A_common)

    # Find focal units by the method of MIS via greedy vertex coloring ---------

    # greedy vertex coloring
    coloring_out <- igraph::greedy_vertex_coloring(graph = G_common)

    # "Maximum color"
    coloring_max <- which.max(table(coloring_out)) %>%
      names() %>%
      as.integer()

    # Unit ID with "maximum color"
    id_focal_unit <- which(coloring_out == coloring_max) %>%
      names() %>%
      as.integer() %>%
      sort()

  }

  # Random selection -----------------------------------------------------------

  if (method == "random") {

    # Set the number of focal units if num_focal_unit = NULL
    if (is.null(num_focal_unit)) {
      num_focal_unit <- floor(num_N_kappa / 2)
    }

    # Find focal units by random selection
    id_focal_unit <- sample(x       = id_N_kappa,
                            size    = num_focal_unit,
                            replace = FALSE)

  }

  # Return ---------------------------------------------------------------------

  # Focal units (n-dimensional logical vector)
  focal_unit <- 1:n %in% id_focal_unit

  return(focal_unit)

}
