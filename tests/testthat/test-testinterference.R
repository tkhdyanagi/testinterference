# Generate artificial data by simulation ---------------------------------------

# Seed
set.seed(1)

# Sample size
n <- 100

# Artificial data
data <- datageneration(n = n)

# Variable definitions ---------------------------------------------------------

Y                 <- data$Y
Z                 <- data$Z
A                 <- data$A
hypothesis        <- "SUTVA"
method            <- "3-net"
design            <- "complete"
prob              <- NULL
focal_unit        <- NULL
num_focal_unit    <- NULL
num_randomization <- 999
strata            <- NULL
zmatrix           <- NULL
kappa             <- NULL
cores             <- 1

# Y, Z, A ----------------------------------------------------------------------

expect_error(testinterference(Y                 = NULL,
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
                              cores             = cores),
             "The elements of Y, Z, and A must be numeric.")

expect_error(testinterference(Y                 = Y,
                              Z                 = rep(1, n / 2),
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
                              cores             = cores),
             "The lengths of Y, Z, and of the row and column of A must be the same.")

expect_error(testinterference(Y                 = Y,
                              Z                 = seq(1, n, by = 1),
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
                              cores             = cores),
             "Z must not contain the values other than 0 and 1.")

# hypothesis -------------------------------------------------------------------

expect_error(testinterference(Y                 = Y,
                              Z                 = Z,
                              A                 = A,
                              hypothesis        = "hypothesis a",
                              method            = method,
                              design            = design,
                              prob              = prob,
                              focal_unit        = focal_unit,
                              num_focal_unit    = num_focal_unit,
                              num_randomization = num_randomization,
                              strata            = strata,
                              zmatrix           = zmatrix,
                              kappa             = kappa,
                              cores             = cores),
             "The argument 'hypothesis' must be one of 'Fisher', 'SUTVA', 'exposure1', and 'exposure2'.")

# method -----------------------------------------------------------------------

expect_error(testinterference(Y                 = Y,
                              Z                 = Z,
                              A                 = A,
                              hypothesis        = hypothesis,
                              method            = "other",
                              design            = design,
                              prob              = prob,
                              focal_unit        = focal_unit,
                              num_focal_unit    = num_focal_unit,
                              num_randomization = num_randomization,
                              strata            = strata,
                              zmatrix           = zmatrix,
                              kappa             = kappa,
                              cores             = cores),
             "The argument 'method' must be one of '3-net', '2-net', 'random', and 'manual'.")

# design -----------------------------------------------------------------------

expect_error(testinterference(Y                 = Y,
                              Z                 = Z,
                              A                 = A,
                              hypothesis        = hypothesis,
                              method            = method,
                              design            = "manual",
                              prob              = prob,
                              focal_unit        = focal_unit,
                              num_focal_unit    = num_focal_unit,
                              num_randomization = num_randomization,
                              strata            = strata,
                              zmatrix           = zmatrix,
                              kappa             = kappa,
                              cores             = cores),
             "The argument 'design' must be one of 'Bernoulli', 'complete', 'stratified', and 'other'.")

# prob -------------------------------------------------------------------------

expect_error(testinterference(Y                 = Y,
                              Z                 = Z,
                              A                 = A,
                              hypothesis        = hypothesis,
                              method            = method,
                              design            = "Bernoulli",
                              prob              = 100,
                              focal_unit        = focal_unit,
                              num_focal_unit    = num_focal_unit,
                              num_randomization = num_randomization,
                              strata            = strata,
                              zmatrix           = zmatrix,
                              kappa             = kappa,
                              cores             = cores),
             "The argument 'prob' must be a probability or an n-dimensional vector of probabilities when `design = 'Bernoulli'`.")

# focal_unit -------------------------------------------------------------------

expect_error(testinterference(Y                 = Y,
                              Z                 = Z,
                              A                 = A,
                              hypothesis        = hypothesis,
                              method            = "manual",
                              design            = design,
                              prob              = prob,
                              focal_unit        = NULL,
                              num_focal_unit    = num_focal_unit,
                              num_randomization = num_randomization,
                              strata            = strata,
                              zmatrix           = zmatrix,
                              kappa             = kappa,
                              cores             = cores),
             "The argument 'focal_unit' must be an n-dimensional logical vector when `method = 'manual'`.")

# num_focal_unit ---------------------------------------------------------------

# num_randomization ------------------------------------------------------------

expect_error(testinterference(Y                 = Y,
                              Z                 = Z,
                              A                 = A,
                              hypothesis        = hypothesis,
                              method            = method,
                              design            = design,
                              prob              = prob,
                              focal_unit        = focal_unit,
                              num_focal_unit    = num_focal_unit,
                              num_randomization = TRUE,
                              strata            = strata,
                              zmatrix           = zmatrix,
                              kappa             = kappa,
                              cores             = cores),
             "The argument 'num_randomization' must be a large positive integer.")

# strata -----------------------------------------------------------------------

expect_error(testinterference(Y                 = Y,
                              Z                 = Z,
                              A                 = A,
                              hypothesis        = hypothesis,
                              method            = method,
                              design            = "stratified",
                              prob              = prob,
                              focal_unit        = focal_unit,
                              num_focal_unit    = num_focal_unit,
                              num_randomization = num_randomization,
                              strata            = NULL,
                              zmatrix           = zmatrix,
                              kappa             = kappa,
                              cores             = cores),
             "The argument 'strata' must be an n-dimensional numerical vector when `design == 'stratified'`.")

# zmatrix ----------------------------------------------------------------------

expect_error(testinterference(Y                 = Y,
                              Z                 = Z,
                              A                 = A,
                              hypothesis        = hypothesis,
                              method            = method,
                              design            = "other",
                              prob              = prob,
                              focal_unit        = focal_unit,
                              num_focal_unit    = num_focal_unit,
                              num_randomization = num_randomization,
                              strata            = strata,
                              zmatrix           = 1:n,
                              kappa             = kappa,
                              cores             = cores),
             "The argument 'zmatrix' must be a large numerical matrix when `design == 'other'`.")

# kappa ------------------------------------------------------------------------

expect_error(testinterference(Y                 = Y,
                              Z                 = Z,
                              A                 = A,
                              hypothesis        = "exposure2",
                              method            = method,
                              design            = design,
                              prob              = prob,
                              focal_unit        = focal_unit,
                              num_focal_unit    = num_focal_unit,
                              num_randomization = num_randomization,
                              strata            = strata,
                              zmatrix           = zmatrix,
                              kappa             = 1,
                              cores             = cores),
             "The argument 'kappa' must be NULL or a positive integer no less than 2.")

# cores ------------------------------------------------------------------------
