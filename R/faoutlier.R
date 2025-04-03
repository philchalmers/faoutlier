#' Influential case detection methods for factor analysis and SEM
#'
#' Implements robust Mahalanobis methods, generalized Cook's distances,
#' likelihood ratio tests, model implied residuals, and various
#' graphical methods to help detect and summarize influential
#' cases that can affect exploratory and confirmatory factor analyses.
#'
#' @name faoutlier
#' @aliases faoutlier-package
#' @title Influential case detection methods for FA and SEM
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#'
#' Chalmers, R. P. & Flora, D. B. (2015). faoutlier: An R Package for Detecting
#'   Influential Cases in Exploratory and Confirmatory Factor Analysis.
#'   \emph{Applied Psychological Measurement, 39}, 573-574. \doi{10.1177/0146621615597894}
#'
#' Flora, D. B., LaBrish, C. & Chalmers, R. P. (2012). Old and new ideas for data
#' screening and assumption testing for exploratory and confirmatory factor analysis.
#'  \emph{Frontiers in Psychology, 3}, 1-21. \doi{10.3389/fpsyg.2012.00055}
#' @import stats MASS parallel lattice mvtnorm graphics sem
#' @importFrom lavaan logLik
#' @importFrom methods is
#' @importFrom pbapply pblapply pbapply
#' @importFrom utils flush.console tail
#' @keywords package
NULL

#' Description of holzinger data
#'
#' A sample of 100 simulated cases from the infamous Holzinger dataset
#' using 9 variables.
#'
#'
#' @name holzinger
#' @references
#'
#' Chalmers, R. P. & Flora, D. B. (2015). faoutlier: An R Package for Detecting
#'   Influential Cases in Exploratory and Confirmatory Factor Analysis.
#'   \emph{Applied Psychological Measurement, 39}, 573-574. \doi{10.1177/0146621615597894}
#'
#' Flora, D. B., LaBrish, C. & Chalmers, R. P. (2012). Old and new ideas for data
#' screening and assumption testing for exploratory and confirmatory factor analysis.
#'  \emph{Frontiers in Psychology, 3}, 1-21. \doi{10.3389/fpsyg.2012.00055}
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords data
NULL

#' Description of holzinger data with 1 outlier
#'
#' A sample of 100 simulated cases from the infamous Holzinger dataset
#' using 9 variables, but with 1 outlier added to the dataset. The first row was replaced
#' by adding 2 to five of the observed variables (odd-numbered items) and subtracting 2 from
#' the other four observed variables (even-numbered items).
#'
#'
#' @name holzinger.outlier
#' @references
#'
#' Chalmers, R. P. & Flora, D. B. (2015). faoutlier: An R Package for Detecting
#'   Influential Cases in Exploratory and Confirmatory Factor Analysis.
#'   \emph{Applied Psychological Measurement, 39}, 573-574. \doi{10.1177/0146621615597894}
#'
#' Flora, D. B., LaBrish, C. & Chalmers, R. P. (2012). Old and new ideas for data
#' screening and assumption testing for exploratory and confirmatory factor analysis.
#'  \emph{Frontiers in Psychology, 3}, 1-21. \doi{10.3389/fpsyg.2012.00055}
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords data
NULL
