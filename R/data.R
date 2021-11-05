#' Information on available distributions
#'
#' A named list containing information (parameters, contraints, notes, warnings, etc.)
#' for all available univariate distributions.
#'
#' @format A named list where each element is itself a list containing:
#' \describe{
#'   \item{parName}{parameters short names}
#'   \item{parLongName}{parameters long names}
#'   \item{parSymbol}{parameters typical symbols}
#'   \item{constraints}{constraints on parameters}
#'   \item{url}{link to more information}
#'   \item{note}{notes}
#'   \item{warning}{warnings: read carefully since this highlights in particular differences with "standard" parameterizations found in e.g. Wikipedia or R.}
#' }
"distInfo"

#' Default estimation options
#'
#' A named list containing the default estimation options.
#' See ?Hydro3_Estimation for more details.
"options_def"

#' Default MCMC options
#'
#' A named list containing the default MCMC options.
#' See ?Hydro3_Estimation for more details.
"mcmcoptions_def"
