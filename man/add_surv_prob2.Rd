% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add-functions.R
\name{add_surv_prob2}
\alias{add_surv_prob2}
\title{Add survival probability estimates}
\usage{
add_surv_prob2(
  newdata,
  object,
  overwrite = FALSE,
  time_var = NULL,
  interval_length = "intlen",
  ...
)
}
\arguments{
\item{newdata}{ A data frame or list containing the values of the model covariates at which predictions
  are required. If this is not provided then predictions corresponding to the
                original data are returned. If \code{newdata} is provided then
                it should contain all the variables needed for prediction: a
                warning is generated if not. See details for use with \code{link{linear.functional.terms}}. }

\item{object}{ a fitted \code{gam} object as produced by \code{gam()}.
                }

\item{overwrite}{Should hazard columns be overwritten if already present in
the data set? Defaults to \code{FALSE}.}

\item{time_var}{Name of the variable used for the baseline hazard. If
not given, defaults to \code{"tend"} for \code{\link[mgcv]{gam}} fits, else
\code{"interval"}. The latter is assumed to be a factor, the former
numeric.}

\item{interval_length}{The variable in newdata containing the interval lengths.
Can be either bare unquoted variable name or character. Defaults to \code{"intlen"}.}

\item{...}{Further arguments passed to \code{\link[xgboost]{predict.xgb.Booster}}.}
}
\description{
Given suitable data (i.e. data with all columns used for estimation of the model),
this functions adds a column \code{surv_prob} containing survival probabilities
for the specified covariate and follow-up information.
}
