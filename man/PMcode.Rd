% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PMcode.R
\name{PMcode}
\alias{PMcode}
\title{Pmetrics GUI Tutor}
\usage{
PMcode(func)
}
\arguments{
\item{func}{Quoted name of a function family used in Pmetrics.  Currently, these are limited to \dQuote{run}, 
for \code{\link{NPrun}}, \code{\link{ITrun}} and \dQuote{plot}.  For the first two, make sure that the model and data files are in your
working directory before calling the function.}
}
\value{
Nothing is returned, but the user interface is launched in the default browser.  Appropriate R code to execute
Pmetrics commands is generated depending on defaults and user-selected input.  For plotting, the resulting plot is previewed
directly in the browser.
}
\description{
Learn Pmetrics R code with user friendly graphical interfaces in the default browser.
}
\details{
PMcode provides a graphical user interface to learn many of the Pmetrics functions and their arguments
using the Shiny package. A graphical user interface will launch in the default browser.  This GUI enables a point
and click approach to generating Pmetrics code (which can be pasted into the R script) and plot previews.
The idea is for users to learn the R code in an intuitive and easier manner.  There are more options available for Pmetrics
functions that are served by the GUI, but it is sufficiently powerful to serve basic needs.  To stop the shiny browser GUI, click
the stop buttton in Rstudio (upper left corner of console window) or ESC or CTRL-C may work when the R window is active.
}
\author{
Michael Neely
}
