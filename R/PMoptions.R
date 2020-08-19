#' Get user options for Pmetrics
#'
#' This function will get user options for Pmetrics.  Current user options are
#' \itemize{
#' \item sep Field separator in data files
#' \item dec Decimal separator in numbers
#' }  
#'
#' @title Get Pmetrics User Options
#' @param opt The option to retrieve.  If omitted, all option values will be returned.
#' @return The user options file will be updated.  This will persist from session to session.
#' @author Michael Neely
#' @export

getPMoptions <- function(opt) {

  #options file name
  PMoptionsFile <- paste(system.file("options", package = "Pmetrics"),"PMoptions.json", sep = "/")
  #if it doesn't exist, create it with defaults
  if (!file.exists(PMoptionsFile)) {
    baseuuid <- paste(sample(c(letters[1:6],0:9),30,replace=TRUE),collapse="")

    installation_code <- paste(
        substr(baseuuid,1,8),
        "-",
        substr(baseuuid,9,12),
        "-",
        "4",
        substr(baseuuid,13,15),
        "-",
        sample(c("8","9","a","b"),1),
        substr(baseuuid,16,18),
        "-",
        substr(baseuuid,19,30),
        sep="",
        collapse=""
    )
    PMopts <- list(sep = ",",
                   dec = ".",
                   server_address = "http://localhost:5000",
                   installation_code = installation_code)
    options(PMopts)
    jsonlite::write_json(PMopts, path = PMoptionsFile, auto_unbox=T)
  }
  #read the options file
  PMopts <- jsonlite::read_json(path = PMoptionsFile)
  if (missing(opt)) {
    return(PMopts)
  } else {
    index <- which(names(PMopts) == opt)
    if (length(index) == 0) { return(NULL) } else { return(PMopts[[index]]) }
  }
}

#' Set user options for Pmetrics
#'
#' This function will set user options for Pmetrics.  
#'
#' @title Set Pmetrics User Options
#' @param sep The field separator character; \dQuote{,} by default, but could be \dQuote{;}
#' @param dec The decimal separator character; \dQuote{.} by default, but could be \dQuote{,}
#' @param server_address Specify address of server for remote runs.  Server must be set up separately.
#' This functionality is coming soon.
#' @return The user preferences file will be updated.  This will persist from session to session.
#' @author Michael Neely
#' @export

setPMoptions <- function(sep, dec, server_address) {
  #read old values first
  PMopts <- getPMoptions()
  #update/add options
  if (!missing(sep)) PMopts$sep <- sep
  if (!missing(dec)) PMopts$dec <- dec
  if (!missing(server_address)) PMopts$server_address <- server_address
  PMopts$installation_code <- installation_code
  #set the options
  options(PMopts)
  #store the options
  PMoptionsFile <- paste(system.file("options", package = "Pmetrics"),"PMoptions.json", sep = "/")
  jsonlite::write_json(PMopts, path = PMoptionsFile, auto_unbox=T)
}

testFunc <- function(){
  return("this is a test\n")
}