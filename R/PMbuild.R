#' \code{PMBuild} will ensure all dependent packages are installed and compile
#' Fortran source code for permanent Pmetrics modules
#'
#' @title Build Pmetrics
#' @author Michael Neely
#' @export



PMbuild <- function(skipRegistration = F) {
  if (.check_and_install_gfortran(skipRegistration)) {

    currwd <- getwd()
    OS <- getOS()

    #load necessary packages
    # packages <- packageDescription("Pmetrics")$Suggests
    # packages <- gsub("\n", "", packages)
    # packages <- unlist(strsplit(packages, ","))
    # cat("\nChecking for required packages...\n")
    # for (i in packages) {
    #   if (system.file(package = i) == "") {
    #     if (getOption("repos")[1] == "") { setRepositories() }
    #     install.packages(i, repos = getOption("repos"), dependencies = T)
    #   }

    # }

    ODEsolver <- c("dvode_v1.f90") # stand alone package, w/module included in f90 file
    NPAGutils <- c("npag_utils.f90") # see detailed useage notes about line 80 below

    compiler <- PMFortranConfig()
    #try again just in case redefined
    compiler <- PMFortranConfig()
    #check if parallel is possible
    if (length(compiler) == 2 & getBits() == "64") {
      parallel <- T
    } else { parallel <- F }
    sourcedir <- system.file("code", package = "Pmetrics")
    destdir <- paste(system.file("", package = "Pmetrics"), "compiledFortran", sep = "/")
    #remove old files if present
    oldfiles <- c(Sys.glob(paste(destdir, "*.o", sep = "/")), Sys.glob(paste(destdir, "*.exe", sep = "/")))

    if (length(oldfiles) > 0) { file.remove(oldfiles) }
    #compile new files
    setwd(sourcedir)
    if (!file.exists(destdir)) dir.create(destdir, showWarnings = F)
    PMfiles <- data.frame(filename = as.character(c("NPprep"
    , "NPeng", "ITprep", "ITeng", "ITerr", "SIMeng"
    , "DOprep", "DOeng", "mb2csv")))
    PMfiles$path <- sapply(PMfiles$filename,
   function(x) shQuote(
    list.files(
      getwd(), pattern = as.character(paste(x, "_[[:digit:]]+\\.f", sep = ""))
    )
  ))

# Create ODEsolver object files (serial and parallel) -----------------
  serialCommand <- sub("<exec>", paste("s", "ODEsolver", ".o -c", sep = ""), compiler[1])
  serialCommand <- sub("<files>", ODEsolver, serialCommand)
  serialFortstatus <- suppressWarnings(system(serialCommand, intern = T, ignore.stderr = F))
  if (parallel) {
    parallelCommand <- sub("<exec>", paste("p", "ODEsolver", ".o -c", sep = ""), compiler[2])
    parallelCommand <- sub("<files>", ODEsolver, parallelCommand)
    parallelFortstatus <- suppressWarnings(system(parallelCommand, intern = T, ignore.stderr = F))
    if (!is.null(attr(parallelFortstatus, "status"))) {
      unlink(switch(OS, "~/.config/Pmetrics",
                    paste(Sys.getenv("APPDATA"), "\\Pmetrics", sep = ""),
                    "~/.config/Pmetrics"), recursive = T)
      stop(paste("\nThere was an error compiling "
        , ODEsolver
        , ".\nDid you select the right fortran compiler?  "
        , "If yes, try reinstalling fortran.\nFor gfortran, "
        , "log into www.lapk.org and access system-specific "
        , "tips on the Pmetrics installation page (step 5).\n", sep = ""))
    }
  }

  # Create C object files (serial OR parallel) -----------------
  # Try to remove the following gcc, it is from the first time that
  # J&W tried to envelope emint in a DLL and we decided to use a 
  # different approach.
  system(paste("gcc -c c_utils.c", sep =  " "), intern = T, ignore.stderr = F)

  # Create NPAGutils module ---------------------------------------------
  serialCommand <- sub("<exec>", paste("s", "npag_utils", ".o -c", sep = ""), compiler[1])
  serialCommand <- sub("<files>", NPAGutils, serialCommand);
  serialCommand
  serialFortstatus <- suppressWarnings(system(serialCommand, intern = T, ignore.stderr = F))
  if (parallel) {
    parallelCommand <- sub("<exec>", paste("p", "npag_utils", ".o -c", sep = ""), compiler[2])
    parallelCommand <- sub("<files>", NPAGutils, parallelCommand);
    parallelCommand
    parallelFortstatus <- suppressWarnings(system(parallelCommand, intern = T, ignore.stderr = F))
    if (!is.null(attr(parallelFortstatus, "status"))) {
      unlink(switch(OS, "~/.config/Pmetrics",
                    paste(Sys.getenv("APPDATA"), "\\Pmetrics", sep = ""),
                    "~/.config/Pmetrics"), recursive = T)
      stop(paste("\nThere was an error compiling "
        , NPAGutils
        , ".\nDid you select the right fortran compiler?  "
        , "If yes, try reinstalling fortran.\nFor gfortran, "
        , "log into www.lapk.org and access system-specific "
        , "tips on the Pmetrics installation page (step 5).\n", sep = ""))
    }
  }

  # Compile dot-o files ---------------------------------------------
  # Note: The above blocks make npag_utils.mod and compiles snpag_utils.o,
  #   pnpag_utils.o, sODEsolver.o and pODEsolver.o.  All five files will
  #   be moved to the compiled fortran directory (which is for most installations
  #   ~/.config/Pmetrics/compiledFortran/) along w/the below compiled engines.
  # Note: To compile the engines below, all you need are the *.mod referenced
  #   by USE statements in the Fortran code  ...
  # Note: ... But to run the code, you will have to link the s- or p- .o files
  #   to the engine. That could be done below, but it is easier to link them
  #   at the same time the model file is linked to the engine, in PMrun().
  # Note: ... with one exception: mb2csv is a standalone program, not 
  #   requiring linking later.  So for that program, linking is done here.

# Note: Transitioning to "use npag_utils" some program do NOT use npag_utils,
#  other do, but fall into two categories: run and prep. Prep programs will
#  not call the GET*, DIFFEQ, JACOB, DVODE, subroutines referenced from 
#  SimConc and USERANAL -- so these programs have dummy routines in them
#  IF they use npag_utils ELSE they will not link npag_utils.
#  TODO: make two npag_utils -- one for prep, and one for run.
#
    for (i in 1:nrow(PMfiles)) {
      cat(paste("\nCompiling ", i, " of ", nrow(PMfiles), ": ", PMfiles$filename[i], "...", sep = ""))
      flush.console()
      if (PMfiles$filename[i] %in% c("mb2csv")) { # "DOprep"
        #list of compiled and linked files
        serialCommand <- sub("<exec>", paste(PMfiles$filename[i], ".exe"
        , sep = ""), compiler[1])
      } else if (PMfiles$filename[i] %in% c("NPeng")) {
      serialCommand <- sub("<exec>", paste("s", PMfiles$filename[i]
         , ".o -c", sep = ""), compiler[1])
      } else {
        serialCommand <- sub("<exec>", paste("s", PMfiles$filename[i]
        , ".o -c", sep = ""), compiler[1])
      }
      serialCommand <- sub("<files>", PMfiles$path[i], serialCommand)
      #
      serialFortstatus <- suppressWarnings(system(serialCommand
      , intern = T, ignore.stderr = F))
      if (!is.null(attr(serialFortstatus, "status"))) {
        unlink(switch(OS, "~/.config/Pmetrics",
                    paste(Sys.getenv("APPDATA"), "\\Pmetrics", sep = ""),
                    "~/.config/Pmetrics"), recursive = T)
        stop(paste("\nThere was an error compiling ", PMfiles$filename[i]
        , ".\nDid you select the right fortran compiler?  If yes, try reinstalling fortran.\nFor gfortran, log into www.lapk.org and access system-specific tips on the Pmetrics installation page (step 5).\n", sep = ""))
      }
      if (i == 2 & parallel) {
        # parallel compilation for NPAG only
        parallelCommand <- sub("<exec>", paste("p", PMfiles$filename[i]
        , ".o -c", sep = ""), compiler[2])
        parallelCommand <- sub("<files>", PMfiles$path[i], parallelCommand)
        parallelFortstatus <- suppressWarnings(system(parallelCommand
        , intern = T, ignore.stderr = F))
        if (!is.null(attr(parallelFortstatus, "status"))) {
          unlink(switch(OS, "~/.config/Pmetrics",
                      paste(Sys.getenv("APPDATA"), "\\Pmetrics", sep = ""),
                      "~/.config/Pmetrics"), recursive = T)
          stop(paste("\nThere was an error compiling ", PMfiles$filename[i]
          , ".\nDid you select the right fortran compiler?  If yes, try reinstalling fortran.\nFor gfortran, log into www.lapk.org and access system-specific tips on the Pmetrics installation page (step 5).\n", sep = ""))
        }
      }

    }


    cat("\nAll packages installed and permanent Fortran modules compiled.\n")
    flush.console()
    invisible(file.copy(from = Sys.glob(c("*.o", "*.exe")), to = destdir))
    invisible(file.remove(Sys.glob(c("*.o", "*.exe"))))
    fort <- paste(system.file("config", package = "Pmetrics"), "newFort.txt", sep = "/")
    writeLines("0", fort) #reset to zero
    setwd(currwd)
  }
}

.check_and_install_gfortran <- function(skipRegistration) {
  #restore user defaults - deprecated
  #if(length(system.file(package="Defaults"))==1){PMreadDefaults()}
  sch_str <- c("which -s gfortran", "where gfortran", "which -s gfortran")
  OS <- getOS()
  env = Sys.getenv("env")
  if (env != "Development") {
    if (!binaries.installed()) {
      cat("Pmetrics cannot find required compiled binary files.\n")
      if (system(sch_str[OS]) != 0) {
        cat("Pmetrics cannot detect gfortran and will attempt to download and install all components.\n")
        input <- tolower(readline(prompt = "Do you agree? (Y/N)"))
        if (substr(input, 1, 1) == "y") {
          if (.installOrUpdateGfortran()) {
            cat("Pmetrics has installed gfortran and will now compile required binary files.\n")
            cat("Pmetrics has anonymously registered your installation of this version.\nLAPKB does not collect or store any personal or identifying information.")
            cat("If the registration time outs, please run PMbuild(skipRegistration=T) ")
            if(skipRegistration == F){
              .PMremote_registerNewInstallation()
            }
            return(T)
          } else {
            cat("ERROR: Pmetrics did not install gfortran automatically.\nPlease install gfortran manually and then run PMbuild().\nGo to http://www.lapk.org/Pmetrics_install.php for help.\n")
            return(F)
          }
        } else {
          cat("You must have gfortran to run Pmetrics.\nPlease install gfortran manually and then run PMbuild().\nGo to http://www.lapk.org/Pmetrics_install.php for help.\n")
          return(F)
        }
      } else {
        cat("Pmetrics has detected gfortran and will compile required binary files.\n")
        cat("Pmetrics has anonymously registered your installation of this version.\nLAPKB does not collect or store any personal or identifying information.")
        cat("If the registration time outs, please run PMbuild(skipRegistration=T) ")
        if(skipRegistration == F){
          .PMremote_registerNewInstallation()
        }
        return(T)
      }
    } else {
      cat("Pmetrics has found required compiled binary files.\n")
      return(F)
    }
  }
  else {
    print("You are inside the development folder, skipping gfortran installation")
    return(F)
  }
}
