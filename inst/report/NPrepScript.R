require(Pmetrics)
cat(commandArgs(), "\n")
wd <- commandArgs()[6]
icen <- commandArgs()[7]
parallel <- as.logical(commandArgs()[8])
cat("NPrepScript thinks this is the wd:\n", wd, "\n")
setwd(wd)
PMreport(wd,icen=icen,type="NPAG",parallel=parallel)



