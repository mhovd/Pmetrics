require(Pmetrics)
cat(commandArgs(), "\n")

# Backup of old 
#wd <- commandArgs()[6]
#icen <- commandArgs()[7]
#parallel <- as.logical(commandArgs()[8])

# New argument indices
wd <- commandArgs()[1]
icen <- commandArgs()[2]
parallel <- as.logical(commandArgs()[3])

cat("NPrepScript thinks this is the wd:\n", wd, "\n")

setwd(wd)
PMreport(wd,icen=icen,type="NPAG",parallel=parallel)



