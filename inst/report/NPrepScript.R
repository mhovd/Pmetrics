require(Pmetrics)
cat(commandArgs(), "\n")

# Backup of old 
#wd <- commandArgs()[6]
#icen <- commandArgs()[7]
#parallel <- as.logical(commandArgs()[8])

# Find the right argument
cat("Argument 1: ", commandArgs()[1], "\n")
cat("Argument 2: ", commandArgs()[2], "\n")
cat("Argument 3: ", commandArgs()[3], "\n")
cat("Argument 4: ", commandArgs()[4], "\n")
cat("Argument 5: ", commandArgs()[5], "\n")
cat("Argument 6: ", commandArgs()[6], "\n")
cat("Argument 7: ", commandArgs()[7], "\n")
cat("Argument 8: ", commandArgs()[8, "\n")
cat("Argument 9: ", commandArgs()[9], "\n")
cat("Argument 10: ", commandArgs()[10], "\n")
                                  
# New argument indices
wd <- commandArgs()[8]
icen <- commandArgs()[9]
parallel <- as.logical(commandArgs()[10])

cat("NPrepScript thinks this is the wd:\n", wd, "\n")

setwd(wd)
PMreport(wd,icen=icen,type="NPAG",parallel=parallel)



