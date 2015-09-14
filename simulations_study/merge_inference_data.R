#############################################
## Bind files together several K
#############################################

# saveresultfile = "../Results/Simulations_Several_K/several_K_estimations"
# datestamp_day <- "2015-03-17"

saveresultfile = "../Results/Simulations_Student/student_estimations"
datestamp_day <- "2015-09-09"
ak <- ""
student <- "_df5"
imax <- 100
favorable <- "favorables"
fav <- "_fav"

simests <- as.name(paste0("simestimations", favorable, ak, student))
assign(paste(simests), NULL)

for (inference.index in 1:imax){
  file <- paste0(saveresultfile, favorable, ak, "-", datestamp_day, "_", inference.index, ".RData")
  if (file.exists(file)) {
    load(file)
    ## simestimations
    siminf <- as.name(paste0("simestimations", fav, ak, "_", inference.index))
    assign(paste(simests), c(eval(simests), 
                             eval(siminf)))
    rm(list = paste0("simestimations", fav, ak, "_", inference.index))
  } else {
    warning(paste0("File number ", inference.index, " do not exists"))
  }
}

save.image(paste0(saveresultfile, favorable, ak, student, "-", datestamp_day, "_all", ".RData"))