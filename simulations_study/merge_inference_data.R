#############################################
## Bind files together several K
#############################################

saveresultfile = "../Results/Simulations_Several_K/several_K_estimations_SUN_rBM"
datestamp_day <- "2016-02-17" # 03-17 11-19 11-24

#simestimations_alpha_known <- NULL
simestimations <- NULL
ak <- "" # alpha_known_

for (inference.index in 1:4){
  file <- paste0(saveresultfile, ak, "-", datestamp_day, "_", inference.index, ".RData")
  if (file.exists(file)) {
    load(file)
    ## simestimations
    siminf <- as.name(paste0("simestimations_", ak, inference.index))
    simestimations <- c(simestimations, 
                                    eval(siminf))
    rm(list = paste0("simestimations_", ak, inference.index))
    datestamp_day <- "2016-02-17" # 03-17 11-19 11-24
  } else {
    warning(paste0("File number ", inference.index, " do not exists"))
  }
}

save.image(paste0(saveresultfile, ak, "-", datestamp_day, "_all", ".RData"))
