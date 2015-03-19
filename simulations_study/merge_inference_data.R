#############################################
## Bind files together several K
#############################################

saveresultfile = "../Results/Simulations_Several_K/several_K_estimations"
datestamp_day <- "2015-03-17"

simestimations_alpha_known <- NULL

for (inference.index in 1:40){
  file <- paste0(saveresultfile, "_alpha_known-", datestamp_day, "_", inference.index, ".RData")
  if (file.exists(file)) {
    load(file)
    ## simestimations
    siminf <- as.name(paste0("simestimations_alpha_known_", inference.index))
    simestimations_alpha_known <- c(simestimations_alpha_known, 
                                    eval(siminf))
    rm(list = paste0("simestimations_alpha_known_", inference.index))
  } else {
    warning(paste0("File number ", inference.index, " do not exists"))
  }
}

save.image(paste0(saveresultfile, "_alpha_known-", datestamp_day, "_all", ".RData"))