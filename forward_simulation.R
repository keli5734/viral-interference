# source("two_pathogen_model.R") # upload the two-pathogen model (RSV-Influenza interaction)
# source("two_pathogen_model_initate_set_up.R") # setting up initial conditions and parameter values
# library(tidyverse)
# library(deSolve)



results <- ode(y = yinit.vector, t = my_times,
               func = two_pathogen_model,
               parms = my_parmset)


results.burned <- results[-c(1:burnN),]
results.burned <- results.burned[c(105:208), ]

results.burned.all <- apply(results.burned, 1, sum)

pop.all <- rowSums(results.burned[,-1])

infected.cols <- results.burned[,c( grep('Xi1s', colnames(results.burned))[1:21],
                                    grep('Xi2s', colnames(results.burned))[1:21],
                                    grep('Xi3s', colnames(results.burned))[1:21],
                                    grep('Xi4s', colnames(results.burned))[1:21],
                                    #grep('Xi1i', colnames(results.burned))[1:21],
                                    #grep('Xi2i', colnames(results.burned))[1:21],
                                    #grep('Xi3i', colnames(results.burned))[1:21],
                                    #grep('Xi4i', colnames(results.burned))[1:21],
                                    grep('Xi1r', colnames(results.burned))[1:21],
                                    grep('Xi2r', colnames(results.burned))[1:21],
                                    grep('Xi3r', colnames(results.burned))[1:21],
                                    grep('Xi4r', colnames(results.burned))[1:21],
                                    grep('Xi1s1', colnames(results.burned))[1:21],
                                    grep('Xi2s1', colnames(results.burned))[1:21],
                                    grep('Xi3s1', colnames(results.burned))[1:21],
                                    grep('Xi4s1', colnames(results.burned))[1:21])]
                                    #grep('Xi1i2', colnames(results.burned))[1:21],
                                    #grep('Xi2i2', colnames(results.burned))[1:21],
                                    #grep('Xi3i2', colnames(results.burned))[1:21],
                                    #grep('Xi4i2', colnames(results.burned))[1:21] )]


infected.cols.coinfection <- results.burned[,c(grep('Xi1i', colnames(results.burned))[1:21],
                                               grep('Xi2i', colnames(results.burned))[1:21],
                                               grep('Xi3i', colnames(results.burned))[1:21],
                                               grep('Xi4i', colnames(results.burned))[1:21],
                                               grep('Xi1i2', colnames(results.burned))[1:21],
                                               grep('Xi2i2', colnames(results.burned))[1:21],
                                               grep('Xi3i2', colnames(results.burned))[1:21],
                                               grep('Xi4i2', colnames(results.burned))[1:21])]




infected.cols.flu <-  results.burned[,c( grep('Xsi1', colnames(results.burned))[1:21],
                                         grep('Xs1i', colnames(results.burned))[1:21],
                                         grep('Xs2i', colnames(results.burned))[1:21],
                                         grep('Xs3i', colnames(results.burned))[1:21],
                                         #grep('Xi1i', colnames(results.burned))[1:21],
                                         #grep('Xi2i', colnames(results.burned))[1:21],
                                         #grep('Xi3i', colnames(results.burned))[1:21],
                                         #grep('Xi4i', colnames(results.burned))[1:21],
                                         grep('Xsi2', colnames(results.burned))[1:21],
                                         grep('Xs1i2', colnames(results.burned))[1:21],
                                         grep('Xs2i2', colnames(results.burned))[1:21],
                                         grep('Xs3i2', colnames(results.burned))[1:21])]
                                         #grep('Xi1i2', colnames(results.burned))[1:21],
                                         #grep('Xi2i2', colnames(results.burned))[1:21],
                                         #grep('Xi3i2', colnames(results.burned))[1:21],
                                         #grep('Xi4i2', colnames(results.burned))[1:21])]


susceptibles <-  results.burned[,c( grep('Xss', colnames(results.burned))[1:21],
                                              grep('Xs1s', colnames(results.burned))[1:21],
                                              grep('Xs2s', colnames(results.burned))[1:21],
                                              grep('Xs3s', colnames(results.burned))[1:21],
                                              
                                              grep('Xsi1', colnames(results.burned))[1:21],
                                              grep('Xs1i', colnames(results.burned))[1:21],
                                              grep('Xs2i', colnames(results.burned))[1:21],
                                              grep('Xs3i', colnames(results.burned))[1:21],
                                              
                                              grep('Xsr', colnames(results.burned))[1:21],
                                              grep('Xs1r', colnames(results.burned))[1:21],
                                              grep('Xs2r', colnames(results.burned))[1:21],
                                              grep('Xs3r', colnames(results.burned))[1:21],
                                         
                                               grep('Xss1', colnames(results.burned))[1:21],
                                               grep('Xs1s1', colnames(results.burned))[1:21],
                                               grep('Xs2s1', colnames(results.burned))[1:21],
                                               grep('Xs3s1', colnames(results.burned))[1:21],
                                               
                                               grep('Xsi2', colnames(results.burned))[1:21],
                                               grep('Xs1i2', colnames(results.burned))[1:21],
                                               grep('Xs2i2', colnames(results.burned))[1:21],
                                               grep('Xs3i2', colnames(results.burned))[1:21]
                                         )]

#grep('Xi1i2', colnames(results.burned))[1:21],
#grep('Xi2i2', colnames(results.burned))[1:21],
#grep('Xi3i2', colnames(results.burned))[1:21],
#grep('Xi4i2', colnames(results.burned))[1:21])]


infected.all <- apply(infected.cols,1, sum)
infected.all.coinfection <- apply(infected.cols.coinfection, 1, sum)
infected.all.flu <- apply(infected.cols.flu,1, sum)
susceptible.all <- apply(susceptibles,1, sum)


cases_prediction <- data.frame(rsv = infected.all,
                               flu = infected.all.flu,
                               coinfected = infected.all.coinfection,
                               susceptible = susceptible.all, 
                               pop = pop.all)
                               
                               

saved_ts[[n]] <- cases_prediction
#saved_ts_sep_rsv[[n]] <- infected.cols
#saved_ts_sep_coinf[[n]] <- infected.cols.coinfection

