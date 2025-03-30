library(doParallel)
source("func_lib.R")
source("mdl.R")

file_path <- list.files("../data/simul", full.names=TRUE,
                        pattern = "scenario1")
                        # pattern="data_mult.*p20_k3.*rep50\\.RData$")
num_iter <- 25000
num_burn <- 10000

for (file in file_path){
  load(file)
  num_rep <- length(data)
  # Parallel computation for 50 replications
  registerDoParallel(50)
  output <- foreach(rep = 1:num_rep,
                    .packages = c("LaplacesDemon",
                                  "tidyverse",
                                  "FactoMineR")
                    ) %dopar%
    {Y <- data[[rep]]$Y
    Lambda0 <- data[[rep]]$Lambda
    eta0 <- data[[rep]]$eta
    Sigma0 <- data[[rep]]$Sigma
    sparse_bayes(Y, num_iter, num_burn, 
                 eta0=eta0, 
                 Lambda0=Lambda0, 
                 Sigma0=Sigma0)}
  stopImplicitCluster()

  # Reorganize the results
  RV_Lambda = RV_eta = theta_est = k_est <- NULL
  for (i in 1:length(output)){
    RV_Lambda <- c(RV_Lambda, output[[i]]$RV_Lambda)
    RV_eta <- c(RV_eta, output[[i]]$RV_eta)
    # theta_est <- c(theta_est, output[[i]]$theta_est)
    k_est <- c(k_est, output[[i]]$k_est)
  }
  output_name <- sprintf("../results/MGPS_%s_iter%d_estK.Rdata", 
                         gsub("^../data/simul_data_|.RData$", "", file), 
                         num_iter)
  save(RV_Lambda, RV_eta, k_est, file=output_name)
}