#####################################################
### Scenario 1: Comparison in the Small K setting ###
#####################################################

####################################################################
# Goal 1: finite-sample performance for the vector estimator	   #
# average of estimation bias (Bias)								   #
# Monte Carlo standard deviation (MCSD)							   #
# average of standard errors (ASE)								   #
# coverage probability (CP) of 95% Wald-type confidence interval   #
#																   #
# Goal 2: evaluation of the var/covariance estimators	  		   #
# average of estimated covariance (AC)			      			   #
# real values in the formula (RC)					               #
# Monte Carlo covariance (MCC)						               #
#																   #
# Goal 3: evaluation of covariance- and precision-matrix estimates #
#         in Frobenius and spectral norms					       #
####################################################################

############################
### Loading the packages ###
############################

REQUIRED_PACKAGES <- c(
	"MASS", 
	### MASS::mvrnorm(), generate multivariate normal random samples
	"POET",		
	### POET, apply POET method for covariance matrix estimation
	"CovTools"	
	### CovTools, apply various methods for covariance/precision matrix estimation
)

CHECK_PACKAGES <- lapply(X = REQUIRED_PACKAGES,
					   FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
    require(x, character.only = TRUE)
  }
)

########################################################
### Loading required functions for Goal 1 and Goal 3 ###
########################################################

CALCULATE_VAR_A_B <- function(A, B, p_vec, n){
	
	K <- length(p_vec)
	VAR_A_MAT <- diag(K)
	VAR_B_MAT <- matrix(0, K, K)
	
	for(k in 1 : K){
		for(kp in 1 : K){
			if(kp == k){
				VAR_A_MAT[k, kp] <- 2 * A[k, k] ^ 2 / ((n - 1) * (p_vec[k] - 1))
				VAR_B_MAT[k, kp] <- 2 * ((A[k, k] + p_vec[k] * B[k, k]) ^ 2 - B[k, k] * (2 * A[k, k] + p_vec[k] * B[k, k])) / ((n - 1) * (p_vec[k] - 1) * p_vec[k])
			} else {
				VAR_A_MAT[k, kp] <- 0
				VAR_B_MAT[k, kp] <- (p_vec[k] * p_vec[kp] * (B[k, kp] ^ 2 + B[kp, k] ^ 2) + 2 * (A[k, k] + p_vec[k] * B[k, k]) * (A[kp, kp] + p_vec[kp] * B[kp, kp])) / (2 * (n - 1) * p_vec[k] * p_vec[kp])
			}
		}
	}
	return(list(VAR_A_MAT = VAR_A_MAT, VAR_B_MAT = VAR_B_MAT))
}

BLOCK_HADAMARD_PRODUCT <- function(A, B, p_vec){
	K <- length(p_vec)
	COL_temp <- c()
	for(k in 1 : K){
		ROW_temp <- c()
		for(kp in 1 : K){
			if(k == kp){
				ROW_temp <- cbind(ROW_temp, A[k, k] * diag(p_vec[k]) + B[k, k] * matrix(1, p_vec[k], p_vec[k]))
			} else{
				ROW_temp <- cbind(ROW_temp, B[k, kp] * matrix(1, p_vec[k], p_vec[kp]))
			}
		}
		COL_temp <- rbind(COL_temp, ROW_temp)
	}
	M <- COL_temp
	return(M)
}

BEST_UNBIASED_ESTIMATOR <- function(S, p_vec){
	## this version can deal with NA values
	K <- length(p_vec)
	A_temp <- matrix(0, K, K)
	B_temp <- matrix(0, K, K)
	for(k in 1 : K){
		for(kp in 1 : K){
			SUB <- S[(sum(p_vec[0 : (k - 1)]) + 1) : sum(p_vec[0 : k]), 
						  (sum(p_vec[0 : (kp - 1)]) + 1) : sum(p_vec[0 : kp])]
			if(kp == k){
				SUB_ON <- mean(diag(SUB), na.rm = TRUE)
				SUB_OF <- (sum(SUB[upper.tri(SUB)], na.rm = TRUE) + sum(SUB[lower.tri(SUB)], na.rm = TRUE)) / (p_vec[k] ^ 2 - p_vec[k] - sum(is.na(SUB[upper.tri(SUB)])) - sum(is.na(SUB[lower.tri(SUB)])))
				A_temp[k, kp] <- SUB_ON - SUB_OF
				B_temp[k, kp] <- SUB_OF
			} else{
				B_temp[k, kp] <- mean(as.vector(SUB), na.rm = TRUE)
			}	
		}
	}
	A <- A_temp
	B <- (B_temp + t(B_temp)) / 2
	return(list(A = A, B = B))
}

#############################################
### Loading required functions for Goal 2 ###
#############################################

CALCULATE_COV_ALPHA_ALPHA <- function(A, B, p_vec, n){
	K <- length(p_vec)
	ALPHA_VAR_MAT <- matrix(0, K, K)
	for(k in seq(K)){
		for(kp in seq(K)){
			if(kp == k){
				ALPHA_VAR_MAT[k, kp] <- 2 * (A[k,k] ^ 2 + 2 * A[k,k] * B[k,k] + p_vec[k] * B[k,k] ^ 2) / ((n - 1) * p_vec[k])
			} else{
				ALPHA_VAR_MAT[k, kp] <- 2 * B[k,kp] * B[kp, k] / (n - 1)
			}
		}
	}
	return(ALPHA_VAR_MAT)
}

CALCULATE_COV_BETA_BETA <- function(A, B, p_vec, n){
	
	Loc <- function(k, kp, K) return((k - 1) * K + kp)
	
	K <- length(p_vec)
	
	BETA_VAR_MAT <- matrix(0, K ^ 2, K ^ 2)
	
	for(k in seq(K)){
		for(kp in seq(K)){
			Loc_temp_1 <- Loc(k, kp, K)
			Loc_temp_2 <- Loc(k, kp, K)
			if(kp == k){
				BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] <- 2 / ((n - 1) * p_vec[k] ^ 2) * (A[k, k] + p_vec[k] * B[k, k]) ^ 2
			} else{
				BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] <- 1 / (2 * (n - 1) * p_vec[k] * p_vec[kp]) * (p_vec[k] * p_vec[kp] * (B[k,kp]^2 + B[kp,k] ^ 2) + 2 * (A[k, k] + p_vec[k] * B[k,k]) * (A[kp,kp] + p_vec[kp] * B[kp, kp]))
			}
		}
	}
	
	for(k in seq(K)){
		for(kp in seq(K)){
			Loc_temp_1 <- Loc(k, k, K)
			Loc_temp_2 <- Loc(kp, kp, K)
			if(kp != k){
				BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] <- 2 / (n - 1) * B[k, kp] * B[kp, k]
				BETA_VAR_MAT[Loc_temp_2, Loc_temp_1] <- BETA_VAR_MAT[Loc_temp_1, Loc_temp_2]
			}
		}
	}

	for(k in seq(K)){
		for(kp in seq(K)){
			for(kpp in seq(K)){
				Loc_temp_1 <- Loc(k, k, K)
				Loc_temp_2 <- Loc(kp, kpp, K)
				
				if(kp != kpp){
					if(k == kp){
						BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] <- 1 / ((n - 1) * p_vec[k]) * (A[k,k] + p_vec[k] * B[k,k]) * (B[kp,kpp] + B[kpp, kp])
						BETA_VAR_MAT[Loc_temp_2, Loc_temp_1] <- BETA_VAR_MAT[Loc_temp_1, Loc_temp_2]
					} else if(k == kpp){
						BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] <- 1 / ((n - 1) * p_vec[k]) * (A[k,k] + p_vec[k] * B[k,k]) * (B[kp,kpp] + B[kpp, kp])
						BETA_VAR_MAT[Loc_temp_2, Loc_temp_1] <- BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] 
					} else{
						BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] <- 1 / (n - 1) * (B[kp, k] * B[k, kpp] + B[kpp, k] * B[k, kp])
						BETA_VAR_MAT[Loc_temp_2, Loc_temp_1] <- BETA_VAR_MAT[Loc_temp_1, Loc_temp_2]
					} 
				
				}
				
			}
		}
	}

	for(k1 in seq(K)){
		for(k2 in seq(K)){
			for(l1 in seq(K)){
				for(l2 in seq(K)){
					Loc_temp_1 <- Loc(k1, k2, K)
					Loc_temp_2 <- Loc(l1, l2, K)
					
					if(k1 != k2 & l1 != l2){
						if(k1 == l1 & k2 == l2){
							### (2-2)
							BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] <- 1 / (2 * (n - 1)) * ((B[l1, l2] ^ 2 + B[l2, l1] ^ 2) + 2 * (A[l2, l2] + p_vec[l2] * B[l2, l2]) * (A[l1, l1] + p_vec[l1] * B[l1, l1]) / (p_vec[l2] * p_vec[l1]))
						} else if(k1 == l2 & k2 == l1){
							### (3-2)
							BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] <- 1 / (2 * (n - 1)) * ((B[l2, l1] ^ 2 + B[l1, l2] ^ 2) + 2 * (A[l2, l2] + p_vec[l2] * B[l2, l2]) * (A[l1, l1] + p_vec[l1] * B[l1, l1]) / (p_vec[l2] * p_vec[l1]))
						} else if(k1 == l2 & k2 != l1 & k2 != l2){
							### (3-1)
							BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] <- 1 / (2 * (n - 1)) * ((A[k1, k1] * B[l1, k2] + A[k1, k1] * B[k2, l1]) / p_vec[l2] + B[l1, k2] * B[k1, l2] + B[l2, k2] * B[k1, l1] + B[l1, k1] * B[k2, l2] + B[l2, k1] * B[k2, l1])
						} else if(k1 == l1 & k2 != l1 & k2 != l2){
							### (2-1)
							BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] <- 1 / (2 * (n - 1)) * ((A[k1, k1] * B[l2, k2] + A[k1, k1] * B[k2, l2]) / p_vec[l1] + B[l1, k2] * B[k1, l2] + B[l2, k2] * B[k1, l1] + B[l1, k1] * B[k2, l2] + B[l2, k1] * B[k2, l1])
						} else if(k2 == l1 & k1 != l1 & k1 != l2){
							### (1-2)
							BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] <- 1 / (2 * (n - 1)) * ((A[k2, k2] * B[l2, k1] + A[k2, k2] * B[k1, l2]) / p_vec[l1] + B[l1, k1] * B[k2, l2] + B[l2, k1] * B[k2, l1] + B[l1, k2] * B[k1, l2] + B[l2, k2] * B[k1, l1])
						} else if(k2 == l2 & k1 != l1 & k1 != l2){
							### (1-3)
							BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] <- 1 / (2 * (n - 1)) * ((A[k2, k2] * B[l1, k1] + A[k2, k2] * B[k1, l1]) / p_vec[l2] + B[l1, k1] * B[k2, l2] + B[l2, k1] * B[k2, l1] + B[l1, k2] * B[k1, l2] + B[l2, k2] * B[k1, l1])
						} else {
							### (1-1)
							BETA_VAR_MAT[Loc_temp_1, Loc_temp_2] <- 1 / (2 * (n - 1)) * (B[l1, k1] * B[k2, l2] + B[l2, k1] * B[k2, l1] + B[l1, k2] * B[k1, l2] + B[l2, k2] * B[k1, l1])
						}
					### End of If
					}
				
				}
			}
		}
	}
	
	return(BETA_VAR_MAT)
}

CALCULATE_COV_ALPHA_BETA <- function(A, B, p_vec, n){
	Loc <- function(k, kp, K) return((k - 1) * K + kp)
	
	K <- length(p_vec)
	
	ALPHA_BETA_VAR_MAT <- matrix(0, K, K ^ 2)
	
	for(k in seq(K)){
		for(kp in seq(K)){
			Loc_temp_2 <- Loc(kp, kp, K)
			if(kp == k){
				ALPHA_BETA_VAR_MAT[k, Loc_temp_2] <- 2 / ((n - 1) * p_vec[k] ^ 2) * (A[k, k] + p_vec[k] * B[k, k]) ^ 2
			} else{
				ALPHA_BETA_VAR_MAT[k, Loc_temp_2] <- 2 / (n - 1) * B[k, kp] * B[kp, k]
			}
		}
	}
	
	for(k in seq(K)){
		for(kp in seq(K)){
			for(kpp in seq(K)){
				Loc_temp_2 <- Loc(kp, kpp, K)
				
				if(kp != kpp){
					if(k == kp){
						ALPHA_BETA_VAR_MAT[k, Loc_temp_2] <- 1 / ((n - 1) * p_vec[kp]) * (B[kp, kpp] + B[kpp, kp]) * (A[k, k] + p_vec[k] * B[k, k])
					} else if(k == kpp){
						ALPHA_BETA_VAR_MAT[k, Loc_temp_2] <- 1 / ((n - 1) * p_vec[kpp]) * (B[kp, kpp] + B[kpp, kp]) * (A[k, k] + p_vec[k] * B[k, k])
					} else{
						ALPHA_BETA_VAR_MAT[k, Loc_temp_2] <- 1 / (n - 1) * (B[kp, k] * B[k, kpp] + B[kpp, k] * B[k, kp])
					}
				}
			}
		}
	}

	return(ALPHA_BETA_VAR_MAT)
}

BEST_UNBIASED_ESTIMATOR_ALPHA_BETA <- function(S, p_vec){
	K <- length(p_vec)
	ALPHA_temp <- matrix(0, K, K)
	BETA_temp <- matrix(0, K, K)
	for(k in 1 : K){
		for(kp in 1 : K){
			SUB <- S[(sum(p_vec[0 : (k - 1)]) + 1) : sum(p_vec[0 : k]), 
						  (sum(p_vec[0 : (kp - 1)]) + 1) : sum(p_vec[0 : kp])]
			if(kp == k){
				SUB_ON <- mean(diag(SUB))
				SUB_ALL <- sum(SUB) / (p_vec[k] ^ 2)
				ALPHA_temp[k, kp] <- SUB_ON
				BETA_temp[k, kp] <- SUB_ALL
			} else{
				BETA_temp[k, kp] <- mean(as.vector(SUB))
			}	
		}
	}
	A <- ALPHA_temp
	B <- (BETA_temp + t(BETA_temp)) / 2
	return(list(A = A, B = B))
}

############################
### Setup for Scenario 1 ###
############################

n <- 50 			###  50, 100, 150
p_ind <- 30 		###  30,  45,  60
p_vec <- rep(p_ind, 5)
K <- length(p_vec)  
p <- sum(p_vec)		### 150, 225, 300
thetaA <- c(0.01595042, 0.21392707, 0.74912381, 0.06771268, 0.10017260)
thetaB <- rep(0, K * K)
thetaB <- c(6.73139386,-1.69034339, 0.69591280,-2.93647430, 1.91315819, 
		   -1.69034339, 5.21462208, 3.81497235,-1.01011751, 0.70298054, 
            0.69591280, 3.81497235, 4.32780351,-3.35737580,-0.26890092,
		   -2.93647430,-1.01011751,-3.35737580, 6.78768893, 0.00018746,
            1.91315819, 0.70298054,-0.26890092, 0.00018746, 3.95418249)
indexB <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 13, 14, 15, 19, 20, 25)
A0 <- diag(thetaA)
B0 <- matrix(thetaB, K, K)
theta0 <- c(thetaA, thetaB[indexB])
mu0 <- rep(0, p)
SIGMA0 <- BLOCK_HADAMARD_PRODUCT(A0, B0, p_vec)
DELTA0 <- A0 + B0 %*% diag(p_vec)
OMEGA0 <- BLOCK_HADAMARD_PRODUCT(solve(A0), - solve(DELTA0) %*% B0 %*% solve(A0), p_vec)

ALPHA <- 0.05
reps <- 1
repsMax <- 1000
set.seed(2021)

### Setting variables for Bias, MCSD, ASE and WCP for Goal 1 ###

theta_EST <- matrix(0, repsMax, length(theta0))
theta_ASE <- matrix(0, repsMax, length(theta0))
theta_WCP <- matrix(0, repsMax, length(theta0))

### Setting variables for AC, RC, MCC for Goal 2 ###

AlphaBetaMC_MAT <- matrix(0, repsMax, K + K ^ 2)
AlphaAlphaMCcov <- matrix(0, K, K)
BetaBetaMCcov <- matrix(0, K ^ 2, K ^ 2)
AlphaBetaMCcov <- matrix(0, K, K ^ 2)
AlphaAlphaAcov <- matrix(0, K, K)
BetaBetaAcov <- matrix(0, K ^ 2, K ^ 2)
AlphaBetaAcov <- matrix(0, K, K ^ 2)

### Setting variables for AC, RC, MCC for Goal 3 ###

NICE_COV <- SOFT_COV <- HARD_COV <- ADAP_COV <- POET_COV <- NICE_PRE <- GLAS_PRE <- BAND_PRE <- BAY4_PRE <- BAY7_PRE <- matrix(0, repsMax, 2) 
# [, 1] = Frobenius, [, 2] = spectral
NICE_COV_TIME <- SOFT_COV_TIME <- HARD_COV_TIME <- ADAP_COV_TIME <- POET_COV_TIME <- NICE_PRE_TIME <- GLAS_PRE_TIME <- BAND_PRE_TIME <- BAY4_PRE_TIME <- BAY7_PRE_TIME <- 0 

###########################
### Computing procedure ###
###########################

while(reps <= repsMax){
	tryCatch({

	DATA <- mvrnorm(n = n, mu = mu0, Sigma = SIGMA0)
	S <- cov(DATA)
	
	### for Goal 1 ###

	RES_theta_tilde <- BEST_UNBIASED_ESTIMATOR(S, p_vec)
	A_tilde <- RES_theta_tilde$A
	B_tilde <- RES_theta_tilde$B
	theta_tilde_temp <- c(diag(A_tilde), (as.vector(B_tilde))[indexB])
	theta_EST[reps, ] <- theta_tilde_temp
	
	RES_variance <- CALCULATE_VAR_A_B(A_tilde, B_tilde, p_vec, n)
	VAR_A_MAT <- RES_variance$VAR_A_MAT
	VAR_B_MAT <- RES_variance$VAR_B_MAT
	theta_VAR_temp <- c(diag(VAR_A_MAT), (as.vector(VAR_B_MAT))[indexB])
	theta_ASE[reps, ] <- sqrt(theta_VAR_temp)
	
	CI_LOWER_temp <- theta_tilde_temp - qnorm(1 - ALPHA / 2) * sqrt(theta_VAR_temp)
	CI_UPPER_temp <- theta_tilde_temp + qnorm(1 - ALPHA / 2) * sqrt(theta_VAR_temp)
	CI_temp <- rbind(CI_LOWER_temp, CI_UPPER_temp)
	theta_WCP_temp <- 1 * (CI_temp[1, ] < theta0 & CI_temp[2, ] > theta0)
	theta_WCP[reps, ] <- theta_WCP_temp
	
	### for Goal 2 ###
	
	RES_S_AlphaBeta_tilde <- BEST_UNBIASED_ESTIMATOR_ALPHA_BETA(S, p_vec)
	AlphaBeta_tilde_temp <- c(diag(RES_S_AlphaBeta_tilde$A), as.vector(RES_S_AlphaBeta_tilde$B))
	AlphaBetaMC_MAT[reps, ] <- AlphaBeta_tilde_temp

	ALPHA_ALPHA_COV_temp <- CALCULATE_COV_ALPHA_ALPHA(A_tilde, B_tilde, p_vec, n)
	BETA_BETA_COV_temp <- CALCULATE_COV_BETA_BETA(A_tilde, B_tilde, p_vec, n)
	ALPHA_BETA_COV_temp <- CALCULATE_COV_ALPHA_BETA(A_tilde, B_tilde, p_vec, n)
	
	AlphaAlphaAcov <- AlphaAlphaAcov + ALPHA_ALPHA_COV_temp
	BetaBetaAcov <- BetaBetaAcov + BETA_BETA_COV_temp
	AlphaBetaAcov <- AlphaBetaAcov + ALPHA_BETA_COV_temp
		
	### for Goal 3 ###

	TIME_TEMP <- Sys.time()
	S_TEMP <- cov(DATA)
	RES_NICE <- BEST_UNBIASED_ESTIMATOR(S_TEMP, p_vec)
	COV_A_NICE <- RES_NICE$A
	COV_B_NICE <- RES_NICE$B
	SIGMA_NICE <- BLOCK_HADAMARD_PRODUCT(COV_A_NICE, COV_B_NICE, p_vec)
	NICE_COV[reps, 1] <- norm(SIGMA_NICE - SIGMA0, "F")
	NICE_COV[reps, 2] <- norm(SIGMA_NICE - SIGMA0, "2")
	NICE_COV_TIME <- NICE_COV_TIME + (Sys.time() - TIME_TEMP)
	
	TIME_TEMP <- Sys.time()
	S_TEMP <- cov(DATA)
	RES_NICE <- BEST_UNBIASED_ESTIMATOR(S_TEMP, p_vec)
	COV_A_NICE <- RES_NICE$A
	COV_B_NICE <- RES_NICE$B
	PRE_A_NICE <- solve(COV_A_NICE)
	PRE_B_NICE <- - solve(COV_A_NICE + COV_B_NICE %*% diag(p_vec)) %*% COV_B_NICE %*% solve(COV_A_NICE)
	OMEGA_NICE <- BLOCK_HADAMARD_PRODUCT(PRE_A_NICE, PRE_B_NICE, p_vec)
	NICE_PRE[reps, 1] <- norm(OMEGA_NICE - OMEGA0, "F")
	NICE_PRE[reps, 2] <- norm(OMEGA_NICE - OMEGA0, "2")
	NICE_PRE_TIME <- NICE_PRE_TIME + (Sys.time() - TIME_TEMP)
		
	TIME_TEMP <- Sys.time()	
	RES_SOFT <- CovEst.soft(X = DATA, thr = exp(seq(from = log(0.1), to = log(10), length.out = 10)))
	SOFT_COV[reps, 1] <- norm(RES_SOFT$S - SIGMA0, "F")
	SOFT_COV[reps, 2] <- norm(RES_SOFT$S - SIGMA0, "2")
	SOFT_COV_TIME <- SOFT_COV_TIME + (Sys.time() - TIME_TEMP)
	
	TIME_TEMP <- Sys.time()
	RES_HARD <- CovEst.hard(X = DATA, thr = exp(seq(from = log(0.1), to = log(10), length.out = 10)))
	HARD_COV[reps, 1] <- norm(RES_HARD$S - SIGMA0, "F")
	HARD_COV[reps, 2] <- norm(RES_HARD$S - SIGMA0, "2")
	HARD_COV_TIME <- HARD_COV_TIME + (Sys.time() - TIME_TEMP)
	
	TIME_TEMP <- Sys.time()
	RES_ADAP <- CovEst.adaptive(X = DATA, thr = exp(seq(from = log(0.1), to = log(10), length.out = 10)))
	ADAP_COV[reps, 1] <- norm(RES_ADAP$S - SIGMA0, "F")
	ADAP_COV[reps, 2] <- norm(RES_ADAP$S - SIGMA0, "2")
	ADAP_COV_TIME <- ADAP_COV_TIME + (Sys.time() - TIME_TEMP)
	
	TIME_TEMP <- Sys.time()
	K_HAT <- POETKhat(t(DATA))$K1HL ### p by n, each row has zero mean
	RES_POET <- POET(Y = t(DATA), K = K_HAT, C = 0.5, thres = "soft", matrix = "vad")
	POET_COV[reps, 1] <- norm(RES_POET$SigmaY - SIGMA0, "F")
	POET_COV[reps, 2] <- norm(RES_POET$SigmaY - SIGMA0, "2")
	POET_COV_TIME <- POET_COV_TIME + (Sys.time() - TIME_TEMP)
	
	TIME_TEMP <- Sys.time()
	RES_GLAS <- PreEst.glasso(X = DATA, method = list(type = "confidence",param = 0.95))
	GLAS_PRE[reps, 1] <- norm(RES_GLAS$C - OMEGA0, "F")
	GLAS_PRE[reps, 2] <- norm(RES_GLAS$C - OMEGA0, "2")
	GLAS_PRE_TIME <- GLAS_PRE_TIME + (Sys.time() - TIME_TEMP)
	
	TIME_TEMP <- Sys.time()
	RES_BAND <- PreEst.2014An(X = DATA, upperK = 2, algorithm = "Bonferroni", alpha = 0.05)
	BAND_PRE[reps, 1] <- norm(RES_BAND$C - OMEGA0, "F")
	BAND_PRE[reps, 2] <- norm(RES_BAND$C - OMEGA0, "2")
	BAND_PRE_TIME <- BAND_PRE_TIME + (Sys.time() - TIME_TEMP)
	
	TIME_TEMP <- Sys.time()
	RES_BAY4 <- PreEst.2014Banerjee(X = DATA, upperK = 2, delta = 10, logpi = function(k) { -k^4 }, loss = "Stein")
	BAY4_PRE[reps, 1] <- norm(RES_BAY4$C - OMEGA0, "F")
	BAY4_PRE[reps, 2] <- norm(RES_BAY4$C - OMEGA0, "2")
	BAY4_PRE_TIME <- BAY4_PRE_TIME + (Sys.time() - TIME_TEMP)
	
	TIME_TEMP <- Sys.time()
	RES_BAY7 <- PreEst.2017Lee(X = DATA, upperK = 2, logpi = function(k) {-k^4})
	BAY7_PRE[reps, 1] <- norm(RES_BAY7$C - OMEGA0, "F")
	BAY7_PRE[reps, 2] <- norm(RES_BAY7$C - OMEGA0, "2")
	BAY7_PRE_TIME <- BAY7_PRE_TIME + (Sys.time() - TIME_TEMP)
	
	cat(" iteration:  ", reps, "\r")
	reps <- reps + 1
	}, error = function(e){})
}

save.image("Scenario_1_SampleSize_50_Dimention_150_PC_1000.RData")

### Computing results for Goal 1 ###

bias1 <- apply(theta_EST, 2, mean) - theta0
bias2 <- abs(apply(theta_EST, 2, mean) - theta0) / theta0
MCSD <- apply(theta_EST, 2, sd)
ASE <- apply(theta_ASE, 2, mean)
WCP <- apply(theta_WCP, 2, mean)
round(cbind(bias1, MCSD, ASE, WCP) * 100, 1)
### view bias
round(cbind(bias2, MCSD, ASE, WCP) * 100, 1)
### view absolute relative bias

### Computing results for Goal 2 ###

for(ak in 1 : K){
	for(akp in seq(K)){
		AlphaAlphaMCcov[ak, akp] <- cov(AlphaBetaMC_MAT[, ak], AlphaBetaMC_MAT[, akp])
	}
}

for(bk in 1 : (K ^ 2)){
	for(bkp in 1 : (K ^ 2)){
		BetaBetaMCcov[bk, bkp] <- cov(AlphaBetaMC_MAT[, (K + bk)], AlphaBetaMC_MAT[, (K + bkp)])
	}
}

for(ak in 1 : K){
	for(bk in 1 : (K ^ 2)){
		AlphaBetaMCcov[ak, bk] <- cov(AlphaBetaMC_MAT[, ak], AlphaBetaMC_MAT[, (K + bk)])
	}
}

AlphaAlphaTcov <- CALCULATE_COV_ALPHA_ALPHA(A0, B0, p_vec, n)
BetaBetaTcov <- CALCULATE_COV_BETA_BETA(A0, B0, p_vec, n)
AlphaBetaTcov <- CALCULATE_COV_ALPHA_BETA(A0, B0, p_vec, n)
AlphaAlphaAcov <- AlphaAlphaAcov / repsMax
BetaBetaAcov <- BetaBetaAcov / repsMax
AlphaBetaAcov <- AlphaBetaAcov / repsMax

### comparison Alpha by True, MC, Averaged ###

round(cbind(AlphaAlphaTcov[, 1], AlphaAlphaMCcov[, 1], AlphaAlphaAcov[, 1]), 4)
round(cbind(AlphaAlphaTcov[, 2], AlphaAlphaMCcov[, 2], AlphaAlphaAcov[, 2]), 4)
round(cbind(AlphaAlphaTcov[, 3], AlphaAlphaMCcov[, 3], AlphaAlphaAcov[, 3]), 4)
round(cbind(AlphaAlphaTcov[, 4], AlphaAlphaMCcov[, 4], AlphaAlphaAcov[, 4]), 4)
round(cbind(AlphaAlphaTcov[, 5], AlphaAlphaMCcov[, 5], AlphaAlphaAcov[, 5]), 4)


### comparison Beta by True, MC, Averaged ###
round(cbind(BetaBetaTcov[indexB, indexB[1]], BetaBetaMCcov[indexB, indexB[1]], BetaBetaAcov[indexB, indexB[1]]), 4)
round(cbind(BetaBetaTcov[indexB, indexB[2]], BetaBetaMCcov[indexB, indexB[2]], BetaBetaAcov[indexB, indexB[2]]), 4)
round(cbind(BetaBetaTcov[indexB, indexB[3]], BetaBetaMCcov[indexB, indexB[3]], BetaBetaAcov[indexB, indexB[3]]), 4)
round(cbind(BetaBetaTcov[indexB, indexB[4]], BetaBetaMCcov[indexB, indexB[4]], BetaBetaAcov[indexB, indexB[4]]), 4)
round(cbind(BetaBetaTcov[indexB, indexB[5]], BetaBetaMCcov[indexB, indexB[5]], BetaBetaAcov[indexB, indexB[5]]), 4)

round(cbind(BetaBetaTcov[indexB, indexB[6]], BetaBetaMCcov[indexB, indexB[6]], BetaBetaAcov[indexB, indexB[6]]), 4)
round(cbind(BetaBetaTcov[indexB, indexB[7]], BetaBetaMCcov[indexB, indexB[7]], BetaBetaAcov[indexB, indexB[7]]), 4)
round(cbind(BetaBetaTcov[indexB, indexB[8]], BetaBetaMCcov[indexB, indexB[8]], BetaBetaAcov[indexB, indexB[8]]), 4)
round(cbind(BetaBetaTcov[indexB, indexB[9]], BetaBetaMCcov[indexB, indexB[9]], BetaBetaAcov[indexB, indexB[9]]), 4)
round(cbind(BetaBetaTcov[indexB, indexB[10]], BetaBetaMCcov[indexB, indexB[10]], BetaBetaAcov[indexB, indexB[10]]), 4)

round(cbind(BetaBetaTcov[indexB, indexB[11]], BetaBetaMCcov[indexB, indexB[11]], BetaBetaAcov[indexB, indexB[11]]), 4)
round(cbind(BetaBetaTcov[indexB, indexB[12]], BetaBetaMCcov[indexB, indexB[12]], BetaBetaAcov[indexB, indexB[12]]), 4)
round(cbind(BetaBetaTcov[indexB, indexB[13]], BetaBetaMCcov[indexB, indexB[13]], BetaBetaAcov[indexB, indexB[13]]), 4)
round(cbind(BetaBetaTcov[indexB, indexB[14]], BetaBetaMCcov[indexB, indexB[14]], BetaBetaAcov[indexB, indexB[14]]), 4)
round(cbind(BetaBetaTcov[indexB, indexB[15]], BetaBetaMCcov[indexB, indexB[15]], BetaBetaAcov[indexB, indexB[15]]), 4)


### comparison Alpha-Beta by True, MC, Averaged ###
round(cbind(AlphaBetaTcov[1, indexB], AlphaBetaMCcov[1, indexB], AlphaBetaAcov[1, indexB]), 4)
round(cbind(AlphaBetaTcov[2, indexB], AlphaBetaMCcov[2, indexB], AlphaBetaAcov[2, indexB]), 4)
round(cbind(AlphaBetaTcov[3, indexB], AlphaBetaMCcov[3, indexB], AlphaBetaAcov[3, indexB]), 4)
round(cbind(AlphaBetaTcov[4, indexB], AlphaBetaMCcov[4, indexB], AlphaBetaAcov[4, indexB]), 4)
round(cbind(AlphaBetaTcov[5, indexB], AlphaBetaMCcov[5, indexB], AlphaBetaAcov[5, indexB]), 4)


### Computing results for Goal 3 ###
rbind(
c(mean(NICE_COV[, 1]), sd(NICE_COV[, 1]), 
  mean(NICE_COV[, 2]), sd(NICE_COV[, 2]), NICE_COV_TIME),
c(mean(SOFT_COV[, 1]), sd(SOFT_COV[, 1]), 
  mean(SOFT_COV[, 2]), sd(SOFT_COV[, 2]), SOFT_COV_TIME),  
c(mean(HARD_COV[, 1]), sd(HARD_COV[, 1]), 
  mean(HARD_COV[, 2]), sd(HARD_COV[, 2]), HARD_COV_TIME),  
c(mean(ADAP_COV[, 1]), sd(ADAP_COV[, 1]), 
  mean(ADAP_COV[, 2]), sd(ADAP_COV[, 2]), ADAP_COV_TIME),
c(mean(POET_COV[, 1]), sd(POET_COV[, 1]), 
  mean(POET_COV[, 2]), sd(POET_COV[, 2]), POET_COV_TIME))

rbind(
c(mean(NICE_PRE[, 1]), sd(NICE_PRE[, 1]), 
  mean(NICE_PRE[, 2]), sd(NICE_PRE[, 2]), NICE_PRE_TIME),
c(mean(GLAS_PRE[, 1]), sd(GLAS_PRE[, 1]), 
  mean(GLAS_PRE[, 2]), sd(GLAS_PRE[, 2]), GLAS_PRE_TIME),  
c(mean(BAND_PRE[, 1]), sd(BAND_PRE[, 1]), 
  mean(BAND_PRE[, 2]), sd(BAND_PRE[, 2]), BAND_PRE_TIME),  
c(mean(BAY4_PRE[, 1]), sd(BAY4_PRE[, 1]), 
  mean(BAY4_PRE[, 2]), sd(BAY4_PRE[, 2]), BAY4_PRE_TIME),
c(mean(BAY7_PRE[, 1]), sd(BAY7_PRE[, 1]), 
  mean(BAY7_PRE[, 2]), sd(BAY7_PRE[, 2]), BAY7_PRE_TIME))

### checking the unit for computational time: mins or secs



