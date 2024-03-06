#############################################################
### Scenario 3: Comparison under misspecification setting ###
#############################################################

##################################################################
# Goal: evaluation of covariance- and precision-matrix estimates #
#         in Frobenius and spectral norms					     #
##################################################################

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

###########################################
### Loading required functions for Goal ###
###########################################

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

############################
### Setup for Scenario 3 ###
############################

ERROR_SIGMA <- 0.1 ### 0.1, 0.5, 0.8
n <- 50
p_ind <- 30
p_vec <- rep(p_ind, 5)
K <- length(p_vec)
p <- sum(p_vec) ### 150
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

set.seed(2021)
reps <- 1
repsMAX <- 1000

### Setting variables ###

NICE_COV <- SOFT_COV <- HARD_COV <- ADAP_COV <- POET_COV <- NICE_PRE <- GLAS_PRE <- BAND_PRE <- BAY4_PRE <- BAY7_PRE <- matrix(0, repsMAX, 2) 
### [, 1] = Frobenius, [, 2] = spectral

###########################
### Computing procedure ###
###########################

while(reps <= repsMAX){
	tryCatch({
	
	ERROR_MAT <- rWishart(1, df = p, Sigma = ERROR_SIGMA * diag(p))[,,1]
	UPSILON0 <- SIGMA0 + ERROR_MAT
	while(!all(eigen(UPSILON0)$values > 0)){
		ERROR_MAT <- rWishart(1, df = p, Sigma = ERROR_SIGMA * diag(p))[,,1]
		UPSILON0 <- SIGMA0 + ERROR_MAT
	}
	
	DATA_MAT <- mvrnorm(n = n, mu = mu0, Sigma = UPSILON0)
	inv_UPSILON0 <- solve(UPSILON0)
	
	S_TEMP <- cov(DATA_MAT)
	RES_NICE <- BEST_UNBIASED_ESTIMATOR(S_TEMP, p_vec)
	COV_A_NICE <- RES_NICE$A
	COV_B_NICE <- RES_NICE$B
	PRE_A_NICE <- solve(COV_A_NICE)
	PRE_B_NICE <- - solve(COV_A_NICE + COV_B_NICE %*% diag(p_vec)) %*% COV_B_NICE %*% solve(COV_A_NICE)
	SIGMA_NICE <- BLOCK_HADAMARD_PRODUCT(COV_A_NICE, COV_B_NICE, p_vec)
	OMEGA_NICE <- BLOCK_HADAMARD_PRODUCT(PRE_A_NICE, PRE_B_NICE, p_vec)
	NICE_COV[reps, 1] <- norm(SIGMA_NICE - UPSILON0, "F")
	NICE_COV[reps, 2] <- norm(SIGMA_NICE - UPSILON0, "2")
	NICE_PRE[reps, 1] <- norm(OMEGA_NICE - inv_UPSILON0, "F")
	NICE_PRE[reps, 2] <- norm(OMEGA_NICE - inv_UPSILON0, "2")
	
	RES_SOFT <- CovEst.soft(X = DATA_MAT, thr = exp(seq(from = log(0.1), to = log(10), length.out = 10)))
	SOFT_COV[reps, 1] <- norm(RES_SOFT$S - UPSILON0, "F")
	SOFT_COV[reps, 2] <- norm(RES_SOFT$S - UPSILON0, "2")
	
	RES_HARD <- CovEst.hard(X = DATA_MAT, thr = exp(seq(from = log(0.1), to = log(10), length.out = 10)))
	HARD_COV[reps, 1] <- norm(RES_HARD$S - UPSILON0, "F")
	HARD_COV[reps, 2] <- norm(RES_HARD$S - UPSILON0, "2")
	
	RES_ADAP <- CovEst.adaptive(X = DATA_MAT, thr = exp(seq(from = log(0.1), to = log(10), length.out = 10)))
	ADAP_COV[reps, 1] <- norm(RES_ADAP$S - UPSILON0, "F")
	ADAP_COV[reps, 2] <- norm(RES_ADAP$S - UPSILON0, "2")
	
	K_HAT <- POETKhat(t(DATA_MAT))$K1HL ### p by n, each row has zero mean
	RES_POET <- POET(Y = t(DATA_MAT), K = K_HAT, C = 0.5, thres = "soft", matrix = "vad")
	POET_COV[reps, 1] <- norm(RES_POET$SigmaY - UPSILON0, "F")
	POET_COV[reps, 2] <- norm(RES_POET$SigmaY - UPSILON0, "2")
	
	RES_GLAS <- PreEst.glasso(X = DATA_MAT, method = list(type = "BIC", param = exp(seq(from = log(0.1), to = log(10), length.out = 10))))
	GLAS_PRE[reps, 1] <- norm(RES_GLAS$C - inv_UPSILON0, "F")
	GLAS_PRE[reps, 2] <- norm(RES_GLAS$C - inv_UPSILON0, "2")
	
	RES_BAND <- PreEst.2014An(X = DATA_MAT, upperK = 2, algorithm = "Bonferroni", alpha = 0.05)
	BAND_PRE[reps, 1] <- norm(RES_BAND$C - inv_UPSILON0, "F")
	BAND_PRE[reps, 2] <- norm(RES_BAND$C - inv_UPSILON0, "2")
	
	RES_BAY4 <- PreEst.2014Banerjee(X = DATA_MAT, upperK = 2, delta = 10, logpi = function(k) { -k^4 }, loss = "Stein")
	BAY4_PRE[reps, 1] <- norm(RES_BAY4$C - inv_UPSILON0, "F")
	BAY4_PRE[reps, 2] <- norm(RES_BAY4$C - inv_UPSILON0, "2")
	
	RES_BAY7 <- PreEst.2017Lee(X = DATA_MAT, upperK = 2, logpi = function(k) {-k^4})
	BAY7_PRE[reps, 1] <- norm(RES_BAY7$C - inv_UPSILON0, "F")
	BAY7_PRE[reps, 2] <- norm(RES_BAY7$C - inv_UPSILON0, "2")
	
	cat(" iteration:  ", reps, "\r")
	reps <- reps + 1
	}, error = function(e){})
}

save.image("Scenario_3_SampleSize_50_error_sigma_010508.RData")

### Computing results ###

rbind(
c(mean(NICE_COV[, 1]), sd(NICE_COV[, 1]), 
  mean(NICE_COV[, 2]), sd(NICE_COV[, 2])),
c(mean(SOFT_COV[, 1]), sd(SOFT_COV[, 1]), 
  mean(SOFT_COV[, 2]), sd(SOFT_COV[, 2])),  
c(mean(HARD_COV[, 1]), sd(HARD_COV[, 1]), 
  mean(HARD_COV[, 2]), sd(HARD_COV[, 2])),  
c(mean(ADAP_COV[, 1]), sd(ADAP_COV[, 1]), 
  mean(ADAP_COV[, 2]), sd(ADAP_COV[, 2])),
c(mean(POET_COV[, 1]), sd(POET_COV[, 1]), 
  mean(POET_COV[, 2]), sd(POET_COV[, 2])))

rbind(
c(mean(NICE_PRE[, 1]), sd(NICE_PRE[, 1]), 
  mean(NICE_PRE[, 2]), sd(NICE_PRE[, 2])),
c(mean(GLAS_PRE[, 1]), sd(GLAS_PRE[, 1]), 
  mean(GLAS_PRE[, 2]), sd(GLAS_PRE[, 2])),  
c(mean(BAND_PRE[, 1]), sd(BAND_PRE[, 1]), 
  mean(BAND_PRE[, 2]), sd(BAND_PRE[, 2])),  
c(mean(BAY4_PRE[, 1]), sd(BAY4_PRE[, 1]), 
  mean(BAY4_PRE[, 2]), sd(BAY4_PRE[, 2])),
c(mean(BAY7_PRE[, 1]), sd(BAY7_PRE[, 1]), 
  mean(BAY7_PRE[, 2]), sd(BAY7_PRE[, 2])))


##################################################
### Additional Goal: create boxplots by MATLAB ###
##################################################

require(R.matlab)
require(MASS)

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
	K <- length(p_vec)
	A_temp <- matrix(0, K, K)
	B_temp <- matrix(0, K, K)
	for(k in 1 : K){
		for(kp in 1 : K){
			SUB <- S[(sum(p_vec[0 : (k - 1)]) + 1) : sum(p_vec[0 : k]), 
						  (sum(p_vec[0 : (kp - 1)]) + 1) : sum(p_vec[0 : kp])]
			if(kp == k){
				SUB_ON <- mean(diag(SUB))
				SUB_OF <- (sum(SUB[upper.tri(SUB)]) + sum(SUB[lower.tri(SUB)])) / (p_vec[k] ^ 2 - p_vec[k])
				A_temp[k, kp] <- SUB_ON - SUB_OF
				B_temp[k, kp] <- SUB_OF
			} else{
				B_temp[k, kp] <- mean(as.vector(SUB))
			}	
		}
	}
	A <- A_temp
	B <- (B_temp + t(B_temp)) / 2
	return(list(A = A, B = B))
}

ERROR_SIGMA <- 0.1 ### 0.1, 0.5, 0.8
n <- 50
p_ind <- 30
p_vec <- rep(p_ind, 5)
K <- length(p_vec)
p <- sum(p_vec) ### 150
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

ERROR_MAT <- rWishart(1, df = p, Sigma = ERROR_SIGMA * diag(p))[,,1]
UPSILON0 <- SIGMA0 + ERROR_MAT
while(!all(eigen(UPSILON0)$values > 0)){
	ERROR_MAT <- rWishart(1, df = p, Sigma = ERROR_SIGMA * diag(p))[,,1]
	UPSILON0 <- SIGMA0 + ERROR_MAT
}
DATA_MAT <- mvrnorm(n = n, mu = mu0, Sigma = UPSILON0)
inv_UPSILON0 <- solve(UPSILON0)
set.seed(2021)	
S_TEMP <- cov(DATA_MAT)
S_CORR <- cov2cor(S_TEMP)
writeMat("Supp_C3_Boxplots_error_sigma_010508.MAT", S_01 = S_CORR, p_vec = p_vec)

### see SECTION SUPPLEMENTARY ###

