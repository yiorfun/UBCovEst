#####################################################
### Scenario 2: Comparison in the Large K setting ###
#####################################################

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
	"CovTools",
	### CovTools, apply various methods for covariance/precision matrix estimation
	"CVTuningCov"
	### CVTuningCov, modify it to perform CV for our method
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

MODIFIED_HARD_THRESHOLDING_RANDOM_CV_REVISED <- function(X, thr, test.size, norm, boot.num, BLOCK_HADAMARD_PRODUCT, BEST_UNBIASED_ESTIMATOR, p_vec){
    n <- nrow(X)
    p <- ncol(X)
	k.grid <- thr
    length.k.grid <- length(k.grid)    
    CV.pre.error <- matrix(0, length.k.grid, 2)
    test.size <- round(test.size)
    for (j in 1 : length.k.grid) {
        k <- k.grid[j]
        CV <- rep(0, 2)
		
        for (cv.v in 1 : boot.num) {
            cv.test <- sample(1 : n, size = test.size, replace = FALSE)
            cv.train <- cv.test
			
			S_RES_test_delt <- BEST_UNBIASED_ESTIMATOR(cov(X[-(cv.test), ]), p_vec)
			S_RES_test_only <- BEST_UNBIASED_ESTIMATOR(cov(X[cv.test, ]), p_vec)
			S_RES_trai_delt <- BEST_UNBIASED_ESTIMATOR(cov(X[-(cv.train), ]), p_vec)
			S_RES_trai_only <- BEST_UNBIASED_ESTIMATOR(cov(X[cv.train, ]), p_vec)
			
			
            if (norm == "L2") {
				
				CV[1] <- CV[1] + L2.norm2(BLOCK_HADAMARD_PRODUCT(hard.thresholding(S_RES_test_delt$A, k), hard.thresholding(S_RES_test_delt$B, k), p_vec) - cov(X[cv.test, ]))
			
            }
            else {
                
				CV[1] <- CV[1] + F.norm2(BLOCK_HADAMARD_PRODUCT(hard.thresholding(S_RES_test_delt$A, k), hard.thresholding(S_RES_test_delt$B, k), p_vec) - cov(X[cv.test, ]))
				
            }
        }
        CV.pre.error[j, ] <- CV / boot.num
    }
    CV.k <- rep(0, 2)
    CV.k[1] <- k.grid[which.min(CV.pre.error[, 1])]
    return(list(CV.k = CV.k, k.grid = k.grid, CV.pre.error = CV.pre.error))
}

############################
### Setup for Scenario 2 ###
############################

n <- 30
K <- 30 ### 30, 40, 50
set.seed(14)
# K    30, 40, 50
# Seed 14,140,290  

p_ind <- 10
p_vec <- rep(p_ind, K)
p <- sum(p_vec)  ### 300, 400, 500


thetaA <- runif(K)
thetaB_VEC <- runif(K * (K + 1) / 2, min = - 0.01, max = 0.01)
A0 <- diag(thetaA)
B0 <- matrix(0, K, K)
B0[upper.tri(B0, diag = TRUE)] <- thetaB_VEC
B0 <- (B0 + t(B0)) / 2
thetaB <- as.vector(B0)
theta0 <- c(thetaA, thetaB_VEC)
mu0 <- rep(0, p)
SIGMA0 <- COMB(A0, B0, p_vec)
DELTA0 <- A0 + B0 %*% diag(p_vec)
OMEGA0 <- COMB(solve(A0), - solve(DELTA0) %*% B0 %*% solve(A0), p_vec)

reps <- 1
repsMAX <- 1000

### Setting variables ###

NICE_COV <- SOFT_COV <- HARD_COV <- ADAP_COV <- POET_COV <- matrix(0, repsMAX, 2)

###########################
### Computing procedure ###
###########################

while(reps <= repsMAX){
	tryCatch({
	
	DATA_MAT <- mvrnorm(n = n, mu = mu0, Sigma = SIGMA0)
	
	opt_thr <- MODIFIED_HARD_THRESHOLDING_RANDOM_CV_REVISED(X = DATA_MAT, thr = seq(from = 1e-3, to = 1, length.out = 20), test.size = floor(n / log(n)), boot.num = 30, norm = "F", BLOCK_HADAMARD_PRODUCT = BLOCK_HADAMARD_PRODUCT, BEST_UNBIASED_ESTIMATOR = BEST_UNBIASED_ESTIMATOR, p_vec = p_vec)
	
	S_TEMP <- cov(DATA_MAT) 
	RES_NICE <- BEST_UNBIASED_ESTIMATOR(S_TEMP, p_vec)
	COV_A_NICE <- hard.thresholding(RES_NICE$A, opt_thr$CV.k[1])
	COV_B_NICE <- hard.thresholding(RES_NICE$B, opt_thr$CV.k[1])
	SIGMA_NICE <- BLOCK_HADAMARD_PRODUCT(COV_A_NICE, COV_B_NICE, p_vec)
	NICE_COV[reps, 1] <- norm(SIGMA_NICE - SIGMA0, "F")
	NICE_COV[reps, 2] <- norm(SIGMA_NICE - SIGMA0, "2")
	
	RES_SOFT <- CovEst.soft(X = DATA_MAT, thr = exp(seq(from = log(1e-3), to = log(1), length.out = 20)), nCV = 30)
	SOFT_COV[reps, 1] <- norm(RES_SOFT$S - SIGMA0, "F")
	SOFT_COV[reps, 2] <- norm(RES_SOFT$S - SIGMA0, "2")
	
	RES_HARD <- CovEst.hard(X = DATA_MAT, thr = exp(seq(from = log(1e-3), to = log(1), length.out = 20)), nCV = 30)
	HARD_COV[reps, 1] <- norm(RES_HARD$S - SIGMA0, "F")
	HARD_COV[reps, 2] <- norm(RES_HARD$S - SIGMA0, "2")
	
	RES_ADAP <- CovEst.adaptive(X = DATA_MAT, thr = exp(seq(from = log(1e-3), to = log(1), length.out = 20)), nCV = 30)
	ADAP_COV[reps, 1] <- norm(RES_ADAP$S - SIGMA0, "F")
	ADAP_COV[reps, 2] <- norm(RES_ADAP$S - SIGMA0, "2")
	
	K_HAT <- POETKhat(t(DATA_MAT))$K1HL ### p by n, each row has zero mean
	RES_POET <- POET(Y = t(DATA_MAT), K = K_HAT, C = 0.5, thres = "soft", matrix = "vad")
	POET_COV[reps, 1] <- norm(RES_POET$SigmaY - SIGMA0, "F")
	POET_COV[reps, 2] <- norm(RES_POET$SigmaY - SIGMA0, "2")
	
	cat(" iteration:  ", reps, "\r")
	reps <- reps + 1
	}, error = function(e){})
}

save.image("Scenario_2_SampleSize_30_NumberBlocks_304050.RData")

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


