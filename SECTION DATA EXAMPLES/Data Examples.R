#####################################################################
### Real Data Analyses for Proteomics Data and Brain Imaging Data ###
#####################################################################

############################
### Loading the packages ###
############################

REQUIRED_PACKAGES <- c(
	"R.matlab", 
	### R.matlab, convert between R files and MATLAB files
	"matrixcalc",		
	### matrixcalc::vech(), create a vector from a symmetric matrix
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

##################################
### Loading required functions ###
##################################

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

#######################################
### Loading the dataset of interest ###
#######################################

DATASET <- "PD"
#DATASET <- "EPSI"

if(DATASET == "PD"){

	EXAMPLES <- readMat("Proteomics_FOUR_MATS.mat")
	WorkingData <- EXAMPLES$PD.107
	FullData <- EXAMPLES$PD.ALL
	RawR <- EXAMPLES$PD.RAW
	n <- 288 
	p_vec <- EXAMPLES$p.PD[,1]
	indexB <- c(1 : 7, 9 : 14, 17 : 21, 25 : 28, 33 : 35, 41 : 42, 49)

} else if (DATASET == "EPSI"){

	EXAMPLES <- readMat("EPSI_FOUR_MATS.mat")
	WorkingData <- EXAMPLES$EPSI.227
	FullData <- EXAMPLES$EPSI.ALL
	RawR <- EXAMPLES$EPSI.RAW
	n <- 78 
	p_vec <- EXAMPLES$p.EPSI[,1]
	indexB <- c(1 : 5, 7 : 10, 13 : 15, 19 : 20, 25)

}

#############################################
### Fit the proposed model to the dataset ###
#############################################

K <- length(p_vec)
p <- sum(p_vec)
ALPHA <- 0.05

EST_RES <- BEST_UNBIASED_ESTIMATOR(WorkingData, p_vec)
A_EST <- EST_RES$A
B_EST <- EST_RES$B
thetaEST <- c(diag(A_EST), as.vector(B_EST)[indexB])
RES_variance <- CALCULATE_VAR_A_B(A_EST, B_EST, p_vec, n)

VAR_A_MAT <- RES_variance$VAR_A_MAT
VAR_B_MAT <- RES_variance$VAR_B_MAT
theta_VAR_temp <- c(diag(VAR_A_MAT), vech(VAR_B_MAT)[, 1])
theta_ASE <- sqrt(theta_VAR_temp)
thetaLB <- thetaEST - qnorm(1 - ALPHA / 2) * theta_ASE
thetaUB <- thetaEST + qnorm(1 - ALPHA / 2) * theta_ASE
round(cbind(thetaEST, theta_ASE, thetaLB, thetaUB), 4)

FullResidual <- FullData
FullResidual[1 : p, 1 : p] <- 0
ResResidual_SOFT <- soft.thresholding(Sigma = FullResidual, c = 0.5)
Sigma_prop <- ResResidual_SOFT
Sigma_prop[1 : p, 1 : p] <- BLOCK_HADAMARD_PRODUCT(A_EST, B_EST, p_vec)

##################################################
### Fit the convenctional model to the dataset ###
##################################################

Sigma_soft <- soft.thresholding(Sigma = FullData, c = 0.5)

############################################
### Convert results to MATLAB data files ###
############################################

if(DATASET == "PD"){
	PD_Sigma_samp <- RawR
	PD_Sigma_full <- FullData
	PD_Sigma_prop <- Sigma_prop
	PD_Sigma_soft <- Sigma_soft
} else if(DATASET == "EPSI"){
	EPSI_Sigma_samp <- RawR
	EPSI_Sigma_full <- FullData
	EPSI_Sigma_prop <- Sigma_prop
	EPSI_Sigma_soft <- Sigma_soft
}

writeMat("PD_EPSI_RESULT.mat", PD_Sigma_prop = PD_Sigma_prop, PD_Sigma_soft = PD_Sigma_soft, PD_Sigma_samp = PD_Sigma_samp, PD_Sigma_full = PD_Sigma_full, EPSI_Sigma_prop = EPSI_Sigma_prop, EPSI_Sigma_soft = EPSI_Sigma_soft, EPSI_Sigma_samp = EPSI_Sigma_samp, EPSI_Sigma_full = EPSI_Sigma_full)

figure;imagesc(PD_Sigma_samp);colormap jet;colorbar;snapnow;caxis([-1,1]);
figure;imagesc(PD_Sigma_prop);colormap jet;colorbar;snapnow;caxis([-1,1]);
figure;imagesc(PD_Sigma_soft);colormap jet;colorbar;snapnow;caxis([-1,1]);
figure;imagesc(PD_Sigma_full);colormap jet;colorbar;snapnow;caxis([-1,1]);

figure;imagesc(EPSI_Sigma_samp);colormap jet;colorbar;snapnow;caxis([-1,1]);
figure;imagesc(EPSI_Sigma_prop);colormap jet;colorbar;snapnow;caxis([-1,1]);
figure;imagesc(EPSI_Sigma_soft);colormap jet;colorbar;snapnow;caxis([-1,1]);
figure;imagesc(EPSI_Sigma_full);colormap jet;colorbar;snapnow;caxis([-1,1]);





