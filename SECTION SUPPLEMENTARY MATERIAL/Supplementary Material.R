#######################################################################
### Supplementary Material Section C: Additional Simulation Studies ###
#######################################################################

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
	"R.matlab", 
	### R.matlab, convert between R files and MATLAB files
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

##################################
### Loading required functions ###
##################################

require(quantreg)
require(lars)
require(dplyr, warn.conflicts = FALSE)

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

FINITE_PERFORMANCE <- function(repsMax, n, theta0, mu0, SIGMA0, p_vec, seed, ALPHA, indexB, CALCULATE_VAR_A_B, BLOCK_HADAMARD_PRODUCT, BEST_UNBIASED_ESTIMATOR){
	
	set.seed(seed)

	theta_EST <- matrix(0, repsMax, length(theta0))
	A_EST <- matrix(0, length(p_vec), length(p_vec))
	B_EST <- matrix(0, length(p_vec), length(p_vec))
	theta_ASE <- matrix(0, repsMax, length(theta0))
	theta_WCP <- matrix(0, repsMax, length(theta0))
	reps <- 1
	
	while(reps <= repsMax){

		DATA <- mvrnorm(n = n, mu = mu0, Sigma = SIGMA0)
		S <- cov(DATA)
	
		RES_theta_tilde <- BEST_UNBIASED_ESTIMATOR(S, p_vec)
		A_tilde <- RES_theta_tilde$A
		B_tilde <- RES_theta_tilde$B
		theta_tilde_temp <- c(diag(A_tilde), (as.vector(B_tilde))[indexB])
		theta_EST[reps, ] <- theta_tilde_temp
		A_EST <- A_EST + A_tilde
		B_EST <- B_EST + B_tilde
		
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
		reps <- reps + 1
	}
	A_EST <- A_EST / repsMax
	B_EST <- B_EST / repsMax
	return(list(theta_EST = theta_EST, theta_ASE = theta_ASE, theta_WCP = theta_WCP, A_EST = A_EST, B_EST = B_EST, S = S))
}

pfa <- function(Z, Sigma, t, Kmax, reg="L1", e=.05, gamma, K, plot="-log") {
  
  Z <- as.vector(Z)
  Sigma <- as.matrix(Sigma)
  p <- length(Z)
  
  # standardize the data

  SD <- sqrt(diag(Sigma))
  Z <- Z/SD
  Sigma <- diag(1/SD) %*% Sigma %*% diag(1/SD)

  # pca 

  pca <- svd(Sigma, nu=0, nv=Kmax)
  lambda <- pca$d
  eigvec <- pca$v
    

  # determine the factor loadings
  
  if (missing(K)) {
      K = 1
      while( K < Kmax & sqrt(sum(lambda[(K+1):length(lambda)]^2)) >= e*sum(lambda)) 
         K = K + 1 
  }    
  sqrt_lambda <- as.matrix(sqrt(lambda[1:K]))
  b <- as.matrix(eigvec[,1:K])
  for (i in 1:K)  {
    b[,i] <- b[,i] * sqrt_lambda[i]  # factor loadings
  }


  # estimate the factors, with 5% largest |Z| eliminated.
   
  
  if (reg == "L1") {
    W.hat <- rq(Z ~ b -1, 0.5)$coef     # L_1 regression (no intercept)
    # W.hat <- W.hat[2:(K+1)] 
  } 
  else if (reg == "L2")  # L_2 regression (no intercept) 
  {    
    #temp<-sort(abs(Z),decreasing=TRUE,index.return=TRUE)
    #len=round(length(Z)*0.05)
    #Ztemp<-temp$x
    #btemp<-as.matrix(b[temp$ix,])
    #Ztemp<-Ztemp[(len+1):length(Z)]
    #btemp<-btemp[(len+1):length(Z),]
    #W.hat<-lm(Ztemp ~ btemp - 1)$coef
    o = order(abs(Z))
    Zperm = Z[o]
    Lperm = as.matrix(b[o,])
    Z.reduce = Zperm[1:(p*0.95)]
    L.reduce = as.matrix(Lperm[1:(p*0.95),]) 
    W.hat = lsfit(x=L.reduce,y=Z.reduce,intercept=F)$coef
  }
  #if (reg == "huber") {
    #W.hat <- rlm(Z ~ b, 0.5)$coef	  # robust/huber regression
  #}

  rs <- rowSums(b^2)
  inv_a <- sqrt( ((1 - rs)+abs(1-rs))/2 )
  bW.est <- b%*%(W.hat)


  # compute p-values & estimate p0
  
  P <- 2*(1-pnorm(abs(Z)))
  sort <- sort(P, index.return=TRUE)
  index <- sort$ix
  P <- sort$x
  if (missing(gamma))  
      gamma <- as.numeric(quantile(P, probs=0.4))
  p0.est <- min(p, sum(P>gamma)/(1-gamma))


  # estimate FDP
  
  t.default <- TRUE
  if (!missing(t)) {
       if (t[1] =="pval")  { ### origonal t == "pval" ###
          t <- P
          t.default <- FALSE
       } 
       if (is.numeric(t)) 
          t.default = ( sum(t>=0)+ sum(t<=1) < 2*length(t))
  }
  if (t.default) {
       logt.l <- max(min(log(P)), log(1e-14))
       logt.u <- max(log(P))
       grid <- (logt.u-logt.l)*seq(from=0.01,to=1,by=0.025)*0.5 + 0.85*logt.l+0.15*logt.u
       t <- exp(grid)
  }
       
  FDPt <- Vt <- Rt<- rep(0, length(t))
  for (l in 1:length(t)) {
  	   P1 <- 2*(1-pnorm(abs(Z)))
       Rt[l] <- sum(P1<=t[l])
       a <- rep(0,p)
       for (j in 1:p)  {
           qtl <- qnorm(t[l]/2)
           if (inv_a[j]>0)  {
                a[j] <- pnorm((qtl + bW.est[j])/inv_a[j])+ pnorm((qtl - bW.est[j])/inv_a[j])
           } else {
                a[j] <- as.numeric(abs(bW.est[j])>abs(qtl))
           }
       }
       Vt[l] <- min(sum(a), Rt[l]) 
       if (Rt[l]==0)   {
           FDPt[l] <- 0
       } else  {
           FDPt[l] <- Vt[l]/Rt[l]
       }
  } 


  # factor adjusted procedure

  adj.P <- as.vector(rep(0,p))
  for (j in 1:p) {
       if (inv_a[j]>0)  {
           adj.P[j] <- 2*( 1- pnorm(abs(Z[j] - bW.est[j])/inv_a[j]) )
       }  else  {
           adj.P[j] <- as.numeric(abs(Z[j]-bW.est[j])==0)
       }
  }
  sort <- sort(adj.P, index.return=TRUE)
  adj.index <- sort$ix
  adj.P <- sort$x

  
  # output
  
  Pvals <- data.frame(p.value = P,  Index =index)
  adjPvals <- data.frame(p.value = adj.P,  Index = adj.index)
  if (t.default)   {
       FDPvals <- data.frame( minus.logt= - log(t), rejects= Rt, false.rejects=Vt, FDP=FDPt )
  } else {
       FDPvals <- data.frame( t= t, rejects= Rt, false.rejects=Vt, FDP=FDPt )
  }
  results <- list("Pvalue"=Pvals, "adjPvalue"=adjPvals, "FDP"=FDPvals, "pi0"=p0.est/p, "K"=K, "sigma"=NULL)
  class(results) <- "FDPresults"
  if (plot=="-log")  {
            par(mfrow=c(2,2))
            hist(P, main = "Histogram of p-values", xlab="p-values")
            plot(-log(t), Rt, xlab="-log(t)", ylab="", main="Number of total rejections", type='o')
            plot(-log(t),Vt, xlab="-log(t)", ylab="", main="Number of estimated false rejections", type='o')
            plot(-log(t),FDPt,xlab="-log(t)", ylab="", main="Estimated FDP", type='o')
   }  else  if (plot=="linear") {
            par(mfrow=c(2,2))
            hist(P, main = "Histogram of p-values", xlab="p-values")
            plot(t, Rt, xlab="t", ylab="", main="Number of total rejections", type='o')
            plot(t,Vt, xlab="t", ylab="", main="Number of estimated false rejections", type='o')
            plot(t,FDPt,xlab="t", ylab="", main="Estimated FDP", type='o')
  }  else if (plot=="log") {
            par(mfrow=c(2,2))
            hist(P, main = "Histogram of p-values", xlab="p-values")
            plot(log(t), Rt, xlab="log(t)", ylab="", main="Number of total rejections", type='o')
            plot(log(t),Vt, xlab="log(t)", ylab="", main="Number of estimated false rejections", type='o')
            plot(log(t),FDPt,xlab="log(t)", ylab="", main="Estimated FDP", type='o')
  }
  return(results)
}

print.FDPresults <- function(x, ...)  {

  cat("P-values:\n")
  print(x$Pvalue)
  
  cat("Factor-ajusted P-values:\n")
  print(x$adjPvalue)

  cat("Estimated false discoveries:\n")
  print(x$FDP)
  
  cat("Estimated true null proportion:\n")
  print(x$pi0)

  cat("Number of factors:\n")
  print(x$K)

  if (! is.null(x$sigma))  {
      cat("Estimated noise standard deviation:\n")
      print(x$sigma)
  }
}

summary.FDPresults <- function(object, ...)  { 
  x <- object

  cat("Method:\n")
  cat(paste("PFA with", x$K, "factors.\n"))
  cat("\n")
  
  p <- length(x$Pvalue$p.value)
  cat("P-values:\n")
  print(x$Pvalue[1:min(20,p),])
  cat("...\n")
  
  cat("Factor-ajusted P-values:\n")
  print(x$adjPvalue[1:min(20,p),])
  cat("...\n")

  T <- length(x$FDP$rejects)
  cat("Estimated false discoveries:\n")
  print(x$FDP[1:min(10,T),])
  cat("...\n")
  
  cat("Estimated true null proportion:\n")
  cat(paste(round(x$pi0, digits=3), "\n"))

  if (! is.null(x$sigma))  {
      cat("Estimated noise standard deviation:\n")
      cat(paste(round(x$sigma,digits=3), "\n"))
  }
}  

pfa.test <- function(X, Y, tval, Sigma, reg="L2", K, e=0.05, gamma, mat_est="poet", plot="-log", pfa) {

	
	# compute the test statistics

	if (missing(Y)) {
		X <- as.matrix(X)
		if (ncol(X)==1)  {
		Z <- as.vector(X)
		}  else  {
		Z <- colMeans(X)*sqrt(nrow(X))
		}
	} else {
		X <- as.matrix(X)
		Y <- as.matrix(Y)
		if (ncol(X)!= ncol(Y))
		stop("Dimensions of X and Y do not match.")
		Z <- (colMeans(X) - colMeans(Y))*sqrt(nrow(X)*nrow(Y)/(nrow(X)+nrow(Y)))
	}
	p <- length(Z)
	Kmax <- p

	# determine number of factors
	if (missing(K)&(ncol(X)>=2)) {
		if(missing(Y)){
		kmax<-floor(0.2*nrow(X))
		if(kmax<=1){
			stop("No enough samples to determine the number of factors!")
		}else{
		sample.cor<-cor(t(X),t(X))
		sample.pca<-eigen(sample.cor)
		sample.eigval<-sample.pca$values
		ratio<-sample.eigval[1:(kmax-1)]/sample.eigval[2:kmax]
		K<-which(ratio==max(ratio))	
		}
		}else{
			n<-nrow(X)
			m<-nrow(Y)
			kmax<-floor(0.2*(n+m))
			if(kmax<=1){
				stop("No enough samples to determine the number of factors!")
			}else{
				Xmean<-apply(X,2,mean)
				Ymean<-apply(Y,2,mean)
				Znew<-as.matrix(rbind(X-Xmean,Y-Ymean))
				sample.cor<-cor(t(Znew),t(Znew))
			sample.pca<-eigen(sample.cor)
			sample.eigval<-sample.pca$values
			ratio<-sample.eigval[1:(kmax-1)]/sample.eigval[2:kmax]
			K<-which(ratio==max(ratio))
		}
		}
	}


	# compute the covariance matrix

	if (missing(Sigma)) {
		if (missing(Y) & nrow(X)==1) {
		stop("No enough samples to estimate the covariance matrix.")
		} 
		if (mat_est =="sample") {
		n <- nrow(X)
		barX <- matrix(colMeans(X), nrow=n, ncol=p, byrow=T)
		SigmaX <- (t(X - barX )%*% (X - barX))/n
		if (missing(Y)) {
			Sigma <- SigmaX
			Kmax <- n
		} else {
			n2 <- nrow(Y)
			barY <- matrix(colMeans(Y), nrow=n2, ncol=p, byrow=T)
			SigmaY <- (t(Y - barY)%*% (Y - barY))/n2
			Sigma <- (SigmaX*n + SigmaY*n2)/(n+n2)
			Kmax <- n + n2
		}  
		} else { 
			if (mat_est=="poet") {
				n <- nrow(X)
				barX <- matrix(colMeans(X), nrow=n, ncol=p, byrow=T)
				if (missing(Y)) {
					# R <- X - barX
					R<-X
				} else {
					n2 <- nrow(Y)
					barY <- matrix(colMeans(Y), nrow=n2, ncol=p, byrow=T)
					# R <- rbind(X-barX, Y-barY)
					R<-rbind(X,Y)
				}
				Sigma <- POET(t(R),K=K,C=0.5,thres='soft',matrix='vad')$SigmaY
				Kmax <- p
			}
		}
	}

	# run the test and compute FDP

	results <- pfa(Z=Z, Sigma=Sigma, t=tval, Kmax=Kmax, reg=reg, e=e, gamma=gamma, K=K, plot=plot)
	return(results)
}

COMPUTER_U_V_T_S <- function(TRUE_NULL_VEC, P_VALUES_MAT, ALPHA){
	
	### TRUE_NULL_VEC: 0 refers null is true  (e.g., mu_i = 0), 
	###			  	   1 refers null is false (e.g., mu_i = 1)
	### DECISION_MAT: 0 refers to remain null
	###				  1 refers to reject null
	### TRUE_NULL_VEC - DECISION_MAT: 0 - 0 remain true null  U
	###							 	  0 - 1 reject true null  V
	###							 	  1 - 0 remain false null T
	###							 	  1 - 1 reject false null S
	
	DECISION_MAT <- 1 * (P_VALUES_MAT <= ALPHA)
	METHOD_MAT <- matrix(0, 4, ncol(DECISION_MAT)) #### U, V, T, S as rows
	
	for(j in 1 : ncol(DECISION_MAT)){
		U <- V <- T <- S <- 0
		COL_VEC <- DECISION_MAT[, j]
		for(i in 1 : nrow(DECISION_MAT)){
			if(TRUE_NULL_VEC[i] == 0 & COL_VEC[i] == 0){
				U <- U + 1
			} else if (TRUE_NULL_VEC[i] == 0 & COL_VEC[i] == 1){
				V <- V + 1
			} else if (TRUE_NULL_VEC[i] == 1 & COL_VEC[i] == 0){
				T <- T + 1
			} else if (TRUE_NULL_VEC[i] == 1 & COL_VEC[i] == 1){
				S <- S + 1
			}
		}
		METHOD_MAT[, j] <- c(U, V, T, S)
	}
	colnames(METHOD_MAT) <- colnames(DECISION_MAT)
	rownames(METHOD_MAT) <- c("U", "V", "T", "S")
	return(METHOD_MAT)
}

COMPUTER_THRESHOLDING_T <- function(RESULT, ALPHA){
	
	RES_FDP <- RESULT$FDP
	RES_FDP <- RES_FDP[order(RES_FDP$FDP), ]
	if(min(RES_FDP$FDP) > ALPHA){
		stop("minimal FDP is greater than ALPHA")
	} else {
		t_ind_temp <- 1
		while(RES_FDP$FDP[t_ind_temp] <= ALPHA){
			t_ind_temp <- t_ind_temp + 1
		}
		### the first t_ind_temp s.t. RES_FDP$FDP[t_ind_temp] > ALPHA
		t_ind <- t_ind_temp - 1
		t_threshold <- RES_FDP$t[t_ind]
	}
	
	RES_P_ADJ <- RESULT$adjPvalue
	P_REJECTION_VEC <- rep(1, length(RES_P_ADJ$Index))
	P_REJECTION_VEC[which(RES_P_ADJ$p.value < t_threshold)] <- 1e-14
	
	return(list(P_REJ = P_REJECTION_VEC, t_th = t_threshold))
}

TABLE_SUMMARY <- function(PFA_OBJ){
    
    summarise_table <- PFA_OBJ %>% group_by(t)  %>% summarise(median(R), sd(R), median(S_adjP), sd(S_adjP), median(FDP), sd(FDP)) %>% as.data.frame()
    
    return(summarise_table)
}

EXTRACT_RAW_ADJ_PDP <- function(PFA_OBJ, CUT_OFFS, m1){

	count_p1 <- count_adj_p1 <- p_val_p1 <- adj_p1 <- adj_total <- adj_p <- c()
  
    for(i in 1 : length(CUT_OFFS)){
      p_val_p1 <-  dim(PFA_OBJ$Pvalue[PFA_OBJ$Pvalue$p.value <= CUT_OFFS[i] &
                                    PFA_OBJ$Pvalue$Index %in% 1:m1, ])[1]
      
      adj_p1   <-  dim(PFA_OBJ$adjPvalue[PFA_OBJ$adjPvalue$p.value <= CUT_OFFS[i] &
                                       PFA_OBJ$adjPvalue$Index %in% 1:m1,])[1]
      
      adj_p <-  dim(PFA_OBJ$adjPvalue[PFA_OBJ$adjPvalue$p.value <= CUT_OFFS[i],])[1]
      
      count_p1 <- rbind(count_p1, p_val_p1)
      
      count_adj_p1 <- rbind(count_adj_p1, adj_p1)
      
      adj_total <- rbind(adj_total, adj_p)
    }
	res <- data.frame('t'        = CUT_OFFS,
					  'R'        = PFA_OBJ$FDP$rejects, 
					  'V'        = PFA_OBJ$FDP$false.rejects,
					  'pi0'      = PFA_OBJ$pi0,
					  'FDP'      = PFA_OBJ$FDP$FDP, 
					  'K'        = PFA_OBJ$K, 
					  'S_unadjP' = count_p1,
					  'R_adjP'   = adj_total,
					  'S_adjP'   = count_adj_p1)
  
	return(res)
}



##################################################################
### Scenario C1: Examination of the accuracy of the covariance ###
### estimator 					            				   ###
##################################################################

### See Goal 2 of Scenario 1: Comparison in the Small K setting

##################################################################
### Scenario C2: Evaluation of the estimated covariance matrix ###
### estimator on multiple testing                              ###
##################################################################

n <- 50 			###  50, 100, 150
p_ind <- 100 		###  30,  45,  60,  100,  200
p_vec <- rep(p_ind, 5)
K <- length(p_vec)  
p <- sum(p_vec)		### 150, 225, 300,  500
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
m1 <- 10
STH <- 1e-1 * 3                     ### the strength 
m0 <- p - m1 		   		        ### true H0: mu_i = 0
mu <- as.vector(c(rep(STH, m1), rep(0, m0)))
SIGMA0 <- BLOCK_HADAMARD_PRODUCT(A0, B0, p_vec)

ALPHA <- 0.05
reps <- 1
repsMax <- 100
RES_ORACLE_SIM <- RES_POET_SIM <- RES_UB_SIM <- c() 
T_THRESHOLD_VEC <- seq(1e-4, 1e-1, length.out = 20)
set.seed(2022)

while(reps <= repsMax){
	tryCatch({
	
	DATA <- mvrnorm(n, mu, SIGMA0)
	X_BAR <- apply(DATA, 2, mean) ### follows N(mu, CORR0 / n)
	
	RES_PFA_ORACLE <- pfa.test(X_BAR, Sigma = SIGMA0 / n, tval = T_THRESHOLD_VEC, reg = "L1", plot = "none", pfa = pfa)
	RES_ORACLE_SIM <- rbind(RES_ORACLE_SIM, EXTRACT_RAW_ADJ_PDP(RES_PFA_ORACLE, T_THRESHOLD_VEC, m1))
	
	
	K_HAT <- POETKhat(t(DATA))$K1HL ### p by n, each row has zero mean
	RES_POET <- POET(Y = t(DATA), K = K_HAT, C = 0.5, thres = "soft", matrix = "vad")	
	RES_PFA_POET <- pfa.test(sqrt(n) * X_BAR, Sigma = RES_POET$SigmaY, tval = T_THRESHOLD_VEC, reg = "L1", plot = "none", pfa = pfa)
	RES_POET_SIM <- rbind(RES_POET_SIM, EXTRACT_RAW_ADJ_PDP(RES_PFA_POET, T_THRESHOLD_VEC, m1))
	
	S <- cov(DATA)
	RES_theta_tilde <- BEST_UNBIASED_ESTIMATOR(S, p_vec)
	A_tilde <- RES_theta_tilde$A
	B_tilde <- RES_theta_tilde$B
	SIGMA_tilde <- BLOCK_HADAMARD_PRODUCT(A_tilde, B_tilde, p_vec)
	RES_PFA_UB <- pfa.test(X_BAR, Sigma = SIGMA_tilde / n, tval = T_THRESHOLD_VEC, reg = "L1", plot = "none", pfa = pfa)
	RES_UB_SIM <- rbind(RES_UB_SIM, EXTRACT_RAW_ADJ_PDP(RES_PFA_UB, T_THRESHOLD_VEC, m1))
	
	cat(" iteration:  ", reps, "\r")
	reps <- reps + 1
	}, error = function(e){})
}

RES_ORACLE_DF <- round(TABLE_SUMMARY(RES_ORACLE_SIM), 4)
RES_POET_DF <- round(TABLE_SUMMARY(RES_POET_SIM), 4)
RES_UB_DF <- round(TABLE_SUMMARY(RES_UB_SIM), 4)


###################################################################
### Scenario C3: Justification of the proposed parameterization ###
### strategy                              					    ###
###################################################################

### See Goal 2 of Scenario 3: Comparison under mis-specification setting

#################################################################
### Scenario C4: Justification of the proposed approach among ###
### alternative covariance structures 					      ###
#################################################################

n <- 50
p_ind <- 10
K <- 5
p_vec <- c(p_ind, p_ind * 1.5, p_ind * 1.5, p_ind * 2, p_ind * 2.5)
p <- sum(p_vec)

thetaAU0 <- c(0.01595042, 0.21392707, 0.74912381, 0.06771268, 0.10017260)
thetaBU0 <- c(6.73139386,-1.69034339, 0.69591280,-2.93647430, 1.91315819, 
		     -1.69034339, 5.21462208, 3.81497235,-1.01011751, 0.70298054, 
              0.69591280, 3.81497235, 4.32780351,-3.35737580,-0.26890092,
		     -2.93647430,-1.01011751,-3.35737580, 6.78768893, 0.00018746,
              1.91315819, 0.70298054,-0.26890092, 0.00018746, 3.95418249)
thetaAH0 <- thetaAU0
thetaBH0 <- c(6.73139386,-1.69034339, 0.53007201, 0.53007201, 0.53007201, 
		     -1.69034339, 5.21462208, 0.53007201, 0.53007201, 0.53007201, 
              0.53007201, 0.53007201, 4.32780351, 0.79039990, 0.79039990,
		      0.53007201, 0.53007201, 0.79039990, 6.78768893, 0.79039990,
              0.53007201, 0.53007201, 0.79039990, 0.79039990, 3.95418249)
thetaAI0 <- thetaAU0
thetaBI0 <- c(6.73139386, 0, 0, 0, 0, 
		      0, 5.21462208, 0, 0, 0, 
              0, 0, 4.32780351, 0, 0,
		      0, 0, 0, 6.78768893, 0,
              0, 0, 0, 0, 3.95418249)
indexB <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 13, 14, 15, 19, 20, 25)
### a symmetric matrix by row along upper triangle
AU0 <- diag(thetaAU0)
BU0 <- matrix(thetaBU0, K, K)
thetaU0 <- c(thetaAU0, thetaBU0[indexB])
AH0 <- diag(thetaAH0)
BH0 <- matrix(thetaBH0, K, K)
thetaH0 <- c(thetaAH0, thetaBH0[indexB])
AI0 <- diag(thetaAI0)
BI0 <- matrix(thetaBI0, K, K)
thetaI0 <- c(thetaAI0, thetaBI0[indexB])

mu0 <- rep(0, p)
SIGMAU0 <- BLOCK_HADAMARD_PRODUCT(AU0, BU0, p_vec)
SIGMAH0 <- BLOCK_HADAMARD_PRODUCT(AH0, BH0, p_vec)
SIGMAI0 <- BLOCK_HADAMARD_PRODUCT(AI0, BI0, p_vec)
repsMax <- 1000
seed <- 2023
ALPHA <- 0.05

RESU <- FINITE_PERFORMANCE(repsMax = repsMax, n = n, theta0 = thetaU0, mu0 = mu0, SIGMA0 = SIGMAU0, p_vec = p_vec, seed = seed, ALPHA = ALPHA, indexB = indexB, CALCULATE_VAR_A_B = CALCULATE_VAR_A_B, BLOCK_HADAMARD_PRODUCT = BLOCK_HADAMARD_PRODUCT, BEST_UNBIASED_ESTIMATOR = BEST_UNBIASED_ESTIMATOR)

RESH <- FINITE_PERFORMANCE(repsMax = repsMax, n = n, theta0 = thetaH0, mu0 = mu0, SIGMA0 = SIGMAH0, p_vec = p_vec, seed = seed, ALPHA = ALPHA, indexB = indexB, CALCULATE_VAR_A_B = CALCULATE_VAR_A_B, BLOCK_HADAMARD_PRODUCT = BLOCK_HADAMARD_PRODUCT, BEST_UNBIASED_ESTIMATOR = BEST_UNBIASED_ESTIMATOR)

RESI <- FINITE_PERFORMANCE(repsMax = repsMax, n = n, theta0 = thetaI0, mu0 = mu0, SIGMA0 = SIGMAI0, p_vec = p_vec, seed = seed, ALPHA = ALPHA, indexB = indexB, CALCULATE_VAR_A_B = CALCULATE_VAR_A_B, BLOCK_HADAMARD_PRODUCT = BLOCK_HADAMARD_PRODUCT, BEST_UNBIASED_ESTIMATOR = BEST_UNBIASED_ESTIMATOR)

save.image("Supp_C4.RData")

biasU <- apply(RESU$theta_EST, 2, mean) - thetaU0
MCSDU <- apply(RESU$theta_EST, 2, sd)
ASEU <- apply(RESU$theta_ASE, 2, mean)
WCPU <- apply(RESU$theta_WCP, 2, mean)
SIGMAU_EST <- BLOCK_HADAMARD_PRODUCT(RESU$A_EST, RESU$B_EST, p_vec)
SU <- RESU$S
round(cbind(biasU, MCSDU, ASEU, WCPU) * 100, 1)

biasH <- apply(RESH$theta_EST, 2, mean) - thetaH0
MCSDH <- apply(RESH$theta_EST, 2, sd)
ASEH <- apply(RESH$theta_ASE, 2, mean)
WCPH <- apply(RESH$theta_WCP, 2, mean)
SIGMAH_EST <- BLOCK_HADAMARD_PRODUCT(RESH$A_EST, RESH$B_EST, p_vec)
SH <- RESH$S
round(cbind(biasH, MCSDH, ASEH, WCPH) * 100, 1)

biasI <- apply(RESI$theta_EST, 2, mean) - thetaI0
MCSDI <- apply(RESI$theta_EST, 2, sd)
ASEI <- apply(RESI$theta_ASE, 2, mean)
WCPI <- apply(RESI$theta_WCP, 2, mean)
SIGMAI_EST <- BLOCK_HADAMARD_PRODUCT(RESI$A_EST, RESI$B_EST, p_vec)
SI <- RESI$S
round(cbind(biasI, MCSDI, ASEI, WCPI) * 100, 1)

writeMat("Supp_C4_Plots.mat", SIGMAU0 = SIGMAU0, SIGMAH0 = SIGMAH0, SIGMAI0 = SIGMAI0, SIGMAU_EST = SIGMAU_EST, SIGMAH_EST = SIGMAH_EST, SIGMAI_EST = SIGMAI_EST, SU = SU, SH = SH, SI = SI)



