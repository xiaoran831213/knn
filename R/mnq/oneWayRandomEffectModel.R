library(MASS)
library(Matrix)
library(rlist)
library(nlme)
library(ape)

# Example: Applying Modified MINQUE to Balanced One-Way Random Effect Model
N  <- 50;                                             # number of samples in each group
m  <- 2;                                              # number of groups
Im <- diag(rep(1,m));                                 # identity matrix of size m
In <- diag(rep(1,N));                                 # identity matrix of size n
Jn <- matrix(1, nrow = N, ncol = N);                  # one matrix of size n
V1 <- Im %x% Jn;                                      # V1 = Im otimes Jn
V2 <- Im %x% In;                                      # V2 = Im otimes In

ViList <- list(V1, V2);



#-------------------Test Correctness of function MINQUE_KKT_Solver------------------------#
A0 <- A <- (1/(m*N))*diag(rep(1, m*N));
t <- 1;


# Correct answer by solving KKT equations directly (time consuming when m*N is large)
V    <- Reduce(f = '+', x = ViList);             # Calculate V matrix
L    <- t(chol(V));                              # Obtain the Cholesky lower triangle
H <- 2*t*(V %x% V) + ginv(A) %x% ginv(A);
Psi <- cbind(c(V1), c(V2))
KKTMat <- cbind(H, Psi)
KKTMat <- rbind(KKTMat, cbind(t(Psi), matrix(0, nrow = 2, ncol = 2)))
g <- 2 * t * (V %x% V) %*% c(A) - c(ginv(A))
bvec <- -c(g,rep(0,2))
xnt <- solve(KKTMat) %*% bvec

# Result from KKT equation 
testOne <- MINQUE_KKT_solver(ViList = ViList, A = A0, t = t, h = FALSE)




#-----------------Test Correctness of function Feasible_MINQUE_Newton----------------------#
testTwo <- Feasible_MINQUE_Newton(ViList = ViList, A = A0, t = t)


#-----------------Test Correctness of function Infeasible_MINQUE_Newton----------------------#
A1        <- diag(rep(1,m*N))
testThree <- Infeasible_MINQUE_Newton(ViList = ViList, A = A0, nu = 0, p = c(1,1), t = 1)
testFour  <- Infeasible_MINQUE_Newton(ViList = ViList, A = A1, nu = 0, p = c(1,1), t = 1)
testFive  <- Infeasible_MINQUE_Newton(ViList = ViList, A = A1, nu = 0, p = c(2,1), t = 1)




#----------------------Compare four different methods for estimating variance components in one-way random effect models----------------------#
sigma_a  <- 1;
sigma    <- 1;
X        <- matrix(1, m*N, 1);
mu       <- 1;

# QR Decomposition of X
qrX      <- qr(X);
rkX      <- rankMatrix(X);
n        <- m*N;
R        <- t(qr.Q(qrX,complete=TRUE)[,(rkX+1):n]);
q        <- n - rkX
V1_M     <- R %*% V1 %*% t(R);
V2_M     <- R %*% V2 %*% t(R);
ViList_M <- list(V1_M, V2_M);
A0       <- diag(rep(1,q));


nS            <- 1000;
mle           <- rep(0, nS);
reml          <- rep(0, nS);
minque        <- rep(0, nS);
M_minque      <- rep(0, nS);
M_minque_time <- rep(0, nS);

mleMat        <- matrix(0, nS, length(ViList));
remlMat       <- matrix(0, nS, length(ViList));
minqueMat     <- matrix(0, nS, length(ViList));
M_minqueMat   <- matrix(0, nS, length(ViList));

p        <- c(2,1);
p1       <- c(1,0);
p2       <- c(0,1);
pm1      <- c(3,1);
pm2      <- c(2,2);
neg_ct_mle <- neg_ct_reml <- neg_ct_minque <- 0;

set.seed(525)
pb <- txtProgressBar(min = 0, max = nS, style = 3)
for(s in 1:nS)
{
  
  # generate one-way random effect model
  a    <- rnorm(m, 0, sigma_a); 
  err  <- rnorm(n,0,sigma);
  Zmat <- diag(rep(1,m)) %x% rep(1,N);
  y    <- X %*% mu + Zmat %*% a + err;
  
  
  # Calculate MLE
  y_bar  <- mean(y);
  yi_bar <- rep(0,m)
  SSE    <- 0;
  SSA    <- 0;
  
  for(i in 1:m)
  {
    yi_vec    <- y[((i-1)*N +1) : (i*N)];
    yi_bar[i] <- mean(yi_vec);
  }
  
  SSE <- crossprod(y - rep(yi_bar, each = N));
  sigma_mle <- MSE <- 1/(m*(N-1)) * SSE;
  
  SSA <- sum(N * (yi_bar - y_bar)^2);
  sigma_a_mle <- (SSA/m - MSE) / N;
  
  if(sigma_a_mle < 0)
  {
    sigma_a_mle <- 0;
    sigma_mle   <- (1/(m*N)) * (SSE+SSA);
    neg_ct_mle  <- neg_ct_mle + 1;
  }
  
  mle[s]     <- crossprod(c(sigma_a_mle, sigma_mle), p);
  mleMat[s,] <- c(sigma_a_mle, sigma_mle);
  
  # Calculate REML
  sigma_reml   <- MSE;
  sigma_a_reml <- (1/N) * (SSA/(m-1) - MSE);
  
  if(sigma_a_reml < 0)
  {
    sigma_a_reml <- 0;
    sigma_reml   <- (1/(m*N-1)) * (SSE+SSA);
    neg_ct_reml  <- neg_ct_reml + 1;
  }
  
  reml[s]     <- crossprod(c(sigma_a_reml, sigma_reml), p);
  remlMat[s,] <- c(sigma_a_reml, sigma_reml);
  
  
  # Calculate MINQUE
  minque[s] <- t(y) %*% LMM_MINQUE_Solver(ViList = ViList, X = X, p = p) %*% y;
  if(minque[s] < 0)
  {
    neg_ct_minque <- neg_ct_minque + 1;
  }
  minque1       <- t(y) %*% LMM_MINQUE_Solver(ViList = ViList, X = X, p = p1) %*% y;
  minque2       <- t(y) %*% LMM_MINQUE_Solver(ViList = ViList, X = X, p = p2) %*% y;
  minqueMat[s,] <- c(minque1, minque2);
  
  # Calculate M_minque
  time1            <- proc.time()
  M_minque_list    <- Infeasible_MINQUE_Newton(ViList = ViList_M, A = A0, nu = 0, p = p, t = 1);
  M_minque[s]      <- t(y) %*% t(R) %*% M_minque_list$A %*% R %*% y;
  M_minque_time[s] <- (proc.time()-time1)["elapsed"];
  
  M_minque_list1   <- Infeasible_MINQUE_Newton(ViList = ViList_M, A = A0, nu = 0, p = pm1, t = 1);
  M_minque_list2   <- Infeasible_MINQUE_Newton(ViList = ViList_M, A = A0, nu = 0, p = pm2, t = 1);
  M_minque1        <- t(y) %*% t(R) %*% M_minque_list1$A %*% R %*% y - t(y) %*% t(R) %*% M_minque_list$A %*% R %*% y;
  M_minque2        <- t(y) %*% t(R) %*% M_minque_list2$A %*% R %*% y - t(y) %*% t(R) %*% M_minque_list$A %*% R %*% y;
  M_minqueMat[s,]  <- c(M_minque1, M_minque2);
  
  
  
  setTxtProgressBar(pb, s)
  
}


mean(mle)
sd(mle)
apply(mleMat, 2, mean)
apply(mleMat, 2, sd)

mean(reml)
sd(reml)
apply(remlMat, 2, mean)
apply(remlMat, 2, sd)

mean(minque)
sd(minque)
apply(minqueMat, 2, mean)
apply(minqueMat, 2, sd)

mean(M_minque)
sd(M_minque)
apply(M_minqueMat, 2, mean)
apply(M_minqueMat, 2, sd)


#----------------------Compare four different methods for estimating variance components in one-way random effect models with Compound Symmetric Covariance matrix----------------------#
rho          <- 0.9;                                 # rho = 0.1, 0.3, 0.7, 0.9
Group        <- rep(1:m, each = N);
group        <- as.factor(Group);
covMat       <- matrix(rho, nrow = m, ncol = m);
diag(covMat) <- 1;

sigma_a  <- 1;
sigma    <- 1;
X        <- matrix(1, m*N, 1);
mu       <- 1;

# QR Decomposition of X
qrX      <- qr(X);
rkX      <- rankMatrix(X);
n        <- m*N;
R        <- t(qr.Q(qrX,complete=TRUE)[,(rkX+1):n]);
q        <- n - rkX;
one_n    <- rep(1, N);
V1       <- (Im %x% one_n) %*% covMat %*% (Im %x% t(one_n));
V2       <- Im %x% In;                                    
ViList   <- list(V1, V2);
V1_M     <- R %*% V1 %*% t(R);
V2_M     <- R %*% V2 %*% t(R);
ViList_M <- list(V1_M, V2_M);
A0       <- diag(rep(1,q));
Zmat     <- diag(rep(1,m)) %x% rep(1,N);

nS            <- 1000;
mle           <- rep(0, nS);
reml          <- rep(0, nS);
minque        <- rep(0, nS);
M_minque      <- rep(0, nS);
M_minque_time <- rep(0, nS);

mleMat        <- matrix(0, nS, length(ViList));
remlMat       <- matrix(0, nS, length(ViList));
minqueMat     <- matrix(0, nS, length(ViList));
M_minqueMat   <- matrix(0, nS, length(ViList));

p        <- c(2,1);
p1       <- c(1,0);
p2       <- c(0,1);
pm1      <- c(3,1);
pm2      <- c(2,2);
neg_ct_mle <- neg_ct_reml <- neg_ct_minque <- 0;

set.seed(525)
pb <- txtProgressBar(min = 0, max = nS, style = 3)
for(s in 1:nS)
{
  a           <- mvrnorm(n = 1, mu = rep(0,m), Sigma = sigma_a*covMat);
  err         <- rnorm(n,0,sigma);
  y           <- X %*% mu + Zmat %*% a + err;
  
  ml_lme      <- lme(y ~ 1, random = ~1|group, method = 'ML', cor = corCompSymm(rho, fixed = TRUE));
  mle[s]      <- crossprod(p, varcomp(ml_lme));
  # mleMat[s,]  <- varcomp(ml_lme); 
  
  reml_lme    <- lme(y ~ 1, random = ~1|group, method = 'REML', cor = corCompSymm(rho, fixed = TRUE));
  reml[s]     <- crossprod(p, varcomp(reml_lme));
  # remlMat[s,] <- varcomp(reml_lme);
  
  minque[s]   <- t(y) %*% LMM_MINQUE_Solver(ViList = ViList, X = X, p = p) %*% y;
  minque1     <- t(y) %*% LMM_MINQUE_Solver(ViList = ViList, X = X, p = p1) %*% y;
  minque2     <- t(y) %*% LMM_MINQUE_Solver(ViList = ViList, X = X, p = p2) %*% y;
  #if(minque1 < 0 || minque2 < 0)
  if(minque[s] < 0)
  {
    neg_ct_minque <- neg_ct_minque + 1;
  }
  # minqueMat[s,] <- c(minque1, minque2)
  
  time1            <- proc.time()
  M_minque_list    <- Infeasible_MINQUE_Newton(ViList = ViList_M, A = A0, nu = 0, p = p, t = 1);
  M_minque[s]      <- t(y) %*% t(R) %*% M_minque_list$A %*% R %*% y;
  M_minque_time[s] <- (proc.time()-time1)["elapsed"];
  
  # M_minque_list1   <- Infeasible_MINQUE_Newton(ViList = ViList_M, A = A0, nu = 0, p = pm1, t = 1);
  # M_minque_list2   <- Infeasible_MINQUE_Newton(ViList = ViList_M, A = A0, nu = 0, p = pm2, t = 1);
  # M_minque1        <- t(y) %*% t(R) %*% M_minque_list1$A %*% R %*% y - t(y) %*% t(R) %*% M_minque_list$A %*% R %*% y;
  # M_minque2        <- t(y) %*% t(R) %*% M_minque_list2$A %*% R %*% y - t(y) %*% t(R) %*% M_minque_list$A %*% R %*% y;
  # M_minqueMat[s,]  <- c(M_minque1, M_minque2);
  
  
  
  setTxtProgressBar(pb, s)
  
  
}




