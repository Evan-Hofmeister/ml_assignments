#-------------------------------------------------------------------------------------#
#                                                                                     #
#                                 Gene Clustering                                     #           #
#                                                                                     #
#-------------------------------------------------------------------------------------#
dat = read.table("genes.txt",header=T)

# converting to numerical type to make coding a bit easier
X = data.matrix(dat) # A=1, C=2, G=3, T=4

# the log-likelihood function
logL = function(X_, p_, pi_, M_){
  N=dim(X_)[1]
  T0=dim(X_)[2]
  K=length(p_)
  sum_n = 0.0
  for (n in 1:N) {
    sum_h = 0.0
    for (h in 1:K) {
      sum_t = 0.0
      for (t in 2:T0) {
        sum_t = sum_t + log(M_[X_[[n,(t-1)]], X_[[n,t]], h]) # prod(A_i) = exp(sum(log(A_i))) to avoid underflow
      }
      sum_h = sum_h + p_[[h]] * pi_[[h,X_[[n,1]]]] * exp(sum_t)
    }
    sum_n = sum_n + log(sum_h)
  }
  return(sum_n)
}

# E-Step
EStep = function(X_, p_, pi_, M_){
  N=dim(X_)[1]
  T0=dim(X_)[2]
  K=length(p_)
  # The variable to store q_(h_n|X_n)
  q_ = matrix(rep(NA, 2*(dim(X_)[[1]])),
              nrow = dim(X_)[[1]], ncol = 2,
              dimnames = list(1:dim(X_)[[1]],
                              c("h=1", "h=2")))
  for (n in 1:N) {
    for (h in 1:K) {
      sum_t = 0.0
      for (t in 2:T0) {
        sum_t = sum_t + log(M_[X_[[n,(t-1)]], X_[[n,t]], h]) # prod(A_i) = exp(sum(log(A_i))) to avoid underflow
      }
      q_[n,h] = p_[[h]] * pi_[[h,X_[[n,1]]]] * exp(sum_t)
    }
    #normalizing q_(h|X_n)
    q_[n,] = q_[n,] / sum(q_[n,]) 
  }
  return(q_)
}

# M_-step
MStep = function(X_,q_){
  # p_ = (P(H=1), P(H=2))
  N=dim(X_)[1]
  T0=dim(X_)[2]
  K=2
  p_ = rep(NA, K)
  names(p_) = c("H=1", "H=2")
  # pi_[h,i] = P(X_{n,1}=i|H=h)
  pi_ = matrix(rep(NA, 8),
               nrow = 2, ncol = 4,
               dimnames = list(c("h=1", "h=2"),
                               c("X1=A", "X1=C", "X1=G", "X1=T"))) # P(X_1=i|H=1)
  # M_[i,j,h] = P(X_{n,t}=j|X_{n,t-1}=i,H=h)
  M_ = array(NA, c(4,4,2),
             dimnames = list(c("X_(t-1)=A", "X_(t-1)=C", "X_(t-1)=G", "X_(t-1)=T"),
                             c("X_(t)=A", "X_(t)=C", "X_(t)=G", "X_(t)=T"),
                             c("h=1", "h=2")))
  for (h in 1:K){
    for (i in 1:4) {
      idx_pi = (X_[,1]==i)
      pi_[h,i] = sum(q_[idx_pi,h])
      for (j in 1:4) {
        M_[i,j,h]=0
        for (n in 1:N) {
          sum_t = sum(((X_[n,]==i)[-N])*((X_[n,]==j)[-1])) # = sum_2^T 1_{X_{n,(t-1)}=i, X_{n,t}=j}
          M_[i,j,h] = M_[i,j,h] + q_[[n,h]] * sum_t  
        }
      }
      #normalizing M_
      M_[i,,h] = M_[i,,h] / sum(M_[i,,h])
    }
    p_[h] = sum(q_[,h])
    #normalizing pi_
    pi_[h,] = pi_[h,] / p_[h]
  }
  # normalizing p_
  p_ = p_ / sum(p_)
  return(list(p=p_, pi=pi_, M=M_))
}

# EM algorithm
EM = function(X, p0, pi0, M0,maxitr){
  p = p0
  pi = pi0
  M = M0
  new_logL = logL(X, p, pi, M)
  itr = 0
  repeat{
    q = EStep(X, p, pi, M)
    param = MStep(X,q)
    p = param$p
    pi = param$pi
    M = param$M
    old_logL = new_logL
    new_logL = logL(X, p, pi, M)
    itr=itr+1
    if((abs((new_logL-old_logL)/new_logL)<0.00001) | (itr>maxitr)){
      break
    }
  }
  return(list(p=p, pi=pi, M=M, q=q, itr=itr, loglik=new_logL))
}

# Initial value of the parameters
# p = (P(H=1), P(H=2))
p = c(0.6, 0.4)
names(p) = c("H=1", "H=2")
# pi[h,i] = P(X_{n,1}=i|H=h)
pi = matrix(c(0.25, 0.25, 0.25, 0.25,
              0.5, 0.3, 0.1, 0.1), # try (0.5, 0.3, 0.1, 0.1) and (0.25, 0.25, 0.25, 0.25)
            nrow = 2, ncol = 4,
            dimnames = list(c("h=1", "h=2"),
                            c("X1=A", "X1=C", "X1=G", "X1=T"))) # P(X_1=i|H=1)
# M[i,j,h] = P(X_{n,t}=j|X_{n,t-1}=i,H=h)
M = array(0.25, c(4,4,2),
          dimnames = list(c("X(t-1)=A", "X(t-1)=C", "X(t-1)=G", "X(t-1)=T"),
                          c("X(t)=A", "X(t)=C", "X(t)=G", "X(t)=T"),
                          c("h=1", "h=2")))

res = EM(X,p,pi,M,1000)
res$itr
res$loglik

totalruns = 100
BestLogL = -Inf
for (run in 1:totalruns) {
  # initializing parameters randomly
  p = runif(2)
  pi = matrix(runif(8),
              nrow = 2, ncol = 4,
              dimnames = list(c("h=1", "h=2"),
                              c("X1=A", "X1=C", "X1=G", "X1=T"))) # P(X_1=i|H=1)
  M = array(runif(32), c(4,4,2),
            dimnames = list(c("X(t-1)=A", "X(t-1)=C", "X(t-1)=G", "X(t-1)=T"),
                            c("X(t)=A", "X(t)=C", "X(t)=G", "X(t)=T"),
                            c("h=1", "h=2")))
  # normalizing p, pi, and M
  p = p / sum(p)
  for (h in 1:2) {
    pi[h,] = pi[h,] / sum(pi[h,])
    for (i in 1:4) {
      M[i,,h] = M[i,,h]/sum(M[i,,h])
    }
  }
  # running the EM algorithm
  res = EM(X,p,pi,M,1000)
  # updating the best estimates, if necessary
  if (res$loglik > BestLogL){
    Bestresult = res
    BestLogL = res$loglik
  }
}
Bestresult$loglik
# The first cluster
dat[(Bestresult$q[,1])>0.5,]

# The second cluster
dat[(Bestresult$q[,1])<=0.5,]
