
rm(list = ls())

# #Question #1) ----------------------------------------------------------


#When there are no missing values/hidden variables
library(naivebayes)
dat = read.table("HW2-Q3.txt",header=T)
dat[,-9] = ifelse(dat[,-9]==0, "n", "y") # transforming from binary to catagorical data
nb = naive_bayes(Class~., data=dat) # training
predict(nb, data.frame(goal="y",
                       football="n",
                       golf="n",
                       defense="y",
                       offense="y",
                       wicket="y",
                       office="y",
                       strategy="n", stringsAsFactors=FALSE),
        type="prob")



# the log-likelihood function
logL = function(theta) log((exp(-(2.75-2*theta)^2) + exp(-(2.75-theta)^2)) / 2/ sqrt(pi))

# The lower bound for log-likelihood
LB = function(theta, q){
  return(-q*log(q) - (1-q)*log(1-q) - log(2*sqrt(pi))
         - q * (2.75-2*theta)^2 - (1-q) * (2.75 - theta)^2)
}

EM = function(th0){
  # Calculates the path of the EM algorithm starting from th0
  th = rep(NaN, 200) # 
  q2 = rep(NaN, 200)
  k=0
  th[[k+1]] = th0
  repeat{
    k = k+1
    # E-step
    q2[[2*k-1]] = (exp(-(2.75-2*th[[2*k-1]])^2) 
                   / (exp(-(2.75-2*th[[2*k-1]])^2) 
                      + exp(-(2.75-th[[2*k-1]])^2)))
    # M-step
    th[[2*k]] = 2.75 * (q2[[2*k-1]]+1) / (3*q2[[2*k-1]]+1)
    q2[[2*k]] = q2[[2*k-1]]
    if(abs(th[[2*k]] - th[[2*k-1]])<0.0001){
      break
    }
    th[[2*k+1]] = th[[2*k]]
  }
  # dropping the NaN elements
  th = th[-which(sapply(th, is.nan))]
  q2 = q2[-which(sapply(q2, is.nan))]
  return(list(th=th, q2=q2))
}

# Finding the MLE
dat = read.table("HW3-4-Q1.txt",header=T)
dat[,] = dat[,8]

thMLE = optimize(logL, dat, maximum=TRUE)[[1]]
EM[dat]
  
  
# #Question #2B) ----------------------------------------------------------
#Part 1)
A = matrix(c(0.5, 0.0, 0.0, 0.3, 0.6, 0.0, 0.2, 0.4, 1.0),
           nrow = 3, ncol = 3, byrow = T)

B = matrix(c(0.7, 0.4, 0.8, 0.3, 0.6, 0.2),
           nrow = 2, ncol = 3, byrow = T)

a = c(0.9, 0.1, 0.0)

v_1_3 = c(1,2,1)

A = t(A)
B = t(B)
a = t(a)

alph_v1 = B[,1]*a

alph_v2 = t(A)%*%t(alph_v1)*(B[,2])

alph_v3 = t(A)%*%(alph_v2)*(B[,1])


#Part 1 answer:
P_v_1_3 <- colSums(alph_v3)
P_v_1_3

#Part 2)

bet_v3 = 1
bet_v2 = A%*%B[,2]
bet_v1 = bet_v2*A%*%B[,1]

bet_v1 = bet_v2*t(a*B[,1])

P_h1_v_1_3 = (1/P_v_1_3)*t(alph_v1)*bet_v1

P_h1_v_1_3 = (1/P_v_1_3)*(alph_v2)*bet_v1

##By hand

#p(x/y)/p(y)

p1 <- .1432/.15380
p2 <- .0106/.15380
p3 <- 0/.15380

#Part 2 Answer:
P_h1_v_1_3 <- c(p1,p2,p3)

#part 3)


# #Question #3A) ----------------------------------------------------------

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


# #Question #3B) ----------------------------------------------------------

#exerciseACGTHMM <- function(){

a <- 1
c <- 2 
g <- 3
t <- 4

#p transition matrix 
p <- matrix(0, nrow = 4,ncol = 4)
p[c,a] <- 1
p[g,c] <- 1
p[t,g] <- 1
p[a,t] <- 1

#create vector of 1 matrix multiplication 
one <- vector(length = 4)
one <- c(.25,.25,.25,.25)

#pnew
pnew <- 0.9*p + 0.2*one 

#q transition matrix 
q <-  matrix(0, nrow = 4,ncol = 4)

q[g,t] <- 1
q[c,g] <- 1
q[a,c] <- 1
q[t,a] <- 1

#qnew
qnew <- 0.9*q+0.1*one

#sequence to check
s <- c(a, a, g, t, a, c, t, t, a, c, c, t, a, c, g, c)

#pnewprobability calc
pnewprob <- 0.25
for (t in 2:length(s)){
  pnewprob <- pnewprob*pnew[s[t],s[t-1]]
}

cat('prob of the sequence s for the given by transition matrix pnew is:', pnewprob)

#qnewprobability calc
qnewprob <- 0.25
for (t in 2:length(s)){
  qnewprob <- qnewprob*qnew[s[t],s[t-1]]
}
cat('prob of the sequence s for the given by transition matrix qnew is:', qnewprob)

# #Question #3C) ----------------------------------------------------------

#generate from pnew
v <- matrix(1, nrow = 100, ncol = 16)
for (n in 1:100){
  v[n, 1] <- sample(1:4, 1, prob=c(1, 1, 1, 1), replace=TRUE)
  for (t in 2:16){
    prob_vec = pnew[,v[n,(t-1)]]
    v[n, t] <- sample(1:4, 1, prob= prob_vec, replace=TRUE)
    
  }
}

#generate from qnew
v2 <- matrix(1, nrow = 100, ncol = 16)
for (n in 1:100){
  v2[n, 1] <- sample(1:4, 1, prob=c(1, 1, 1, 1), replace=TRUE)
  for (t in 2:16){
    prob_vec2 = qnew[,v2[n,(t-1)]]
    v2[n, t] <- sample(1:4, 1, prob= prob_vec2, replace=TRUE)
    
  }
}

#create 1 matrix as instructed with first 100 row from pnew , second 100 from qnew 
v <- rbind(v,v2)

# converting to numerical type to make coding a bit easier
dat = data.matrix(v)
X = (dat) # A=1, C=2, G=3, T=4

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

first_cluster_Q3C <- dat[(Bestresult$q[,1])>0.5,]

# The second cluster
dat[(Bestresult$q[,1])<=0.5,]

second_cluster_Q3C <- dat[(Bestresult$q[,1])>0.5,]

#We can see that the first cluster is from the first 100 rows of v (made up from pnew)
first_cluster_Q3C - v[1:100,]

#We can see that the second cluster is from the first 100 rows of v (made up from qnew)
second_cluster_Q3C - v[1:100,]
