context("Checking if the cross-validation is correct")

Sim <- function(xm, ngroup, sigma=rep(1,p)){
  ## xm = mean of the group,
  ## ngroup = number of person in the group,
  ## sigma vector or variances
  ## nb of lines = nb of groups

  p <- ncol(xm)
  q <- nrow(xm)

  if(!inherits(sigma, c("matrix", "Matrix")))
      sigma <- matrix(sigma,nrow=nrow(xm), ncol=p)

  cl <- rep(1:q, ngroup)

  Y <- do.call(rbind, lapply(1:q, function(i) {
    matrix(t(replicate(ngroup[i], xm[i, ] + rnorm(p,0,sigma))),ncol=p)
  }))

  return(list(Y=as.matrix(Y), class=cl))
}

# function that enable to check that our CV code is OK
# doing manually the job with a dirty for loop. 
cv.check <- function(X,class,folds,weights,gamma=0,min.ratio=1e-8,nlambda=100,W=matrix(nrow=0, ncol=0)){

  err <- numeric(0)
  class <-as.factor(class)

  firstrun   <- fusedanova(X, class, weights=weights, gamma=gamma,W=W)
  lambdamax  <- max(unlist(lapply(firstrun@result, function(x) {max(x$table[, "lambda"])})))
  lambdalist <- 10^seq(log10(min.ratio * lambdamax), log10(lambdamax), len = nlambda)
  
  for (i in 1:length(folds)){

    foldtest <- folds[[i]]

    X <- normalize.cv(X,class,foldtest)
    xtest  <- X[foldtest,1]
    ytest  <- class[foldtest]
    ytrain <- class[-foldtest]
    xtrain <- X[-foldtest,1]

    xm <- tapply(xtrain,ytrain,mean)
    xmtest <- tapply(xtest,ytest,mean)
    ngrouptest <-tapply(ytest,ytest,length)

    fa <- fusedanova(as.matrix(xtrain),ytrain,weights=weights,gamma=gamma, W=W,
                     checkargs=FALSE)
    pred  <- toMat(predict(fa,lambda=lambdalist))
    err <- rbind(err,colSums(sweep(pred[ytest, ], 1L, xtest)^2))
  }

  error <- colMeans(err)
  error2 <- colMeans(err^2)
  sd <- sqrt((error2 -error^2)/(nrow(err)-1))
  return (list(cv.error = error, cv.sd =sd))

}

toMat <- function(x){
  nlambda <- length(unique(x[,"lambda"]))
  x <- x[with(x, order(x[,"lambda"], x[,"class"])), "beta"]
  x <- matrix(x, ncol=nlambda)
  return(x)
}

normalize.cv <- function(x,group,omit){
  xtrain <- x[-omit]
  grouptrain <- group[-omit]
  Pooled <- (nlevels(grouptrain)!= length(xtrain))
  if (Pooled){
    ntrain <- length(xtrain)
    m <- mean(xtrain)
    ngrouptrain <- tabulate(grouptrain)
    s <- rowsum(xtrain^2,grouptrain) - (1/ngrouptrain)*(rowsum(xtrain,grouptrain))^2
    s[which(ngrouptrain==1)] <- 0
    if (sum(s)==0){
      res <- (x-mean(xtrain))/as.numeric((sqrt(var(xtrain))))
    } else {
      s = sqrt(sum(s)/(ntrain-length(ngrouptrain)))
      res <- (x-m)/s
    }
  } else{
    res <- (x-mean(xtrain))/as.numeric((sqrt(var(xtrain))))
  }
  return(res)
}

test_that("Consistency between fusedanova cross-validation and a cv by hand", {

  weights <- "laplace" 
  gamma <- 0.5
  xm <- as.matrix(runif(10,min=0,max=20))
  ng <- sample(10:100,size=10,replace=TRUE)
  sigma <- 1
  data <- Sim(xm,ng,sigma)
  class <- data$class
  Y <- data$Y
  
  n <- sum(ng)
  K <- 10
  folds <- split(sample(1:n), rep(1:K, length = n))
  
  # no split
  cvfa  <- crossval(Y,class,folds=folds,weights=weights,gamma=gamma)
  check <- cv.check(Y,class,folds=folds,weights=weights,gamma=gamma)
  
  # with possible splits
  cvfa2  <- crossval(Y,class,folds=folds,weights="naivettest")
  check2 <- cv.check(Y,class,folds=folds,weights="naivettest") 
  
  # test if same error and sd for each lambda
  expect_that(as.vector(cvfa@global$cv.error$err), is_equivalent_to(as.vector(check$cv.error)))
  expect_that(as.vector(cvfa@global$cv.error$sd), is_equivalent_to(as.vector(check$cv.sd)))
  
  expect_that(as.vector(cvfa2@global$cv.error$err), is_equivalent_to(as.vector(check2$cv.error)))
  expect_that(as.vector(cvfa2@global$cv.error$sd), is_equivalent_to(as.vector(check2$cv.sd))) 
  
})