# a function that finds the optimal penality parameter lambda
obtainOptimalLambda <- function(Y, lambda.set = seq(0, 1, 0.01), k=6, settings.for.spcov){
  lambda.set[lambda.set==0] <- 0.001
  n = nrow(Y)
  k.folds <- create.K.folds(n, k)
  k <- length(k.folds)

  ll.lambda <- t(sapply(lambda.set, function(lbd){
    #print(lbd)
    ll.lbd.fold <- sapply(seq_len(k), function(i){
      Yi  <- Y[k.folds[[i]],]
      Yic <- Y[!(seq_len(n) %in% k.folds[[i]]),]


      Z1 <- spcov::spcov(Sigma = diag(apply(Yic, 2, var)),
                  S = var(Yic)+diag(rep(settings.for.spcov$eps, ncol(Yic))),
                  lambda = lbd,
                  step.size = settings.for.spcov$step.size,
                  trace = 0)

      Z2 <- var(Yi)+diag(rep(settings.for.spcov$eps, ncol(Yi)))

      ll.lbd <- -log(det(Z1$Sigma)) - trace.mat(Z2 %*% solve(Z1$Sigma))
      ll.lbd
    })
    mean(ll.lbd.fold)
  }))
  lambda.chosen <- lambda.set[which.max(ll.lambda)]
  if(length(lambda.chosen)>=1){
    lambda.chosen <- sample(lambda.chosen, 1)
  }
}

# A function to make k-folds for the cross validation
create.K.folds <- function(n, k){
  # a minimum of 6 plates are required
  if(n<6) stop("The number of plates must be 6 or more!")
  if(n<=k){
    k <- floor(n/3)
  }else if(n>k & n<= k*2){
    k <- 3
  }else if(n>k*2 & n<= k*4){
    k <- 5
  }
  suppressWarnings(split(seq_len(n), f = seq_len(k)))
}


# A function th at calculates the trace of a function
trace.mat <- function(X){sum(diag(X))}
