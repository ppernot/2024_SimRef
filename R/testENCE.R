# ENCE and ZMSE for synthetic datasets

ENCEfun = function(x, data, intrv) {
  # !!!! REQUIRES ORDERED uE !!!!!
  E  = data[x,1]
  uE = data[x,2]

  RMV = RMSE = c()
  for(i in 1:intrv$nbr) {
    sel     = intrv$lwindx[i]:intrv$upindx[i]
    RMV[i]  = sqrt(mean(uE[sel]^2))
    RMSE[i] = sqrt(mean(E[sel]^2))
  }
  ENCE  = mean( abs(RMV - RMSE) / RMV )
  return(ENCE)
}
ZMSEfun = function(x, data, intrv) {
  # !!!! REQUIRES ORDERED uE !!!!!
  E  = data[x,1]
  uE = data[x,2]
  Z  = E / uE

  AE = c()
  for(i in 1:intrv$nbr) {
    sel   = intrv$lwindx[i]:intrv$upindx[i]
    AE[i] = abs(log(mean(Z[sel]^2)))
  }
  ZMSE = mean(AE)
  return(ZMSE)
}

# Parameters #####
nMC     = 5000 # Nb of Monte Carlo samples
nuSeq   = c(3:6,12,24) # Shape of uncertainties IG distribution
mSeq    = c(2000,4000,8000,12000,16000) # Dataset sizes
nBinSeq = 10*(1:5) # Binning sizes(for equal-size bins)

# D = N(0,1) ####
smct =  matrix(NA, nrow = nMC, ncol = 2)
scorest = uscorest =
  array(
    NA,
    dim = c(length(mSeq), length(nuSeq), length(nBinSeq), 2),
    dimnames = list(mSeq, nuSeq, nBinSeq, 1:2)
  )
for(i in seq_along(mSeq)) {
  Mt = mSeq[i]
  for(j in seq_along(nuSeq)) {
    nut = nuSeq[j]
    for(k in seq_along(nBinSeq)) {
      nBint = nBinSeq[k]
      print(c(Mt,nut,nBint))
      intrvt = ErrViewLib::genIntervals(1:Mt, nBint)
      for(kt in 1:nMC) {
        uEt = sort(sqrt(MCMCpack::rinvgamma(Mt,nut,nut)))
        Et  = uEt * rnorm(Mt)
        smct[kt,1:2] = ENCEfun(1:Mt, cbind(Et,uEt), intrvt)
      }
      scorest[i,j,k,]  = apply(smct,2,mean)
      uscorest[i,j,k,] = apply(smct,2,sd) / sqrt(nMC)
    }
  }
}
save(nMC, nuSeq, mSeq, nBinSeq, scorest, uscorest, ENCEfun,
     file = file.path(tmpDir,"testENCE.Rda"))


smct =  matrix(NA, nrow = nMC, ncol = 2)
scorest = uscorest =
  array(
    NA,
    dim = c(length(mSeq), length(nuSeq), length(nBinSeq), 2),
    dimnames = list(mSeq, nuSeq, nBinSeq, 1:2)
  )
for(i in seq_along(mSeq)) {
  Mt = mSeq[i]
  for(j in seq_along(nuSeq)) {
    nut = nuSeq[j]
    for(k in seq_along(nBinSeq)) {
      nBint = nBinSeq[k]
      print(c(Mt,nut,nBint))
      intrvt = ErrViewLib::genIntervals(1:Mt, nBint)
      for(kt in 1:nMC) {
        uEt = sort(sqrt(MCMCpack::rinvgamma(Mt,nut,nut)))
        Et  = uEt * rnorm(Mt)
        smct[kt,1:2] = ZMSEfun(1:Mt, cbind(Et,uEt), intrvt)
      }
      scorest[i,j,k,]  = apply(smct,2,mean)
      uscorest[i,j,k,] = apply(smct,2,sd) / sqrt(nMC)
    }
  }
}
save(nMC, nuSeq, mSeq, nBinSeq, scorest, uscorest, ZMSEfun,
     file = file.path(tmpDir,"testZMSE.Rda"))


# D = ts(6) ####
smct =  matrix(NA, nrow = nMC, ncol = 2)
scorest = uscorest =
  array(
    NA,
    dim = c(length(mSeq), length(nuSeq), length(nBinSeq), 2),
    dimnames = list(mSeq, nuSeq, nBinSeq, 1:2)
  )
for(i in seq_along(mSeq)) {
  Mt = mSeq[i]
  for(j in seq_along(nuSeq)) {
    nut = nuSeq[j]
    for(k in seq_along(nBinSeq)) {
      nBint = nBinSeq[k]
      print(c(Mt,nut,nBint))
      intrvt = ErrViewLib::genIntervals(1:Mt, nBint)
      for(kt in 1:nMC) {
        uEt = sort(sqrt(MCMCpack::rinvgamma(Mt,nut,nut)))
        Et  = uEt * rT4(Mt, df = 6)
        smct[kt,1:2] = ENCEfun(1:Mt, cbind(Et,uEt), intrvt)
      }
      scorest[i,j,k,]  = apply(smct,2,mean)
      uscorest[i,j,k,] = apply(smct,2,sd) / sqrt(nMC)
    }
  }
}
save(nMC, nuSeq, mSeq, nBinSeq, scorest, uscorest, ENCEfun,
     file = file.path(tmpDir,"testENCE_T6.Rda"))


smct =  matrix(NA, nrow = nMC, ncol = 2)
scorest = uscorest =
  array(
    NA,
    dim = c(length(mSeq), length(nuSeq), length(nBinSeq), 2),
    dimnames = list(mSeq, nuSeq, nBinSeq, 1:2)
  )
for(i in seq_along(mSeq)) {
  Mt = mSeq[i]
  for(j in seq_along(nuSeq)) {
    nut = nuSeq[j]
    for(k in seq_along(nBinSeq)) {
      nBint = nBinSeq[k]
      print(c(Mt,nut,nBint))
      intrvt = ErrViewLib::genIntervals(1:Mt, nBint)
      for(kt in 1:nMC) {
        uEt = sort(sqrt(MCMCpack::rinvgamma(Mt,nut,nut)))
        Et  = uEt * rT4(Mt, df = 6)
        smct[kt,1:2] = ZMSEfun(1:Mt, cbind(Et,uEt), intrvt)
      }
      scorest[i,j,k,]  = apply(smct,2,mean)
      uscorest[i,j,k,] = apply(smct,2,sd) / sqrt(nMC)
    }
  }
}
save(nMC, nuSeq, mSeq, nBinSeq, scorest, uscorest, ZMSEfun,
     file = file.path(tmpDir,"testZMSE_T6.Rda"))

