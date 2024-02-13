# Functions ####

# Non-stantardized Student's-t
dt_ls <- function(x, df=1, mu=0, sigma=1) 1/sigma * dt((x - mu)/sigma, df)
pt_ls <- function(q, df=1, mu=0, sigma=1)  pt((q - mu)/sigma, df)
qt_ls <- function(p, df=1, mu=0, sigma=1)  qt(p, df)*sigma + mu
rt_ls <- function(n, df=1, mu=0, sigma=1)  rt(n,df)*sigma + mu

# Unit-variance Student's-t
rT4 = function(N, df = 4)
  rt(N, df = df) / sqrt(df/(df-2))

calScoresBS1 = function(x, data) {
  # Statistics to be bootstrapped
  E = data[x,1]; uE = data[x,2]
  Z = E / uE
  c(
    ZMS = mean(Z^2),
    CC  = cor(abs(E), uE, method="spearman")
  )
}

calScoresBS2 = function(x, data, intrv, intrvJack) {

  # Binning-dependent statistics with unknown targets
  uE = data[x,2]
  io = order(uE) # order for uE-binning
  uE = uE[io]
  E  = data[x,1][io]
  Z  = E / uE

  int = intrv
  if(length(uE) < NROW(data))
    int = intrvJack

  AE = RMV = RMSE = c()
  # MV = MSE = c()
  for(i in 1:int$nbr) {
    sel     = int$lwindx[i]:int$upindx[i]
    AE[i]   = abs(log(mean(Z[sel]^2)))
    RMV[i]  = sqrt(mean(uE[sel]^2))
    RMSE[i] = sqrt(mean(E[sel]^2))
  }
  ZMSE = mean(AE)
  ENCE = mean( abs(RMV - RMSE) / RMV )

  c(
    ZMS  = mean(Z^2),
    CC   = cor(abs(E),uE, method = 'spearman'),
    ENCE = ENCE,
    ZMSE = ZMSE
  )
}

fPredBS = function(X, statistic, cl = NULL, nBoot = 5000,
                   method = 'bca', ...) {
  cl0 = cl
  if (is.null(cl0))
    cl <- makeCluster(detectCores())
  bs = nptest::np.boot(
    x=1:NROW(X), data = X,
    statistic = statistic, R = nBoot,
    level = 0.95, method = method,
    parallel = TRUE, cl = cl,
    boot.dist = FALSE, ...)
  if(is.null(cl0))
    stopCluster(cl)

  return(bs)
}

fZetaBS = function(bs, target, Utarget=0, method = 'bca') {
  # Estimate zeta-scores
  score = bs$t0
  delta = score - target
  ## Pick half-CI closest to the target
  ## (to deal with non-symmetrical CIs)
  lim  = (delta >= 0) * bs[[method]][1,] +
         (delta <  0) * bs[[method]][2,]
  # Estimate length of half-CI, including target expanded uncertainty
  width = sqrt((score - lim)^2 + Utarget^2)
  return(
    delta/width
  )
}
