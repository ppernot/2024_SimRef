# Functions ####

# Non-stantardized Student's-t
dt_ls <- function(x, df=1, mu=0, sigma=1) 1/sigma * dt((x - mu)/sigma, df)
pt_ls <- function(q, df=1, mu=0, sigma=1)  pt((q - mu)/sigma, df)
qt_ls <- function(p, df=1, mu=0, sigma=1)  qt(p, df)*sigma + mu
rt_ls <- function(n, df=1, mu=0, sigma=1)  rt(n,df)*sigma + mu

fractHetero = function(S) {
  fHetero = c()
  for(i in seq_along(S)) {
    compo = CHNOSZ::count.elements(S[i])
    tot   = sum(compo, na.rm = TRUE)
    CandH = sum(compo['C'], compo['H'], na.rm = TRUE)
    fHetero[i] = 1 - CandH / tot
  }
  fHetero
}

smile2Formula = function(S) {
  molecule <- parse.smiles(S)[[1]]
  convert.implicit.to.explicit(molecule)
  get.mol2formula(molecule,charge=0)@string
}

smile2Fingerprint = function(S) {
  molecule <- parse.smiles(S)
  fp = lapply(molecule, get.fingerprint, type = 'standard')
  fingerprint::fp.to.matrix(fp)
}

rT4 = function(N, df = 4)
  rt(N, df = df) / sqrt(df/(df-2))


calScoresBS1 = function(x, data) {
  # Statistics with known targets
  E = data[x,1]; uE = data[x,2]
  Z = E / uE
  RMV  = sqrt(mean(uE^2))
  RMSE = sqrt(mean(E^2))
  
  c(
    ZMS = mean(Z^2),
    RCE = (RMV - RMSE) / RMV,
    CC  = cor(abs(E), uE, method="spearman")
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
    boot.dist = TRUE, ...)
  if(is.null(cl0))
    stopCluster(cl)
  
  return(bs)
}

fZetaBS = function(bs, target, Utarget=0, method = 'bca') {

  score = bs$t0
  
  delta = score - target
  
  # Pick half-CI closest to the target 
  # (to deal with non-symmetrical CIs)
  lim  = (delta >= 0) * bs[[method]][1,] + 
         (delta <  0) * bs[[method]][2,]
  # Estimate length of half-CI, including target expanded uncertainty
  width = sqrt((score - lim)^2 + Utarget^2)
  return( 
    delta/width
  )
}