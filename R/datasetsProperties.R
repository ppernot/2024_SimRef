# Calculations ####
cl <- makeCluster(detectCores())

size = meanE = sigE = meanZ = sigZ =  biasE = biasZ = c()
EdistPars = ZdistPars = matrix(NA,ncol=3,nrow=length(setList))
colnames(EdistPars) = colnames(ZdistPars) = c("df","mu","sigma")

for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]; print(setList[i])
  E  = D2$E
  uE = D2$uE
  Z  = E / uE
  size[i] = M = length(uE)

  # E
  mu  = mean(E)
  sig = sd(E)/sqrt(length(E))
  meanE[i] = ErrViewLib::prettyUnc(mu,sig)

  bs = nptest::np.boot(
    x = E, statistic = sd,
    R = nBoot, level = 0.95, method = "bca",
    parallel = TRUE, cl = cl)
  mus  = bs$t0
  sigs = bs$se
  sigE[i] = ErrViewLib::prettyUnc(mus,sigs)
  biasE[i] = abs(mu/mus)

  ## Fit by Stutent's t to get effective df
  fit.t = fitdistrplus::fitdist(
    E, "t_ls", start = list(df=20,mu=mean(E),sigma=sd(E)),
    keepdata = FALSE)
  pars = summary(fit.t)$estimate
  EdistPars[i,] = summary(fit.t)$estimate

  # Z
  mu  = mean(Z)
  sig = sd(Z)/sqrt(length(Z))
  meanZ[i] = ErrViewLib::prettyUnc(mu,sig)

  bs = nptest::np.boot(
    x = Z, statistic = sd,
    R = nBoot, level = 0.95, method = "bca",
    parallel = TRUE, cl = cl)
  mus  = bs$t0
  sigs = bs$se
  sigZ[i] = ErrViewLib::prettyUnc(mus,sigs)
  biasZ[i] = abs(mu/mus)

  ## Fit by Stutent's t to get effective df
  fit.t<-fitdistrplus::fitdist(
    Z, "t_ls", start = list(df=20,mu=mean(Z),sigma=sd(Z)),
    keepdata = FALSE)
  pars = summary(fit.t)$estimate
  ZdistPars[i,] = summary(fit.t)$estimate
}
stopCluster(cl)


