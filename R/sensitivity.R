set.seed(123)

setListCal = setList[7:8]
nMC  = 10000
nBin = 50

# Get stats names
intrv = ErrViewLib::genIntervals(1:10, nBin)
stats = names(calScoresBS2(1:10,cbind(1:10,1:10),intrv))

# Sensitivity to D ####
nuSeq = c(3,4,5,7,10,15,20)
smc     = matrix(NA, nrow = nMC, ncol = length(stats))
scores = uscores =
  array(
    NA,
    dim = c(length(nuSeq), length(setListCal), length(stats)),
    dimnames = list(nuSeq,paste0('Set',1:length(setListCal)),stats)
  )

for(i in seq_along(setListCal)) {
  D2 = dataList[[paste0(setListCal[i],'_cal')]]
  print(setListCal[i])
  uE = D2$uE
  M  = length(uE)

  intrv = ErrViewLib::genIntervals(1:M, nBin)

  for(j in seq_along(nuSeq)) {
    df = nuSeq[j]; print(df)
    for(k in 1:nMC) {
      Ep = uE * rT4(M, df = df)
      smc[k,] = calScoresBS2(1:M, cbind(Ep,uE), intrv)
    }
    scores[j,i,]  = apply(smc, 2, mean, na.rm = TRUE)
    uscores[j,i,] = apply(smc, 2, sd, na.rm = TRUE)/sqrt(nMC)
  }
}


save(nuSeq, stats, setListCal, scores, uscores,
     file ='sensitivity.Rda')
