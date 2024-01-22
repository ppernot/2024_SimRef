# Get data ####
dataList = list()
setList = c()
isl = 0

# ## Synthetic data ###
#
# set.seed(123)
#
# unit = ''
# M  = 1e4
# nu = 6
# uE = 0.1 * sqrt(MCMCpack::rinvgamma(M, nu/2, nu/2 ))
#
# ### D = Normal
# E  = uE * rnorm(M)
# set = 'Synth_Norm'
# isl = isl + 1
# setList[isl] = set
# dataList[[paste0(set,'_cal')]] =
#   list(
#     unit = unit,
#     X1   = NA,
#     X2   = NA,
#     E    = E,
#     uE   = uE
#   )
# sdE = sd(E)
# # calScores(E,uE)
#
# ### D = Normal + shift 20%
# bias = 0.2 * sdE
# E  = bias + uE * rnorm(M)
# set = 'Synth_Norm_Shift20pc'
# isl = isl + 1
# setList[isl] = set
# dataList[[paste0(set,'_cal')]] =
#   list(
#     unit = unit,
#     X1   = NA,
#     X2   = NA,
#     E    = E,
#     uE   = uE
#   )
# # calScores(E,uE)
#
# ### D = Normal + shift 20% + adj var
# uEc = sqrt(bias^2 + uE^2)
# E   = bias + uE * rnorm(M)
# set = 'Synth_Norm_Shift20pcAdj'
# isl = isl + 1
# setList[isl] = set
# dataList[[paste0(set,'_cal')]] =
#   list(
#     unit = unit,
#     X1   = NA,
#     X2   = NA,
#     E    = E,
#     uE   = uEc
#   )
# # calScores(E,uEc)
#
# ### D = t4
# E  = uE * rT4(M)
# set = 'Synth_T4'
# isl = isl + 1
# setList[isl] = set
# dataList[[paste0(set,'_cal')]] =
#   list(
#     unit = unit,
#     X1   = NA,
#     X2   = NA,
#     E    = E,
#     uE   = uE
#   )

## PAL2022 ####

unit = '[eV]'
sets = c('Diffusion', 'Perovskite')[1:2]
methods = c('RF', 'LR', 'GPR_Bayesian')[1:3]
cases = c('Test_uncal', 'Test_cal')

for (method in methods) {
  for (set in sets) {
    fRoot = paste0(set, '_', method)
    isl = isl + 1
    setList[isl] = paste0(fRoot,'_Test')
    for (case in cases) {
      dataName = paste0(fRoot, '_', case)
      fName = file.path('..', 'Data', 'PAL2022',
                        paste0(fRoot, '_', case, '.csv'))
      if (!file.exists(fName))
        next
      data = read.csv(fName, header = TRUE)
      sel = data$uE > 1e-6 * sd(data$E)
      if (sum(!sel) > 0) {
        print(c(case, ' non-pos unc:', sum(!sel)))
        data = data[sel,]
      }
      dataList[[dataName]] =
        list(
          unit = unit,
          X1   = NA,
          X2   = NA,
          E    = data$E,
          uE   = data$uE
        )
    }
  }
}
## BUS2022 ####

unit = '[eV]'
sets = c('qm9_qm9_E_uncalibrated_test',
         'qm9_qm9_E_calibrated_isotonic_test')
aliases = c('QM9_E_Test_uncal', 'QM9_E_Test_cal')

isl = isl + 1
setList[isl] = 'QM9_E_Test'
for (i in seq_along(sets)) {
  fName = file.path('..', 'Data', 'BUS2022', paste0(sets[i],'.csv'))
  D = read.table(fName,
                 sep = ',',
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE
  )

  S  = D[, "formula"]
  X1 = CHNOSZ::mass(S)
  X2 = fractHetero(S)
  R  = D[, "E_"]
  C  = D[, "prediction"]
  uC = D[, "uncertainty_total"] ^ 0.5
  E  = R - C
  uE = uC

  dataList[[aliases[i]]] =
    list(
      unit = unit,
      X1   = X1,
      X2   = X2,
      E    = E,
      uE   = uE
    )
}

## RAS2023 ####
### logP ####
unit = ''
sets = c(
  'logP_10k_a_LS-GCN_test',
  # 'logP_10k_a_LS-NN_test',
  'logP_150k_LS-GCN_test'#,
  # 'logP_150k_LS-GCN_test_flex'
)

for (i in seq_along(sets)) {
  isl = isl + 1
  setList[isl] = sets[i]
  fName = file.path('..', 'Data', 'RAS2023', paste0(sets[i],'.csv'))
  D = read.table(fName,
                 sep = ',',
                 colClasses = c('integer','character','numeric',
                                'numeric', 'numeric', 'numeric'),
                 header = TRUE,
                 comment.char = "",
                 check.names = FALSE,
                 stringsAsFactors = FALSE
  )

  R  = D[, "logP"]
  uC = D[, "uq"]
  C  = D[, "y_pred"]
  E  = R - C
  uE = uC
  S  = D[, "smiles"]
  S  = unname(sapply(S,smile2Formula))
  X1 = CHNOSZ::mass(S)
  X2 = fractHetero(S)

  dataList[[paste0(setList[isl],'_cal')]] =
    list(
      unit = unit,
      X1   = X1,
      X2   = X2,
      E    = E,
      uE   = uE
    )
}

# ### IP ####
# unit = '[eV]'
# seed = 'seed42'
# sets = c(
#   paste0('vertIP_evidential_UQ_',seed),
#   paste0('vertIP_evidential_UQ_',seed,'_lam01'),
#   paste0('vertIP_LS_UQ_',seed),
#   paste0('vertIP_LS_UQ_',seed,'_flex')
# )
# reps = c(11, 11, 6, 6)
#
# for (i in seq_along(sets)) {
#   isl = isl + 1
#   setList[isl] = sets[i]
#   fName = file.path('..', 'Data', 'RAS2023', paste0(sets[i],'.csv'))
#   D = read.table(fName,
#                  sep = ',',
#                  colClasses = c('integer','character',rep('numeric',reps[i])),
#                  header = TRUE,
#                  comment.char = "",
#                  check.names = FALSE,
#                  stringsAsFactors = FALSE
#   )
#
#
#   E  = D[, "error"]
#
#   if(!"uq" %in% colnames(D))
#     uE = sqrt(D[, "uq_ale"]^2 + D[, "uq_epi"]^2)
#   else
#     uE = D[, "uq"]
#
#   dataList[[paste0(setList[isl],'_cal')]] =
#     list(
#       unit = unit,
#       X1   = NA,
#       X2   = NA,
#       E    = E,
#       uE   = uE
#     )
# }
save(setList, dataList, file = 'data.Rda')
