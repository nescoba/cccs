#' Computes Canonical correlation p-values from regression coefficients
#'
#' @param data A data frame containing the variables to be used in the canonical correlation analysis
#' @param vars1 A vector of variable names for the first set of variables
#' @param vars2 A vector of variable names for the second set of variables
#' @param finwgt Name of the column that contains the survey weights
#' @param n_cc Number of canonical correlations to be computed
#' @param rawcoef_var1 A matrix containing the raw coefficients for the first set of variables
#' @param rawcoef_var2 A matrix containing the raw coefficients for the second set of variables
#' @param ccorr A matrix containing the canonical correlations
#' @param howmany Number of canonical correlations to be computed
#'
#' @return A data frame with the p-values for the canonical correlations
#' @export

csdcanon <- function(data, vars1, vars2, finwgt, n_cc, rawcoef_var1, rawcoef_var2, howmany = NA) {
  # determine which of vars1 and vars2 is the longest
  if (length(vars1) >= length(vars2)) {
    longvars <- vars1
    shortvars <- vars2
  } else {
    longvars <- vars2
    shortvars <- vars1
  }

  # Define OgX as the columns of data corresponding to longvars
  OgX <- data[, longvars]
  # Define OgY as the columns of data corresponding to shortvars
  OgY <- data[, shortvars]

  # Define finwgt as the column of data corresponding to finwgt
  finwgt <- data[, finwgt]
  diag_W <- diag(finwgt)

  # Standarize the colums of OgX and store the result in X
  X <- scale(OgX)
  # Standarize the colums of OgY and store the result in Y
  Y <- scale(OgY)

  OgX <- as.matrix(OgX)
  OgY <- as.matrix(OgY)
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  pvals <- calcpval(OgX, OgY, X, Y, diag_W, finwgt, n_cc, rawcoef_var1, rawcoef_var2)

  Results <- pvals$Results

  colnames(Results) <- c("Statistic", "df1", "df2", "Chi-Sq/F", "p-val", "Index")

  OGStataUVW <- pvals$OGStataUVW

  # calculate weightindex and set survey design using svydesign function from survey package
  weightindex <- (2 * n_cc) + 1

  des <- survey::svydesign(weights = OGStataUVW[[weightindex]], ids = ~1, data = OGStataUVW)

  names_OGStataUVW <- names(OGStataUVW)
  # loop through 1 to n_cc
  for (i in 1:n_cc) {
    secondindex <- i + n_cc
    formula_str_1 <- paste(names_OGStataUVW[i], "~", names_OGStataUVW[secondindex])
    # run svyglm function from survey package
    model1 <- survey::svyglm(as.formula(formula_str_1), design = des)
    # save coefficients to Results matrix
    aux <- summary(model1)$coefficients[, 1]
    Results[(8 * (i - 1)) + 7, 1] <- aux[1]
    Results[(8 * (i - 1)) + 7, 2] <- summary(model1)$df.residual
    Results[(8 * (i - 1)) + 7, 3] <- summary(model1)$df.null
    Results[(8 * (i - 1)) + 7, 4] <- summary(model1)$coefficients[2, 3]
    Results[(8 * (i - 1)) + 7, 5] <- summary(model1)$coefficients[2, 3]
    Results[(8 * (i - 1)) + 7, 6] <- i
  }

  # TODO: complex survey design simple regression
  for (i in 1:n_cc) {
    # secondindex <- i + n_cc
    # formula_str_1 <- paste(names_OGStataUVW[i], "~", names_OGStataUVW[secondindex])
    # # run svyglm function from survey package
    # model1 <- svyglm(as.formula(formula_str_1), design = des)
    # save coefficients to Results matrix
    aux <- NA
    Results[(8 * (i - 1)) + 8, 1] <- aux[1]
    Results[(8 * (i - 1)) + 8, 2] <- NA
    Results[(8 * (i - 1)) + 8, 3] <- NA
    Results[(8 * (i - 1)) + 8, 4] <- NA
    Results[(8 * (i - 1)) + 8, 5] <- NA
    Results[(8 * (i - 1)) + 8, 6] <- i
  }

  if (is.na(howmany)) {
    howmany <- n_cc
  }

  ResultsAux <- Results[, 1:5]

  rownames(ResultsAux) <- rep(NA, nrow(ResultsAux))

  for (i in 1:howmany) {
    # Selecting rows in information for Canonical correlation i
    aux <- Results[(8 * (i - 1)) + 1, 1]
    mag <- round(aux, 0.0001)
    ResultsAux[((8 * (i - 1)) + 2):((8 * (i - 1)) + 8), 1:5] <- Results[((8 * (i - 1)) + 2):((8 * (i - 1)) + 8), 1:5]
    rownames(ResultsAux)[((8 * (i - 1)) + 2):((8 * (i - 1)) + 8)] <- c("Wilks' Lambda", "Pillai's Trace", "Hotelling-Lawley Trace", "Roy's Greatest Root", "Wilks Lambda FREQ", "Weighted Reg", "Complex Survey Design Reg")
    # colnames(ResultsAux) <- c("Estimate", "DF Residuals", "DF Model", "F-statistic", "p-value")
    # print(matlist(ResultsAux, title=paste("Statistics for Canonical Correlation: ", i), twidth=30, format="%10.4f", rowtitle=paste("Canonical Correlation=", mag), border="top bottom"))
  }

  return(ResultsAux)
}


#' Computes the p-values for the canonical correlations
#'
#' @param OgX A matrix containing the first set of variables
#' @param OgY A matrix containing the second set of variables
#' @param X A matrix containing the standardized first set of variables
#' @param Y A matrix containing the standardized second set of variables
#' @param diag_W A matrix containing the diagonal of the survey weights
#' @param finwgt A vector containing the survey weights
#' @param n_cc Number of canonical correlations to be computed
#' @param rawcoef_var1 A matrix containing the raw coefficients for the first set of variables
#' @param rawcoef_var2 A matrix containing the raw coefficients for the second set of variables
#'
#' @return A list containing the p-values for the canonical correlations
calcpval <- function(OgX, OgY, X, Y, diag_W, finwgt, n_cc, rawcoef_var1, rawcoef_var2) {
  output <- list()
  myResults <- matrix(0, nrow = 8 * (n_cc - 1) + 8, ncol = 6)

  myX <- X
  myY <- Y
  mydiag_W <- diag_W
  myOgX <- OgX
  myOgY <- OgY
  myrawcoef_var1 <- rawcoef_var1
  myrawcoef_var2 <- rawcoef_var2
  myW <- finwgt


  # calculate U and V in R
  U <- myOgX %*% myrawcoef_var1
  V <- myOgY %*% myrawcoef_var2

  # store U, V, and UVW matrices in R
  OGStataU <- U
  OGStataV <- V
  UV <- cbind(U, V)
  UVW <- cbind(UV, myW)
  output$OGStataUVW <- as.data.frame(UVW)

  # calculate weighted variance and correlation coefficient
  weighted_varU <- 97 * rep(1, ncol(U))
  weighted_varV <- 98 * rep(1, ncol(V))
  weighted_rho <- 99 * rep(1, ncol(U))

  for (i in 1:ncol(U)) {
    meanU <- mean(U[, i])
    meanV <- mean(V[, i])
    U[, i] <- U[, i] - meanU
    V[, i] <- V[, i] - meanV
    meanU <- mean(U[, i])
    meanV <- mean(V[, i])

    weighted_varU[i] <- t(U[, i]) %*% mydiag_W %*% U[, i] / sum(diag(mydiag_W))
    weighted_varV[i] <- t(V[, i]) %*% mydiag_W %*% V[, i] / sum(diag(mydiag_W))
    weighted_rho[i] <- (t(U[, i]) %*% mydiag_W %*% V[, i] / sum(diag(mydiag_W))) / sqrt(weighted_varU[i] * weighted_varV[i])

    myResults[(8 * (i - 1)) + 1, 1] <- weighted_rho[i] # row 1, 9, ...
    myResults[(8 * (i - 1)) + 1, 2] <- NA
    myResults[(8 * (i - 1)) + 1, 3] <- NA
    myResults[(8 * (i - 1)) + 1, 4] <- NA
    myResults[(8 * (i - 1)) + 1, 5] <- NA
  }

  sumFREQ <- sum(trunc(myW))
  Lambda <- numeric(length(weighted_rho))
  df <- numeric(length(weighted_rho))
  ChiSq_Wilks_Lambda <- numeric(length(weighted_rho))
  ChiSq_Wilks_LambdaFREQ <- numeric(length(weighted_rho))
  p_val_Wilks_Lambda <- numeric(length(weighted_rho))
  p_val_Wilks_LambdaFREQ <- numeric(length(weighted_rho))

  for (i in 1:length(weighted_rho)) {
    Lambda[i] <- (1 - (weighted_rho[i]^2))
    for (j in (i + 1):length(weighted_rho)) {
      Lambda[i] <- Lambda[i] * (1 - (weighted_rho[j]^2))
    }

    df[i] <- (ncol(myX) + 1 - i) * (ncol(myY) + 1 - i)
    ChiSq_Wilks_Lambda[i] <- (((nrow(myX) - 1) - 0.5 * (ncol(myX) + ncol(myY) + 1)) * log(Lambda[i])) * -1
    ChiSq_Wilks_LambdaFREQ[i] <- (((sumFREQ - 1) - 0.5 * (ncol(myX) + ncol(myY) + 1)) * log(Lambda[i])) * -1
    p_val_Wilks_Lambda[i] <- pchisq(ChiSq_Wilks_Lambda[i], df[i], lower.tail = FALSE)
    p_val_Wilks_LambdaFREQ[i] <- pchisq(ChiSq_Wilks_LambdaFREQ[i], df[i], lower.tail = FALSE)
    myResults[(8 * (i - 1)) + 2, 1] <- Lambda[i] # row 2, 10, ...
    myResults[(8 * (i - 1)) + 6, 1] <- Lambda[i] # row 6, 14, ...
    myResults[(8 * (i - 1)) + 2, 2] <- df[i] # row 2, 10, ...
    myResults[(8 * (i - 1)) + 6, 2] <- df[i] # row 6, 14, ...
    myResults[(8 * (i - 1)) + 2, 3] <- NA
    myResults[(8 * (i - 1)) + 6, 3] <- NA
    myResults[(8 * (i - 1)) + 2, 4] <- ChiSq_Wilks_Lambda[i] # row 2, 10, ...
    myResults[(8 * (i - 1)) + 6, 4] <- ChiSq_Wilks_LambdaFREQ[i] # row 6, 14, ...
    myResults[(8 * (i - 1)) + 2, 5] <- p_val_Wilks_Lambda[i] # row 2, 10, ...
    myResults[(8 * (i - 1)) + 6, 5] <- p_val_Wilks_LambdaFREQ[i] # row 6, 14, ...
  }

  Lawley <- numeric(length(weighted_rho))

  for (j in 1:length(weighted_rho)) {
    for (i in 1:j) {
      Lawley[j] <- Lawley[j] + (1 / (weighted_rho[i]^2))
    }
  }

  V <- numeric(length(weighted_rho))
  Pillais_Trace_stat <- numeric(length(weighted_rho))
  Pillais_Trace_p_value <- numeric(length(weighted_rho))

  for (j in 1:length(weighted_rho)) {
    for (i in j:length(weighted_rho)) {
      V[j] <- V[j] + (weighted_rho[i]^2)
    }

    # degrees of freedom
    df <- nrow(myX) - 1 - (2 * j) + Lawley[j]

    # Pillai's Trace statistic
    Pillais_Trace_stat[j] <- df * V[j]

    # Pillai's Trace p-value
    Pillais_Trace_p_value[j] <- pchisq(Pillais_Trace_stat[j],
                                       df = (ncol(myX) + 1 - j) * (ncol(myY) + 1 - j),
                                       lower.tail = FALSE
    )

    myResults[(8 * (j - 1)) + 3, 1] <- V[j]
    myResults[(8 * (j - 1)) + 3, 2] <- (ncol(myX) + 1 - j) * (ncol(myY) + 1 - j)
    myResults[(8 * (j - 1)) + 3, 3] <- NA
    myResults[(8 * (j - 1)) + 3, 4] <- Pillais_Trace_stat[j]
    myResults[(8 * (j - 1)) + 3, 5] <- Pillais_Trace_p_value[j]
  }


  U <- numeric(length(weighted_rho))
  Hotelling_Lawley_Trace_stat <- numeric(length(weighted_rho))
  Hotelling_Lawley_Trace_p_value <- numeric(length(weighted_rho))

  for (j in 1:length(weighted_rho)) {
    for (i in j:length(weighted_rho)) {
      U[j] <- U[j] + (weighted_rho[i]^2) / (1 - (weighted_rho[i]^2))
    }

    # degrees of freedom
    df <- nrow(myX) - ncol(myX) - ncol(myY) - 2 + Lawley[j]

    # Hotelling-Lawley Trace statistic
    Hotelling_Lawley_Trace_stat[j] <- df * U[j]

    # Hotelling-Lawley Trace p-value
    Hotelling_Lawley_Trace_p_value[j] <- pchisq(Hotelling_Lawley_Trace_stat[j],
                                                df = (ncol(myX) + 1 - j) * (ncol(myY) + 1 - j),
                                                lower.tail = FALSE
    )

    myResults[(8 * (j - 1)) + 4, 1] <- U[j]
    myResults[(8 * (j - 1)) + 4, 2] <- (ncol(myX) + 1 - j) * (ncol(myY) + 1 - j)
    myResults[(8 * (j - 1)) + 4, 3] <- NA
    myResults[(8 * (j - 1)) + 4, 4] <- Hotelling_Lawley_Trace_stat[j]
    myResults[(8 * (j - 1)) + 4, 5] <- Hotelling_Lawley_Trace_p_value[j]
  }

  p <- ncol(myX)
  q <- ncol(myY)

  pq <- matrix(1, nrow = 1, ncol = 2)
  pq[1, 1] <- p
  pq[1, 2] <- q

  largest_root <- numeric(length(weighted_rho))
  v1 <- numeric(length(weighted_rho))
  v2 <- numeric(length(weighted_rho))

  Roys_Greatest_Root_stat <- numeric(length(weighted_rho))
  Roys_Greatest_Root_p_value <- numeric(length(weighted_rho))

  # Overall test
  s <- min(pq)
  m <- (abs(p - q) - 1) / 2
  n <- (nrow(myX) - p - q - 2) / 2
  v1[1] <- s + 2 * m + 1
  v2[1] <- s + 2 * n + 1
  largest_root[1] <- weighted_rho[1]^2
  Roys_Greatest_Root_stat[1] <- (v2[1] * largest_root[1]) / (v1[1] * (1 - largest_root[1]))
  Roys_Greatest_Root_p_value[1] <- pf(Roys_Greatest_Root_stat[1], df1 = v1[1], df2 = v2[1], lower.tail = FALSE)

  myResults[5, 1] <- largest_root[1]
  myResults[5, 2] <- v1[1]
  myResults[5, 3] <- v2[1]
  myResults[5, 4] <- Roys_Greatest_Root_stat[1]
  myResults[5, 5] <- Roys_Greatest_Root_p_value[1]

  # Roots k >= 2
  for (k in 2:length(weighted_rho)) {
    s <- min(pq) - k - 1
    m <- (abs(p - q - k) - 1) / 2
    n <- (nrow(myX) - p - q - k) / 2
    v1[k] <- s + 2 * m + 1
    v2[k] <- s + 2 * n + 1
    largest_root[k] <- weighted_rho[k]^2
    Roys_Greatest_Root_stat[k] <- (v2[k] * largest_root[k]) / (v1[k] * (1 - largest_root[k]))
    Roys_Greatest_Root_p_value[k] <- pf(Roys_Greatest_Root_stat[k], df1 = v1[k], df2 = v2[k], lower.tail = FALSE)

    myResults[(8 * (k - 1)) + 5, 1] <- largest_root[k]
    myResults[(8 * (k - 1)) + 5, 2] <- v1[k]
    myResults[(8 * (k - 1)) + 5, 3] <- v2[k]
    myResults[(8 * (k - 1)) + 5, 4] <- Roys_Greatest_Root_stat[k]
    myResults[(8 * (k - 1)) + 5, 5] <- Roys_Greatest_Root_p_value[k]
  }

  for (i in 1:length(weighted_rho)) {
    myResults[(8 * (i - 1)) + 1, 6] <- i # row 1, 9,...
    myResults[(8 * (i - 1)) + 2, 6] <- i # row 2, 10,...
    myResults[(8 * (i - 1)) + 3, 6] <- i # row 3, 11,...
    myResults[(8 * (i - 1)) + 4, 6] <- i # row 4, 12,...
    myResults[(8 * (i - 1)) + 5, 6] <- i # row 5, 13,...
    myResults[(8 * (i - 1)) + 6, 6] <- i # row 6, 14,...
  }



  output$Results <- myResults

  return(output)
}

# data <- read.csv("BCC/toy/fit.csv")
# cc <- cancor(cbind(Chins, Situps, Jumps) ~ Weight + Waist + Pulse, data = data, weights = Weight)

# Results <- csdcanon(data, c("Chins", "Situps", "Jumps"), c("Weight", "Waist", "Pulse"), "Weight", cc$ndim, cc$coef$X, cc$coef$Y)

# stargazer(Results,
#     title = "Complete List of Canonical Correlation p-values", type = "text"
#     # column.labels = c("CC Index", "CC Magnitude", "Wilk's Lambda", "Wilk's Lambda (n)", "Roy's Greatest Root", "Pillai's Trace", "Hotelling-Lawley Trace", "Weighted Regression", "CSD Regression"),
#     # align = TRUE, digits = 4, column.sep.width = "-",
#     # row.sep = "-",
#     # float = FALSE
# )
