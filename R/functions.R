library("MASS")
library("glmnet")
library("clusterSim")
library("ncvreg")
library("mvtnorm")
library("caret")
library("nnet")
#' @name get_DB_index
#' @rdname get_DB_index
#'
#' @title Discriminative power weights (DB) function
#'
#' @description Discriminative power weights (DB) function
#'
#' @details Caclulates the weight of DB
#'
#'
#' @param x: matrix with continuous covariates
#' @param y: response vector as factor variable
#'
#' @return DP weight DB
#'
#' @references Fuetterer, C., Nalenz, M., and Augustin, T. (2021).
#' Discriminative Power Lasso – Incorporating Discriminative Power of Genes into
#' Regularization-Based Variable Selection. Technical Report. Available under:
#' \url{https://epub.ub.uni-muenchen.de/91666/1/DPL_TR_2022_03.pdf}
#' @export
#'
get_DB_index = function(x, y){
  cluster_index = c()
  x_sc = scale(x)
  for(k in 1:ncol(x)){
    cluster_index[k] = index.DB(x_sc[,k], y)$DB
  }
  cluster_index
}
#' @rdname get_silhouette_index
#' @export
#'
#' @title Discriminative power weights (Si) function
#'
#' @description Discriminative power weights (Si) function
#'
#' @details Caclulates the weight of DB
#'
#'
#' @param x: matrix with continuous covariates
#' @param y: response vector as factor variable
#'
#' @return DP weight Si
#'
#' @references Fuetterer, C., Nalenz, M., and Augustin, T. (2021).
#' Discriminative Power Lasso – Incorporating Discriminative Power of Genes into
#' Regularization-Based Variable Selection. Technical Report. Available under:
#' \url{https://epub.ub.uni-muenchen.de/91666/1/DPL_TR_2022_03.pdf}
#' @export
#'
get_silhouette_index = function(x, y){
  cluster_index = c()
  x_sc = scale(x)
  for(k in 1:ncol(x)){
    cluster_index[k] = summary(silhouette(y, dist(x_sc[,k])))$si.summary[4]
  }
  sil_ind <- abs(1/cluster_index)
  sil_ind
}
#' @rdname get_ANOVA_index
#' @export
#'
#' @title Discriminative power weights (AN) function
#'
#' @description Discriminative power weights (AN) function
#' @details Caclulates the weight of AN
#'
#'
#' @param x: matrix with continuous covariates
#' @param y: response vector as factor variable
#'
#' @return DP weight AN
#'
#' @references Fuetterer, C., Nalenz, M., and Augustin, T. (2021).
#' Discriminative Power Lasso – Incorporating Discriminative Power of Genes into
#' Regularization-Based Variable Selection. Technical Report. Available under:
#' \url{https://epub.ub.uni-muenchen.de/91666/1/DPL_TR_2022_03.pdf}
#' @export
get_ANOVA_index = function(x, y){
  cluster_index <- c()
  x_sc = scale(x)
  for(k in 1:ncol(x)){

    data <-data.frame(gene= x_sc[,k], y =y)
    colnames(data) <- c("gene", "y")

    anov<- summary(aov(gene ~ y, data = data))
    cluster_index[k] <- anov[[1]]$`F value`[1]
  }
  1/cluster_index
}
#' @rdname simdata_mult
#' @export
#'
#' @title Creation of simulation data
#'
#' @description Function for creation of simulation data
#' @details Constructs matrix of continuous covariates and response vector
#'
#'
#' @param p: number of covariates
#' @param non_zero: number of covariates that were not simulated as associated with the outcome
#' @param sigma: value of variance parameter for the construction of x
#' @param coef_val1: value of predictor that is simulated as truly associated with the first category of the response
#' @param coef_val2: value of predictor that is simulated as truly associated with the second and fifth category of the response
#' @param coef_val3: value of predictor that is simulated as truly associated with the third and sixth category of the response
#' @param coef_val4: value of predictor that is simulated as truly associated with the fourth and seventh category of the response
#' @param struc_diff: simulated structure of truly associated predictors
#' @param corr: simulated correlation between covariates
#' @param N: number of simulated number of observations
#'
#'
#' @return Simulated data.frame for K=7 categories with x and y
#'
#' @references Fuetterer, C., Nalenz, M., and Augustin, T. (2021).
#' Discriminative Power Lasso – Incorporating Discriminative Power of Genes into
#' Regularization-Based Variable Selection. Technical Report. Available under:
#' \url{https://epub.ub.uni-muenchen.de/91666/1/DPL_TR_2022_03.pdf}
#' @export
#'
simdata_mult <- function(p, non_zero, sigma, coef_val1, coef_val2, coef_val3, coef_val4, struc_diff, corr, N){
  sig_mat<- matrix(rep(0, p*p), ncol=p)
  diag(sig_mat) <- sigma
  #set.seed(1)




  S = rep(sqrt(sigma), times = p)
  corr_matrix = diag(rep(1, times = p))

  for(k in 1:(p/20)){
    corr_matrix[(k*20-19):(k*20), (k*20-19):(k*20)] <- corr
  }
  diag(corr_matrix) <- 1
  cov_matrix = diag(S) %*% corr_matrix %*% diag(S)

  x <- rmvnorm(n = N, mean=rep(0, times=p), sigma=cov_matrix)
  x <- as.matrix(x)
  x

  #Ziel: 10 Kovariablen, 5 mit outcome assoziierte, rest =0
  not_ass <- p-non_zero
  coef1 <- c(rep(coef_val1, times=non_zero), rep(0, times=not_ass))
  coef2 <- c(rep(coef_val2, times=non_zero), rep(0, times=not_ass))
  if(struc_diff==1){
    coef3 <- c(rep(0, times=10), rep(coef_val3, times=non_zero), rep(0, times=not_ass-10))
    coef3
  }
  if(struc_diff==0){
    coef3 <- c(rep(coef_val3, times=non_zero), rep(0, times=not_ass))
    coef3
  }
  coef4 <- c(rep(coef_val4, times=non_zero), rep(0, times=not_ass))

  coef5 <- coef2
  coef6 <- coef3
  coef7 <- coef4

  lp1 <- x%*%coef1
  lp2 <- x%*%coef2
  lp3 <- x%*%coef3
  lp4 <- x%*%coef4
  lp5 <- x%*%coef1
  lp6 <- x%*%coef2
  lp7 <- x%*%coef3



  #Daraus die Klassenzugehörigkeit generieren.
  #Wir wollen die Wahrscheinlichkeit von Klasse 1 direkt berücksichtigen und mit K=3 und nicht K=3-1 simulaieren
  # denominator, ensures probabilities sum to 1
  den <- (exp(lp1) + exp(lp2) + exp(lp3)+ exp(lp4) + exp(lp5) + exp(lp6) + exp(lp7))
  p1 <- exp(lp1)/den
  p2 <- exp(lp2)/den
  p3 <- exp(lp3)/den
  p4 <- exp(lp4)/den
  p5 <- exp(lp5)/den
  p6 <- exp(lp6)/den
  p7 <- exp(lp7)/den

  P <- cbind(p1, p2, p3, p4, p5, p6, p7)
  P
  dim(P)

  #Random sample for assigning 4 categories based on
  y<- c()
  for(i in 1:dim(x)[1]){
    tmp <- rmultinom(1, size = 1, prob = P[i,])
    y[i] <- which(tmp==1)
  }

  prop.table(table(y))

  d <- data.frame(x, y = factor(y))
  d

}
