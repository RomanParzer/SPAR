
#' Sparse Projected Averaged Regression
#'
#' Apply Sparse Projected Averaged Regression to High-dimensional Data (see Parzer, Vana-Guer and Filzmoser 2023).
#'
#' @param x n x p matrix of predictor variables.
#' @param y quantitative response vector of length n.
#' @param family 'family'-objected used for glm (except the quasi), default gaussian("identity").
#' @param nscreen number of variables kept after screening in each marginal model, multiples of n are suggested; defaults to 2n.
#' @param nfolds number of folds to use for cross-validation >2, defaults to 10.
#' @param nlambda number of different lambdas to consider for thresholding; ignored when lambdas are given; defaults to 20.
#' @param lambdas optional vector of lambdas to consider for thresholding; if not provided, nlam values ranging from 0 to the maximum ablsolute marginal coefficient are used.
#' @param nummods vector of numbers of marginal models to consider for validation; defaults to c(20).
#' @param split_data logical to indicate whether data for calculation of scr_coef and fitting of mar mods should be split 1/4 to 3/4 to avoid overfitting; default FALSE
#' @param type.measure loss to use for validation; defaults to "deviance" available for all families. Other options are "mse" or "mae" (between responses and predicted means, for all families),
#' "class" (misclassification error) and "1-auc" (One minus area under the ROC curve) both just for "binomial" family.
#' @param type.rpm  type of random projection matrix to be employed; one of "cwdatadriven", "cw", "gaussian", "sparse"; defaults to "cwdatadriven".
#' @param mslow lower bound for unifrom random goal dimensions in marginal models; defaults to log(p).
#' @param msup upper bound for unifrom random goal dimensions in marginal models; defaults to n/2.
#' @param control a list optional arguments to be passed to functions creating the random projection matrices.
#' @returns object of class "spar" with elements
#' \itemize{
#'  \item betas p x max(nummods) matrix of standardized coefficients from each marginal model
#'  \item scr_coef p-vector of HOLP coefficient used for screening
#'  \item inds list of index-vectors corresponding to variables kept after screening in each marginal model of length max(nummods)
#'  \item RPMs list of sparse CW projection matrices used in each marginal model of length max(nummods)
#'  \item val_sum data.frame with CV results (mean and sd validation measure and mean number of active variables) for each element of lambdas and nummods
#'  \item lambdas vector of lambdas considered for thresholding
#'  \item nummods vector of numbers of marginal models considered for validation
#'  \item ycenter empirical mean of initial response vector
#'  \item yscale empirical standard deviation of initial response vector
#'  \item xcenter p-vector of empirical means of initial predictor variables
#'  \item xscale p-vector of empirical standard deviations of initial predictor variables
#' }
#' @examples
#' \dontrun{
#' data("example_data")
#' spar_res <- spar.cv(example_data$x,example_data$y,nummods=c(5,10,15,20,25,30))
#' spar_res
#' coefs <- coef(spar_res)
#' pred <- predict(spar_res,example_data$x)
#' plot(spar_res)
#' plot(spar_res,"Val_Meas","nummod")
#' plot(spar_res,"Val_numAct","lambda")
#' plot(spar_res,"coefs",prange=c(1,400))}
#' @seealso [spar],[coef.spar.cv],[predict.spar.cv],[plot.spar.cv],[print.spar.cv]
#' @export

spar.cv <- function(x,
                    y,
                    family = gaussian("identity"),
                    nscreen = 2*nrow(x),
                    nfolds = 10,
                    nlambda = 20,
                    lambdas = NULL,
                    nummods = c(20),
                    split_data = FALSE,
                    type.measure = c("deviance","mse","mae","class","1-auc"),
                    type.rpm = c("cwdatadriven", "cw", "gaussian", "sparse"),
                    mslow = ceiling(log(ncol(x))),
                    msup = ceiling(nrow(x)/2),
                    control = list(rpm = NULL)) {
  stopifnot("matrix" %in% class(x) |"data.frame" %in% class(x))
  x <- as.matrix(x)
  if (!class(x[1,1])%in%c("numeric","integer")) {
    stop("There are non-numeric data entries, numerical matrix needed!")
  }
  p <- ncol(x)
  n <- nrow(x)

  SPARres <- spar(x,y,family = family, nscreen = nscreen,nlambda = nlambda,
                  mslow=mslow,msup=msup,nummods=nummods,split_data=split_data,
                  type.measure = type.measure, type.rpm = type.rpm,
                  control = control)

  val_res <- SPARres$val_res
  folds <- sample(cut(1:n,breaks=nfolds,labels=FALSE))
  for (k in 1:nfolds) {
    fold_ind <- which(folds==k)
    foldSPARres <- spar(x[-fold_ind,SPARres$xscale>0],y[-fold_ind],family = family,
                        xval = x[fold_ind,SPARres$xscale>0], yval = y[fold_ind],
                        nscreen = nscreen, lambdas = SPARres$lambdas,
                        mslow = mslow, msup = msup,
                        inds = SPARres$inds, RPMs = SPARres$RPMs,
                        nummods = nummods, split_data = split_data,
                        type.measure = type.measure, type.rpm = type.rpm)
    val_res <- rbind(val_res,foldSPARres$val_res)
  }

  val_sum <- dplyr::group_by(val_res, nlam, lam, nummod)
  suppressMessages(
    val_sum <- dplyr::summarise(val_sum, mMeas = mean(Meas,na.rm=TRUE),
                                sdMeas = sd(Meas,na.rm=TRUE),
                                mNumAct = mean(numAct,na.rm=TRUE))
  )

  res <- list(betas = SPARres$betas, intercepts = SPARres$intercepts,
              scr_coef = SPARres$scr_coef, inds = SPARres$inds,
              RPMs = SPARres$RPMs,
              val_sum = val_sum, lambdas = SPARres$lambdas, nummods=nummods,
              family = family, type.measure = type.measure, type.rpm = type.rpm,
              ycenter = SPARres$ycenter, yscale = SPARres$yscale,
              xcenter = SPARres$xcenter, xscale = SPARres$xscale)
  attr(res,"class") <- "spar.cv"
  return(res)
}

#' coef.spar.cv
#'
#' Extract coefficients from spar object
#' @param spar_res result of spar.cv function of class "spar.cv".
#' @param opt_par one of c("1se","best"), chooses whether to select the best pair of lambdas and nummods according to CV-Meas, or the sparsest solution within one sd of that optimal CV-Meas;
#' ignored when nummod and lambda are given
#' @param nummod optional number of models used to form coefficients
#' @param lambda optional threshold level used to form coefficients
#' @return List of coefficients with elements
#' \itemize{
#'  \item intercept
#'  \item beta
#'  \item nummod
#'  \item lambda
#' }
#' @export

coef.spar.cv <- function(spar_res,
                         opt_par = c("best","1se"),
                         nummod = NULL,
                         lambda = NULL) {
  opt_lamnum <- match.arg(opt_par)
  if (is.null(nummod) & is.null(lambda)) {
    best_ind <- which.min(spar_res$val_sum$mMeas)
    if (opt_lamnum=="1se") {
      allowed_ind <- spar_res$val_sum$mMeas<spar_res$val_sum$mMeas[best_ind]+spar_res$val_sum$sdMeas[best_ind]
      ind_1cv <- which.min(spar_res$val_sum$mNumAct[allowed_ind])
      par <- spar_res$val_sum[allowed_ind,][ind_1cv,]
    } else {
      par <- spar_res$val_sum[best_ind,]
    }
    nummod <- par$nummod
    lambda <- par$lam
  } else if (is.null(nummod)) {
    if (!lambda %in% spar_res$val_sum$lam) {
      stop("Lambda needs to be among the previously fitted values when nummod is not provided!")
    }
    tmp_val_sum <- spar_res$val_sum[spar_res$val_sum$lam==lambda,]
    if (opt_lamnum=="1se") {
      allowed_ind <- tmp_val_sum$mMeas<tmp_val_sum$mMeas[best_ind]+tmp_val_sum$sdMeas[best_ind]
      ind_1cv <- which.min(tmp_val_sum$mNumAct[allowed_ind])
      par <- tmp_val_sum[allowed_ind,][ind_1cv,]
    } else {
      par <- tmp_val_sum[which.min(tmp_val_sum$mMeas),]
    }
    nummod <- par$nummod
  } else if (is.null(lambda)) {
    if (!nummod %in% spar_res$val_res$nummod) {
      stop("Number of models needs to be among the previously fitted values when lambda is not provided!")
    }
    tmp_val_sum <- spar_res$val_sum[spar_res$val_sum$nummod==nummod,]
    if (opt_lamnum=="1se") {
      allowed_ind <- tmp_val_sum$mMeas<tmp_val_sum$mMeas[best_ind]+tmp_val_sum$sdMeas[best_ind]
      ind_1cv <- which.min(tmp_val_sum$mNumAct[allowed_ind])
      par <- tmp_val_sum[allowed_ind,][ind_1cv,]
    } else {
      par <- tmp_val_sum[which.min(tmp_val_sum$mMeas),]
    }
    lambda <- par$lam
  } else {
    if (length(nummod)!=1 | length(lambda)!=1) {
      stop("Length of nummod and lambda must be 1!")
    }
  }

  if (nummod > ncol(spar_res$betas)) {
    warning("Number of models is too high, maximum of fitted is used instead!")
    nummod <- ncol(spar_res$betas)
  }

  # calc for chosen parameters
  final_coef <- spar_res$betas[spar_res$xscale>0,1:nummod,drop=FALSE]
  final_coef[abs(final_coef)<lambda] <- 0
  p <- length(spar_res$xscale)
  beta <- numeric(p)
  beta[spar_res$xscale>0] <- spar_res$yscale*Matrix::rowMeans(final_coef)/(spar_res$xscale[spar_res$xscale>0])
  intercept <- spar_res$ycenter + mean(spar_res$intercepts[1:nummod]) - sum(spar_res$xcenter*beta)
  return(list(intercept=intercept,beta=beta,nummod=nummod,lambda=lambda))
}

#' predict.spar.cv
#'
#' Predict responses for new predictors from spar object
#' @param spar_res result of spar function of class "spar".
#' @param xnew matrix of new predictor variables; must have same number of columns as x.
#' @param type the type of required predictions; either on response level (default) or on link level
#' @param avg_type type of averaging the marginal models; either on link (default) or on response level
#' @param opt_par one of c("best","1se"), chooses whether to select the best pair of lambdas and nummods according to CV-Meas, or the sparsest solution within one sd of that optimal CV-Meas;
#' ignored when nummod and lambda, or coef are given
#' @param nummod number of models used to form coefficients; value with minimal validation Meas is used if not provided.
#' @param lambda threshold level used to form coefficients; value with minimal validation Meas is used if not provided.
#' @param coef optional; result of \code{\link{coef.spar.cv}}, can be used if \code{\link{coef.spar.cv}} has already been called.
#' @return Vector of predictions
#' @export
predict.spar.cv <- function(spar_res,
                            xnew,
                            type = c("response","link"),
                            avg_type = c("link","response"),
                            opt_par = c("best","1se"),
                            nummod = NULL,
                            lambda = NULL,
                            coef = NULL) {
  if (ncol(xnew)!=length(spar_res$xscale)) {
    stop("xnew must have same number of columns as initial x!")
  }
  type <- match.arg(type)
  avg_type <- match.arg(avg_type)
  if (is.null(coef)) {
    coef <- coef(spar_res,opt_par,nummod,lambda)
  }
  if (avg_type=="link") {
    if (type=="link") {
      res <- as.numeric(xnew%*%coef$beta + coef$intercept)
    } else {
      eta <- as.numeric(xnew%*%coef$beta + coef$intercept)
      res <- spar_res$family$linkinv(eta)
    }
  } else {
    if (type=="link") {
      res <- as.numeric(xnew%*%coef$beta + coef$intercept)
    } else {
      # do diff averaging
      final_coef <- spar_res$betas[,1:coef$nummod,drop=FALSE]
      final_coef[abs(final_coef)<coef$lambda] <- 0

      preds <- sapply(1:coef$nummod,function(j){
        tmp_coef <- final_coef[spar_res$xscale>0,j]
        beta <- numeric(length(spar_res$xscale))
        beta[spar_res$xscale>0] <- spar_res$yscale*tmp_coef/(spar_res$xscale[spar_res$xscale>0])
        intercept <- spar_res$ycenter + spar_res$intercepts[j]  - sum(spar_res$xcenter*beta)
        eta <- as.numeric(xnew%*%beta + coef$intercept)
        spar_res$family$linkinv(eta)
      })
      res <- rowMeans(preds)
    }
  }
  return(res)
}

#' plot.spar.cv
#'
#' Plot errors or number of active variables over different thresholds or number of models of spar.cv result, or residuals vs fitted
#' @param spar_res result of spar.cv function of class "spar.cv".
#' @param plot_type one of c("Val_Measure","Val_numAct","res-vs-fitted","coefs").
#' @param plot_along one of c("lambda","nummod"); ignored when plot_type="res-vs-fitted".
#' @param opt_par one of c("1se","best"), chooses whether to select the best pair of lambdas and nummods according to CV-Meas, or the sparsest solution within one sd of that optimal CV-Meas;
#' ignored when nummod and lambda, or coef are given
#' @param nummod fixed value for nummod when plot_along="lambda" for plot_type="Val_Measure" or "Val_numAct"; same as for \code{\link{predict.spar.cv}} when plot_type="res-vs-fitted".
#' @param lambda fixed value for lambda when plot_along="nummod" for plot_type="Val_Measure" or "Val_numAct"; same as for \code{\link{predict.spar.cv}} when plot_type="res-vs-fitted".
#' @param xfit data used for predictions in "res-vs-fitted".
#' @param yfit data used for predictions in "res-vs-fitted".
#' @param opt_par one of c("best","1se"), only needed for plot_type="res-vs-fitted" to set type of predictions, see \code{\link{predict.spar.cv}}.
#' @param prange optional vector of length 2 for "coefs"-plot to give the limits of the predictors' plot range; defaults to c(1,p).
#' @param coef_order optional index vector of length p for "coefs"-plot to give the order of the predictors; defaults to 1:p.
#' @return ggplot2 object
#' @import ggplot2
#' @export
plot.spar.cv <- function(spar_res,
                         plot_type = c("Val_Measure","Val_numAct","res-vs-fitted","coefs"),
                         plot_along = c("lambda","nummod"),
                         nummod = NULL,
                         lambda = NULL,
                         xfit = NULL,
                         yfit = NULL,
                         opt_par = c("best","1-se"),
                         prange = NULL,
                         coef_order = NULL) {
  plot_type <- match.arg(plot_type)
  plot_along <- match.arg(plot_along)
  opt_par <- match.arg(opt_par)
  mynummod <- nummod
  my_val_sum <- dplyr::rename(spar_res$val_sum, Meas="mMeas",numAct="mNumAct")

  if (plot_type=="res-vs-fitted") {
    if (is.null(xfit) | is.null(yfit)) {
      stop("xfit and yfit need to be provided for res-vs-fitted plot!")
    }
    pred <- predict(spar_res,xfit,opt_par=opt_par,nummod=nummod,lambda=lambda)
    res <- ggplot2::ggplot(data = data.frame(fitted=pred,residuals=yfit-pred),ggplot2::aes(x=fitted,y=residuals)) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0,linetype=2,linewidth=0.5)
  } else if (plot_type=="Val_Measure") {
    if (plot_along=="lambda") {
      if (is.null(nummod)) {
        mynummod <- my_val_sum$nummod[which.min(my_val_sum$Meas)]
        tmp_title <- "Fixed optimal nummod="
      } else {
        tmp_title <- "Fixed given nummod="
      }
      tmp_df <- dplyr::filter(my_val_sum,nummod==mynummod)
      ind_min <- which.min(tmp_df$Meas)

      allowed_ind <- tmp_df$Meas<tmp_df$Meas[ind_min]+tmp_df$sdMeas[ind_min]
      ind_1se <- which.min(tmp_df$numAct[allowed_ind])

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=nlam,y=Meas)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks=seq(1,nrow(my_val_sum),1),labels=round(my_val_sum$lam,3)) +
        ggplot2::labs(x=expression(lambda),y=spar_res$type.measure) +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nlam[ind_min],y=tmp_df$Meas[ind_min]),ggplot2::aes(x=x,y=y),col="red") +
        ggplot2::ggtitle(paste0(tmp_title,mynummod)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=Meas-sdMeas,ymax=Meas+sdMeas),alpha=0.2,linetype=2,show.legend = FALSE) +
        ggplot2::geom_point(ggplot2::aes(x = x, y = y),
                   color=2,show.legend = FALSE,
                   data=data.frame(x = c(tmp_df$nlam[ind_min],tmp_df$nlam[allowed_ind][ind_1se]),
                                   y = c(tmp_df$Meas[ind_min],tmp_df$Meas[allowed_ind][ind_1se]))) +
        ggplot2::annotate("segment",x = tmp_df$nlam[ind_min], y = tmp_df$Meas[ind_min] + tmp_df$sdMeas[ind_min],
                          xend = tmp_df$nlam[allowed_ind][ind_1se]+1, yend = tmp_df$Meas[ind_min] + tmp_df$sdMeas[ind_min],
                          color=2,linetype=2)
    } else {
      if (is.null(lambda)) {
        lambda <- my_val_sum$lam[which.min(my_val_sum$Meas)]
        tmp_title <- "Fixed optimal "
      } else {
        tmp_title <- "Fixed given "
      }
      tmp_df <- dplyr::filter(my_val_sum,lam==lambda)
      ind_min <- which.min(tmp_df$Meas)

      allowed_ind <- tmp_df$Meas<tmp_df$Meas[ind_min]+tmp_df$sdMeas[ind_min]
      ind_1se <- which.min(tmp_df$numAct[allowed_ind])

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=nummod,y=Meas)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::labs(y=spar_res$type.measure) +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nummod[ind_min],y=tmp_df$Meas[ind_min]),ggplot2::aes(x=x,y=y),col="red")+
        ggplot2::ggtitle(substitute(paste(txt,lambda,"=",v),list(txt=tmp_title,v=round(lambda,3)))) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=Meas-sdMeas,ymax=Meas+sdMeas),alpha=0.2,linetype=2,show.legend = FALSE)+
        ggplot2::geom_point(ggplot2::aes(x = x, y = y),
                            color=2,show.legend = FALSE,
                            data=data.frame(x = c(tmp_df$nummod[ind_min],tmp_df$nummod[allowed_ind][ind_1se]),
                                            y = c(tmp_df$Meas[ind_min],tmp_df$Meas[allowed_ind][ind_1se]))) +
        ggplot2::annotate("segment",x = tmp_df$nummod[ind_min], y = tmp_df$Meas[ind_min] + tmp_df$sdMeas[ind_min],
                          xend = tmp_df$nummod[allowed_ind][ind_1se], yend = tmp_df$Meas[ind_min] + tmp_df$sdMeas[ind_min],
                          color=2,linetype=2)
    }
  } else if (plot_type=="Val_numAct") {
    if (plot_along=="lambda") {
      if (is.null(nummod)) {
        mynummod <- my_val_sum$nummod[which.min(my_val_sum$Meas)]
        tmp_title <- "Fixed optimal nummod="
      } else {
        tmp_title <- "Fixed given nummod="
      }
      tmp_df <- dplyr::filter(my_val_sum,nummod==mynummod)
      ind_min <- which.min(tmp_df$Meas)

      allowed_ind <- tmp_df$Meas<tmp_df$Meas[ind_min]+tmp_df$sdMeas[ind_min]
      ind_1se <- which.min(tmp_df$numAct[allowed_ind])

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=nlam,y=numAct)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks=seq(1,nrow(my_val_sum),1),labels=round(my_val_sum$lam,3)) +
        ggplot2::labs(x=expression(lambda)) +
        ggplot2::geom_point(ggplot2::aes(x = x, y = y),
                            color=2,show.legend = FALSE,
                            data=data.frame(x = c(tmp_df$nlam[ind_min],tmp_df$nlam[allowed_ind][ind_1se]),
                                            y = c(tmp_df$numAct[ind_min],tmp_df$numAct[allowed_ind][ind_1se]))) +
        ggplot2::ggtitle(paste0(tmp_title,mynummod))
    } else {
      if (is.null(lambda)) {
        lambda <- my_val_sum$lam[which.min(my_val_sum$Meas)]
        tmp_title <- "Fixed optimal "
      } else {
        tmp_title <- "Fixed given "
      }
      tmp_df <- dplyr::filter(my_val_sum,lam==lambda)
      ind_min <- which.min(tmp_df$Meas)

      allowed_ind <- tmp_df$Meas<tmp_df$Meas[ind_min]+tmp_df$sdMeas[ind_min]
      ind_1se <- which.min(tmp_df$numAct[allowed_ind])

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=nummod,y=numAct)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::geom_point(ggplot2::aes(x = x, y = y),
                            color=2,show.legend = FALSE,
                            data=data.frame(x = c(tmp_df$nummod[ind_min],tmp_df$nummod[allowed_ind][ind_1se]),
                                            y = c(tmp_df$numAct[ind_min],tmp_df$numAct[allowed_ind][ind_1se]))) +
        ggplot2::ggtitle(substitute(paste(txt,lambda,"=",v),list(txt=tmp_title,v=round(lambda,3))))

    }
  } else if (plot_type=="coefs") {
    p <- nrow(spar_res$betas)
    nummod <- ncol(spar_res$betas)
    if (is.null(prange)) {
      prange <- c(1,p)
    }
    if (is.null(coef_order)) {
      coef_order <- 1:p
    }


    tmp_mat <- data.frame(t(apply(as.matrix(spar_res$betas)[coef_order,],1,function(row)row[order(abs(row),decreasing = TRUE)])),
                          predictor=1:p)
    colnames(tmp_mat) <- c(1:nummod,"predictor")
    tmp_df <- tidyr::pivot_longer(tmp_mat,tidyselect::all_of(1:nummod),names_to = "marginal model",values_to = "value")
    tmp_df$`marginal model` <- as.numeric(tmp_df$`marginal model`)

    mrange <- max(apply(spar_res$betas,1,function(row)sum(row!=0)))
    res <- ggplot2::ggplot(tmp_df,ggplot2::aes(x=predictor,y=`marginal model`,fill=value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2() +
      ggplot2::coord_cartesian(xlim=prange,ylim=c(1,mrange)) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.border = ggplot2::element_blank())

  } else {
    res <- NULL
  }
  return(res)
}


#' print.spar.cv
#'
#' Print summary of spar.cv result
#' @param spar_res result of spar.cv function of class "spar.cv".
#' @return text summary
#' @export
print.spar.cv <- function(spar_res) {
  mycoef_best <- coef(spar_res,opt_par = "best")
  mycoef_1se <- coef(spar_res,opt_par = "1se")
  cat(sprintf("SPAR.cv object:\nSmallest CV-Meas %.1f reached for nummod=%d, lambda=%.3f leading to %d / %d active predictors.\n",
              min(spar_res$val_sum$mMeas),mycoef_best$nummod,mycoef_best$lambda,sum(mycoef_best$beta!=0),length(mycoef_best$beta)))
  cat("Summary of those non-zero coefficients:\n")
  print(summary(mycoef_best$beta[mycoef_best$beta!=0]))
  cat(sprintf("\nSparsest coefficient within one standard error of best CV-Meas reached for nummod=%d, lambda=%.3f \nleading to %d / %d active predictors with CV-Meas %.1f.\n",
              mycoef_1se$nummod,mycoef_1se$lambda,sum(mycoef_1se$beta!=0),length(mycoef_1se$beta),
              spar_res$val_sum$mMeas[spar_res$val_sum$nummod==mycoef_1se$nummod & spar_res$val_sum$lam==mycoef_1se$lambda]))
  cat("Summary of those non-zero coefficients:\n")
  print(summary(mycoef_1se$beta[mycoef_1se$beta!=0]))
}

