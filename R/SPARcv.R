
#' Sparse Projected Averaged Regression
#'
#' Apply Sparse Projected Averaged Regression to High-dimensional Data (see Parzer, Vana-Guer and Filzmoser 2023).
#'
#' @param x n x p matrix of predictor variables
#' @param y response vector of length n
#' @param family 'family'-objected used for glm (excwept the quasi), default gaussian("identity")
#' @param nscreen number of variables kept after screening in each marginal model, multiples of n are suggested; defaults to 2n.
#' @param nfolds number of folds to use for cross-validation >2, defaults to 10.
#' @param nlambda number of different lambdas to consider for thresholding; ignored when lambdas are given; defaults to 20.
#' @param lambdas optional vector of lambdas to consider for thresholding; if not provided, nlam values ranging from 0 to the maximum ablsolute marginal coefficient are used.
#' @param nummods vector of numbers of marginal models to consider for validation; defaults to c(20).
#' @param split_data logical to indicate whether data for calculation of scr_coef and fitting of mar mods should be split 1/4 to 3/4 to avoid overfitting; default FALSE
#' @param mslow lower bound for unifrom random goal dimensions in marginal models; defaults to log(p).
#' @param msup upper bound for unifrom random goal dimensions in marginal models; defaults to n/2.
#' @returns object of class "spar" with elements
#' \itemize{
#'  \item betas p x max(nummods) matrix of standardized coefficients from each marginal model
#'  \item scr_coef p-vector of HOLP coefficient used for screening
#'  \item inds list of index-vectors corresponding to variables kept after screening in each marginal model of length max(nummods)
#'  \item RPMs list of sparse CW projection matrices used in each marginal model of length max(nummods)
#'  \item val_sum data.frame with CV results (mean and sd Dev and mean number of active variables) for each element of lambdas and nummods
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
#' plot(spar_res,"Dev","nummod")
#' plot(spar_res,"numAct","lambda")}
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
                    mslow = ceiling(log(ncol(x))),
                    msup = ceiling(nrow(x)/2)) {
  stopifnot("matrix" %in% class(x) |"data.frame" %in% class(x))
  x <- as.matrix(x)
  if (!class(x[1,1])%in%c("numeric","integer")) {
    stop("There are non-numeric data entries, numerical matrix needed!")
  }
  p <- ncol(x)
  n <- nrow(x)

  SPARres <- spar(x,y,family = family, nscreen = nscreen,nlambda = nlambda,mslow=mslow,msup=msup,nummods=nummods,split_data=split_data)

  val_res <- SPARres$val_res
  folds <- sample(cut(1:n,breaks=nfolds,labels=FALSE))
  for (k in 1:nfolds) {
    fold_ind <- which(folds==k)
    foldSPARres <- spar(x[-fold_ind,],y[-fold_ind],family = family,
                        xval = x[fold_ind,], yval=y[fold_ind],
                               nscreen = nscreen, lambdas = SPARres$lambdas,mslow=mslow,msup=msup,
                               inds = SPARres$inds, RPMs = SPARres$RPMs,nummods=nummods,split_data=split_data)
    val_res <- rbind(val_res,foldSPARres$val_res)
  }

  val_sum <- dplyr::group_by(val_res,nlam,lam,nummod)
  suppressMessages(
    val_sum <- dplyr::summarise(val_sum,mDev=mean(Dev),sdDev=sd(Dev),mNumAct=mean(numAct))
  )

  res <- list(betas = SPARres$betas, intercepts = SPARres$intercepts, scr_coef = SPARres$scr_coef, inds = SPARres$inds, RPMs = SPARres$RPMs,
              val_sum = val_sum, lambdas = SPARres$lambdas, nummods=nummods, family = family,
              ycenter = SPARres$ycenter, yscale = SPARres$yscale, xcenter = SPARres$xcenter, xscale = SPARres$xscale)
  attr(res,"class") <- "spar.cv"
  return(res)
}

#' coef.spar.cv
#'
#' Extract coefficients from spar object
#' @param spar_res result of spar.cv function of class "spar.cv".
#' @param opt_par one of c("1se","best"), chooses whether to select the best pair of lambdas and nummods according to CV-Dev, or the sparsest solution within one sd of that optimal CV-Dev;
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
    best_ind <- which.min(spar_res$val_sum$mDev)
    if (opt_lamnum=="1se") {
      allowed_ind <- spar_res$val_sum$mDev<spar_res$val_sum$mDev[best_ind]+spar_res$val_sum$sdDev[best_ind]
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
      allowed_ind <- tmp_val_sum$mDev<tmp_val_sum$mDev[best_ind]+tmp_val_sum$sdDev[best_ind]
      ind_1cv <- which.min(tmp_val_sum$mNumAct[allowed_ind])
      par <- tmp_val_sum[allowed_ind,][ind_1cv,]
    } else {
      par <- tmp_val_sum[which.min(tmp_val_sum$mDev),]
    }
    nummod <- par$nummod
  } else if (is.null(lambda)) {
    if (!nummod %in% spar_res$val_res$nummod) {
      stop("Number of models needs to be among the previously fitted values when lambda is not provided!")
    }
    tmp_val_sum <- spar_res$val_sum[spar_res$val_sum$nummod==nummod,]
    if (opt_lamnum=="1se") {
      allowed_ind <- tmp_val_sum$mDev<tmp_val_sum$mDev[best_ind]+tmp_val_sum$sdDev[best_ind]
      ind_1cv <- which.min(tmp_val_sum$mNumAct[allowed_ind])
      par <- tmp_val_sum[allowed_ind,][ind_1cv,]
    } else {
      par <- tmp_val_sum[which.min(tmp_val_sum$mDev),]
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
  final_coef <- spar_res$betas[,1:nummod,drop=FALSE]
  final_coef[abs(final_coef)<lambda] <- 0
  beta <- spar_res$yscale*Matrix::rowMeans(final_coef)/spar_res$xscale
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
#' @param opt_par one of c("1se","best"), chooses whether to select the best pair of lambdas and nummods according to CV-Dev, or the sparsest solution within one sd of that optimal CV-Dev;
#' ignored when nummod and lambda, or coef are given
#' @param nummod number of models used to form coefficients; value with minimal validation Dev is used if not provided.
#' @param lambda threshold level used to form coefficients; value with minimal validation Dev is used if not provided.
#' @param coef optional; result of coef.spar.cv, can be used if coef.spar.cv has already been called.
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
  if (ncol(xnew)!=nrow(spar_res$betas)) {
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
        tmp_coef <- final_coef[,j]
        beta <- spar_res$yscale*tmp_coef/spar_res$xscale
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
#' @param plot_type one of c("Dev","numAct","res-vs-fitted").
#' @param plot_along one of c("lambda","nummod"); ignored when plot_type="res-vs-fitted".
#' @param opt_par one of c("1se","best"), chooses whether to select the best pair of lambdas and nummods according to CV-Dev, or the sparsest solution within one sd of that optimal CV-Dev;
#' ignored when nummod and lambda, or coef are given
#' @param nummod fixed value for nummod when plot_along="lambda" for plot_type="Dev" or "numAct"; same as for predict.spar when plot_type="res-vs-fitted".
#' @param lambda fixed value for lambda when plot_along="nummod" for plot_type="Dev" or "numAct"; same as for predict.spar when plot_type="res-vs-fitted".
#' @param xfit data used for predictions in "res-vs-fitted".
#' @param yfit data used for predictions in "res-vs-fitted".
#' @return ggplot2 object
#' @import ggplot2
#' @export
plot.spar.cv <- function(spar_res,
                         plot_type = c("Dev","numAct","res-vs-fitted"),
                         plot_along = c("lambda","nummod"),
                         opt_par = c("1se","best"),
                         nummod = NULL,
                         lambda = NULL,
                         xfit = NULL,
                         yfit = NULL) {
  plot_type <- match.arg(plot_type)
  plot_along <- match.arg(plot_along)
  mynummod <- nummod
  my_val_sum <- dplyr::rename(spar_res$val_sum, Dev="mDev",numAct="mNumAct")

  if (plot_type=="res-vs-fitted") {
    if (is.null(xfit) | is.null(yfit)) {
      stop("xfit and yfit need to be provided for res-vs-fitted plot!")
    }
    pred <- predict(spar_res,xfit,opt_par,nummod,lambda)
    res <- ggplot2::ggplot(data = data.frame(fitted=pred,residuals=yfit-pred),ggplot2::aes(x=fitted,y=residuals)) +
      ggplot2::geom_point() +
      # ggplot2::geom_smooth(size=0.5,alpha=0.2,method = 'loess',formula='y ~ x') +
      ggplot2::geom_hline(yintercept = 0,linetype=2,size=0.5)
  } else if (plot_type=="Dev") {
    if (plot_along=="lambda") {
      if (is.null(nummod)) {
        mynummod <- my_val_sum$nummod[which.min(my_val_sum$Dev)]
        tmp_title <- "Fixed optimal nummod="
      } else {
        tmp_title <- "Fixed given nummod="
      }
      tmp_df <- dplyr::filter(my_val_sum,nummod==mynummod)
      ind_min <- which.min(tmp_df$Dev)

      allowed_ind <- tmp_df$Dev<tmp_df$Dev[ind_min]+tmp_df$sdDev[ind_min]
      ind_1se <- which.min(tmp_df$numAct[allowed_ind])

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=nlam,y=Dev)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks=seq(1,nrow(my_val_sum),1),labels=round(my_val_sum$lam,3)) +
        ggplot2::labs(x=expression(lambda)) +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nlam[ind_min],y=tmp_df$Dev[ind_min]),ggplot2::aes(x=x,y=y),col="red") +
        ggplot2::ggtitle(paste0(tmp_title,mynummod)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=Dev-sdDev,ymax=Dev+sdDev),alpha=0.2,linetype=2,show.legend = FALSE)+
        ggplot2::scale_y_log10() +
        ggplot2::geom_point(ggplot2::aes(x = x, y = y),
                   color=2,show.legend = FALSE,
                   data=data.frame(x = c(tmp_df$nlam[ind_min],tmp_df$nlam[allowed_ind][ind_1se]),
                                   y = c(tmp_df$Dev[ind_min],tmp_df$Dev[allowed_ind][ind_1se]))) +
        ggplot2::geom_segment(ggplot2::aes(x = tmp_df$nlam[ind_min], y = tmp_df$Dev[ind_min] + tmp_df$sdDev[ind_min],
                                           xend = tmp_df$nlam[allowed_ind][ind_1se]+1, yend = tmp_df$Dev[ind_min] + tmp_df$sdDev[ind_min]),
                              color=2,show.legend = FALSE,linetype=2)
    } else {
      if (is.null(lambda)) {
        lambda <- my_val_sum$lam[which.min(my_val_sum$Dev)]
        tmp_title <- "Fixed optimal "
      } else {
        tmp_title <- "Fixed given "
      }
      tmp_df <- dplyr::filter(my_val_sum,lam==lambda)
      ind_min <- which.min(tmp_df$Dev)

      allowed_ind <- tmp_df$Dev<tmp_df$Dev[ind_min]+tmp_df$sdDev[ind_min]
      ind_1se <- which.min(tmp_df$numAct[allowed_ind])

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=nummod,y=Dev)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nummod[ind_min],y=tmp_df$Dev[ind_min]),ggplot2::aes(x=x,y=y),col="red")+
        ggplot2::ggtitle(substitute(paste(txt,lambda,"=",v),list(txt=tmp_title,v=round(lambda,3)))) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=Dev-sdDev,ymax=Dev+sdDev),alpha=0.2,linetype=2,show.legend = FALSE)+
        ggplot2::scale_y_log10() +
        ggplot2::geom_point(ggplot2::aes(x = x, y = y),
                            color=2,show.legend = FALSE,
                            data=data.frame(x = c(tmp_df$nummod[ind_min],tmp_df$nummod[allowed_ind][ind_1se]),
                                            y = c(tmp_df$Dev[ind_min],tmp_df$Dev[allowed_ind][ind_1se]))) +
        ggplot2::geom_segment(ggplot2::aes(x = tmp_df$nummod[ind_min], y = tmp_df$Dev[ind_min] + tmp_df$sdDev[ind_min],
                                           xend = tmp_df$nummod[allowed_ind][ind_1se], yend = tmp_df$Dev[ind_min] + tmp_df$sdDev[ind_min]),
                              color=2,show.legend = FALSE,linetype=2)
    }
  } else if (plot_type=="numAct") {
    if (plot_along=="lambda") {
      if (is.null(nummod)) {
        mynummod <- my_val_sum$nummod[which.min(my_val_sum$Dev)]
        tmp_title <- "Fixed optimal nummod="
      } else {
        tmp_title <- "Fixed given nummod="
      }
      tmp_df <- dplyr::filter(my_val_sum,nummod==mynummod)
      ind_min <- which.min(tmp_df$Dev)

      allowed_ind <- tmp_df$Dev<tmp_df$Dev[ind_min]+tmp_df$sdDev[ind_min]
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
        lambda <- my_val_sum$lam[which.min(my_val_sum$Dev)]
        tmp_title <- "Fixed optimal "
      } else {
        tmp_title <- "Fixed given "
      }
      tmp_df <- dplyr::filter(my_val_sum,lam==lambda)
      ind_min <- which.min(tmp_df$Dev)

      allowed_ind <- tmp_df$Dev<tmp_df$Dev[ind_min]+tmp_df$sdDev[ind_min]
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
  cat(sprintf("SPAR.cv object:\nSmallest CV-Dev %.1f reached for nummod=%d, lambda=%.3f leading to %d / %d active predictors.\n",
              min(spar_res$val_sum$mDev),mycoef_best$nummod,mycoef_best$lambda,sum(mycoef_best$beta!=0),length(mycoef_best$beta)))
  cat("Summary of those non-zero coefficients:\n")
  print(summary(mycoef_best$beta[mycoef_best$beta!=0]))
  cat(sprintf("\nSparsest coefficient within one standard error of best CV-Dev reached for nummod=%d, lambda=%.3f \nleading to %d / %d active predictors with CV-Dev %.1f.\n",
              mycoef_1se$nummod,mycoef_1se$lambda,sum(mycoef_1se$beta!=0),length(mycoef_1se$beta),
              spar_res$val_sum$mDev[spar_res$val_sum$nummod==mycoef_1se$nummod & spar_res$val_sum$lam==mycoef_1se$lambda]))
  cat("Summary of those non-zero coefficients:\n")
  print(summary(mycoef_1se$beta[mycoef_1se$beta!=0]))
}
