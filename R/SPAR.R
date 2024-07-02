###########################################
### Main implementation of sparse projected averaged regression (SPAR)
##########################################

# fix errors for auc evaluation

#' Sparse Projected Averaged Regression
#'
#' Apply Sparse Projected Averaged Regression to High-dimensional Data (see Parzer, Vana-Guer and Filzmoser 2023).
#' This function performs the procedure for given thresholds lambda and numbers of marginal models, and acts as a help-function for the full cross-validated procedure [spar.cv].
#'
#' @param x n x p numeric matrix of predictor variables.
#' @param y quantitative response vector of length n.
#' @param family 'family'-objected used for glm (except the quasi), default gaussian("identity").
#' @param xval optional matrix of predictor variables observations used for validation of threshold lambda and number of models; x is used if not provided.
#' @param yval optional response observations used for validation of threshold lambda and number of models; y is used if not provided.
#' @param nscreen number of variables kept after screening in each marginal model, multiples of n are suggested; defaults to 2n.
#' @param nlambda number of different lambdas to consider for thresholding; ignored when lambdas are given; defaults to 20.
#' @param lambdas optional vector of lambdas to consider for thresholding; if not provided, nlam values ranging from 0 to the maximum ablsolute marginal coefficient are used.
#' @param nummods vector of numbers of marginal models to consider for validation; defaults to c(20).
#' @param split_data logical to indicate whether data for calculation of scr_coef and fitting of mar mods should be split 1/4 to 3/4 to avoid overfitting; default FALSE
#' @param type.measure loss to use for validation; defaults to "deviance" available for all families. Other options are "mse" or "mae" (between responses and predicted means, for all families),
#' "class" (misclassification error) and "1-auc" (one minus area under the ROC curve) both just for "binomial" family.
#' @param type.rpm  type of random projection matrix to be employed; one of "cwdatadriven", "cw", "gaussian", "sparse"; defaults to "cwdatadriven".
#' @param mslow lower bound for unifrom random goal dimensions in marginal models; defaults to log(p).
#' @param msup upper bound for unifrom random goal dimensions in marginal models; defaults to n/2.
#' @param inds optional list of index-vectors corresponding to variables kept after screening in each marginal model of length max(nummods),dimensions need to fit those of RPMs.
#' @param RPMs optional list of sparse CW projection matrices used in each marginal model of length max(nummods), diagonal elements will be overwritten with a coefficient only depending on the given x and y.
#' @returns object of class "spar" with elements
#' \itemize{
#'  \item betas p x max(nummods) matrix of standardized coefficients from each marginal model
#'  \item scr_coef p-vector of HOLP coefficient used for screening
#'  \item inds list of index-vectors corresponding to variables kept after screening in each marginal model of length max(nummods)
#'  \item RPMs list of sparse CW projection matrices used in each marginal model of length max(nummods)
#'  \item val_res data.frame with validation results (validation measure and number of active variables) for each element of lambdas and nummods
#'  \item val_set logical flag, whether validation data were provided; if FALSE, training data were used for validation
#'  \item lambdas vector of lambdas considered for thresholding
#'  \item nummods vector of numbers of marginal models considered for validation
#'  \item xcenter p-vector of empirical means of initial predictor variables
#'  \item xscale p-vector of empirical standard deviations of initial predictor variables
#' }
#' @examples
#' \dontrun{
#' data("example_data")
#' spar_res <- spar(example_data$x,example_data$y,
#' xval=example_data$xtest,yval=example_data$ytest,nummods=c(5,10,15,20,25,30))
#' spar_res
#' coefs <- coef(spar_res)
#' pred <- predict(spar_res,xnew=example_data$x)
#' plot(spar_res)
#' plot(spar_res,"Val_Meas","nummod")
#' plot(spar_res,"Val_numAct","lambda")
#' plot(spar_res,"coefs",prange=c(1,400))}
#' @seealso [spar.cv],[coef.spar],[predict.spar],[plot.spar],[print.spar]
#' @export
#' @importFrom stats rnorm
spar <- function(x,
                 y,
                 family = gaussian("identity"),
                 xval = NULL,
                 yval = NULL,
                 nscreen = 2*nrow(x),
                 nlambda = 20,
                 lambdas = NULL,
                 nummods = c(20),
                 split_data = FALSE,
                 type.measure = c("deviance","mse","mae","class","1-auc"),
                 type.rpm = c("cwdatadriven", "cw", "gaussian", "sparse"),
                 mslow = ceiling(log(ncol(x))),
                 msup = ceiling(nrow(x)/2),
                 inds = NULL,
                 RPMs = NULL,
                 control = list()) {

  stopifnot(mslow <= msup)
  stopifnot(msup <= nscreen)

  stopifnot(is.numeric(y))
  p <- ncol(x)
  n <- nrow(x)
  stopifnot(length(y)==n)

  type.measure <- match.arg(type.measure)
  type.rpm <- match.arg(type.rpm)

  if (split_data==TRUE) {
    scr_inds <- sample(1:n,n%/%4)
    mar_inds <- (1:n)[-scr_inds]
  } else {
    mar_inds <- scr_inds <- 1:n
  }

  xcenter <- apply(x,2,mean)
  xscale <- apply(x,2,sd)

  if (is.null(inds) | is.null(RPMs)) {
    actual_p <- sum(xscale>0)
    z <- scale(x[,xscale>0],center = xcenter[xscale>0],scale = xscale[xscale>0])
  } else {
    actual_p <- p
    xscale[xscale==0] <- 1
    z <- scale(x,center = xcenter,scale = xscale)
  }


  if (family$family=="gaussian" & family$link=="identity") {
    fit_family <- "gaussian"
    ycenter <- mean(y)
    yscale <- sd(y)
    yz <- scale(y,center = ycenter,scale = yscale)
    if (actual_p < n/2) {
      scr_coef <- tryCatch( solve(crossprod(z[scr_inds,]),crossprod(z[scr_inds,],yz[scr_inds])),
                            error=function(error_message) {
                              return(solve(crossprod(z[scr_inds,])+(sqrt(actual_p)+sqrt(n))*diag(actual_p),crossprod(z[scr_inds,],yz[scr_inds])))
                            })
    } else if (actual_p < 2*n) {
      scr_coef <- crossprod(z[scr_inds,],solve(tcrossprod(z[scr_inds,])+(sqrt(actual_p)+sqrt(n))*diag(n),yz[scr_inds]))
    } else {
      solve_res <- tryCatch( solve(tcrossprod(z[scr_inds,]),yz[scr_inds]),
                             error=function(error_message) {
                               return(solve(tcrossprod(z[scr_inds,])+(sqrt(actual_p)+sqrt(n))*diag(n),yz[scr_inds]))
                             })
      scr_coef <- crossprod(z[scr_inds,],solve_res)
    }
  } else {
    if (family$family=="binomial" & family$link=="logit") {
      fit_family <- "binomial"
    } else if (family$family=="poisson" & family$link=="log") {
      fit_family <- "poisson"
    } else {
      fit_family <- family
    }

    ycenter <- 0
    yscale <- 1
    glmnet_res <- glmnet::glmnet(x=z[scr_inds,],y=y[scr_inds],family = fit_family,alpha=0)
    lam <- min(glmnet_res$lambda)
    scr_coef <- coef(glmnet_res,s=lam)[-1]
  }

  inc_probs <- abs(scr_coef)
  max_inc_probs <- max(inc_probs)
  inc_probs <- inc_probs/max_inc_probs

  max_num_mod <- max(nummods)
  intercepts <- numeric(max_num_mod)
  betas_std <- Matrix::Matrix(data=c(0),actual_p,max_num_mod,sparse = TRUE)

  drawRPMs <- FALSE
  if (is.null(RPMs)) {
    RPMs <- vector("list",length=max_num_mod)
    drawRPMs <- TRUE
    ms <- sample(seq(floor(mslow),ceiling(msup)),max_num_mod,replace=TRUE)
  }
  drawinds <- FALSE
  if (is.null(inds)) {
    inds <- vector("list",length=max_num_mod)
    drawinds <- TRUE
  }

  for (i in 1:max_num_mod) {
    if (drawinds) {
      if (nscreen<p) {
        ind_use <- sample(1:actual_p,nscreen,prob=inc_probs)
      } else {
        ind_use <- 1:actual_p
      }
      inds[[i]] <- ind_use
    } else {
      ind_use <- inds[[i]]
    }
    p_use <- length(ind_use)

    if (drawRPMs) {
      m <- ms[i]
      if (p_use < m) {
        m <- p_use
        RPM <- Matrix::Matrix(diag(1,m),sparse=TRUE)
        RPMs[[i]] <- RPM
      } else {
        RPM <- switch(
            type.rpm,
            "cwdatadriven" = generate_cw_rp(m = m, p = p_use, coef = scr_coef[ind_use]/max_inc_probs),
            "cw"           = generate_cw_rp(m = m, p = p_use,
                                            coef = sample(c(-1,1), p_use, replace = TRUE)),
            "gaussian"     = generate_gaussian_rp(m = m, p = p_use),
            "sparse"       = generate_sparse_rp(m = m, p = p_use, psi = control$rpm$psi))
        # RPM <- generate_RPM(m,p_use,coef=scr_coef[ind_use]/max_inc_probs)
        RPMs[[i]] <- RPM
      }
    } else {
      RPM <- RPMs[[i]]
      if (type.rpm == "cwdatadriven")  RPM@x <- scr_coef[ind_use]/max_inc_probs
      ## TODO: think about this
    }

    znew <- Matrix::tcrossprod(z[mar_inds,ind_use],RPM)
    if (family$family=="gaussian" & family$link=="identity") {
      mar_coef <- tryCatch( Matrix::solve(Matrix::crossprod(znew),Matrix::crossprod(znew,yz[mar_inds])),
                            error=function(error_message) {
                              return(Matrix::solve(Matrix::crossprod(znew)+0.01*diag(ncol(znew)),
                                                   Matrix::crossprod(znew,yz[mar_inds])))
                            })
      intercepts[i] <- 0
      betas_std[ind_use,i] <- Matrix::crossprod(RPM,mar_coef)
    } else {
      glmnet_res <- glmnet::glmnet(znew,y[mar_inds],family = fit_family,alpha=0)
      mar_coef <- coef(glmnet_res,s=min(glmnet_res$lambda))
      intercepts[i] <- mar_coef[1]
      betas_std[ind_use,i] <- Matrix::crossprod(RPM,mar_coef[-1])
    }
  }

  if (is.null(lambdas)) {
    if (nlambda>1) {
      lambdas <- c(0,exp(quantile(log(abs(betas_std@x)),probs=1:(nlambda-1)/(nlambda-1))))
    } else {
      lambdas <- c(0)
    }
  } else {
    nlambda <- length(lambdas)
  }

  val_res <- data.frame(nlam=NULL,lam=NULL,nummod=NULL,numAct=NULL,Meas=NULL)
  if (!is.null(yval) & !is.null(xval)) {
    val_set <- TRUE
  } else {
    val_set <- FALSE
    yval <- y
    xval <- x
  }

  if (type.measure=="deviance") {
    val.meas <- function(yval,eta_hat) {
      return(sum(family$dev.resids(yval,family$linkinv(eta_hat),1)))
    }
  } else if (type.measure=="mse") {
    val.meas <- function(yval,eta_hat) {
      return(mean((yval-family$linkinv(eta_hat))^2))
    }
  } else if (type.measure=="mae") {
    val.meas <- function(yval,eta_hat) {
      return(mean(abs(yval-family$linkinv(eta_hat))))
    }
  } else if (type.measure=="class") {
    stopifnot(family$family=="binomial")
    val.meas <- function(yval,eta_hat) {
      return(mean(yval!=round(family$linkinv(eta_hat))))
    }
  } else if (type.measure=="1-auc") {
    stopifnot(family$family=="binomial")
    val.meas <- function(yval,eta_hat) {
      if (var(yval)==0) {
        res <- NA
      } else {
        res <- 1-ROCR::performance(ROCR::prediction(family$linkinv(eta_hat),yval),measure="auc")@y.values[[1]]
      }
      return(res)
    }
  }

  for (nummod in nummods) {
    coef <- betas_std[,1:nummod,drop=FALSE]
    abscoef <- abs(coef)
    tabres <- sapply(1:nlambda, function(l){
      thresh <- lambdas[l]
      tmp_coef <- coef
      tmp_coef[abscoef<thresh] <- 0

      avg_coef <- Matrix::rowMeans(tmp_coef)
      tmp_beta <- numeric(p)
      tmp_beta[xscale>0] <- yscale*avg_coef/(xscale[xscale>0])
      tmp_intercept <- mean(intercepts[1:nummod]) + as.numeric(ycenter - sum(xcenter*tmp_beta) )
      eta_hat <- xval%*%tmp_beta + tmp_intercept
      c(l,
        thresh,
        nummod,
        sum(tmp_beta!=0),
        val.meas(yval,eta_hat)
      )
    })
    rownames(tabres) <- c("nlam","lam","nummod","numAct","Meas")
    val_res <- rbind(val_res,data.frame(t(tabres)))
  }
  betas <- Matrix::Matrix(data=c(0),p,max_num_mod,sparse = TRUE)
  betas[xscale>0,] <- betas_std

  res <- list(betas = betas, intercepts = intercepts, scr_coef = scr_coef, inds = inds, RPMs = RPMs,
              val_res = val_res, val_set = val_set, lambdas = lambdas, nummods = nummods,
              ycenter = ycenter, yscale = yscale, xcenter = xcenter, xscale = xscale,
              family = family, type.measure = type.measure, type.rpm = type.rpm)
  attr(res,"class") <- "spar"

  return(res)
}

#' coef.spar
#'
#' Extract coefficients from spar object
#' @param spar_res result of spar function of class "spar".
#' @param nummod number of models used to form coefficients; value with minimal validation Meas is used if not provided.
#' @param lambda threshold level used to form coefficients; value with minimal validation Meas is used if not provided.
#' @return List of coefficients with elements
#' \itemize{
#'  \item intercept
#'  \item beta
#'  \item nummod
#'  \item lambda
#' }
#' @export

coef.spar <- function(spar_res,
                      nummod = NULL,
                      lambda = NULL) {
  if (is.null(nummod) & is.null(lambda)) {
    best_ind <- which.min(spar_res$val_res$Meas)
    par <- spar_res$val_res[best_ind,]
    nummod <- par$nummod
    lambda <- par$lam
  } else if (is.null(nummod)) {
    if (!lambda %in% spar_res$val_res$lam) {
      stop("Lambda needs to be among the previously fitted values when nummod is not provided!")
    }
    tmp_val_res <- spar_res$val_res[spar_res$val_res$lam==lambda,]
    nummod <- tmp_val_res$nummod[which.min(tmp_val_res$Meas)]
  } else if (is.null(lambda)) {
    if (!nummod %in% spar_res$val_res$nummod) {
      stop("Number of models needs to be among the previously fitted values when lambda is not provided!")
    }
    tmp_val_res <- spar_res$val_res[spar_res$val_res$nummod==nummod,]
    lambda <- tmp_val_res$lam[which.min(tmp_val_res$Meas)]
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

#' predict.spar
#'
#' Predict responses for new predictors from spar object
#' @param spar_res result of spar function of class "spar".
#' @param xnew matrix of new predictor variables; must have same number of columns as x.
#' @param type the type of required predictions; either on response level (default) or on link level
#' @param avg_type type of averaging the marginal models; either on link (default) or on response level
#' @param nummod number of models used to form coefficients; value with minimal validation Meas is used if not provided.
#' @param lambda threshold level used to form coefficients; value with minimal validation Meas is used if not provided.
#' @param coef optional; result of coef.spar, can be used if coef.spar has already been called.
#' @return Vector of predictions
#' @export

predict.spar <- function(spar_res,
                         xnew,
                         type = c("response","link"),
                         avg_type = c("link","response"),
                         nummod = NULL,
                         lambda = NULL,
                         coef = NULL) {
  if (ncol(xnew)!=length(spar_res$xscale)) {
    stop("xnew must have same number of columns as initial x!")
  }
  type <- match.arg(type)
  avg_type <- match.arg(avg_type)
  if (is.null(coef)) {
    coef <- coef(spar_res,nummod,lambda)
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
      final_coef <- spar_res$betas[spar_res$xscale>0,1:coef$nummod,drop=FALSE]
      final_coef[abs(final_coef)<coef$lambda] <- 0

      preds <- sapply(1:coef$nummod,function(j){
        tmp_coef <- final_coef[,j]
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

#' plot.spar
#'
#' Plot errors or number of active variables over different thresholds or number of models of spar result, or residuals vs fitted
#' @param spar_res result of spar function of class "spar".
#' @param plot_type one of c("Val_Measure","Val_numAct","res-vs-fitted","coefs").
#' @param plot_along one of c("lambda","nummod"); ignored when plot_type="res-vs-fitted".
#' @param nummod fixed value for nummod when plot_along="lambda" for plot_type="Val_Measure" or "Val_numAct"; same as for \code{\link{predict.spar}} when plot_type="res-vs-fitted".
#' @param lambda fixed value for lambda when plot_along="nummod" for plot_type="Val_Measure" or "Val_numAct"; same as for \code{\link{predict.spar}} when plot_type="res-vs-fitted".
#' @param xfit data used for predictions in "res-vs-fitted".
#' @param yfit data used for predictions in "res-vs-fitted".
#' @param prange optional vector of length 2 for "coefs"-plot to give the limits of the predictors' plot range; defaults to c(1,p).
#' @param coef_order optional index vector of length p for "coefs"-plot to give the order of the predictors; defaults to 1:p.
#' @return ggplot2 object
#' @import ggplot2
#' @export

plot.spar <- function(spar_res,
                      plot_type = c("Val_Measure","Val_numAct","res-vs-fitted","coefs"),
                      plot_along = c("lambda","nummod"),
                      nummod = NULL,
                      lambda = NULL,
                      xfit = NULL,
                      yfit = NULL,
                      prange = NULL,
                      coef_order = NULL) {
  plot_type <- match.arg(plot_type)
  plot_along <- match.arg(plot_along)
  mynummod <- nummod
  if (plot_type=="res-vs-fitted") {
    if (is.null(xfit) | is.null(yfit)) {
      stop("xfit and yfit need to be provided for res-vs-fitted plot!")
    }
    pred <- predict(spar_res,xfit,nummod,lambda)
    res <- ggplot2::ggplot(data = data.frame(fitted=pred,residuals=yfit-pred),ggplot2::aes(x=fitted,y=residuals)) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0,linetype=2,linewidth=0.5)
  } else if (plot_type=="Val_Measure") {
    if (plot_along=="lambda") {
      if (is.null(nummod)) {
        mynummod <- spar_res$val_res$nummod[which.min(spar_res$val_res$Meas)]
        tmp_title <- "Fixed optimal nummod="
      } else {
        tmp_title <- "Fixed given nummod="
      }
      tmp_df <- dplyr::filter(spar_res$val_res,nummod==mynummod)
      ind_min <- which.min(tmp_df$Meas)

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=nlam,y=Meas)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks=seq(1,nrow(spar_res$val_res),1),labels=round(spar_res$val_res$lam,3)) +
        ggplot2::labs(x=expression(lambda),y=spar_res$type.measure) +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nlam[ind_min],y=tmp_df$Meas[ind_min]),ggplot2::aes(x=x,y=y),col="red") +
        ggplot2::ggtitle(paste0(tmp_title,mynummod))
    } else {
      if (is.null(lambda)) {
        lambda <- spar_res$val_res$lam[which.min(spar_res$val_res$Meas)]
        tmp_title <- "Fixed optimal "
      } else {
        tmp_title <- "Fixed given "
      }
      tmp_df <- dplyr::filter(spar_res$val_res,lam==lambda)
      ind_min <- which.min(tmp_df$Meas)

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=nummod,y=Meas)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::labs(y=spar_res$type.measure) +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nummod[ind_min],y=tmp_df$Meas[ind_min]),ggplot2::aes(x=x,y=y),col="red")+
        ggplot2::ggtitle(substitute(paste(txt,lambda,"=",v),list(txt=tmp_title,v=round(lambda,3))))
    }
  } else if (plot_type=="Val_numAct") {
    if (plot_along=="lambda") {
      if (is.null(nummod)) {
        mynummod <- spar_res$val_res$nummod[which.min(spar_res$val_res$Meas)]
        tmp_title <- "Fixed optimal nummod="
      } else {
        tmp_title <- "Fixed given nummod="
      }
      tmp_df <- dplyr::filter(spar_res$val_res,nummod==mynummod)
      ind_min <- which.min(tmp_df$Meas)

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=nlam,y=numAct)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks=seq(1,nrow(spar_res$val_res),1),labels=round(spar_res$val_res$lam,3)) +
        ggplot2::labs(x=expression(lambda)) +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nlam[ind_min],y=tmp_df$numAct[ind_min]),ggplot2::aes(x=x,y=y),col="red")+
        ggplot2::ggtitle(paste0(tmp_title,mynummod))
    } else {
      if (is.null(lambda)) {
        lambda <- spar_res$val_res$lam[which.min(spar_res$val_res$Meas)]
        tmp_title <- "Fixed optimal "
      } else {
        tmp_title <- "Fixed given "
      }
      tmp_df <- dplyr::filter(spar_res$val_res,lam==lambda)
      ind_min <- which.min(tmp_df$Meas)

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=nummod,y=numAct)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nummod[ind_min],y=tmp_df$numAct[ind_min]),ggplot2::aes(x=x,y=y),col="red")+
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

#' print.spar
#'
#' Print summary of spar result
#' @param spar_res result of spar function of class "spar".
#' @return text summary
#' @export
print.spar <- function(spar_res) {
  mycoef <- coef(spar_res)
  beta <- mycoef$beta
  cat(sprintf("SPAR object:\nSmallest Validation Measure reached for nummod=%d, lambda=%.3f leading to %d / %d active predictors.\n",mycoef$nummod,mycoef$lambda,sum(beta!=0),length(beta)))
  cat("Summary of those non-zero coefficients:\n")
  print(summary(beta[beta!=0]))
}


