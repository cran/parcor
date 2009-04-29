
mylars<-function (X, y, k = 10, fraction = seq(from = 0, to = 1, length = 1000),use.Gram=TRUE,normalize=TRUE) 
{
    x<-X
    all.folds <- cv.folds(length(y), k)
    residmat <- matrix(0, length(fraction), k)
    l1.ols.cv<-vector(length=k) # length of ols estimate on cv splits
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        fit <- lars(x[-omit, ,drop=FALSE], y[-omit],use.Gram=use.Gram,normalize=normalize)
      ols<-predict(fit, type="coefficients", mode = "fraction", 
            s = 1)$coefficients
    l1.ols.cv[i]<-sum(abs(ols))
        fit <- predict(fit, x[omit, , drop = FALSE], mode = "fraction", 
            s = fraction)$fit
        if (length(omit) == 1) 
            fit <- matrix(fit, nrow = 1)
        residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
    }
    l1.cv<-mean(l1.ols.cv) # mean length of ols estimate on cv splits
    cv <- apply(residmat, 1, mean)
    cv.lasso<-min(cv)
    cv.error <- sqrt(apply(residmat, 1, var)/k)
    s.old<-fraction[which.min(cv)]
    fit<-lars(x,y,use.Gram=use.Gram,normalize=normalize)
    ols<-predict(fit,type="coefficients",mode="fraction",s=1)$coefficients
    l1.ols<-sum(abs(ols))
    normalization=l1.cv/l1.ols
    s.opt=min(1,s.old*normalization)
    coefficients=predict(fit,type="coefficients",mode="fraction",s=s.opt)$coefficients
    
    object <- list(coefficients=coefficients,cv.lasso=cv.lasso)
    invisible(object)
}
