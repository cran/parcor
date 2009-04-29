adalasso<-function(X, y, k=10,use.Gram=TRUE,both=TRUE){
    colnames(X)=1:ncol(X)
    fraction = seq(from = 0, to = 1, length = 1000)
    cv.adalasso<-NULL
    if (both==TRUE){ # cross-validation for adaptive lasso
        all.folds <- cv.folds(length(y), k)
        residmat <- matrix(0, length(fraction), k)
        l1.ols.cv<-vector(length=k) # length of ols estimate on cv splits
        for (i in seq(k)) {
            omit <- all.folds[[i]]
            Xtrain<-X[-omit,,drop=FALSE]
            ytrain<-y[-omit]
            Xtest<-X[omit,,drop=FALSE]
            ytest<-y[omit]
            coef.lasso<-mylars(Xtrain,ytrain,k=k,normalize=TRUE,use.Gram=use.Gram)$coefficients
            weights <- 1/abs(coef.lasso[ abs(coef.lasso)>0 ])
            if (length(weights)==0){
                residmat[,i]<-mean((mean(ytrain)-ytest)^2)
            }
            if (length(weights)>0){
                XXtrain <- Xtrain[ , names(weights), drop=FALSE]
                XXtest<-Xtest[ , names(weights), drop=FALSE]
                if ( length(weights)==1 ) {
                    XXtrain <- XXtrain/weights
                    XXtest<-XXtest/weights
                }
                if (length(weights)>1){        
                    XXtrain <- scale(XXtrain, center=FALSE, scale=weights)
                    XXtest<-scale(XXtest, center=FALSE, scale=weights)
                }
                fit<-lars(XXtrain,ytrain,use.Gram=use.Gram,normalize=FALSE)
                ols<-predict(fit, type="coefficients", mode = "fraction", s = 1)$coefficients
                l1.ols.cv[i]<-sum(abs(ols))
                pred<-predict(fit,XXtest,mode="fraction",s=fraction)$fit
                if (length(omit) == 1){
                    pred <- matrix(pred, nrow = 1)
                }
                residmat[, i] <- apply((ytest - pred)^2, 2, mean)
            }
    }
    l1.cv<-mean(l1.ols.cv) # mean length of ols estimate on cv splits
    cv <- apply(residmat, 1, mean)
    s.old<-fraction[which.min(cv)]
    cv.adalasso<-min(cv)
    }
    fit<-mylars(X,y,k=k,fraction=fraction,use.Gram=use.Gram)
    coefficients.lasso=fit$coefficients
    cv.lasso<-fit$cv.lasso
    coefficients.adalasso=NULL
    if (both==TRUE){
        weights <- 1/abs(coefficients.lasso[ abs(coefficients.lasso)>0 ])
        coefficients.adalasso<-rep(0,ncol(X))
        names(coefficients.adalasso)<-1:ncol(X)
        if (length(weights)>0){
            XX <- X[ , names(weights), drop=FALSE]
        if ( length(weights)==1 )  XX <- XX/weights        
        else  XX <- scale(XX, center=FALSE, scale=weights)
        fit<-lars(XX,y,normalize=FALSE,use.Gram=use.Gram)
        ols<-predict(fit,type="coefficients",mode="fraction",s=1)$coefficients
        l1.ols<-sum(abs(ols))
        normalization=l1.cv/l1.ols
        s.opt=min(1,s.old*normalization)
        coefficients=predict(fit,type="coefficients",mode="fraction",s=s.opt)$coefficients
        coefficients.adalasso[names(weights)]<-coefficients/weights
        }
    
    }
    return(list(coefficients.lasso=coefficients.lasso,coefficients.adalasso=coefficients.adalasso,cv.lasso=cv.lasso,cv.adalasso=cv.adalasso))
}
