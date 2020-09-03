tsdiagnose.tsvets.estimate <- function(object, ...)
{
    Amat <- t(object$spec$vets_env$Amat[,object$spec$vets_env$Amat_index[1]:object$spec$vets_env$Amat_index[2]])
    D <- object$spec$vets_env$Fmat  - object$spec$vets_env$Gmat %*% Amat %*% object$spec$vets_env$Hmat
    E <- eigen(D)
    e <- abs(Re(E$values))
    cat("Real Eigenvalues (D):", round(e, 3))
    cat("\n")
    colnames(object$Error) <- object$spec$target$y_names
    mvn_test <- mvn(as.data.frame(object$Error), univariatePlot = "none", multivariatePlot = "none", 
                    mvnTest = "dh", univariateTest = "SF", multivariateOutlierMethod = "adj", showOutliers = TRUE)
    y_index <- object$spec$target$index
    if (!is.null(mvn_test$multivariateOutliers$Observation) | NROW(mvn_test$multivariateOutliers$Observation) > 0) {
        has_outliers <- TRUE
        mvn_test$multivariateOutliers$Observation <- y_index[as.integer(mvn_test$multivariateOutliers$Observation)]
        mvn_test$multivariateOutliers <- mvn_test$multivariateOutliers[order(mvn_test$multivariateOutliers$Observation),]
    } else {
        has_outliers <- FALSE
    }
    cat("\nMultivariate Normality Tests\n")
    print(mvn_test$multivariateNormality)
    cat("\nUnivariate Normality Tests\n")
    print(mvn_test$univariateNormality)
    if (has_outliers) {
        cat("\nMultivariate Outliers (Mahalanobis Distance)\n")
        print(mvn_test$multivariateOutliers)
    }
    L <- list(D.eigenvalues = e)
    L <- c(L, mvn_test)
    return(invisible(L))
}
