

#' @export
coloc.fast <- function(beta1,
                       se1,
                       beta2,
                       se2,
                       priorsd1 = 1,
                       priorsd2 = 1,
                       priorc1 = 1e-4,
                       priorc2 = 1e-4,
                       priorc12 = 1e-5,
                       rounded = 6) {
    abf1 <- abf.Wakefield(beta1, se1, priorsd1, log = TRUE)
    abf2 <- abf.Wakefield(beta2, se2, priorsd2, log = TRUE)
    w <- which(!is.na(abf1) & !is.na(abf2))
    ## Fill in zeros not filter out, thanks Karl and Karsten!!!!
    nv <- length(w)
    abf1 <- norm1(c(0, abf1[w]), log = TRUE)
    abf2 <- norm1(c(0, abf2[w]), log = TRUE)
    res <- data.frame(hypothesis = paste0("H", 0:4),
    label = c("No association",
                "One variant associated with phenotype 1 only",
                "One variant associated with phenotype 2 only",
                "Two variants separately associated with phenotypes 1 and 2",
                "One variant associated with phenotypes 1 and 2"),
    prior = norm1(c(1, priorc1*nv, priorc2*nv, priorc1*priorc2*nv*(nv - 1), priorc12*nv)),
                bf = if (nv > 0) c(abf1[1]*abf2[1],
                sum(abf1[-1])*abf2[1]/nv,
                abf1[1]*sum(abf2[-1])/nv,
                (sum(abf1[-1])*sum(abf2[-1]) - sum(abf1[-1]*abf2[-1]))/(nv*(nv - 1)),
                sum(abf1[-1]*abf2[-1])/nv) else rep(NA, 5))
    res$bf <- res$bf/max(res$bf) # was: res$bf/res$bf[1] # caused division by zero for strongly colocalized signals
    res$posterior <- norm1(res$prior*res$bf)
    ## compute model averaged effect size ratios
    mw <- abf1[-1]*abf2[-1] # model weights
    alpha12 <- sum(beta1[w]/beta2[w]*mw)/sum(mw)
    alpha21 <- sum(beta2[w]/beta1[w]*mw)/sum(mw)
    if (is.finite(rounded)) {
            res$posterior = round(res$posterior, rounded)
    }
    return(list(results = res, nvariants = length(w), alpha12 = alpha12, alpha21 = alpha21))
}


abf.Wakefield <- function(beta, se, priorsd, log = FALSE) {
    if (log) {
        return(log(sqrt(se^2/(se^2 + priorsd^2))) +
            (beta/se)^2/2 * priorsd^2/(se^2 + priorsd^2))
    } else {
        return(sqrt(se^2/(se^2 + priorsd^2)) *
            exp((beta/se)^2/2 * priorsd^2/(se^2 + priorsd^2)))
    }
}


norm1 <- function(x, log = FALSE) {
    if (all(is.na(x))) return(x)
    if (log) {
        x <- x - max(x, na.rm = TRUE)
        x <- exp(x)
    } else {
        ## This does not work if x contains NaNs or +Infs
        stopifnot(all(x >= 0, na.rm = TRUE))
        x <- x / max(x, na.rm = TRUE)
    }
    return(x / sum(x, na.rm = TRUE))
}