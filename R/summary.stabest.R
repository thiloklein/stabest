#' @method summary stabest
#' @export
summary.stabest <- function(object, ...) {
  
  # students' college selection
  student <- data.frame(coef=object$beta.hat, se=sqrt(diag(object$beta.vcov)))
  student$t <- student$coef/student$se
  student$p <- 2*(1 - pnorm(abs(student$t)))
  student[['2.5%']]  <- object$beta.ci95[1,]
  student[['97.5%']] <- object$beta.ci95[2,]
  
  cat("Students' selection of colleges:\n")
  print(student)
  
  # colleges' student selection
  college <- data.frame(coef=object$gamma.hat, se=sqrt(diag(object$gamma.vcov)))
  college$t <- college$coef/college$se
  college$p <- 2*(1 - pnorm(abs(college$t)))
  college[['2.5%']]  <- object$gamma.ci95[1,]
  college[['97.5%']] <- object$gamma.ci95[2,]
  
  cat("colleges' selection of students:\n")
  print(college)
}
