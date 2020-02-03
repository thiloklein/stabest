#' print summary stats
#' 
#' @param x an estimation object of type "stabest"
#' 
#' @return nothing
#' 
#' @export
summary.stabit <- function(x) {
  
  # students' college selection
  student <- data.frame(coef=x$beta.hat, se=sqrt(diag(x$beta.vcov)))
  student$t <- student$coef/student$se
  student$p <- 2*(1 - pnorm(abs(student$t)))
  student[['2.5%']]  <- x$beta.ci95[1,]
  student[['97.5%']] <- x$beta.ci95[2,]
  
  cat("Students' selection of colleges:\n")
  print(student)
  
  # colleges' student selection
  college <- data.frame(coef=x$gamma.hat, se=sqrt(diag(x$gamma.vcov)))
  college$t <- college$coef/college$se
  college$p <- 2*(1 - pnorm(abs(college$t)))
  college[['2.5%']]  <- x$gamma.ci95[1,]
  college[['97.5%']] <- x$gamma.ci95[2,]
  
  cat("colleges' selection of students:\n")
  print(college)
}
