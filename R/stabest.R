#' Estimation of student and college preferences in a m:1 school market. 
#' Identifying restrictions are imposed by (i) stability conditions of the observed matching, and 
#' (ii) the submitted rank order lists (which may be incomplete or contain unacceptable students / schools).
#' 
#' WARNING: The method 
#' 
#' @param student.prefs formula object, rank ~ explanatory variables. Encoding: lower rank = better, -1 = unacceptable / last, NA = unobserved. 
#' Careful when including factor variables with mutually exclusive levels - one level must be missing so as to not introduce an intercept through the back door!
#' @param college.prefs formula object, rank ~ explanatory variables. Encoding: lower rank = better, -1 = unacceptable, NA = unobserved
#' #' Careful when including factor variables with mutually exclusive levels - one level must be missing so as to not introduce an intercept through the back door!
#' @param student.id variable containing the IDs of the students
#' @param college.id variable containing the IDs of the schools
#' @param match.id variable indicating a match (1=match, 0=no match). If NULL, no stability conditions will be imposed.
#' @param data a data.frame
#' @param nSeats numeric vector indicating the number of seats for each school (order by college.id)
#' @param demean whether to de-mean variables by decision making unit (beta)
#' @param niter number of data augmentation iterations
#' @param thin thinning parameter
#' @param burnin number of discarded samples (in the thinned series)
#' @param initparm if not NULL, this must be a list with elements beta (coefficients driving students' prefs)
#' and gamma (coefficients driving schools' decisions). The number of coefficients must match the number of parameters
#' supplied via student.prefs and school.prefs, respectively (including factors).
#' @param nCores number of parallel cores to use (currently only 1 is permitted)
#' 
#' @return An object of type "stabest" 
#' 
#' @examples put example here
#' 
#' @import Matrix
#' @import Rcpp
#' @useDynLib stabest
#' 
#' 
#' @export
stabest <- function(student.prefs, college.prefs, 
                    student.id='s.id', college.id='c.id', match.id=NULL,
                    data, nSeats=NULL, demean=FALSE, 
                    niter=500, thin=10, burnin=50, initparm=NULL, nCores=1) {
  
  
  
  # parse student preference equation
  
  student.prefs.terms <- terms(student.prefs)
  stopifnot(attr(student.prefs.terms, 'intercept')==0) 														# intercept is not identified
  if  (all.vars(student.prefs)[1]!='.') {
    student.rankvar <- all.vars(student.prefs)[1]	
  }  else {    
    student.rankvar <- NULL
  }
  student.varnames <- attr(student.prefs.terms, "term.labels")
  student.regform <- update(student.prefs, Vc~.) 															# regression formula for latent valuations over colleges
  
  
  
  
  # parse college preference equation
  
  college.prefs.terms <- terms(college.prefs)																# parse equation object
  stopifnot(attr(college.prefs.terms, 'intercept')==0)														# intercept is not identified
  if (all.vars(college.prefs)[1]!='.') {
    college.rankvar <- all.vars(college.prefs)[1]				# check if rank information is provided or not
  } else {
    college.rankvar <- NULL
  }
  college.varnames <- attr(college.prefs.terms, "term.labels")
  college.regform <- update(college.prefs, Vs~.) 															# regression formula for latent valuations over colleges
  
  
  
  
  # bring data set in correct order 
  
  data <- data[order(data[[college.id]], data[[student.id]]),]
  data$ID_ <- 1:nrow(data)  # R style index
  data$sid_ <- as.numeric(factor(data[[student.id]])) 														# generate generic consecutive student and college IDs
  #data <- data[order(data[[college.id]], data[[student.id]]),]
  data$cid_ <- as.numeric(factor(data[[college.id]]))
  #data <- with(data, data[order(cid_, sid_),])
  nColleges <- length(unique(data$cid_))
  nStudents <- length(unique(data$sid_))
  
  
  
  # find ID_ of next college for each student/college pair
  # since data will be sorted cid_ - sid_, the next sid_ is sid_+1 so we don't need a lookup for the next student record
  
  data <- data %>% group_by(sid_) %>% mutate(student_obs = n()) %>% 
    mutate(ID_nextCollege=if_else(student_obs>1, c(ID_[2:length(ID_)], ID_[1]), ID_)) %>% select(-student_obs)
  stopifnot(!is.na(data$ID_nextCollege))
  
  
  
  
  # helper variables
  
  nSamples <- floor(niter/thin)
  stopifnot(burnin<=nSamples)

  
  
  
  # create the ranks of better and worse alternatives for students
  
  if (!is.null(student.rankvar)) {
    
    cat("Process students' rank variables ...\n")
	  # create consecutive rank variable, crank: how student ranks the College
	  # coding: NA - unknown; -1 - school is unacceptable; 1,2,3,... - ranks
    data$student.rankvar <- ifelse(data[[student.rankvar]]>=0, data[[student.rankvar]], NA)
    data <- data %>% group_by(sid_) %>% mutate(cRank_ = rank(student.rankvar, na.last='keep', ties.method='random'))
    data$cRank_[data[[student.rankvar]]==-1] <- -1 # recode -1 to -1
    
	  # find better and worse ranks
    data <- data %>% group_by(sid_) %>% mutate(cRank_max = ifelse(all(is.na(cRank_)), NA, max(cRank_, na.rm=T)))	# max rank among all ranked alternatives, NA if student did not rank anything
    data$cRank_better <- with(data, ifelse(cRank_==1 , NA,               # if top-ranked - there is no better school
                                    ifelse(cRank_==-1, cRank_max,        # else if unacceptable - last ranked school is better
                                    cRank_-1)))                          # else the next lower ranked school is better (or NA if cRank_ == NA)
    data$cRank_worse <- with(data, ifelse(cRank_==cRank_max, NA,         # if last ranked school - all unacceptable schools are worse, need to deal with that separately
                                   ifelse(cRank_==-1, NA,                # else if unacceptable - don't know which schools are worse
                                   cRank_+1)))                           # else - higher ranked school is better (or NA if unknown ranking)
    
    # get the R-style row indices (ID_) of better and worse alternatives
    # these are needed to bound valuations using the submitted rank order lists
    data <- left_join(x=data, y=data[data$cRank_>0,c('sid_','cRank_','ID_')], 
                      by=c('sid_'='sid_','cRank_better'='cRank_'),
                      suffix=c('','cBetter'))
    data <- left_join(x=data, y=data[data$cRank_>0,c('sid_','cRank_','ID_')], 
                  by=c('sid_'='sid_','cRank_worse'='cRank_'),
                  suffix=c('','cWorse'))
    
    # encode some more information in these indices
    # 0 -- NA, no better/worse is determined
    # >=1 -- ID_, R-style index of better / worse alternative
    # for ID_cWorse: 
    #   -1 indicates that obs is the last ranked, so lower bound is max of unacceptable IDs
    #   -2 indicates that the current option is unacceptable (so there is no lower bound)            
    data$ID_cBetter[is.na(data$ID_cBetter) | is.na(data$cRank_)] <- 0
    data$ID_cWorse[is.na(data$ID_cWorse) | is.na(data$cRank_)] <- 0
    data$ID_cWorse[data$cRank_==data$cRank_max & !is.na(data$cRank_max)] <- -1
    data$ID_cWorse[data$cRank_==-1] <- -2
    data$cRank_[is.na(data$cRank_)] <- 0 # encode NA as 0 too
    
    data$cRank_better <- data$cRank_worse <- data$student.rankvar <- NULL
    
  } else {
    cat("Students' rank variables not provided ...\n")
    data$cRank_ <- data$ID_cBetter <- data$ID_cWorse <- data$cRank_max <- 0
  }
  
  
  
  
  # create the ranks of better and worse alternatives for colleges
  
  if (!is.null(college.rankvar)) {
    
    cat("Process colleges' rank variables ...\n")
    # create consecutive rank variable
    data$college.rankvar <- ifelse(data[[college.rankvar]]>=0, data[[college.rankvar]], NA)
    data <- data %>% group_by(cid_) %>% mutate(sRank_ = rank(college.rankvar, na.last="keep", ties.method = "random"))
    data$college.rankvar <- NULL

    # find better and worse rank
    data$sRank_[data[[college.rankvar]]==-1] <- -1
    data <- data %>% group_by(cid_) %>% mutate(sRank_max = ifelse(all(is.na(sRank_)), NA, max(sRank_, na.rm=T)) )
    data$sRank_better <- with(data, ifelse(sRank_==1 , NA,               # if top-ranked - there is no better student
                                    ifelse(sRank_==-1, sRank_max,        # else if unacceptable - last ranked student is better
                                    sRank_-1)))                          # else the next lower ranked student is better (or NA if sRank_ == NA)
    data$sRank_worse <- with(data, ifelse(sRank_==sRank_max, NA,         # if last ranked student - all unacceptable students are worse, need to deal with that separately
                                   ifelse(sRank_==-1, NA,                # else if unacceptable - don't know which students are worse
                                   sRank_+1)))                           # else - higher ranked student is better (or NA if unknown ranking)

    # get the R-style row indices IDs of better and worse alternatives
    # these are needed to bound valuations using the submitted rank order lists
    data <- left_join(x=data, y=data[data$sRank_>0, c('cid_','sRank_','ID_')], 
                  by=c('cid_'='cid_','sRank_better'='sRank_'),
                  suffix=c('','sBetter'))
    data <- left_join(x=data, y=data[data$sRank_>0, c('cid_','sRank_','ID_')], 
                  by=c('cid_'='cid_','sRank_worse'='sRank_'),
                  suffix=c('','sWorse'))

    # encode some more information in these indices
    # 0 -- NA, no better/worse is determined
    # >=1 -- ID_, R-style index of better / worse alternative
    # for ID_sWorse: 
    #   -1 indicates that obs is the last ranked, so lower bound is max of unacceptable IDs
    #   -2 indicates that the current option is unacceptable (so there is no lower bound)            
    data$ID_sBetter[is.na(data$ID_sBetter) | is.na(data$sRank_)] <- 0
    data$ID_sWorse[is.na(data$ID_sWorse) | is.na(data$sRank_)] <- 0
    data$ID_sWorse[data$sRank_==data$sRank_max & !is.na(data$sRank_max)] <- -1
    data$ID_sWorse[data$sRank_==-1] <- -2
    data$sRank_[is.na(data$sRank_)] <- 0 # encode NA as 0 too
    
    data$sRank_better <- data$sRank_worse <- NULL
    
  } else {
    cat("Colleges' rank variables not provided ...\n")
    data$sRank_ <- data$ID_sBetter <- data$ID_sWorse <- data$sRank_max <- 0
  }
  
  
  
  
  # get equilibrium ID of assigned college for each student id
  
  StudentStats <- data.frame(sid_=unique(data$sid_))
  if (!is.null(match.id)) {
    cat('Process assignment information ...\n')
    StudentStats <- merge(StudentStats, data[data[[match.id]]==1, c('sid_', 'cid_', 'ID_')], by='sid_',all.x=T)
    # check that at most one college is matched to every student
    tmp2 <- as.numeric(table(StudentStats$sid_))
    stopifnot( all(tmp2<=1) )
    rm(tmp2)
    match.info <- T
  } else { # if match.id is unspecified, we cannot compute stability bounds
    StudentStats$cid_ <- -1
    StudentStats$ID_ <- NA
    data$D <- 0
    match.id <- 'D'
    match.info <- F
    cat("Assignment not provided.\n")
  }
  StudentStats$ID_ <- ifelse(is.na(StudentStats$ID_), 0, StudentStats$ID_) # ID_=0 denotes not matched
  StudentStats <- with(StudentStats, StudentStats[order(sid_),])
  
  
  
  # get number and ID_ of assigned students per college, amongst others
  
  SchoolStats <- data.frame(cid_=unique(data$cid_))
  StudentStats$nAssigned <- ifelse(StudentStats$ID_>0, 1, 0)
  tmp1 <- aggregate(nAssigned~cid_, data=StudentStats, FUN=sum)
  StudentStats$nAssigned <- NULL
  SchoolStats <- merge(SchoolStats, tmp1, by='cid_',all.x=T)
  SchoolStats$nAssigned <- with(SchoolStats, ifelse(is.na(nAssigned), 0, nAssigned))
  # create capacity vector for colleges if it is not supplied
  if (!is.null(nSeats)) SchoolStats$nSeats <- nSeats else SchoolStats$nSeats <- SchoolStats$nAssigned
  SchoolStats$VacantSeats <- with(SchoolStats, nSeats>nAssigned)
  if (!match.info) SchoolStats$VacantSeats <- T # if no match info exists, we assume that all schools are unconstrained, so no eq. bounds can be computed
  # tmp1 <- aggregate(ID_~cid_, data, FUN=min)
  # SchoolStats <- merge(SchoolStats, tmp1, by='cid_', suffix=c('','start'), all.x=T)
  # tmp1 <- aggregate(ID_~cid_, data, FUN=max)
  # SchoolStats <- merge(SchoolStats, tmp1, by='cid_', suffix=c('','end'), all.x=T)
  tmp1 <- data %>% group_by(cid_) %>% summarise(ones_=n())
  SchoolStats <- merge(SchoolStats, tmp1, by='cid_', suffix=c('','end'), all.x=T)
  rm(tmp1)
  
  
  
  
  # generate design matrices
  
  cat('Generate design matrices ...\n')
  # attention: Xc are variables that drive student's preferences over colleges (i.e., cRank_, ID_sBetter and ID_sWorse)
  data <- data[order(data$ID_),]
  Xc <- sparse.model.matrix(update(student.prefs,ID_~.), model.frame(~., data, na.action='na.pass')) # ID_ is a workaround incase the depvar is '.'
  data <- data[order(data$ID_),]
  Xs <- sparse.model.matrix(update(college.prefs,ID_~.), model.frame(~., data, na.action='na.pass'))
  students.nCoef <- ncol(Xc)
  college.nCoef <- ncol(Xs)
  # check that we have the same number of obs on both sides - if not true, check for missings on both sides
  stopifnot(nrow(Xs)==nrow(Xc))
  # generate indicator which data to use (Q: why would data not be used??)
  data$use <- F; data$use[match(rownames(Xs), rownames(data))] <- T
  gc()
  if (demean) { # if required, de-mean variables by decision making unit  
    cat('Demeaning design matrices. Could cause memory issues with many dummy variables ...\n')
    # It could cause trouble, because this transformation assigns non-zero values to virtually all
    # entries of the design matrix, so that the benefit of using sparse matrices is void, and
    # memory is likely to be a constraint when that matrix includes many dummy variables.
    # Xs: variables driving college preferences over students, so need to de-mean by college
    # for (i in 1:ncol(Xs))  Xs[,i] <- Xs[,i] - ave(Xs[,i], data$cid_[data$use], FUN=mean) # too slow
    groups <- sparse.model.matrix(~factor(cid_)-1, data[data$use,])
    Xs <- Xs - groups %*% ( t(groups) %*% Xs / colSums(groups) )
    # Xc: variables driving student prefs over colleges, so demean by sid_
    groups <- sparse.model.matrix(~factor(sid_)-1, data[data$use,])
    #for (i in 1:students.nCoef) Xc[,i] <- Xc[,i] - groups %*% ( t(groups) %*% Xc[,i] / colSums(groups) )
    Xc <- Xc - groups %*% ( t(groups) %*% Xc / colSums(groups) )
    rm(groups)
  }
  XsXsInv <- solve(t(Xs)%*%Xs) # dealing with sparse matrices in Cpp requries extra library, so we do it here
  XcXcInv <- solve(t(Xc)%*%Xc) # dealing with sparse matrices in Cpp requries extra library, so we do it here
  
  
  
  
  # initial parameter
  
  if (is.null(initparm)) {
    beta <- matrix(0, nrow=students.nCoef)
    gamma <- matrix(0, nrow=college.nCoef)
  } else {
    stopifnot(length(initparm$beta)==students.nCoef)
    stopifnot(length(initparm$gamma)==college.nCoef)
    beta <- matrix(initparm, ncol=1, nrow=students.nCoef)
    gamma <- matrix(initparm, ncol=1, nrow=college.nCoef)
  }
  
  
  
  # initialize valuations - this saves a lot of iterations
  cat('Initialize latent valuations ...\n')
  
  # if ranks are observed, we init them to the quantiles of the N(0,1) distribution
  data$tmp1 <- with(data, ifelse(cRank_>0, cRank_, ifelse(cRank_==-1, cRank_max+cid_, NA)))
  data <- data %>% group_by(sid_) %>% mutate(tmp2 = rank(tmp1, na.last='keep'))
  data <- data %>% group_by(sid_) %>% mutate(tmp3 = sum(!is.na(tmp1)))
  data$Vc_ <- with(data, ifelse(!is.na(tmp2), qnorm(1 - (2*tmp2 - 1) / (2*tmp3)), 0))
  
  data$tmp1 <- with(data, ifelse(sRank_>0, sRank_, ifelse(sRank_==-1, sRank_max+sid_, NA)))
  data <- data %>% group_by(cid_) %>% mutate(tmp2 = rank(tmp1, na.last='keep'))
  data <- data %>% group_by(cid_) %>% mutate(tmp3 = sum(!is.na(tmp1)))
  data$Vs_ <- with(data, ifelse(!is.na(tmp2), qnorm(1 - (2*tmp2 - 1) / (2*tmp3)), 0))
  data$tmp1 <- data$tmp2 <- data$tmp3 <- NULL
  
  
  # sort the data cid_ - sid_
  
  data <- data[order(data$ID_),]
  data <- as.data.frame(data)
  gc() # clean up garbage
  
  
  # run the estimator
  
  cat('Call the estimation routine ...\n')
  res = stabest_internal(Vc=data$Vc_,
                         Vs=data$Vs_,
                         Xc=Xc, 
                         XcXcInv=XcXcInv, 
                         Xs=Xs, 
                         XsXsInv=XsXsInv, 
                         betaR=beta, 
                         gammaR=gamma, 
                         VacantSeats = SchoolStats$VacantSeats,
                         nObs=SchoolStats$ones_, 
                         ID_eqCollege = StudentStats$ID_,
                         match=data[[match.id]],
                         sid=data$sid_,
                         ID_cBetter=data$ID_cBetter,
                         ID_cWorse=data$ID_cWorse,
                         ID_sBetter=data$ID_sBetter,
                         ID_sWorse=data$ID_sWorse,
                         ID_nextCollege = data$ID_nextCollege,
                         niter=niter,
                         thin=thin,
                         demean=demean)
  
  
  
  
  
  
  
  
  
  # get parameter estimates, compute statistics
  
  colnames(res$betadraws) <- colnames(Xc)
  colnames(res$gammadraws) <- colnames(Xs)
  if (students.nCoef>1) {
    beta.hat <- colMeans(res$betadraws[burnin:nSamples,])
    beta.vcov <- cov(res$betadraws[burnin:nSamples,])
    beta.ci95 <- apply(res$betadraws, 2, FUN=function(z) quantile(z[burnin:nSamples], c(.025,.975)))
  } else {
    beta.hat <- as.vector(mean(res$betadraws[burnin:nSamples,]))
    names(beta.hat) <- colnames(Xc)
    beta.vcov <- matrix(var(res$betadraws[burnin:nSamples,]), ncol=1, nrow=1)
    colnames(beta.vcov) <- colnames(Xc)
    rownames(beta.vcov) <- colnames(Xc)
    beta.ci95 <- matrix(quantile(res$betadraws[burnin:nSamples], c(.025,.975)), ncol=1)
    colnames(beta.ci95) <- colnames(Xc)
    rownames(beta.ci95) <- c('2.5%','97.5%')
  }
  if (college.nCoef>1) {
    gamma.hat <- colMeans(res$gammadraws[burnin:nSamples,])
    gamma.vcov <- cov(res$gammadraws[burnin:nSamples,])
    gamma.ci95 <- apply(res$gammadraws, 2, FUN=function(z) quantile(z[burnin:nSamples], c(.025,.975)))
  } else {
    gamma.hat <- as.vector(mean(res$gammadraws[burnin:nSamples,]))
    names(gamma.hat) <- colnames(Xs)
    gamma.vcov <- matrix(var(res$gammadraws[burnin:nSamples,]), ncol=1, nrow=1)
    colnames(gamma.vcov) <- colnames(Xs)
    rownames(gamma.vcov) <- colnames(Xs)
    gamma.ci95 <- matrix(quantile(res$gammadraws[burnin:nSamples], c(.025,.975)), ncol=1)
    colnames(gamma.ci95) <- colnames(Xs)
    rownames(gamma.ci95) <- c('2.5%','97.5%')
  }
  
  data$Vc_hat[data$use==T] <- res$Vc_hat
  data$Vs_hat[data$use==T] <- res$Vs_hat
  est <- list(beta.hat=beta.hat, beta.vcov=beta.vcov, beta.ci95=beta.ci95, betadraws=res$betadraws, 
              gamma.hat=gamma.hat, gamma.vcov=gamma.vcov, gamma.ci95=gamma.ci95, gammadraws=res$gammadraws,
              niter=niter, burnin=burnin, thin=thin,
              data=data[,c(student.id, college.id, match.id, 'Vc_','Vs_','Vc_hat','Vs_hat')])
  class(est) <- "stabest"
  # for debugging: 'sid_','cid_','ID_','sRank_','ID_sBetter','ID_sWorse','cRank_','ID_cBetter','ID_cWorse',

  return(est)
}

