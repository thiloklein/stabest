library(stabest)

niter <- 100
thin <- 2
burnin <- 5
dat <- schoolmarket200
nSeats <- dat$sch_capacity[dat$stu_id==1]


#--------------------------------------------------------------------------#
#--- STABILTIY -------------------------------------------------------------
#--------------------------------------------------------------------------#

fit1 <- stabest(.~-1 + distance + stu_score:sch_mscore + sch_FE, .~stu_score-1,
                data=dat, nSeats=nSeats, 
                student.id='stu_id',college.id='sch_id',match.id='sch_assignment',
                niter=niter, burnin=burnin, thin=thin # burnin is with respect to the thinned series of estimates
)
summary(fit1)


#--------------------------------------------------------------------------#
#--- STABILITY + UNDOM (full pref matrices) ----------------------------------
#--------------------------------------------------------------------------#

fit2 <- stabest(true_rk~-1 + distance + stu_score:sch_mscore + sch_FE, sch_rk_true~stu_score-1,
                data=dat, nSeats=nSeats, 
                student.id='stu_id',college.id='sch_id',match.id='sch_assignment',
                niter=niter, burnin=burnin, thin=thin # burnin is with respect to the thinned series of estimates
)
summary(fit2)


#--------------------------------------------------------------------------#
#--- STABILITY+UNDOM w/ partial pref matrices  -----------------------------
#--------------------------------------------------------------------------#

fit3 <- stabest(choice_rk~-1 + distance + stu_score:sch_mscore + sch_FE, sch_rk_observed~stu_score-1,
                data=dat, nSeats=nSeats, 
                student.id='stu_id',college.id='sch_id',match.id='sch_assignment',
                niter=niter, burnin=burnin, thin=thin # burnin is with respect to the thinned series of estimates
)
summary(fit3)


#--------------------------------------------------------------------------#
#--- UNDOM w/ full pref matrices -------------------------------------------
#--------------------------------------------------------------------------#


fit4 <- stabest(true_rk~-1 + distance + stu_score:sch_mscore + sch_FE, sch_rk_true~stu_score-1,
                data=dat, nSeats=nSeats,
                student.id='stu_id',college.id='sch_id',match.id=NULL, 
                niter=niter, burnin=burnin, thin=thin # burnin is with respect to the thinned series of estimates
)
summary(fit4)


#--------------------------------------------------------------------------#
#--- UNDOM w/ partial pref matrices  ---------------------------------------
#--------------------------------------------------------------------------#

fit5 <- stabest(choice_rk~-1 + distance + stu_score:sch_mscore + sch_FE, sch_rk_observed~stu_score-1,
                data=dat, nSeats=nSeats, 
                student.id='stu_id',college.id='sch_id',match.id=NULL,
                niter=niter, burnin=burnin, thin=thin # burnin is with respect to the thinned series of estimates
)
summary(fit5)


#--------------------------------------------------------------------------#
#--- Weak Truth Telling (WTT) ----------------------------------------------
#--------------------------------------------------------------------------#

fit6 <- stabest(choice_rk.WTT~-1 + distance + stu_score:sch_mscore + sch_FE, sch_rk_observed~stu_score-1,
                data=dat, nSeats=nSeats, 
                student.id='stu_id',college.id='sch_id',match.id=NULL,
                niter=niter, burnin=burnin, thin=thin # burnin is with respect to the thinned series of estimates
)
summary(fit6)



#--------------------------------------------------------------------------#
#--- combine and compare all results ---------------------------------------
#--------------------------------------------------------------------------#

res_ <- rbind(c(fit1$beta.hat, fit1$gamma.hat),
              c(fit2$beta.hat, fit3$gamma.hat),
              c(fit3$beta.hat, fit3$gamma.hat),
              c(fit4$beta.hat, fit4$gamma.hat),
              c(fit5$beta.hat, fit5$gamma.hat),
              c(fit6$beta.hat, fit6$gamma.hat),
              c(-1, 1, 3, 1)) # true parameters
res_ <- as.data.frame(round(res_,2))
res_$method <- c('stability',
                 'stability + UNDOM (full prefs)',
                 'stability + UNDOM (short prefs)',
                 'UNDOM (full prefs)',
                 'UNDOM (short prefs)',
                 'WTT',
                 'true')
res_ <- res_[,c(5,1:4)]
res_

