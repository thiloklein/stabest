# load raw data, generated using the DGP from Fack et al. (2019)'s Appendix

fpath <- paste(Sys.getenv('dropbox'), '\\8.IntegrationBudapest\\MCsims\\data\\SIM_DATA_200students_normal.csv', sep='')

# select first sample for a trial
dta <- read.table(fpath, stringsAsFactors = F, header=T, sep=',')
dat <- dta[dta$Sample_id==1,]
dat <- dat[,c('stu_id','sch_id','sch_assignment',
              'true_rk','choice_rk','sch_capacity','stu_priority',
              'stu_score','sch_mscore','distance','sch_FE')]

# create true and observed school's ranking over students
dat$sch_rk_true <- with(dat, ave(-stu_priority, sch_id, FUN=rank)) # lower rank = better
dat$stu_priority <- NULL
dat$sch_rk_observed <- with(dat, ifelse(choice_rk>0, sch_rk_true, NA))
dat$sch_rk_observed <- with(dat, ave(sch_rk_observed, sch_id, FUN = function(z) rank(z, na.last="keep")))

# 0 means NA here
dat$choice_rk <- with(dat, ifelse(choice_rk>0, choice_rk, NA))

# WTT: unranked are unacceptable, so encode NA as -1
dat$choice_rk.WTT <- with(dat, ifelse(is.na(choice_rk), -1, choice_rk))

table(dat$true_rk, dat$choice_rk)
schoolmarket200 <- dat
#setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
save(schoolmarket200, file='data/schoolmarket200.RData')
