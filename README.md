# stabest
This is an **R package** developed with [Thilo Klein](https://github.com/thiloklein), and described in our [joint paper](https://sites.google.com/view/robertaue/startseite#h.p_RZ2BD-FWtmTy) with [Josue Ortega](https://pure.qub.ac.uk/en/persons/josue-ortega). Its aim is to estimate students' preferences over schools, and vice versa, in a school choice market when rank order lists are submitted strategically. The method provides a way to estimate preferences under different identifying assumptions. It can be installed from Github (requires R build tools):
```r
library(devtools)

devtools::install_github("robertaue/stabest", build_vignettes = TRUE, upgrade=FALSE)
```
Then, you can use it to estimate students' and schools' preferences jointly, under the identifying assumptions of stability and undominated strategies:
```r
library(stabest)
fit <- stabest(choice_rk~-1 + distance + stu_score:sch_mscore + sch_FE,

               sch_rk_observed~stu_score-1,

               data=schoolmarket200, nSeats=schoolmarket200$sch_capacity[schoolmarket200$stu_id==1],

               student.id='stu_id',college.id='sch_id',match.id='sch_assignment',

               niter=20000, burnin=1000, thin=10

)
summary(fit)
```
For more details, read our joint paper, or see the package vignette: `vignette("stabest-vignette")`
