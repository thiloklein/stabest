#' Simulated data from a market with 200 students and six schools
#' 
#' This dataset was generate based on the procedure (and code) described in
#' the appendix to Fack et al. (2019). That code simulates a market in which
#' students can choose to rank any market, but face a small application cost.
#' Thus, students face an optimal portoflio choice problem that depends
#' on their expectations about other students' behvaiour and priorities.
#' 
#' We changed the code slightly along two dimensions:
#' 
#' 1. we used normally distributed errors instead of T1EV distributed errors,
#' in line with our empirical model, and
#' 2. we also made the schools' choices depend on an unobservable component
#' 
#' @format A data frame with 1200 rows and 11 variables. Each row refers to the 
#' combination of one student (stu_id) and one school (sch_id).
#' \describe{
#'     \item{stu_id}{Student ID}
#'     \item{sch_id}{School ID}
#'     \item{sch_assignment}{=1 if student is assigned to school}
#'     \item{true_rk}{Student's true ranking of the school (lower=more preferred)}
#'     \item{choice_rk}{Student's optimal ranking in the Bayesian equilibrium, thus subject to the 'skipping problem' and 'truncation problem'}
#'     \item{sch_capacity}{Capacity at the school}
#'     \item{stu_score}{Student's test global score}
#'     \item{sch_mscore}{Average of all students' test scores at that school}
#'     \item{distance}{Geographical distance between the school and the student}
#'     \item{sch_FE}{School quality}
#'     \item{choice_rk.WTT}{Same as choice_rk, but unranked are coded -1 (unacceptable)}
#'     \item{sch_rk_observed}{Priority rank of student at the school, considering only observable applications}
#'     \item{sch_rk_true}{True priority rank of student at the school, across all students}
#' }
#' 
#' @source{https://www.aeaweb.org/articles?id=10.1257/aer.20151422}
"schoolmarket200"
