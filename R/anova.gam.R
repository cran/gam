"anova.Gam" <-
  function(object, ..., test = c("Chisq", "F", "Cp"))
{
  test=match.arg(test)
  margs <- function(...)
    nargs()
  if(margs(...))
    anova(structure(list(object, ...),class="glmlist"), test = test)
  else summary.Gam(object)$anova
}
