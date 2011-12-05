"anova.gam" <-
  function(object, ..., test = c("Chisq", "F", "Cp"))
{
  test=match.arg(test)
  margs <- function(...)
    nargs()
  if(margs(...))
    anova.glmlist(list(object, ...), test = test)
  else summary.gam(object)$anova
}
