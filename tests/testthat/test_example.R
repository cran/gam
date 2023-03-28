set.seed(20101)
data(kyphosis)
fit1=gam(Kyphosis ~ s(Age,4) + Number, family = binomial, data=kyphosis)
data(airquality)
fit2=gam(Ozone^(1/3) ~ lo(Solar.R) + lo(Wind, Temp), data=airquality, na=na.gam.replace)
fit3=gam(Kyphosis ~ poly(Age,2) + s(Start), data=kyphosis, family=binomial, subset=Number>2)
data(gam.data)
fit4=Gam.object <- gam(y ~ s(x,6) + z,data=gam.data)
sum1=summary(Gam.object)
an1=anova(Gam.object)
fit5 <- update(fit4, ~.-z)
an2=anova(fit4, fit5, test="Chisq")
fit6 <- gam(y~x+z, data=gam.data)
step1 <-step.Gam(Gam.object, scope=list("x"=~1+x+s(x,4)+s(x,6)+s(x,12),"z"=~1+z+s(z,4)))
data(gam.newdata)
pred1=predict(Gam.object,type="terms",newdata=gam.newdata)
if (getRversion() >= numeric_version("4.3.0")) {
    ## Fix to account for changes in family objects for R >= 4.3.0
    fit1$family$dispersion <- NULL
    fit2$family$dispersion <- NULL
    fit3$family$dispersion <- NULL
    fit4$family$dispersion <- NULL
    fit5$family$dispersion <- NULL
    fit6$family$dispersion <- NULL
    step1$family$dispersion <- NULL
}


objects  <- list(
    fit1=fit1,
    fit2=fit2,
    fit3=fit3,
    fit4=fit4,
    fit5=fit5,
    fit6=fit6,
    sum1=sum1,
    an1=an1,
    an2=an2,
    step1=step1)
##saveRDS(objects, "test_results/gam-1.20-results.RDS")

expected  <- readRDS("test_results/gam-1.20-results.RDS")
for (x in names(objects)) {
    cat(sprintf("Testing %s\n", x))
    expect_equal(objects[[x]], expected[[x]])
}


