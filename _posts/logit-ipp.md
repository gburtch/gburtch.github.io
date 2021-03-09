---
title: Non-linear Regression and the Incidental Parameters Problem
output:
  md_document:
    variant: gfm
    preserve_yaml: TRUE
knit: (function(inputFile, encoding) {
  rmarkdown::render(inputFile, encoding = encoding, output_dir = "../_posts") })
author: "steve"
date: '2021-03-08'
excerpt: "Here is a simulation-based demonstration of the incidental parameters problem with logistic regression."
layout: post
permalink: /posts/2021/03/logit-ipp/
tags:
  - econometrics
  - logistic
  - endogeneity
---

This is the first in what I hope will be a series of basic blog posts in
which I’ll discuss basic issues in econometric analysis, and common
pitfalls. My plan is to demonstrate the issues intuitively, and
concretely, via simulation. Hopefully by providing this code, anyone who
comes across the material can download and play with the code themselves
to better understand the issues at play. All simulations will be
conducted in R.

## Incidental Parameters Problem

The incidental parameters problem (henceforth IPP) is a well-known
problem in econometrics, dating back more than 80 years. The seminal
paper introducing the problem was published by Jerzy Neyman and Betty
Scott in
[Econometrica](https://www.jstor.org/stable/1914288?seq=1 "Econometrica")
in 1948. I won’t rehash the details of the issue here, but I’ll give you
the punch line.

1.  With linear regression, we generally address time invariant
    confounds in panel data via fixed effects (either panel dummies, or
    the mathematically equivalent within-transformation).

2.  Unfortunately, those methods do not ‘work’ in non-linear regression,
    because the fixed effects enter the specification multiplicatively
    (the only exception here is Poisson regression, where the math still
    works out as a matter of convenience).

3.  Fortunately, Gary Chamberlain derived [conditional
    logit](https://www.jstor.org/stable/2297110?seq=1 "conditional logit"),
    wherein the panel fixed effects are instead addressed via a
    transformation to the likelihood function.

While I encourage anyone to go and read these papers, my goal here is to
demonstrate concretely, via simulation, just how poorly a dummy fixed
effect estimator performs in the case of panel logistic regression. I
will first simulate a simple scenario of static confounding, wherein a
standard logistic regression model fails to recover the parameter of
interest. Subsequently, I will compare its performance against
alternative correction procedures, namely i) a dummy fixed effect
estimator and ii) the Chamberlain conditional logit estimator. Notably a
more recent literature pursues bias correction methods for general
non-linear panel data models, e.g., Fernandez-Val (2009) in the [Journal
of
Econometrics](https://www.sciencedirect.com/science/article/pii/S0304407609000463 "Journal of Econometrics").
Accordingly, iii) I will also explore the efficacy of Fernandez-Val’s
bias correction following dummy fixed effect logit.

# Simulation Setup

The simulation setup is as follows. I begin with 25 panel unit and
simulate a time-invariant characteristic that serves as the confound. I
then combine the units into a data frame.

``` r
# I begin by setting some simulation parameters.
set.seed(12345)
n <- 100
sims <- 100
panel_max <- 15

# Lastly, for each panel unit
confound <- rnorm(n, mean = 0, sd=1)
id <- seq(n)
treatment_grp <- rbinom(n=n, size=1, prob=confound > median(confound))
df <- data.frame(confound=confound,id = id, treatment_grp = treatment_grp)
```

Next, in each run of the simulation, I create a new panel data set,
extending off the basic setup. I expand the original data-frame into
repeated observations of each panel, to varying lengths, *T* = 2..50.
With each new panel dataset, I turn the treatment ‘on’ for the treatment
group at the panel midpoint (rounded down).

After I create a dataset, I estimate four alternative logistic
regression models:

-   A naive regression (ignoring the confound)
-   A dummy fixed effect regression (the incidental parameters problem
    applies here)
-   A dummy fixed effect regression with analytical correction (using
    the ‘bife’ package)
-   A conditional logit regression (this is the idea approach)
-   An oracle regression (standard logit, conditioning on the confound)

``` r
# I create an empty data frame to store all the simulation results, based on the simulation paramters 
sim_result <- expand.grid(periods = seq(from=2,to=panel_max,by=1),run = seq(sims))

# Simulate each panel length 50 times and re-run all models. 
for (r in 1:sims){
  for (p in seq(from=2,to=panel_max,by=1)){
    
      # Simulate a panel of p periods.
      # Make the data with whole panel length
      df_tmp <- expandRows(df, count=p, count.is.col=FALSE)
      
      # Treatment begins midway through the panel for treated units.
      df_tmp <- df_tmp %>% group_by(id) %>% mutate(period = seq(from=1,to=p,by=1)) %>% mutate(Treat = (period > floor((p+1)/2))*treatment_grp) %>% arrange(treatment_grp,id,period) %>% select(-c(treatment_grp))
      df_tmp$Treat_center <- df_tmp$Treat - mean(df_tmp$Treat)
      df_tmp$xb <- 2*df_tmp$Treat_center + 4*df_tmp$confound
      df_tmp$p <- inv.logit(df_tmp$xb)
      df_tmp$Y <- rbinom(nrow(df_tmp),size=1, prob=df_tmp$p)
      
      # Estimate our four logit models.
      logit_mod_control <- glm(data=df_tmp,Y ~ Treat_center + confound, family = "binomial")
      logit_mod_omit <- glm(data=df_tmp,Y ~ Treat_center, family = "binomial")
      logit_mod_fe <- clogit(data=df_tmp,Y ~ Treat_center + strata(id))
      logit_mod_ipp <- glm(data=df_tmp,Y ~ Treat_center + factor(id), family = "binomial")
      
      # Sometimes the analytical correction fails, and it doesn't handle it well.
      # Note: usually when bife fails, clogit will also fail to converge and reports NA.
      tryCatch(
        {
          # Let's set a placeholder result to NA in case the estimation fails.
          logit_mod_bife <- data.frame(coefficients=NA)
          suppressWarnings(logit_mod_bife <- bias_corr(bife(Y~Treat_center | id,data=df_tmp)))
        },error=function(cond){
          # BIFE crapped out.
          return(logit_mod_bife)
        },warning=function(cond){
          # BIFE generated a warning.
          return(logit_mod_bife)
        }
      )
      
      sim_result_row <- data.frame(cbind(periods=p,run=r,b_control=logit_mod_control$coefficients[2],b_omit=logit_mod_omit$coefficients[2],b_fe=logit_mod_fe$coefficients[1],b_ipp=logit_mod_ipp$coefficients[2],b_bife=logit_mod_bife$coefficients[1]))
      
      sim_result <- sim_result %>% left_join(sim_result_row,by=c("periods","run"))
      sim_result$b_control.y[is.na(sim_result$b_control.y) ] <- sim_result$b_control.x[ is.na(sim_result$b_control.y)]
      sim_result$b_omit.y[is.na(sim_result$b_omit.y) ] <- sim_result$b_omit.x[ is.na(sim_result$b_omit.y)]
      sim_result$b_fe.y[is.na(sim_result$b_fe.y) ] <- sim_result$b_fe.x[ is.na(sim_result$b_fe.y)]
      sim_result$b_ipp.y[is.na(sim_result$b_ipp.y) ] <- sim_result$b_ipp.x[ is.na(sim_result$b_ipp.y)]
      sim_result$b_bife.y[is.na(sim_result$b_bife.y) ] <- sim_result$b_bife.x[ is.na(sim_result$b_bife.y)]
      
      if("b_control.y" %in% colnames(sim_result)){
        sim_result <- sim_result %>% select(-c(b_control.x, b_omit.x, b_fe.x, b_ipp.x,b_bife.x)) %>% rename(b_control=b_control.y, b_omit = b_omit.y, b_fe = b_fe.y, b_ipp = b_ipp.y, b_bife = b_bife.y)
      }
   }
}
```

# Regression Output

Let’s take a quick look at the outpout from the final set of
regressions, in the last iteration. Note that I am omitting the
intercept and dummy coefficients in the output.

``` r
screenreg(list(logit_mod_omit,logit_mod_ipp,logit_mod_bife,logit_mod_fe,logit_mod_control),digits=3,omit.coef=c("factor|Intercept|confound"),custom.model.names=c("logit","logit + dummies","bife","clogit","logit + control"))
```

    ## 
    ## ========================================================================================
    ##                 logit         logit + dummies  bife         clogit       logit + control
    ## ----------------------------------------------------------------------------------------
    ## Treat_center       4.684 ***     1.791 **         1.683 **     1.685 **     1.627 **    
    ##                   (0.506)       (0.581)          (0.565)      (0.567)      (0.536)      
    ## ----------------------------------------------------------------------------------------
    ## AIC             1627.689       715.132                       413.710      607.953       
    ## BIC             1638.316      1251.768                                    623.892       
    ## Log Likelihood  -811.845      -256.566         -256.584                  -300.976       
    ## Deviance        1623.689       513.132          513.168                   601.953       
    ## Num. obs.       1500          1500              615         1500         1500           
    ## R^2                                                            0.008                    
    ## Max. R^2                                                       0.246                    
    ## Num. events                                                  857                        
    ## Missings                                                       0                        
    ## ========================================================================================
    ## *** p < 0.001; ** p < 0.01; * p < 0.05

# Visualizing the Distribution of Estimates

Okay, but how consistently did each alternative estimator perform,
across the 50 runs? Let’s visualize this. One thing to note here is that
sometimes the estimators fail to converge, and thus spit out
**woefully** inaccurate estimates (off by orders of magnitude). So, we
will first look at the general plots, so you can see how often these
happen across the different estimators. Then, we will zoom in on the
“valid” estimates!

![](/Users/gburtch/Documents/GitHub/gburtch.github.io/_posts/logit-ipp_files/figure-gfm/plots-1.png)<!-- -->

Now we will zoom in around 2, to the estimates that were generated by
models that converged.

![](/Users/gburtch/Documents/GitHub/gburtch.github.io/_posts/logit-ipp_files/figure-gfm/plots_trunc-1.png)<!-- -->
