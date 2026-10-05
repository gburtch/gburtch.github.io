---
title: "Why Your Logit (and Cox) Coefficients Move When You Add Controls, Even Without Confounding"
date: '2026-09-08'
permalink: /posts/2026/09/noncollapsibility/
excerpt: "Odds ratios and hazard ratios are non-collapsible: adding a covariate that predicts the outcome changes the treatment coefficient even when the covariate is unrelated to treatment. So coefficient stability across specifications tells you nothing about confounding. Marginal effects do not have this problem."
tags:
  - econometrics
  - logistic
  - survival
  - marginal effects
  - simulation
---

There is a ritual in empirical papers that goes like this. Estimate the effect of $D$ on $Y$. Add some controls. Add some more. If the coefficient on $D$ barely moves, declare it “robust” and move on; if it moves a lot, worry about confounding. [Oster (2019)](https://doi.org/10.1080/07350015.2016.1227711) and Altonji, Elder and Taber (2005) formalised a version of this logic for linear models, and it is a reasonable heuristic there.

It is not a reasonable heuristic for logit, probit, or Cox proportional hazards models, and this post is a simulation to show why. The short version: the coefficients from those models are **non-collapsible**. Add a covariate that predicts the outcome but has *nothing to do with treatment*, and the treatment coefficient changes anyway, always in the direction of larger magnitude. So coefficient movement is not evidence of confounding, and coefficient stability is not evidence of its absence. If you want a quantity you can compare across specifications, look at marginal effects (risk differences, average marginal effects, differences in survival probabilities), which are collapsible.

None of this is new. Gail, Wieand and Piantadosi ([1984](https://doi.org/10.1093/biomet/71.3.431)) showed it for randomized experiments; Greenland, Robins and Pearl ([1999](https://doi.org/10.1214/ss/1009211805)) untangled collapsibility from confounding; Mood ([2010](https://doi.org/10.1093/esr/jcp006)) wrote the paper every social scientist using logit should read; Hernán ([2010](https://doi.org/10.1097/EDE.0b013e3181c1ea43)) made the case for survival models. But it is still routinely ignored, so here is the phenomenon in a form you can run yourself.

## The setup

I simulate a *randomized* binary treatment $D$, so confounding is impossible by construction, and one covariate $X \sim N(0,1)$ that is independent of $D$. Outcomes are generated from a single latent index

$$z = -0.5 + 1 \cdot D + \gamma X,$$

and I vary $\gamma$, the strength of the covariate, from 0 to 3. From that index I generate four outcomes: a continuous one ($z$ plus normal noise), a binary one through a logit link, a binary one through a probit link, and an event time with hazard $\exp(z)$ (so proportional hazards hold *conditional on* $X$), randomly censored. For each outcome I fit the matching model twice, with and without $X$, and record two things: the coefficient on $D$, and a marginal effect of $D$.

The marginal effect is the same object in every case: the average difference in the predicted outcome when everyone is set to $D = 1$ versus $D = 0$ (g-computation, or “standardization”). For the binary outcomes that is the average risk difference; for the survival outcome it is the difference in the probability of surviving past $t = 1$; for OLS it is just the coefficient. `marginaleffects::avg_comparisons()` computes all of these with proper standard errors, but they are two lines of `predict()` each, and I want the mechanics visible.

``` r
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(survival) })
set.seed(20260908)

sim_data <- function(n = 50000, gamma = 1, confounded = FALSE, beta_d = 1, alpha = -0.5) {
  X <- rnorm(n)
  D <- if (confounded) rbinom(n, 1, plogis(1.2 * X)) else rbinom(n, 1, 0.5)
  z <- alpha + beta_d * D + gamma * X
  t_event <- rexp(n, rate = exp(z)); t_cens <- rexp(n, rate = 0.3)
  data.frame(D, X,
             y_ols    = z + rnorm(n),
             y_logit  = rbinom(n, 1, plogis(z)),
             y_probit = rbinom(n, 1, pnorm(z)),
             time = pmin(t_event, t_cens), event = as.integer(t_event <= t_cens))
}

# Average marginal effect of D by g-computation: E[Y | do(D=1)] - E[Y | do(D=0)].
ame <- function(m, d) {
  d1 <- d; d1$D <- 1; d0 <- d; d0$D <- 0
  mean(predict(m, d1, type = "response")) - mean(predict(m, d0, type = "response"))
}
# Same idea for a Cox model: difference in standardized survival probability at t0.
surv_diff <- function(m, d, t0 = 1) {
  H0 <- basehaz(m, centered = FALSE)
  H0t <- approx(H0$time, H0$hazard, xout = t0, rule = 2)$y
  d1 <- d; d1$D <- 1; d0 <- d; d0$D <- 0
  lp1 <- predict(m, d1, type = "lp", reference = "zero")
  lp0 <- predict(m, d0, type = "lp", reference = "zero")
  mean(exp(-H0t * exp(lp1))) - mean(exp(-H0t * exp(lp0)))
}

fit_all <- function(d, with_X) {
  rhs <- if (with_X) "D + X" else "D"
  f <- function(lhs) as.formula(paste(lhs, "~", rhs))
  m <- list(OLS    = lm(f("y_ols"), d),
            Logit  = glm(f("y_logit"), d, family = binomial()),
            Probit = glm(f("y_probit"), d, family = binomial(link = "probit")),
            Cox    = coxph(f("Surv(time, event)"), d))
  data.frame(model = names(m),
             coef  = sapply(m, function(x) unname(coef(x)["D"])),
             me    = c(unname(coef(m$OLS)["D"]), ame(m$Logit, d), ame(m$Probit, d), surv_diff(m$Cox, d)),
             spec  = if (with_X) "With X" else "Without X")
}

# Model-free benchmarks for the marginal effects (D is randomized, so these are unbiased).
model_free <- function(d) {
  km <- survfit(Surv(time, event) ~ D, data = d)
  s1 <- summary(km, times = 1)
  data.frame(model = c("OLS", "Logit", "Probit", "Cox"),
             me = c(diff(tapply(d$y_ols, d$D, mean)), diff(tapply(d$y_logit, d$D, mean)),
                    diff(tapply(d$y_probit, d$D, mean)), diff(s1$surv)),
             spec = "Model-free (diff. in means / Kaplan-Meier)")
}
```

## Simulation 1: an independent covariate of increasing strength

``` r
gammas <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
res1 <- bind_rows(lapply(gammas, function(g) {
  d <- sim_data(gamma = g)
  bind_rows(fit_all(d, TRUE), fit_all(d, FALSE), model_free(d)) %>% mutate(gamma = g)
})) %>% mutate(model = factor(model, levels = c("OLS", "Logit", "Probit", "Cox")))
```

First, the coefficients on $D$. The true conditional coefficient is 1 in every model, and every model recovers it when $X$ is included.

![](/images/noncollapsibility_files/figure-gfm/sim1_coef-1.png)<!-- -->

For OLS, dropping $X$ changes nothing, because $X$ is independent of $D$: omitted variable bias requires the omitted variable to be correlated with the regressor. For logit, probit and Cox, dropping $X$ pulls the coefficient toward zero, and the more $X$ matters for the outcome, the further it is pulled. At $\gamma = 3$ the logit coefficient without $X$ is 0.45, versus 1.00 with it. Nothing about confounding changed between those two regressions; $D$ was a coin flip both times.

The mechanism is that the coefficients in these models are scaled by the residual variation in the index. In a probit, the estimated coefficient on $D$ when you omit $X$ is approximately $\beta_D / \sqrt{1 + \gamma^2}$, because the omitted $\gamma X$ gets absorbed into an error term whose standard deviation is no longer 1. Logit behaves the same way with a different constant. Mood (2010) calls this the “unobserved heterogeneity” problem; it is the same phenomenon as non-collapsibility of the odds ratio. The Cox model has an additional wrinkle: once you marginalise over $X$, the hazards are no longer proportional, and the coefficient you get is a weighted average of a hazard ratio that shrinks toward 1 over time as the high-$X$ subjects are selected out of the risk set (Struthers and Kalbfleisch, [1986](https://doi.org/10.1093/biomet/73.2.363); Hernán, 2010).

Now the marginal effects from the same fitted models.

![](/images/noncollapsibility_files/figure-gfm/sim1_me-1.png)<!-- -->

The marginal effects agree with each other and with the model-free benchmark at every $\gamma$, whether or not $X$ is in the model. Note that the marginal effect *itself* changes with $\gamma$: a stronger independent covariate spreads the latent index out, which flattens the link function, which makes a one-unit shift in $D$ move fewer people across the threshold. That is a real feature of the population, not an artefact of the specification. The coefficient conflates that population fact with the model’s scaling.

## Simulation 2: a specification search

Here is what the ritual looks like when you actually run it on a logit. Same randomized $D$, now with eight independent covariates $X_1, \dots, X_8$ with coefficients from 1.5 down to 0.2. I fit 200 logits, each with a random subset of the controls, and record the coefficient on $D$ and the average risk difference.

``` r
n <- 20000; K <- 8; gam <- c(1.5, 1.2, 1.0, 0.8, 0.6, 0.4, 0.3, 0.2)
Xmat <- matrix(rnorm(n * K), n, K, dimnames = list(NULL, paste0("X", 1:K)))
D <- rbinom(n, 1, 0.5)
z <- -0.5 + 1 * D + drop(Xmat %*% gam)
d2 <- data.frame(y = rbinom(n, 1, plogis(z)), D = D, Xmat)

subsets <- c(list(character(0), colnames(Xmat)),
             replicate(198, sample(colnames(Xmat), sample(0:K, 1)), simplify = FALSE))
res2 <- bind_rows(lapply(subsets, function(s) {
  m <- glm(reformulate(c("D", s), "y"), d2, family = binomial())
  data.frame(n_controls = length(s), coef = unname(coef(m)["D"]), me = ame(m, d2))
}))
```

![](/images/noncollapsibility_files/figure-gfm/sim2_plot-1.png)<!-- -->

The coefficient ranges from 0.55 to 1.01 across specifications that differ only in which irrelevant-to-$D$ controls are included. It climbs steadily with the explanatory power of the control set, exactly as the theory says. If you had been told this was an observational study, you would have read that as a treatment whose effect “grows once you account for” the controls, and you would have been wrong. The risk difference, meanwhile, sits between 0.131 and 0.139.

## Simulation 3: what actual confounding looks like

To be clear about what the marginal effect buys you, here is the contrast. Same design with $\gamma = 1.5$, but now in one version $D$ depends on $X$ (probability of treatment $= \Lambda(1.2X)$), so $X$ is a genuine confounder.

``` r
res3 <- bind_rows(lapply(c(FALSE, TRUE), function(cf) {
  d <- sim_data(gamma = 1.5, confounded = cf)
  bind_rows(fit_all(d, TRUE), fit_all(d, FALSE)) %>%
    filter(model == "Logit") %>%
    mutate(design = if (cf) "D depends on X (confounded)" else "D randomized")
}))
```

| Design | Coef, with X | Coef, without X | Risk diff., with X | Risk diff., without X |
|:---|:---|:---|:---|:---|
| D randomized | 1.01 | 0.68 | 0.176 | 0.169 |
| D depends on X (confounded) | 1.02 | 1.82 | 0.180 | 0.426 |

Logit, gamma = 1.5, n = 50,000. True conditional coefficient on D is 1.

In both rows the coefficient moves when $X$ is added. Only in the second row does the risk difference move. The marginal effect changing under adjustment is informative about confounding; the coefficient changing is not, because it changes either way.

## What to do

1.  **Report marginal effects for logit, probit, Cox and friends, and compare *those* across specifications.** Average risk differences, average marginal effects, standardized survival differences at a horizon of interest, restricted mean survival time. These are collapsible in the sense that matters here: adding a covariate that is independent of treatment leaves them, in expectation, where they were.
2.  **Do not read coefficient movement in a non-linear model as evidence about confounding**, in either direction. The Oster-style coefficient-stability logic is derived for linear models; applying it to logit coefficients mixes confounding with rescaling. Karlson, Holm and Breen (2012) propose a decomposition that separates the two if you really need to compare logit coefficients across nested models, but the marginal effect route is simpler and answers the question people usually mean to ask.
3.  **Likewise for hazard ratios.** A “marginal” hazard ratio from a Cox model without covariates is not the average of the conditional ones, is not constant over time even when the conditional one is, and depends on the follow-up window. If the question is causal, contrast survival probabilities or restricted mean survival times.
4.  **If you want an odds ratio or hazard ratio that is stable across control sets, target the marginal one directly**, for instance with inverse probability of treatment weighting rather than conditioning. That gets you a fixed estimand (the marginal odds ratio), at the cost of the coefficient no longer being comparable to a conditional one from another paper. Pick your estimand, then pick your specification; not the reverse.

*References.* Altonji, Elder and Taber (2005), *JPE*; Daniel, Zhang and Farewell (2021), *Biometrical Journal*; Gail, Wieand and Piantadosi (1984), *Biometrika*; Greenland, Robins and Pearl (1999), *Statistical Science*; Hernán (2010), *Epidemiology*; Karlson, Holm and Breen (2012), *Sociological Methodology*; Martinussen, Vansteelandt and Andersen (2020), *Lifetime Data Analysis*; Mood (2010), *European Sociological Review*; Oster (2019), *JBES*; Struthers and Kalbfleisch (1986), *Biometrika*.
