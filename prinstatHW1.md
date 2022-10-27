---
title: 'Prinstat: Home work 1'
author: "Leah Jimmiy Shitindi, Thi Lan Anh Tran"
output:
  pdf_document:
    number_sections: yes
  html_document:
    df_print: paged
editor_options: 
  markdown: 
    wrap: 72
---

# Question 1

The distribution of the incubation times for both Alpha and Delta
variants are not normally distributed. The means of the distribution for
all variants are greater than their respective medians. The box plots
(figure1) and the histograms(figure 2) show the positive skewness for
both variants which confirm that the incubation time is not normally
distributed. The incubation period is between one day to 16 days for the
alpha variant and 1day to 12 days for the delta variant while the median
median is 4 days for both variant again suggesting the right skewness of
the Covid-19 Alpha and Delta variants.

```{r prep, echo=FALSE, include=FALSE}
library(tidyverse)
library(kableExtra)
library(stats) 
library(gridExtra)
load("Incubations.RData")
```

```{r data explorarion1, echo=FALSE, out.width = "60%"}
# sample size for both variants
N <- table(Incubations$Variant)

summary_statistics <- Incubations %>% 
  group_by(Variant) %>% summarise(Means=mean(Time),
                                  Median=median(Time),
                                  Minimum=min(Time),
                                  Maximum = max(Time),
                                  Std = sd(Time))

summary_statistics %>%
  kbl(caption = "Summary Statistics of the Incubation Period for Alpha and Delta Variants  ") %>%
  kable_classic(full_width = F, html_font = "Cambria")

g1 <- Incubations %>% ggplot(mapping=aes(x= Variant,y=Time))+
  geom_boxplot()+
  labs(title=paste("Boxplots of Incubation Time for Alpha and Delta Covid-19 Variant"),tag="Figure 1: Box plot")+
  theme_bw() +
  theme(plot.title = element_text(size = 7, face = "bold")) 


g2 <- Incubations %>% ggplot()+
  geom_bar(aes(x=Time),stat="count")+
 geom_density(aes(x=Time,after_stat(count)))+
  facet_wrap(~Variant)+ 
  theme_classic()
grid.arrange(g1, g2, nrow=1)
```

# Question 2

Here we have the distributions of Alpha and Delta variants are gamma
distributed with the shape parameter $\alpha$ and the rate parameter
$\beta$. The density function for gamma distribution is as follows. $$
\begin{equation}
\tag{1}
f(x)=\frac{\beta^\alpha}{\Gamma(\alpha)} x^{\alpha-1} \exp (-\beta x)
\end{equation}
$$

The likelihood function can then be derived from equation (1) as
follows: $$
\begin{equation}
\tag{2}
L(\alpha, \beta)=\prod_{i=1}^n \frac{\beta^\alpha x_i^{\alpha-1}}{\Gamma(\alpha)} \exp (-\beta x_i)
\end{equation}
$$ and hence the log-likelihood for for the rate parameter $\beta$ is
given by: $$
\log (L(\alpha, \beta))=\log \left [\prod_{i=1}^n \frac{\beta^\alpha x_i^{\alpha-1}}{\Gamma(\alpha)} \exp (-\beta x_i)\right]\\
$$ 
$$\begin{equation}
\tag{3}
\ = n \alpha \log \beta+(\alpha-1) \log \Sigma x_i-\beta \Sigma x_i-n \log (\Gamma(\alpha))
\end{equation} $$ 


From equation 3, we derive the log-likelihood
functions of Alpha variant ($n_Alpha$ = 193 and $\alpha_Alpha$ =3.08)
and delta variant ($n_Delta$ =89 and $\alpha_Delta$ =4.43) as follows:

$$
\log (L(\beta\__Alpha))=594.44\log \beta+2.08\log \Sigma x_i-\beta \Sigma x_i-193 \log (\Gamma(3.08))\\
\log (L(\beta\__Delta))=394.27\log \beta+3.43\log \Sigma x_i-\beta \Sigma x_i-89 \log (\Gamma(4.43))
$$

# Question 3

Using the log-likelihood function obtained from equation 3, we plotted
the log-likelihood functions for Alpha and Delta variants with the range
of the $\beta$ values from 0 to 5, with the increment of 0.1.

```{r logLH, echo=FALSE, out.width = "50%"}
shape_Alpha <- 3.08
n_Alpha <- 193
beta <- seq(0,5, 0.1)
x_a <- Incubations$Time[Incubations$Variant=="Alpha"]
loglikelihood_Alpha <- n_Alpha*(shape_Alpha)*log(beta)-beta*sum(x_a)+
  (1-shape_Alpha)*sum(log(x_a)) - n_Alpha*log(gamma(shape_Alpha))

# Delta 
shape_Delta <- 4.43
n_Delta <- 89
beta <- seq(0,5, 0.1)
x_D<- Incubations$Time[Incubations$Variant=="Delta"]
loglikelihood_Delta <- n_Delta*(shape_Delta)*log(beta)-beta*sum(x_D)+
  (1-shape_Delta)*sum(log(x_D)) - n_Delta*log(gamma(shape_Delta))

par(mfrow=c(1,2))
plot(x=beta,y=loglikelihood_Alpha, 
     main = expression(paste("Log-likelihood function of the rate parameter ", (beta), " for Alpha variant")))
 plot(x=beta,y=loglikelihood_Delta, 
     main = expression(paste("Log-likelihood function of the rate parameter ", (beta), " for Delta variant")))

```

We can see that the optimal value for $\beta_Apha$ and $\beta_Delta$ (
which maximizes the log-likelihood function) is between 0.5-0.7 and
0.9-1, respectively. The exact value we will be calculated in the next
question.

# Question 4

The formula for the maximum likelihood estimate parameter $\beta_Apha$
and $\beta_Delta$, knowing that $\alpha_Alpha$ = 3.08, and
$\alpha_Delta$ = 4.43 is derived by maximizing the first order
differential equation of equation 3 as follows: $$
\log (L(\alpha, \beta))=n \alpha \log \beta+(\alpha-1) \log \Sigma x_i-\beta \Sigma x_i-n \log (\Gamma(\alpha) \\
\frac{\partial log(L(\alpha, \beta))}{\partial\beta} = \frac {n\alpha}{\beta} - \sum x_i \\
0 = \frac {n\alpha}{\beta} - \sum x_i \\
\sum x_i = \frac {n\alpha}{\beta} \\
\bar{x} = \frac {\alpha}{\beta}\\
\beta = \frac {\alpha}{\bar{x}}\\
Hence, 
\beta_Alpha = \frac {3.08}{\bar{x}_Alpha}\\
\beta_Delta = \frac {4.43}{\bar{x}_Delta}
$$ It can be show that the maximum likelihood of the rate parameter is
in fact the ratio of the shape parameter and the sample mean.

# Question 5

```{r optim, echo = FALSE}
lh_func_a <- function(beta, data=Incubations){
shape_Alpha = 3.08
n_Alpha <- 193
x_a <-data$Time[data$Variant=="Alpha"]
loglikelihood_Alpha <- n_Alpha*(shape_Alpha)*log(beta)-beta*sum(x_a)+
  (1-shape_Alpha)*sum(log(x_a)) - n_Alpha*log(gamma(shape_Alpha))
}
beta.estimation_a <-optimize(lh_func_a, c(0,1),maximum = TRUE)
beta.max_alpha <- beta.estimation_a$maximum

lh_func_d <- function(beta, data=Incubations){
shape_Alpha = 4.43
n_Alpha <- 89
x_a <-data$Time[data$Variant=="Delta"]
loglikelihood_Alpha <- n_Alpha*(shape_Alpha)*log(beta)-beta*sum(x_a)+
  (1-shape_Alpha)*sum(log(x_a)) - n_Alpha*log(gamma(shape_Alpha))
}
beta.estimation_d <-optimize(lh_func_d, c(0,1),maximum = TRUE)
beta.max_delta <- beta.estimation_d$maximum
b<-mean(Incubations$Time[Incubations$Variant=='Alpha'])

```

Using the formula derived in question 4, maximum likelihood estimates
for $\beta_Apha$ and $\beta_Delta$ are
$$ \beta_Alpha = `r 3.08/mean(Incubations$Time[Incubations$Variant=='Alpha'])`\\
\beta_Delta =`r 4.43/mean(Incubations$Time[Incubations$Variant!='Alpha'])`$$

Using the function optimize in R. The arguments are the mle function we
constructed in the previous question together with the interval which we
set at [0, 1]. The optimal value for $\beta_Alpha$ is `r beta.max_alpha`
and for $\beta_Delta$ is `r beta.max_delta`. Both methods give the same
estimates for the shape parameters $\beta_Alpha$ and $\beta_Delta$ of
Alpha and Delta variants.

# Question 6

```{r CI beta, echo=FALSE, out.width = "50%"}
alpha <- 0.05
c.alpha <- qchisq(1-alpha,df=1)

#max LH for beta of delta variant
max.loglik<- lh_func_d(beta.estimation_d$maximum,data=Incubations)

alpha.optim_d<-function(a,data) {
  (max.loglik-0.5*c.alpha-(lh_func_d(a,data=data)))^2
}


o.up_d<-optimize(alpha.optim_d,
               data=Incubations,
               upper=15,
               lower=(beta.estimation_d$maximum)*0.999)

o.low_d<-optimize(alpha.optim_d,
                data=Incubations,
                lower=0,
                upper=(beta.estimation_d$maximum)*1.1)


gamma.seq<-seq(0,2,0.01)
ll_d<- -sapply(gamma.seq,FUN=lh_func_d,data=Incubations)

#max LH for beta alpha variant
max.loglik<- lh_func_a(beta.estimation_a$maximum,data=Incubations)

alpha.optim_a<-function(a,data) {
  (max.loglik-0.5*c.alpha-(lh_func_a(a,data=data)))^2
}


o.up_a<-optimize(alpha.optim_a,
               data=Incubations,
               upper=15,
               lower=(beta.estimation_a$maximum)*0.999)

o.low_a<-optimize(alpha.optim_a,
                data=Incubations,
                lower=0,
                upper=(beta.estimation_a$maximum)*1.1)


gamma.seq<-seq(0,2,0.01)
ll_a<- -sapply(gamma.seq,FUN=lh_func_a,data=Incubations)

#plot
par(mfrow=c(1,2))
plot(gamma.seq,ll_a,type="l",
     xlab="beta parameter for Alpha variant",
     ylab="log-likelihood(beta)")
abline(v=beta.estimation_a$maximum,lty=1,col=2)
abline(v=o.low_a$minimum,lty=2,col=4)
abline(v=o.up_a$minimum,lty=2,col=4)
abline(h=max.loglik-0.5*c.alpha,lty=2,col=1)


plot(gamma.seq,ll_d,type="l",
     xlab="beta parameter for Delta variant",
     ylab="log-likelihood(beta)")
abline(v=beta.estimation_d$maximum,lty=1,col=2)
abline(v=o.low_d$minimum,lty=2,col=4)
abline(v=o.up_d$minimum,lty=2,col=4)
abline(h=max.loglik-0.5*c.alpha,lty=2,col=1)
```

The confidence interval (CI) for beta parameter of Alpha and Delta
variants are [`r o.low_a$minimum`, `r o.up_a$minimum`] and
[`r o.low_d$minimum`, `r o.up_d$minimum`].s

# Question 7


```{r, echo=FALSE}
a= 3.0
p_alpha<- sum(sapply(1:6, function(i) (beta.max_alpha^a*(i^(a-1))*exp(-beta.max_alpha*i))/gamma(a)))
d <- 4.43
p_delta<- sum(sapply(1:6, function(i) (beta.max_delta^d*(i^(d-1))*exp(-beta.max_delta*i))/gamma(d)))
```

$$P(x < 7|Alpha) =\sum_{x=1}^6\frac{\beta^\alpha}{\Gamma(\alpha)} x^{\alpha-1} \exp (-\beta x)\ where\ \alpha\ =\ 3.08\ and\ \beta\ =\ 0.609 \\
P(x < 7|Alpha) = `r p_alpha`$$

$$P(x < 7|Delta) =\sum_{x=1}^6\frac{\beta^\alpha}{\Gamma(\alpha)} x^{\alpha-1} \exp (-\beta x)\ where\ \alpha\ =\ 4.43\ and\ \beta\ =\ 0.97 \\
P(x < 7|Delta) = `r p_delta`$$
The estimated probability that the incubation period is less than one week (7day) for Alpha and Delta variants us 0.76 and 0.83 respectively.

# Question 8

```{r, echo=FALSE}
p_14_a <- beta.max_alpha^a*(14^(a-1))*exp(-beta.max_alpha*14)/gamma(a)
p_14_d <- beta.max_delta^d*(14^(d-1))*exp(-beta.max_delta*14)/gamma(d)
p_alpha_14 <- (p_14_a*0.2)/(p_14_a+p_14_d)
```

We calculated the probability, $P(Alpha | x\ =\ 14\ days)$ given the
prevalence of Alpha and Delta variants are 0.2 and 0.8. 
$$ \text {Using bayes rule;}\\
P(Alpha|x\ =\ 14) = \frac{P(x\ =\ 14|Alpha) * P(Alpha)}{P(x\ =\ 14)}\\
\text{having}\ P(x\ =\ 14) = P( x\ =\ 14|Alpha) + P(x\ =\ 14|Delta)\\
= `r p_14_a` + `r p_14_d`\\
\text{then}\space P(Alpha|x\ =\ 14) = `r p_alpha_14`$$

# Question 9
Using the estimates obtained above and the formula of the
mean and variance of gamma distribution, we obtained the confidence
interval for the mean incubation time for the Alpha variant using the
formula. $$ 95\% \space CI(\mu) = \mu \pm t(n-1;\alpha/2)SE\\ \text{where}\space SE= \frac{S}{\sqrt{n}}$$

```{r, echo=FALSE}
sd_a <- 2.97
mean_a <- 3.08/beta.max_alpha
sd_a <-sqrt (3.08/(beta.max_alpha)^2)
ci <- (mean_a+c(-1,1))*qt(0.975, n_Alpha-1)*sd_a/sqrt(n_Alpha)
```

The 95% confidence interval indicate that if the study is repeated many times 95% of the study will cover the true population mean of Incubation period between 1.6596515 and 2.4778407.
