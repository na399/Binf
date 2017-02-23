load("~/Downloads/exam.RData")


### Problem 2

n <- 10000

## Simulate n independent uniform numbers
u_list <- replicate(n, runif(1))

## Sort the stimultaed numbers
sorted_u_list <- sort(u_list)
sorted_u_list <- c(sorted_u_list,1)   #add U(n+1)


## Calculate the n+1 intervals of u(i) - u(i-1)
s_list = c()
for(i in c(seq(1,n+1))){
  if (i == 1){
    s <- sorted_u_list[i] - 0
  }
  else{
    s <- sorted_u_list[i] - sorted_u_list[i-1]
  }
  s_list <- c(s_list, s)
}

## Maximum likelihood estimation for the exponential distribution.
estimator = 1/mean(s_list)

plot(sort(s_list))

# Kernel plot
plot(density(s_list))
curve(dexp(x, estimator), col="red", add=TRUE, from=0, to=1e-3)
abline(v=0, col="blue")

# Kernel plot (adjusted bandwidth to minimise the graph to the left of 0)
plot(density(s_list, bw=2e-6))
curve(dexp(x, estimator), col="red", add=TRUE, from=0, to=1e-3)
abline(v=0, col="blue")


# Histogram
hist(s_list, breaks=12)
curve(dexp(x, estimator), col="red", add=TRUE, from=0, to=1e-3)

# QQ Plot
n = length(s_list)
tq_exp   = qexp((1:n-0.5)/n, estimator)

plot(tq_exp, sort(s_list), col="red")
abline(0,1)


### Problem 3

## Q3.1

linear_fit = glm(organisms ~ concentration, data = Cerio, family=poisson)
plot(log(organisms) ~ concentration, data = Cerio)
abline(linear_fit, col="red")


polynomial_fit = glm(organisms ~ poly(concentration,2), data = Cerio, family=poisson)

newdata = data.frame(concentration = seq(0,2,length.out=100))
lines(newdata$concentration, predict(polynomial_fit, newdata), col="blue")

legend("bottomleft", lty=c(1,1) ,col = c("red","blue"),
       c("linear model", "polynomial model"))


par(mfrow=c(2,2))
plot(linear_fit); plot(polynomial_fit)


anova(linear_fit, polynomial_fit, test = "LRT")

## Q3.2

linear_log = glm(dead ~ log(conc), data = fly.death, family = binomial("logit"))

polynomial_log = glm(dead ~ poly(log(conc),2), data = fly.death, family = binomial("logit"))

anova(linear_log, polynomial_log, test = "LRT")



### Problem 4

## Q4.4
dlaplace = function(x, mu, b) ifelse(x<mu,1/2*exp((x-mu)/b),1-1/2*exp(-(x-mu)/b))

plaplace = function(x, mu, b) 1/(2*b)*exp(-abs(x-mu)/b)
  
qlaplace = function(p, mu, b) ifelse(0<p & p<.5, mu + b*log(2*p), mu - b*(log(2*(1-p))))
  
## Q4.5
rlaplace = function(n, mu, b) qlaplace(runif(n), mu, b)

## Q4.6
curve(dlaplace(x, 0, 2), xlim=c(-15,15), main="CDF of the Laplace distribution")

curve(plaplace(x, 0, 2), xlim=c(-15,15), main="PDF of the Laplace distribution")

curve(qlaplace(x, 0, 2), xlim=c(1e-2,1), main="Quantile Function of the Laplace distribution", xlab="p", ylab="qlaplace(p, 0, 2)")
# x=0 and thus log(0) is undefined 

## Q4.7
mu   = integrate( function(x) x*plaplace(x, 1, 2), lower=-Inf, upper=Inf )$value
sig2 = integrate( function(x) (x-mu)^2*plaplace(x, 1, 2), lower=-Inf, upper=Inf )$value

# mu = 1 and sig2 = 8


## Q4.10

mle_laplace = function(x) {
  mu = median(x)
  b = 1/length(x)*(sum(abs(x-mu)))
  c(mu, b)
}

## Q4.11

muHat = mle_laplace(data)[1]
bHat = mle_laplace(data)[2]

# muHat = 1.214393 
# bHat = 2.050347


plot(density(data), ylim=c(0,0.25))
curve(plaplace(x, muHat, bHat), xlim=c(-15,15), add=TRUE, col="red")

# Q-Q plot
tq_laplace   = qlaplace((1:length(data)-0.5)/length(data), muHat, bHat)

plot(tq_laplace, sort(data), col="red")
abline(0,1)

# simulate bootstrap data
B = 10000; n = length(data)
bootstrap_data = matrix(
  rlaplace(B*n, muHat, bHat),
  ncol = B
)

# fit the model for each bootstrap sample
bootstrap_dist = apply(bootstrap_data, 2, mle_laplace)

# obtain the confidence interval from the bootstrap empirical distribution
conf_int1 = quantile(bootstrap_dist[1,], c(0.025, 0.975))
#     2.5%     97.5% 
# 0.8377695 1.6058310 

conf_int2 = quantile(bootstrap_dist[2,], c(0.025, 0.975))
#    2.5%    97.5% 
# 1.700070 2.409023



## 4.15
# Least absolute residual estimation (LARE)
# x is a vector for regressors
# y is a vector for observed response data
lare = function(x, y) {
  # search for the LARE numerically by optim
  opt = optim(
    par = c(1,1), # initial guess for parameters
    fn = function(params){
      sum(abs(y - (params[1] + params[2] * x)))
    }
  )
  # return the results
  opt$par
}

# simulate data for a given Laplace model with beta0 and beta1
simLaplaceRegression = function(beta0, beta1){
  x = seq(-10, 10, 0.01)
  y = beta0 + beta1 * x + rlaplace(length(x), 0, 1)
  sim = data.frame(regressor = x, dependent = y)
}

# test lare() by using simulated data
lare(simLaplaceRegression(0, 1)$regressor, simLaplaceRegression(0, 1)$dependent)
# -0.009056626  1.000992811

lare(simLaplaceRegression(2, 3)$regressor, simLaplaceRegression(2, 3)$dependent)
# 2.016784 3.002886

lare(simLaplaceRegression(-4, -5)$regressor, simLaplaceRegression(-4, -5)$dependent)
# -3.977439 -4.997884

## 4.16
plot(reg.data)
plot(reg.data.outlier)

linear_fit = lm(y ~ x, data=reg.data)
abline(reg=linear_fit, col="blue")

linear_fit_out = lm(y ~ x, data=reg.data.outlier)
abline(reg=linear_fit_out, col="blue", lty=2)

laplace_fit = lare(reg.data$x, reg.data$y)
abline(laplace_fit[1], laplace_fit[2], col="red")

laplace_fit_out = lare(reg.data.outlier$x, reg.data.outlier$y)
abline(laplace_fit_out[1], laplace_fit_out[2], col="red", lty=2)

legend("bottomright", lty=c(1,2,1,2), col = c("blue","blue","red","red"),
       c("linear", "linear with outliers", "Laplace", "Laplace with outliers"))

# a modified function of plot.lm() to plot Residuals vs Fitted (rvsf) of novel Laplace regression model
rvsf = function(fit, x, y) {
  beta0 = fit[1]
  beta1 = fit[2]
  
  fitted_values = beta0 + beta1 * x
  residuals = y - fitted_values
  
  ylim <- range(residuals, na.rm = TRUE)
  ylim <- extendrange(r = ylim, f = 0.08)
  
  label.pos = c(4, 2)  
  cex.id = 0.75
  labels.id <- paste(1:length(residuals))
                       
  text.id <- function(x, y, ind, adj.x = TRUE) {
    labpos <- if (adj.x) 
      label.pos[1 + as.numeric(x > mean(range(x)))]
    else 3
    text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE, 
         pos = labpos, offset = 0.25)
  }
  
  iid <- 1:3
  show.r <- sort.list(abs(residuals), decreasing = TRUE)[iid]
  
  plot(fitted_values, residuals, ylim=ylim, ylab="Residuals", 
       xlab="Fitted values", main="Residuals vs Fitted", type="n")
  panel.smooth(fitted_values, residuals)
  abline(h = 0, lty = 3, col = "gray")
  
  y.id <- residuals[show.r]
  y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
  text.id(fitted_values[show.r], y.id, show.r)

  
}

# Test new function, each pair of plots should be exactly the same
par(mfrow=c(2,2))
plot(linear_fit, which=1)
plot(linear_fit_out, which=1)
rvsf(linear_fit$coefficients, reg.data$x, reg.data$y)
rvsf(linear_fit_out$coefficients, reg.data.outlier$x, reg.data.outlier$y)
# They are the same, so my function works.

# Plot Residuals vs Fitted of the linear regression models and Laplace regresion models
rvsf(linear_fit$coefficients, reg.data$x, reg.data$y)
rvsf(linear_fit_out$coefficients, reg.data.outlier$x, reg.data.outlier$y)
rvsf(laplace_fit, reg.data$x, reg.data$y)
rvsf(laplace_fit_out, reg.data.outlier$x, reg.data.outlier$y)
