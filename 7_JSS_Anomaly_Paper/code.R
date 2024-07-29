library("anomaly")
set.seed(0)
x <- rnorm(5000)
x[401:500] <- rnorm(100, 4, 1)
x[1601:1800] <- rnorm(200, 0, 0.01)
x[3201:3500] <- rnorm(300, 0, 10)
x[c(1000, 2000, 3000, 4000)] <- rnorm(4, 0, 100)
x <- (x - median(x)) / mad(x)
res <- capa(x)
summary(res)

plot(res)

res <- capa(x, type = "mean")
collective_anomalies(res)

res <- capa(1 + 2 * x, type = "mean")
nrow(collective_anomalies(res))

data("machinetemp")
attach(machinetemp)
x <- (temperature - median(temperature)) / mad(temperature)
res <- capa(x, type = "mean")
canoms <- collective_anomalies(res)
dim(canoms)[1]

library("robustbase")
n <- length(x)
x.lagged <- matrix(c(x[1:(n - 1)], x[2:n]), n - 1, 2)
rho_hat <- covMcd(x.lagged, cor = TRUE)$cor[1,2]


inflated_penalty <- 3 * (1 + phi) / (1 - phi) * log(n)
res <- capa(x, type = "mean", beta = inflated_penalty, beta_tilde = inflated_penalty)
summary(res)


data("simulated")
res <- capa(sim.data, type = "mean", min_seg_len = 2)
plot(res, subset = 1:20)

beta <- 2 * log(ncol(sim.data):1)
beta[1] <- beta[1] + 3 * log(nrow(sim.data))
res <- capa(sim.data, type= "mean", min_seg_len = 2,beta = beta)
plot(res, subset = 1:20)


set.seed(0)
x1 <- rnorm(500)
x2 <- rnorm(500)
x3 <- rnorm(500)
x4 <- rnorm(500)
x1[151:200] <- x1[151:200] + 2
x2[171:200] <- x2[171:200] + 2
x3[161:190] <- x3[161:190] - 3
x1[351:390] <- x1[371:390] + 2
x3[351:400] <- x3[351:400] - 3
x4[371:400] <- x4[371:400] + 2
x4[451] <- x4[451] * max(1, abs(1 / x4[451])) * 6
x4[100] <- x4[100] * max(1, abs(1 / x4[100])) * 6
x2[050] <- x2[050] * max(1, abs(1 / x2[050])) * 6
x1 <- (x1 - median(x1)) / mad(x1)
x2 <- (x2 - median(x2)) / mad(x2)
x3 <- (x3 - median(x3)) / mad(x3)
x4 <- (x4 - median(x4)) / mad(x4)
 x <- cbind(x1, x2, x3, x4)
 res <- capa(x, max_lag = 20, type = "mean")
plot(res)

library("anomaly")
data("simulated")
res <- pass(sim.data, max_seg_len = 20, alpha = 3)
collective_anomalies(res)

library("anomaly")
data("simulated")
bard.res <- bard(sim.data)
sampler.res <- sampler(bard.res, gamma = 1/3, num_draws = 1000)
show(sampler.res)

plot(sampler.res, marginals = TRUE)

library("ecp")
data("ACGH")
acgh <- ACGH[[1]][,1:20]
library("dplyr")
ac_corrected <- function(X){
  n <- length(X)
  rcor <- covMcd(matrix(c(X[2:n], X[1:(n-1)]), ncol = 2), cor = TRUE)
  psi <- rcor$cor[1,2]
  correction_factor <- sqrt((1 - psi) / (1 + psi))
  return(correction_factor * (X - median(X)) / mad(X))
}
acgh %<>% data.frame %>% mutate_all(ac_corrected)
