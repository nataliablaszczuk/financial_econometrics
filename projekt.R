###########################################
# Ekonometria Finansowa II - projekt
###########################################

# Wczytanie pakietów
library(zoo)
library(xts)
library(moments)
library(forecast)
require(MASS)
library(tseries)
library(rugarch)
library(knitr)
library(ggplot2)
source("Block2Functions.R")

rm(list = ls())

# 1. Ładowanie danych i obliczanie zwrotów portfela
####################################

filename <- "wig20.Rdata"
load(filename)

dates     <- index(data)
startDate <- as.Date("2005-01-01")
endDate   <- as.Date("2050-01-01")
y         <- window(data, start = startDate, end = endDate)
y         <- na.omit(y[, c("lpp", "kru")])

dy <- 100 * diff(log(y))  # Logarytmiczne zwroty - mierzenie cen aktywów w czasie
dates <- index(dy)
w <- c(0.5, 0.5)           # Wagi portfela
r <- zoo(dy %*% w, dates)  # Zwroty portfela - średnia ważona
R <- as.numeric(coredata(r)) # Wartości zwrotów

# 2. Wykresy szeregów czasowych i statystyki
###########################################

# Wykres szeregu czasowego stóp zwrotu portfela
jpeg("~/Desktop/SGH/Ekonometria finansowa II/zwroty_portfela.jpg", width = 1200, height = 400, quality = 100)
plot(r, main = "Logarytmiczne stopy zwrotu portfela", col = "green", lwd = 2)
abline(h = 0, col = "grey", lty = 2)  # Linia oznaczająca zero
dev.off()


# Statystyki opisowe
Nyear <- 252  # Liczba dni handlowych w roku
mu <- mean(R) * Nyear
sig <- sd(R) * sqrt(Nyear)
mom <- data.frame(
  Stat = c("Mean (annualized)", "Std Dev (annualized)", "Min", "Max", "Skewness", "Kurtosis", "JB Test Statistic"),
  Value = c(mu, sig, min(R), max(R), skewness(R), kurtosis(R), jarque.bera.test(R)$statistic)
)
kable(mom, digits = 3)

# QQ plot
jpeg("~/Desktop/SGH/Ekonometria finansowa II/qq_plot.jpg", width = 1200, height = 400, quality = 300)
Rstar <- (R - mean(R)) / sd(R)  # Standaryzowane zwroty

qqplot(qnorm(ppoints(length(Rstar))), Rstar, 
       main = "QQ Wykres Standaryzowanych Zwrotów", 
       xlab = "Teoretyczne Kwantyle", ylab = "Kwantyle Próbki", 
       col = "green", pch = 19, cex = 0.6, bg = "grey")

abline(0, 1, col = "grey", lwd = 2)
dev.off()

jpeg("~/Desktop/SGH/Ekonometria finansowa II/ACF_zwroty_dzienne.jpg", width = 1500, height = 200, quality = 100)
# ACF dla zwrotów
Acf(R, main = "ACF dla dziennych zwrotów")
dev.off()

jpeg("~/Desktop/SGH/Ekonometria finansowa II/ACF_zwroty_kwadratowe.jpg", width = 1500, height = 200, quality = 100)
# ACF dla kwadratowych zwrotów
Acf(R^2, main = "ACF dla kwadratowych dziennych zwrotów")
dev.off


# 3. Estymacja najlepszego modelu GARCH
########################################

# Funkcja do wyboru najlepszych opóźnień GARCH(p, q)
LagSel <- function(x, Pmax = 4, Qmax = 4, crit = "SIC") {
  IC <- matrix(NA, Pmax, Qmax + 1)
  for (p in 1:Pmax) {
    for (q in 0:Qmax) {
      spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(p, q)),
                         mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                         distribution.model = "std")
      fit <- tryCatch(ugarchfit(data = x, spec = spec), error = function(e) NA)
      if (!is.na(fit)) {
        IC[p, q + 1] <- infocriteria(fit)[crit == c("AIC", "SIC", "HQ")]
      }
    }
  }
  rownames(IC) <- paste0("p=", 1:Pmax)
  colnames(IC) <- paste0("q=", 0:Qmax)
  return(IC)
}

# Wybór najlepszego modelu na podstawie SIC
# niższa wartość SIC - lepszy model
ICtab <- LagSel(R, Pmax = 4, Qmax = 4, crit = "SIC")
kable(ICtab, digits = 3)

# Model 1 - klasyczny GARCH
# Wybrane opóźnienia
# W tym przypadku najniższa wartość wynosi 4.047 dla p=1 i q=1
pq   = c(1,1)
PQ   = c(0,0)
dist = "std"

spec1 = ugarchspec(variance.model=list(model="sGARCH", garchOrder=pq), 
                   mean.model=list(armaOrder=PQ, include.mean=TRUE),  
                   distribution.model=dist)
fit1 = ugarchfit(data=r, spec=spec1)
plot(fit1, which=12)

# Model 2 - eGARCH
spec.e   = ugarchspec(variance.model=list(model="eGARCH", garchOrder=pq), 
                      mean.model=list(armaOrder=PQ, include.mean=TRUE),  
                      distribution.model=dist)
fit.e   = ugarchfit(data=r, spec=spec.e)

# Model 3 - gjrGARCH
spec.gjr = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=pq), 
                      mean.model=list(armaOrder=PQ, include.mean=TRUE),  
                      distribution.model=dist)
fit.gjr = ugarchfit(data=r, spec=spec.gjr)
mod = "gjrGARCH"

# Model 4 - GARCH-in-Mean 
spec.m = ugarchspec(variance.model=list(model=mod, garchOrder=pq), 
                    mean.model=list(armaOrder=PQ, include.mean=TRUE, archm = TRUE), distribution.model=dist)
fit.m  = ugarchfit(data=r, spec=spec.m)

IC <- cbind(infocriteria(fit1), infocriteria(fit.e), infocriteria(fit.gjr), infocriteria(fit.m))
colnames(IC) <- c("GARCH", "eGARCH", "gjrGARCH","GARCH-in-mean")
kable(IC,digits=3)

# Wybór najlepszego modelu na podstawie Akaike i Bayes
# Najniższe wartości wskazują na model eGARCH jako najlepszy

# Wyświetlenie estymacji parametrów dla najlepszego modelu
param_table <- as.data.frame(fit.e@fit$matcoef)
colnames(param_table) <- c("Oszacowanie", "Błąd std.", "t-Value", "Pr(>|t|)")
kable(param_table, digits = 4, caption = "Estymacje parametrów dla modelu eGARCH")

# Wykres warunkowego odchylenia standardowego
jpeg("~/Desktop/SGH/Ekonometria finansowa II/wykres_odch_stnd_poprawiony.jpg", 
     width = 1500, height = 800, quality = 100, res = 300)

# Zmniejszenie marginesów
par(mar = c(4, 4, 3, 2))  # Ustawienia marginesów: dół, lewo, góra, prawo

# Tworzenie wykresu
plot(sigma_best, 
     type = "l", 
     col = "green", 
     lwd = 2, 
     main = "Warunkowe odchylenie standardowe (Conditional SD)", # Tytuł
     xlab = "Data", # Etykieta osi X
     ylab = "Odchylenie standardowe", # Etykieta osi Y
     xaxt = "n")  # Wyłączamy automatyczne oznaczenia osi X

# Dodanie osi czasu z rzadszymi oznaczeniami
axis(1, at = seq(from = index(sigma_best)[1], 
                 to = index(sigma_best)[length(index(sigma_best))], 
                 by = "1 year"), 
     labels = format(seq(from = index(sigma_best)[1], 
                         to = index(sigma_best)[length(index(sigma_best))], 
                         by = "1 year"), "%Y"))

# Dodanie siatki
grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")

# Kończymy zapis do pliku
dev.off()

### Model DCC-GARCH

require(rmgarch)
mod  = "eGARCH"
pq   = c(1,1)
PQ   = c(0,0)
dist = "std"
specU = ugarchspec(variance.model=list(model=mod, garchOrder=pq), 
                   mean.model=list(armaOrder=PQ, include.mean=TRUE), 
                   distribution.model=dist)
mspec = multispec(c(specU, specU))
specDCC <- dccspec(mspec, dccOrder = c(1,1), model = "DCC")
fitDCC  <- dccfit(specDCC,dy)  # start.pars = list())
fitDCC

### Symulacje dla VaR i ES
Nsim <- 10000     # Liczba symulacji
w <- c(0.5, 0.5)  # Wagi portfela

# Obliczanie statystyk portfela
m <- mean(R)  # Średnia zwrotów
s <- sd(R)  # Odchylenie standardowe zwrotów
v <- 4 + 6 / (kurtosis(R) - 3)  # Stopnie swobody dla t-Studenta

# Metody dla VaR i ES

### H = 1, p = 1%

H <- 1
p <- 0.01

# Normal
VaRnorm_1_1 <- qnorm(p)*s
VaRnorm_1_1
qf          <- function(x) qdist("norm", p=x)
ESnorm_1_1 <- m + s*(1/p * integrate(qf, 0, p)$value)
ESnorm_1_1

# T-student
VaRt11  <- m + s*qdist("std",shape=v,p=p)
VaRt11
qf    <- function(x) qdist("std", p=x, shape=v)
ESt11   <- m + s*(1/p * integrate(qf, 0, p)$value) 
ESt11

# Historical simulation
R0    <- sort(R)              
N0    <- floor(length(R0)*p)                            
VaRhsHH11 <- R0[N0]
VaRhsHH11
EShsHH11  <- mean(R0[1:N0])
EShsHH11

# EMWA
lambda <- 0.94
mu <- mean(R)
N <- length(R)

sigma2_EWMA <- rep(0, N)
sigma2_EWMA[1] <- var(R)  # Wariancja początkowa
for (t in 2:N) {
  sigma2_EWMA[t] <- lambda * sigma2_EWMA[t - 1] + (1 - lambda) * R[t - 1]^2
}
sigma_last <- sqrt(sigma2_EWMA[N])  # Ostatnie odchylenie standardowe
VaREWMA_1_1 <- qnorm(p) * sigma_last
ES_EWMA_1_1 <- mu - (sigma_last * dnorm(qnorm(p)) / p)

VaREWMA_1_1
ES_EWMA_1_1

# GARCH
fct = ugarchforecast(fit.e,data=r, n.ahead = 1)
sig <- sigma(fct)
mu  <- fitted(fct)
v <- fit.e@fit$coef["shape"]
qf   <- function(x) qdist("std", p=x, shape=v)

VaRgarch <- mu + sig*qdist("std", p, shape=v)
ESgarch  <- mu + sig*(1/p * integrate(qf, 0, p)$value)  
VaRgarch
ESgarch

# DCC-GARCH
simDCC <- dccsim(fitDCC, n.sim = H, m.sim = Nsim, startMethod = c("sample"), rseed = 7)
draws  <- simDCC@msim$simX # draws for returns
Rdraws <- numeric(Nsim)
for (m in 1:Nsim) {
  Rdraws[m] <- draws[[m]] %*% w  # Portfolio returns based on weights
}

VaRHdcc <- quantile(Rdraws, probs = p)  # Value at Risk
ESHdcc <- mean(Rdraws[Rdraws <= VaRHdcc])
VaRHdcc
ESHdcc

### H = 1, p = 5%

H <- 1
p <- 0.05
Nsim <- 10000     # Liczba symulacji
w <- c(0.5, 0.5)  # Wagi portfela
m <- mean(R)  # Średnia zwrotów
s <- sd(R)  # Odchylenie standardowe zwrotów
v <- 4 + 6 / (kurtosis(R) - 3)

# Normal
VaRnorm_1_1 <- qnorm(p)*s
VaRnorm_1_1
qf          <- function(x) qdist("norm", p=x)
ESnorm_1_1 <- m + s*(1/p * integrate(qf, 0, p)$value)
ESnorm_1_1

# T-student
VaRt11  <- m + s*qdist("std",shape=v,p=p)
VaRt11
qf    <- function(x) qdist("std", p=x, shape=v)
ESt11   <- m + s*(1/p * integrate(qf, 0, p)$value) 
ESt11

# Historical simulation
R0    <- sort(R)              
N0    <- floor(length(R0)*p)                            
VaRhsHH11 <- R0[N0]
VaRhsHH11
EShsHH11  <- mean(R0[1:N0])
EShsHH11

# EMWA
lambda <- 0.94
mu <- mean(R)
N <- length(R)

sigma2_EWMA <- rep(0, N)
sigma2_EWMA[1] <- var(R)  # Wariancja początkowa
for (t in 2:N) {
  sigma2_EWMA[t] <- lambda * sigma2_EWMA[t - 1] + (1 - lambda) * R[t - 1]^2
}
sigma_last <- sqrt(sigma2_EWMA[N])  # Ostatnie odchylenie standardowe
VaREWMA_1_1 <- qnorm(p) * sigma_last
ES_EWMA_1_1 <- mu - (sigma_last * dnorm(qnorm(p)) / p)

VaREWMA_1_1
ES_EWMA_1_1

# GARCH
fct = ugarchforecast(fit.e,data=r, n.ahead = 1)
sig <- sigma(fct)
mu  <- fitted(fct)
v <- fit.e@fit$coef["shape"]
qf   <- function(x) qdist("std", p=x, shape=v)

VaRgarch <- mu + sig*qdist("std", p, shape=v)
ESgarch  <- mu + sig*(1/p * integrate(qf, 0, p)$value)  
VaRgarch
ESgarch

# DCC-GARCH
simDCC <- dccsim(fitDCC, n.sim = H, m.sim = Nsim, startMethod = c("sample"), rseed = 7)
draws  <- simDCC@msim$simX # draws for returns
Rdraws <- numeric(Nsim)
for (m in 1:Nsim) {
  Rdraws[m] <- draws[[m]] %*% w  # Portfolio returns based on weights
}

VaRHdcc <- quantile(Rdraws, probs = p)  # Value at Risk
ESHdcc <- mean(Rdraws[Rdraws <= VaRHdcc])
VaRHdcc
ESHdcc

### H = 10, p = 5%

H <- 10
p <- 0.05
Nsim <- 10000     # Liczba symulacji
w <- c(0.5, 0.5)  # Wagi portfela
m <- mean(R)  # Średnia zwrotów
s <- sd(R)  # Odchylenie standardowe zwrotów
v <- 4 + 6 / (kurtosis(R) - 3)

# Normalny
VaRHnorm <- sqrt(1:H)*qnorm(p)*s + (1:H)*m    
ESHnorm  <- (1:H)*m - sqrt(1:H)*s*dnorm(qnorm(p))/p
VaRHnorm
ESHnorm

# T-studenta
Nsim    <- 10000
Rdraws  <- matrix(rdist(distribution="std", Nsim*H, mu = m, sigma = s, shape = v),H,Nsim) 
RdrawsC <- apply(Rdraws,2,cumsum)
VaRHt   <- ESHt <- rep(NaN,H)
M0  <- floor(Nsim*p)     # observation for p-th quantile 
for(h in 1:H){
  temp   = sort(RdrawsC[h,])
  VaRHt[h] = temp[M0]
  ESHt[h]  = mean(temp[1:M0])
}
VaRHt
ESHt

# Historical simulation - Bootstrap resampling zwrotów dziennych dla H = 10
Rdraws <- matrix(sample(R, Nsim * H, replace = TRUE), H, Nsim)  # Symulacja zwrotów dziennych
RdrawsC <- apply(Rdraws, 2, cumsum)  # Skumulowane zwroty dla H
VaRHhs <- quantile(RdrawsC[H, ], probs = p)  # Kwantyl dla skumulowanych zwrotów
EShs <- mean(RdrawsC[H, ][RdrawsC[H, ] <= VaRHhs])  # Średnia poniżej VaR

VaRHhs
EShs

# EMWA
sigma_ewma <- numeric(length(R))
sigma_ewma[1] <- var(R)
for (t in 2:length(R)) {
  sigma_ewma[t] <- lambda * sigma_ewma[t - 1] + (1 - lambda) * R[t - 1]^2
}
sigma_ewma_t <- sqrt(sigma_ewma[length(R)])
Rdraws <- matrix(rnorm(Nsim * H, mean = 0, sd = sigma_ewma_t), H, Nsim)  # Generowanie zwrotów dziennych
RdrawsC <- apply(Rdraws, 2, cumsum)  # Skumulowane zwroty dla H
VaREWMA <- quantile(RdrawsC[H, ], probs = p)  # Kwantyl dla skumulowanych zwrotów
ESEWMA <- mean(RdrawsC[H, ][RdrawsC[H, ] <= VaREWMA])  # Średnia poniżej VaR

VaREWMA
ESEWMA

# eGARCH
sim2 <- ugarchsim(fit.e, n.sim = H, n.start = 0, m.sim = Nsim, startMethod = "sample")
Rdraws  <- fitted(sim2)
temp   <- RdrawsToVaRES(RDraws,p)
VaRHgarch <- temp$VaR
ESHgarch <- temp$ES
VaRHgarch
ESHgarch

#DCC-GARCH
simDCC <- dccsim(fitDCC, n.sim = H, m.sim = Nsim, startMethod = c("sample"), rseed = 7)
draws  <- simDCC@msim$simX # draws for returns
Rdraws <- matrix(NA,H,Nsim) # return from investing in our portfolio of 2 assets
for (m in 1:Nsim){
  Rdraws[,m]  = draws[[m]]%*%w    # w should be a vector of weights
}
temp      <- RdrawsToVaRES(RDraws,p)
VaRHdcc   <- temp$VaR
ESHdcc    <- temp$ES
VaRHdcc
ESHdcc

### H = 10, p = 1%

H <- 10
p <- 0.01
Nsim <- 10000     # Liczba symulacji
w <- c(0.5, 0.5)  # Wagi portfela
m <- mean(R)  # Średnia zwrotów
s <- sd(R)  # Odchylenie standardowe zwrotów
v <- 4 + 6 / (kurtosis(R) - 3)

# Normalny
VaRHnorm <- sqrt(1:H)*qnorm(p)*s + (1:H)*m    
ESHnorm  <- (1:H)*m - sqrt(1:H)*s*dnorm(qnorm(p))/p
VaRHnorm
ESHnorm

# T-studenta
Nsim    <- 10000
Rdraws  <- matrix(rdist(distribution="std", Nsim*H, mu = m, sigma = s, shape = v),H,Nsim) 
RdrawsC <- apply(Rdraws,2,cumsum)
VaRHt   <- ESHt <- rep(NaN,H)
M0  <- floor(Nsim*p)     # observation for p-th quantile 
for(h in 1:H){
  temp   = sort(RdrawsC[h,])
  VaRHt[h] = temp[M0]
  ESHt[h]  = mean(temp[1:M0])
}
VaRHt
ESHt

# Historical simulation - Bootstrap resampling zwrotów dziennych dla H = 10
Rdraws <- matrix(sample(R, Nsim * H, replace = TRUE), H, Nsim)  # Symulacja zwrotów dziennych
RdrawsC <- apply(Rdraws, 2, cumsum)  # Skumulowane zwroty dla H
VaRHhs <- quantile(RdrawsC[H, ], probs = p)  # Kwantyl dla skumulowanych zwrotów
EShs <- mean(RdrawsC[H, ][RdrawsC[H, ] <= VaRHhs])  # Średnia poniżej VaR

VaRHhs
EShs

# EMWA
sigma_ewma <- numeric(length(R))
sigma_ewma[1] <- var(R)
for (t in 2:length(R)) {
  sigma_ewma[t] <- lambda * sigma_ewma[t - 1] + (1 - lambda) * R[t - 1]^2
}
sigma_ewma_t <- sqrt(sigma_ewma[length(R)])
Rdraws <- matrix(rnorm(Nsim * H, mean = 0, sd = sigma_ewma_t), H, Nsim)  # Generowanie zwrotów dziennych
RdrawsC <- apply(Rdraws, 2, cumsum)  # Skumulowane zwroty dla H
VaREWMA <- quantile(RdrawsC[H, ], probs = p)  # Kwantyl dla skumulowanych zwrotów
ESEWMA <- mean(RdrawsC[H, ][RdrawsC[H, ] <= VaREWMA])  # Średnia poniżej VaR

VaREWMA
ESEWMA

# eGARCH
sim2 <- ugarchsim(fit.e, n.sim = H, n.start = 0, m.sim = Nsim, startMethod = "sample")
Rdraws  <- fitted(sim2)
temp   <- RdrawsToVaRES(RDraws,p)
VaRHgarch <- temp$VaR
ESHgarch <- temp$ES
VaRHgarch
ESHgarch

#DCC-GARCH
simDCC <- dccsim(fitDCC, n.sim = H, m.sim = Nsim, startMethod = c("sample"), rseed = 7)
draws  <- simDCC@msim$simX # draws for returns
Rdraws <- matrix(NA,H,Nsim) # return from investing in our portfolio of 2 assets
for (m in 1:Nsim){
  Rdraws[,m]  = draws[[m]]%*%w    # w should be a vector of weights
}
temp      <- RdrawsToVaRES(RDraws,p)
VaRHdcc   <- temp$VaR
ESHdcc    <- temp$ES
VaRHdcc
ESHdcc


#### BACKTESTING
# Ustawienia początkowe
p <- 0.05  # Poziom tolerancji dla VaR/ES
M <- 250   # Liczba obserwacji do backtestingu
R <- 500   # Długość okna kroczącego
realized <- tail(r, M) # realizacje
# Wyniki rolling VaR/ES dla różnych metod

### Historical simulation
temp   <- HSfct(r, M, R, p)
ESrhs  <- temp$ES
VaRrhs <- temp$VaR

# Backtesting dla HS
HR_HS     <- as.numeric(realized < VaRrhs)
a_HS <- mean(HR_HS)

# Kupiec test dla HS
temp_HS <- VaRTest(alpha = p, coredata(realized), coredata(VaR_HS))
n_exceed_HS <- temp_HS$actual.exceed
decision_kupiec_HS <- temp_HS$uc.Decision

n  = length(HR_HS)
n1 = sum(HR_HS) #Liczba przekroczeń
n1
n0 = n - n1
LR_uc   = (p/a_HS)^n1 *((1-p)/(1-a_HS))^n0
stat_uc = -2*log(LR_uc)
prob_uc = 1 - pchisq(stat_uc,1, lower.tail=TRUE)
prob_uc
# p-value = 10,16 nie odrzucamy H0

# Christoffersen Independence dla HS
HR0 <- HR_HS[2:M]
HR1 <- HR_HS[1:(M-1)]
n00 = sum(!HR1 & !HR0) # no exceedance after no exceedance
n01 = sum(!HR1 & HR0)  # exceedance after no exceedance
n10 = sum(HR1 & !HR0)  # no exceedance after exceedance
n11 = sum(HR1 & HR0)   # exceedance after exceedance

n0  = n00 + n10
n1  = n01 + n11
n1

pi01 = n01 / (n00+n01) # prob. of exceedance after no exceedance 
pi11 = n11 / (n10+n11) # prob. of exceedance after exceedance 
pi   = (n01+n11) / (n00+n01+n10+n11) # podobna wartosc do alpha

LR_ind   = pi^n1*(1-pi)^n0 / (pi01^n01 * pi11^n11 * (1-pi01)^n00  * (1-pi11)^n10)
stat_ind = -2*log(LR_ind)
prob_ind = 1 - pchisq(stat_ind,1)
prob_ind

# Christoffersen conditional coverage test
pi01 = n01 / (n00+n01) # prob. of exceedance after no exceedance 
pi11 = n11 / (n10+n11) # prob. of exceedance after exceedance 

# LR_cc   = (p/pi01)^n01 * (p/pi11)^n11 * ((1-p)/(1-pi01))^n00 * ((1-p)/(1-pi11))^n10
LR_cc=LR_ind*LR_uc
stat_cc = -2*log(LR_cc)
prob_cc = 1 - pchisq(stat_cc,2)
prob_cc

# McNeil-Frey dla HS
temp_MF_HS <- ESTest(alpha = p, realized, ESrhs, VaRrhs)
temp_MF_HS
decision_mf_HS <- temp_MF_HS$p.value


#### EWMA ###
temp <- EWMAfct(r, M, R, p, lambda = 0.94)
VaR_EWMA <- temp$VaR
ES_EWMA <- temp$ES

# Backtesting dla EWMA
HR_EWMA <- as.numeric(realized < VaR_EWMA)
a_EWMA <- mean(HR_EWMA)  # Hit ratio


# Kupiec test dla EWMA
# 1 sposob
temp_EWMA <- VaRTest(alpha = p, coredata(realized), coredata(VaR_EWMA))
n_exceed_EWMA <- temp_EWMA$actual.exceed
decision_kupiec_EWMA <- temp_EWMA$uc.Decision
temp_EWMA
# 2 sposob
n  = length(HR_EWMA)
n1 = sum(HR_EWMA) #Liczba przekroczeń
n1
n0 = n - n1
LR_uc   = (p/a_HS)^n1 *((1-p)/(1-a_HS))^n0
stat_uc = -2*log(LR_uc)
prob_uc = 1 - pchisq(stat_uc,1, lower.tail=TRUE)
prob_uc

# Christoffersen Independence test
HR0 <- HR_EWMA[2:M]
HR1 <- HR_EWMA[1:(M-1)]
n00 = sum(!HR1 & !HR0) # no exceedance after no exceedance
n01 = sum(!HR1 & HR0)  # exceedance after no exceedance
n10 = sum(HR1 & !HR0)  # no exceedance after exceedance
n11 = sum(HR1 & HR0)   # exceedance after exceedance

n0  = n00 + n10
n1  = n01 + n11
n1

pi01 = n01 / (n00+n01) # prob. of exceedance after no exceedance 
pi11 = n11 / (n10+n11) # prob. of exceedance after exceedance 
pi   = (n01+n11) / (n00+n01+n10+n11) # podobna wartosc do alpha

LR_ind   = pi^n1*(1-pi)^n0 / (pi01^n01 * pi11^n11 * (1-pi01)^n00  * (1-pi11)^n10)
stat_ind = -2*log(LR_ind)
prob_ind = 1 - pchisq(stat_ind,1)
prob_ind

# Christoffersen Conditional Coverage test
pi01 = n01 / (n00+n01) # prob. of exceedance after no exceedance 
pi11 = n11 / (n10+n11) # prob. of exceedance after exceedance 

# LR_cc   = (p/pi01)^n01 * (p/pi11)^n11 * ((1-p)/(1-pi01))^n00 * ((1-p)/(1-pi11))^n10
LR_cc=LR_ind*LR_uc
stat_cc = -2*log(LR_cc)
prob_cc = 1 - pchisq(stat_cc,2)
prob_cc

# McNeil-Frey dla EWMA
temp_MF_HS <- ESTest(alpha = p, realized, ES_EWMA, VaR_EWMA)
temp_MF_HS

### Best GARCH (eGARCH) ###
temp <- GARCHfct(r, M, R, p, v = 5)
VaR_GARCH <- temp$VaR
ES_GARCH <- temp$ES

# Backtesting dla GARCH
HR_GARCH <- as.numeric(realized < VaR_GARCH)
a_GARCH <- mean(HR_GARCH)

# Kupiec test dla GARCH
temp_GARCH <- VaRTest(alpha = p, coredata(realized), coredata(VaR_GARCH))
temp_GARCH
n_exceed_GARCH <- temp_GARCH$actual.exceed
decision_kupiec_GARCH <- temp_GARCH$uc.Decision

# Christoffersen Independence test
HR0 <- HR_GARCH[2:M]
HR1 <- HR_GARCH[1:(M-1)]
n00 = sum(!HR1 & !HR0) # no exceedance after no exceedance
n01 = sum(!HR1 & HR0)  # exceedance after no exceedance
n10 = sum(HR1 & !HR0)  # no exceedance after exceedance
n11 = sum(HR1 & HR0)   # exceedance after exceedance

n0  = n00 + n10
n1  = n01 + n11
n1

pi01 = n01 / (n00+n01) # prob. of exceedance after no exceedance 
pi11 = n11 / (n10+n11) # prob. of exceedance after exceedance 
pi   = (n01+n11) / (n00+n01+n10+n11) # podobna wartosc do alpha

LR_ind   = pi^n1*(1-pi)^n0 / (pi01^n01 * pi11^n11 * (1-pi01)^n00  * (1-pi11)^n10)
stat_ind = -2*log(LR_ind)
prob_ind = 1 - pchisq(stat_ind,1)
prob_ind

# Christoffersen conditional coverage
pi01 = n01 / (n00+n01) # prob. of exceedance after no exceedance 
pi11 = n11 / (n10+n11) # prob. of exceedance after exceedance 
LR_cc=LR_ind*LR_uc
stat_cc = -2*log(LR_cc)
prob_cc = 1 - pchisq(stat_cc,2)
prob_cc

# McNeil-Frey dla GARCH
temp_MF_GARCH <- ESTest(alpha = p, realized, ES_GARCH, VaR_GARCH)
decision_mf_GARCH <- temp_MF_GARCH$p.value
temp_MF_GARCH

### DCC-GARCH ###
temp <- DCCGARCHfct(dy, w, M, R, p, v = 5)
temp
VaR_DCC <- temp$VaR
ES_DCC <- temp$ES
VaR_DCC

# Backtesting dla DCC-GARCH
HR_DCC <- as.numeric(realized < VaR_DCC)
a_DCC <- mean(HR_DCC)

# Kupiec test dla DCC-GARCH
temp_DCC <- VaRTest(alpha = p, coredata(realized), coredata(VaR_DCC))
n_exceed_DCC <- temp_DCC$actual.exceed
decision_kupiec_DCC <- temp_DCC$uc.Decision
temp_DCC

# Christoffersen Independence Test
HR0 <- HR_DCC[2:M]
HR1 <- HR_DCC[1:(M-1)]
n00 = sum(!HR1 & !HR0) # no exceedance after no exceedance
n01 = sum(!HR1 & HR0)  # exceedance after no exceedance
n10 = sum(HR1 & !HR0)  # no exceedance after exceedance
n11 = sum(HR1 & HR0)   # exceedance after exceedance

n0  = n00 + n10
n1  = n01 + n11
n1

pi01 = n01 / (n00+n01) # prob. of exceedance after no exceedance 
pi11 = n11 / (n10+n11) # prob. of exceedance after exceedance 
pi   = (n01+n11) / (n00+n01+n10+n11) # podobna wartosc do alpha

LR_ind   = pi^n1*(1-pi)^n0 / (pi01^n01 * pi11^n11 * (1-pi01)^n00  * (1-pi11)^n10)
stat_ind = -2*log(LR_ind)
prob_ind = 1 - pchisq(stat_ind,1)
prob_ind

# Christoffersen Conditional Coverage
pi01 = n01 / (n00+n01) # prob. of exceedance after no exceedance 
pi11 = n11 / (n10+n11) # prob. of exceedance after exceedance 
LR_cc=LR_ind*LR_uc
stat_cc = -2*log(LR_cc)
prob_cc = 1 - pchisq(stat_cc,2)
prob_cc


# McNeil-Frey dla DCC-GARCH
temp_MF_DCC <- ESTest(alpha = p, realized, ES_DCC, VaR_DCC)
decision_mf_DCC <- temp_MF_DCC$p.value
temp_MF_DCC
