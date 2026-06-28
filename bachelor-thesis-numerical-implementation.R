
# (4.2-4.4) Numerical sensitivity of backward-looking caplets/floorlets in the 
# Discrete Hull-White model


# Base case 

Delta   <- 1/360          # ACT/360, one step = one calendar day
a0      <- 0.75           # annualized mean reversion 
sigHW0  <- 0.01           # annualized HW short-rate vol (100 bp)
sig0    <- sigHW0 * Delta # thesis per-period sigma = sigma_HW * Delta
m0      <- 90             # base pre-accrual length (3M, daily steps)
L0      <- 90             # base accrual length     (3M, daily steps)
Kstar0  <- 0.035          # ATM strike 
y0_lvl  <- 0.03           # base yield level for 1e curves


# S_ell = (1-c^ell)/(1-c),  c = 1 - a*Delta
S_ell <- function(ell, a, Delta) {
  c <- 1 - a * Delta
  (1 - c^ell) / (1 - c)
}

# in-accrual variance nu_{m,N},  L = N - m 
nu_in <- function(L, a, sigma, Delta) {
  if (L < 2) return(0)
  ell <- 1:(L - 1)
  sigma^2 * Delta * sum(S_ell(ell, a, Delta)^2)
}

# pre-accrual variance nu_W at n=0   
nu_pre <- function(m, L, a, sigma, Delta) {
  if (m < 1) return(0)
  c    <- 1 - a * Delta
  SL   <- S_ell(L, a, Delta)
  geom <- (1 - c^(2 * m)) / (1 - c^2) 
  sigma^2 * Delta * SL^2 * geom
}

# total variance
nu_X <- function(m, L, a, sigma, Delta){
  nu_pre(m, L, a, sigma, Delta) + nu_in(L, a, sigma, Delta)
}

# normalized caplet/floorlet (p(0,m)=1), fixed maturities 
# depends only on (nu, K, Kstar, alpha)
# ATM => depends only on nu

caplet_norm <- function(nu, K, Kstar, alpha) {
  ratio <- log((1 + alpha * Kstar) / (1 + alpha * K))   
  d1 <- (ratio + 0.5 * nu) / sqrt(nu)
  d2 <- d1 - sqrt(nu)
  pnorm(d1) - (1 + alpha * K) / (1 + alpha * Kstar) * pnorm(d2)
}

floorlet_norm <- function(nu, K, Kstar, alpha) {
  ratio <- log((1 + alpha * Kstar) / (1 + alpha * K))
  d1 <- (ratio + 0.5 * nu) / sqrt(nu)
  d2 <- d1 - sqrt(nu)
  (1 + alpha * K) / (1 + alpha * Kstar) * pnorm(-d2) - pnorm(-d1)
}

# absolute caplet at n=0 on a linear-yield curve (for 1e)
p0 <- function(N, y0, s, Delta) exp(-(y0 + s * N * Delta) * N * Delta)
caplet_n0 <- function(m, L, a, sigma, Delta, K, y0, s) {
  N  <- m + L
  pm <- p0(m, y0, s, Delta)
  pN <- p0(N, y0, s, Delta)
  alpha <- L * Delta; kj <- 1 + alpha * K
  nu <- nu_X(m, L, a, sigma, Delta)
  d1 <- (log(pm / (kj * pN)) + 0.5 * nu) / sqrt(nu)
  d2 <- d1 - sqrt(nu)
  pm * pnorm(d1) - kj * pN * pnorm(d2)
}


#  BASE-CASE DIAGNOSTICS

nW  <- nu_pre(m0, L0, a0, sig0, Delta)
nI  <- nu_in(L0, a0, sig0, Delta)
nX  <- nW + nI
alpha0 <- L0 * Delta
capATM <- caplet_norm(nX, Kstar0, Kstar0, alpha0)

cat("================ BASE CASE ================\n")
cat(sprintf("Delta = 1/360,  a = %.2f,  sigma_HW = %.4f,  sigma = %.3e\n", a0, sigHW0, sig0))
cat(sprintf("m = %d, L = %d  (c = %.6f, c^L = %.4f, S_L = %.2f vs L = %d)\n",
            m0, L0, 1 - a0*Delta, (1-a0*Delta)^L0, S_ell(L0,a0,Delta), L0))
cat(sprintf("nu_W      = %.4e   sqrt = %.4e\n", nW, sqrt(nW)))
cat(sprintf("nu_{m,N}  = %.4e   sqrt = %.4e\n", nI, sqrt(nI)))
cat(sprintf("nu_X      = %.4e   sqrt = %.4e\n", nX, sqrt(nX)))
cat(sprintf("pre/total share = %.1f%%,  in/total share = %.1f%%\n", 100*nW/nX, 100*nI/nX))
cat(sprintf("total vol sqrt(nu_X) = %.4e\n", sqrt(nX)))
cat(sprintf("ATM normalized caplet = %.4e  (= %.2f bp of p(0,m))\n", capATM, 1e4*capATM))


col_W <- "blue"
col_I <- "red"
col_X <- "black"
col_up <- "darkgreen"
col_fl <- "gray"
col_in <- "purple"


# (4.2a) 
draw_1a <- function() {
  par(mfrow = c(1,2))
  
  # panel (i): nu_W, nu_in, nu_X vs m  (fixed L=L0)
  mm <- 1:720
  vW <- sapply(mm, nu_pre, L=L0, a=a0, sigma=sig0, Delta=Delta) * 1e6
  vI <- rep(nu_in(L0,a0,sig0,Delta), length(mm)) * 1e6
  vX <- vW + vI
  plot(mm, vX, type="l", lwd=2, col=col_X, ylim=c(0, max(vX)),
       xlab="pre-accrual length  m  (days)",
       ylab=expression(variance ~~ "(10^-6)"), main="(i)  fixed accrual L = 90")
  lines(mm, vW, lwd=2, col=col_W)
  lines(mm, vI, lwd=2, col=col_I)
  legend("topleft", lwd=2, col=c(col_X,col_W,col_I),
         legend=c(expression(nu[X]), expression(nu[W]~"(pre-accrual)"),
                  expression(nu[m*","*N]~"(in-accrual)")),cex = 0.7)
  
  # panel (ii): nu_in vs L
  LL <- 2:360
  vIL <- sapply(LL, nu_in, a=a0, sigma=sig0, Delta=Delta) * 1e6
  plot(LL, vIL, type="l", lwd=2, col=col_I,
       xlab="accrual length  L = N - m  (days)",
       ylab=expression(nu[m*","*N] ~~ "(10^-6)"),
       main="(ii)  in-accrual variance")
  par(mfrow=c(1,1))
}


# (4.2b) 

draw_1b <- function() {
  sHW <- seq(0.001, 0.025, length.out=200)
  cap <- sapply(sHW, function(s) {
    nu <- nu_X(m0, L0, a0, s*Delta, Delta)
    1e4 * caplet_norm(nu, Kstar0, Kstar0, alpha0)
  })
  plot(100*sHW, cap, type="l", lwd=2, col=col_X,
       xlab=expression("annualized short-rate volatility  "*sigma[HW]*"  (%)"),
       ylab="ATM caplet  (bp of p(0,m))",
       main=expression("caplet vs "*sigma*"  (ATM, fixed maturities)"))
  abline(v=100*sigHW0, lty=3, col=col_I)
  text(100*sigHW0, max(cap)*0.1, "base", col=col_I, pos=4)
}

## (4.2c)

draw_1c <- function() {
  par(mfrow=c(1,2))
  aa <- seq(0.05, 3, length.out=200)
  cap <- sapply(aa, function(a) 1e4*caplet_norm(nu_X(m0,L0,a,sig0,Delta), Kstar0, Kstar0, alpha0))
  plot(aa, cap, type="l", lwd=2, col=col_X,
       xlab="mean reversion  a", ylab="ATM caplet  (bp of p(0,m))",
       main="caplet vs a")
  abline(v=a0, lty=3, col=col_I)
  text(a0, min(cap)+0.05*diff(range(cap)), "base", col=col_I, pos=4)
  vW <- sapply(aa, nu_pre, m=m0, L=L0, sigma=sig0, Delta=Delta)*1e6
  vI <- sapply(aa, nu_in,  L=L0, sigma=sig0, Delta=Delta)*1e6
  plot(aa, vW+vI, type="l", lwd=2, col=col_X, ylim=c(0,max(vW+vI)),
       xlab="mean reversion  a", ylab=expression(variance~~"(10^-6)"),
       main="variance components vs a")
  lines(aa, vW, lwd=2, col=col_W); lines(aa, vI, lwd=2, col=col_I)
  abline(v=a0, lty=3, col="gray")
  legend("topright", lwd=2, col=c(col_X,col_W,col_I),
         legend=c(expression(nu[X]),expression(nu[W]),expression(nu[m*","*N])))
  par(mfrow=c(1,1))
}

# (4.2d)

draw_1d <- function() {
  KK <- seq(Kstar0-0.02, Kstar0+0.02, length.out=200)
  cap <- sapply(KK, function(K) 1e4*caplet_norm(nX, K, Kstar0, alpha0))
  flr <- sapply(KK, function(K) 1e4*floorlet_norm(nX, K, Kstar0, alpha0))
  plot(100*KK, cap, type="l", lwd=2, col=col_X, ylim=range(c(cap,flr)),
       xlab="strike  K  (%)", ylab="price  (bp of p(0,m))",
       main="caplet & floorlet vs strike")
  lines(100*KK, flr, lwd=2, col=col_I)
  abline(v=100*Kstar0, lty=3, col="gray")
  legend("top", lwd=2, col=c(col_X,col_I),
         legend=c("caplet","floorlet"), horiz=TRUE)
  text(100*Kstar0, max(cap), expression(K^"*"), col="gray", pos=3, xpd=TRUE)
}

## (4.2e)
draw_1e <- function() {
  
  ## fixed external strike = base ATM on the FLAT curve at base maturities,
  ## so all three regimes start near ATM and tilt symmetrically
  sl   <- 0.004   # +-40 bp of yield per year of maturity
  Kfix <- (1/alpha0) * (p0(m0,y0_lvl, 0,Delta)/p0(m0+L0,y0_lvl,0,Delta) - 1)
  par(mfrow=c(1,2))
  
  # pre-accrual: vary m, fixed L
  
  mm <- seq(1, 720, by=3)
  pe <- function(s) sapply(mm, caplet_n0, L=L0, a=a0, sigma=sig0, Delta=Delta, K=Kfix, y0=y0_lvl, s=s)
  up <- 1e4*pe(sl)
  fl <- 1e4*pe(0); iv <- 1e4*pe(-sl)
  xmo <- mm/30
  plot(xmo, up, type="l", lwd=2, col=col_up, ylim=range(c(up,fl,iv)),
       xlab="pre-accrual length (months)", ylab="caplet  (bp)",
       main="pre-accrual, fixed L")
  lines(xmo, fl, lwd=2, col=col_fl, lty=2); lines(xmo, iv, lwd=2, col=col_in)
  legend("topleft", lwd=2, lty=c(1,2,1), col=c(col_up,col_fl,col_in),
         legend=c("upward","flat","inverted"))
  
  # in-accrual: vary L, fixed m
  
  LL <- seq(2, 720, by=3)
  pin <- function(s) sapply(LL, function(L) caplet_n0(m0, L, a0, sig0, Delta, Kfix, y0_lvl, s))
  up2 <- 1e4*pin(sl)
  fl2 <- 1e4*pin(0)
  iv2 <- 1e4*pin(-sl)
  xmo2 <- LL/30
  plot(xmo2, up2, type="l", lwd=2, col=col_up, ylim=range(c(up2,fl2,iv2)),
       xlab="accrual length (months)", ylab="caplet  (bp)",
       main="in-accrual, fixed m")
  lines(xmo2, fl2, lwd=2, col=col_fl, lty=2); lines(xmo2, iv2, lwd=2, col=col_in)
  legend("topleft", lwd=2, lty=c(1,2,1), col=c(col_up,col_fl,col_in),
         legend=c("upward","flat","inverted"))
  par(mfrow=c(1,1))
}

draw_1a()
draw_1b()
draw_1c()
draw_1d()
draw_1e()


#  (4.4)  Numerical tests of the backward-forward gap
#  Valuation at n = 0; gaps normalized per unit p(0,m), in bp.

gap_norm <- function(nuW, nuX, K, Kstar, alpha) {
  caplet_norm(nuX, K, Kstar, alpha) - caplet_norm(nuW, K, Kstar, alpha)
}

# ATM leading-order closed form (normalized): (sqrt(nuX)-sqrt(nuW))/sqrt(2*pi)
gap_atm_cf <- function(nuW, nuX) (sqrt(nuX) - sqrt(nuW)) / sqrt(2 * pi)

# ATM relative gap = sqrt(1 + nu_in/nuW) - 1  
gap_rel <- function(nuW, nuI) sqrt(1 + nuI / nuW) - 1

# base-case variances
nW0 <- nu_pre(m0, L0, a0, sig0, Delta)
nI0 <- nu_in(L0, a0, sig0, Delta)
nX0 <- nW0 + nI0

# (4.4a)
draw_4a <- function() {
  par(mfrow = c(1,2))
  sHW <- seq(0.001, 0.025, length.out = 200)
  ex <- sapply(sHW, function(s) { sg <- s*Delta
  w <- nu_pre(m0,L0,a0,sg,Delta)
  x <- w + nu_in(L0,a0,sg,Delta)
  1e4 * gap_norm(w, x, Kstar0, Kstar0, alpha0) })
  cf <- sapply(sHW, function(s) { sg <- s*Delta
  w <- nu_pre(m0,L0,a0,sg,Delta)
  x <- w + nu_in(L0,a0,sg,Delta)
  1e4 * gap_atm_cf(w, x) })
  plot(100*sHW, ex, type = "l", lwd = 2, col = col_X,
       xlab = expression("annualized short-rate volatility  "*sigma[HW]*"  (%)"),
       ylab = "ATM gap  (bp of p(0,m))", main = "gap vs volatility (ATM)")
  points(100*sHW[seq(1,200,by=20)], cf[seq(1,200,by=20)], pch = 1, col = col_I)
  abline(v = 100*sigHW0, lty = 3, col = "gray")
  legend("topleft", lwd = c(2,NA), pch = c(NA,1), col = c(col_X,col_I),
         legend = c("exact gap", "leading-order formula"),cex=0.8)
  
  #relative gap vs volatility
  
  rel <- sapply(sHW, function(s) { sg <- s*Delta
  100 * gap_rel(nu_pre(m0,L0,a0,sg,Delta), nu_in(L0,a0,sg,Delta)) })
  plot(100*sHW, rel, type = "l", lwd = 2, col = col_X,
       xlab = expression(sigma[HW]*"  (%)"), ylab = "relative gap  (%)",
       main = "relative gap vs volatility")
  abline(v = 100*sigHW0, lty = 3, col = "gray")
  par(mfrow = c(1,1))
}

# (4.4b) 
draw_4b <- function() {
  par(mfrow = c(1,2))
  LL <- seq(5, 720, by = 3)
  ex <- sapply(LL, function(L) { 
  w <- nu_pre(m0,L,a0,sig0,Delta)
  x <- w + nu_in(L,a0,sig0,Delta)
  1e4 * gap_norm(w, x, Kstar0, Kstar0, L*Delta) 
  })
  plot(LL/30, ex, type = "l", lwd = 2, col = col_X,
       xlab = "accrual length  L  (months)", ylab = "ATM gap  (bp of p(0,m))",
       main = "gap vs accrual length")
  abline(v = L0/30, lty = 3, col = col_I); text(L0/30, max(ex)*0.08, "base", col = col_I, pos = 4)
  rel <- sapply(LL, function(L) 100 * gap_rel(nu_pre(m0,L,a0,sig0,Delta), nu_in(L,a0,sig0,Delta)))
  plot(LL/30, rel, type = "l", lwd = 2, col = col_X,
       xlab = "accrual length  L  (months)", ylab = "relative gap  (%)",
       main = "relative gap vs accrual length")
  abline(v = L0/30, lty = 3, col = col_I); abline(h = 100, lty = 2, col = "gray")
  par(mfrow = c(1,1))
}

# (4.4c) 

draw_4c <- function() {
  par(mfrow = c(1,2))
  aa <- seq(0.05, 3, length.out = 200)
  ex <- sapply(aa, function(a) { 
  w <- nu_pre(m0,L0,a,sig0,Delta)
  x <- w + nu_in(L0,a,sig0,Delta)
  1e4 * gap_norm(w, x, Kstar0, Kstar0, alpha0) 
  })
  plot(aa, ex, type = "l", lwd = 2, col = col_X,
       xlab = "mean reversion  a", ylab = "ATM gap  (bp of p(0,m))",
       main = "gap vs mean reversion")
  abline(v = a0, lty = 3, col = col_I); text(a0, min(ex), "base", col = col_I, pos = 4)
  rel <- sapply(aa, function(a) 100 * gap_rel(nu_pre(m0,L0,a,sig0,Delta), nu_in(L0,a,sig0,Delta)))
  plot(aa, rel, type = "l", lwd = 2, col = col_X,
       xlab = "mean reversion  a", ylab = "relative gap  (%)",
       main = "relative gap vs mean reversion")
  abline(v = a0, lty = 3, col = col_I)
  par(mfrow = c(1,1))
}

# (4.4d) 
draw_4d <- function() {
  KK <- seq(Kstar0 - 0.02, Kstar0 + 0.02, length.out = 200)
  ex <- sapply(KK, function(K) 1e4 * gap_norm(nW0, nX0, K, Kstar0, alpha0))
  plot(100*KK, ex, type = "l", lwd = 2, col = col_X,
       xlab = "strike  K  (%)", ylab = "gap  (bp of p(0,m))", main = "gap vs strike")
  abline(v = 100*Kstar0, lty = 3, col = "gray")
  text(100*Kstar0, max(ex), expression(K^"*"), col = "gray", pos = 3, xpd = TRUE)
}

draw_4a()
draw_4b()
draw_4c()
draw_4d()

