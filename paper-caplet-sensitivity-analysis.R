# =========================================================
# Numerical analysis: discrete Hull-White backward-looking
# caplets and the backward-forward gap
# Reproduces all figures and numbers of Section 4.
# =========================================================

# ---------- Base case ----------
Delta   <- 1/360          # ACT/360, one step = one calendar day
a0      <- 0.75           # mean reversion
sigHW0  <- 0.01           # annualized short-rate vol (100 bp)
sig0    <- sigHW0 * Delta # per-period sigma = sigma_ann * Delta
m0      <- 90             # base pre-accrual length (3M)
L0      <- 90             # base accrual length     (3M)
Katm0   <- 0.035          # ATM strike K^ATM for normalized plots

out_dir <- "figures"
if (!dir.exists(out_dir)) dir.create(out_dir)

save_png <- function(name, draw, width = 2000, height = 1000, res = 200) {
  png(file.path(out_dir, paste0(name, ".png")),
      width = width, height = height, res = res)
  draw()
  dev.off()
}

# ---------- Variance components ----------

S_ell <- function(ell, a, Delta) {
  c <- 1 - a * Delta
  (1 - c^ell) / (1 - c)
}

nu_in <- function(L, a, sigma, Delta) {           #computationally efficient closed form of the summation
  if (L <= 1) return(0)
  c <- 1 - a * Delta 
  s <- ((L - 1) - 2*c*(1 - c^(L-1))/(1 - c) + c^2*(1 - c^(2*(L-1)))/(1 - c^2)) / (1 - c)^2 
  sigma^2 * Delta * s 
}

nu_pre <- function(m, L, a, sigma, Delta) {       
  if (m < 1) return(0)
  c    <- 1 - a * Delta
  SL   <- S_ell(L, a, Delta)
  geom <- (1 - c^(2 * m)) / (1 - c^2)
  sigma^2 * Delta * SL^2 * geom
}

nu_X <- function(m, L, a, sigma, Delta) {
  nu_pre(m, L, a, sigma, Delta) + nu_in(L, a, sigma, Delta)
}

nu_pre_inf <- function(m, a, sigma, Delta) {      
  c <- 1 - a * Delta
  sigma^2 * Delta * (1 - c^(2 * m)) / ((a * Delta)^2 * (1 - c^2))
}

# ---------- Prices (normalized by p(0,m)) ----------

caplet_norm <- function(nu, K, Katm, alpha) {
  ratio <- log((1 + alpha * Katm) / (1 + alpha * K))
  d1 <- (ratio + 0.5 * nu) / sqrt(nu)
  d2 <- d1 - sqrt(nu)
  pnorm(d1) - (1 + alpha * K) / (1 + alpha * Katm) * pnorm(d2)
}

caplet_atm <- function(nu) 2 * pnorm(0.5 * sqrt(nu)) - 1

# ---------- Gap: EXACT closed forms (used in all plots) ----------

gap_norm <- function(nuW, nuX, K, Katm, alpha) {
  caplet_norm(nuX, K, Katm, alpha) - caplet_norm(nuW, K, Katm, alpha)
}

gap_atm_exact <- function(nuW, nuX) {            
  2 * (pnorm(0.5 * sqrt(nuX)) - pnorm(0.5 * sqrt(nuW)))
}

gap_rel_exact <- function(nuW, nuX) {             
  caplet_atm(nuX) / caplet_atm(nuW) - 1
}

M_max <- function(nuW, nuX) -0.5 * sqrt(nuX * nuW)
K_max <- function(nuW, nuX, Katm, alpha) {
  ((1 + alpha * Katm) * exp(-M_max(nuW, nuX)) - 1) / alpha
}

# ---------- Gap: leading-order forms ----------

gap_atm_lo <- function(nuW, nuX) (sqrt(nuX) - sqrt(nuW)) / sqrt(2 * pi)
gap_rel_lo <- function(nuW, nuI) sqrt(1 + nuI / nuW) - 1

# ---------- Base-case values ----------

nW0 <- nu_pre(m0, L0, a0, sig0, Delta)
nI0 <- nu_in(L0, a0, sig0, Delta)
nX0 <- nW0 + nI0
alpha0 <- L0 * Delta

# ---------- Colors ----------

col_W <- "blue"
col_I <- "red"
col_X <- "black"


# =================
# Section 4: plots
# =================

draw_4a <- function() {   
  par(mfrow = c(1,2))
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
                  expression(nu[m*","*N]~"(in-accrual)")), cex=0.7)
  LL <- 2:360
  vIL <- sapply(LL, nu_in, a=a0, sigma=sig0, Delta=Delta) * 1e6
  plot(LL, vIL, type="l", lwd=2, col=col_I,
       xlab="accrual length  L = N - m  (days)",
       ylab=expression(nu[m*","*N] ~~ "(10^-6)"),
       main="(ii)  in-accrual variance")
  par(mfrow=c(1,1))
}

draw_4b <- function() {   
  par(mfrow = c(1,2))
  LL <- seq(5, 720, by = 1)
  vw <- sapply(LL, function(L) nu_pre(m0,L,a0,sig0,Delta))
  vx <- vw + sapply(LL, function(L) nu_in(L,a0,sig0,Delta))
  ex <- 1e4 * gap_atm_exact(vw, vx)
  plot(LL/30, ex, type = "l", lwd = 2, col = col_X,
       xlab = "accrual length  L  (months)", ylab = "ATM gap  (bp of p(0,m))",
       main = "gap vs accrual length")
  abline(v = L0/30, lty = 3, col = col_I)
  text(L0/30, max(ex)*0.08, "base", col = col_I, pos = 4)
  rel <- 100 * gap_rel_exact(vw, vx)
  plot(LL/30, rel, type = "l", lwd = 2, col = col_X,
       xlab = "accrual length  L  (months)", ylab = "relative gap  (%)",
       main = "relative gap vs accrual length")
  abline(v = L0/30, lty = 3, col = col_I)
  abline(h = 100, lty = 2, col = "gray")
  par(mfrow = c(1,1))
}

draw_4c <- function() {   
  par(mfrow = c(1,2))
  aa <- seq(0.05, 3, length.out = 2000)
  vw <- sapply(aa, function(a) nu_pre(m0,L0,a,sig0,Delta))
  vx <- vw + sapply(aa, function(a) nu_in(L0,a,sig0,Delta))
  ex <- 1e4 * gap_atm_exact(vw, vx)
  plot(aa, ex, type = "l", lwd = 2, col = col_X,
       xlab = "mean reversion  a", ylab = "ATM gap  (bp of p(0,m))",
       main = "gap vs mean reversion")
  abline(v = a0, lty = 3, col = col_I)
  text(a0, min(ex), "base", col = col_I, pos = 4)
  rel <- 100 * gap_rel_exact(vw, vx)
  plot(aa, rel, type = "l", lwd = 2, col = col_X,
       xlab = "mean reversion  a", ylab = "relative gap  (%)",
       main = "relative gap vs mean reversion")
  abline(v = a0, lty = 3, col = col_I)
  par(mfrow = c(1,1))
}

draw_4d <- function() {   
  KK <- seq(Katm0 - 0.02, Katm0 + 0.02, length.out = 4000)
  ex <- sapply(KK, function(K) 1e4 * gap_norm(nW0, nX0, K, Katm0, alpha0))
  Kmx <- K_max(nW0, nX0, Katm0, alpha0)
  plot(100*KK, ex, type = "l", lwd = 2, col = col_X,
       xlab = "strike  K  (%)", ylab = "gap  (bp of p(0,m))",
       main = "gap vs strike")
  abline(v = 100*Katm0, lty = 3, col = "gray")
  abline(v = 100*Kmx,   lty = 3, col = col_I)   # coincides visually with K^ATM
  text(100*Katm0, max(ex), expression(K^"ATM"), col = "gray", pos = 3, xpd = TRUE)
  cat(sprintf("check 4d: grid argmax K = %.6f%%, closed-form K^max = %.6f%%\n",
              100*KK[which.max(ex)], 100*Kmx))
}

# ---------- Render everything ----------
save_png("fig_4a", draw_4a)
save_png("fig_4b", draw_4b)
save_png("fig_4c", draw_4c)
save_png("fig_4d", draw_4d)

# ===================================================
#  PAPER CLAIM CHECKS
#  Every number quoted in Section 4, recomputed here.
# ===================================================


cat("\n---- Section 4.1 text ----\n")
cat(sprintf("a*Delta = %.1e   (paper: 2.1e-3)\n", a0*Delta))
cat(sprintf("c = %.4f   (paper: 0.9979)\n", 1 - a0*Delta))
cat(sprintf("sigma = %.2e   (paper: 2.78e-5)\n", sig0))

cat("\n---- Table 1, base case ----\n")
cat(sprintf("c^L = %.3f   (paper: 0.829)\n", (1-a0*Delta)^L0))
cat(sprintf("S_L = %.1f   (paper: 82.1, vs L = 90)\n", S_ell(L0,a0,Delta)))
cat(sprintf("nu_W = %.2e, sqrt = %.2e   (paper: 1.09e-6, 1.04e-3)\n", nW0, sqrt(nW0)))
cat(sprintf("nu_{m,N} = %.2e, sqrt = %.2e   (paper: 4.47e-7, 6.69e-4)\n", nI0, sqrt(nI0)))
cat(sprintf("nu_X = %.2e, sqrt = %.2e   (paper: 1.53e-6, 1.24e-3)\n", nX0, sqrt(nX0)))
cat(sprintf("variance shares: pre %.1f%%, in %.1f%%   (paper: 70.9%%, 29.1%%)\n",
            100*nW0/nX0, 100*nI0/nX0))

cat("\n---- Remark on the leading-order form ----\n")
cat(sprintf("rel. error of LO vs exact ATM gap = %.1e   (paper: 1.6e-7)\n",
            abs(gap_atm_lo(nW0,nX0)/gap_atm_exact(nW0,nX0) - 1)))

cat("\n---- Remark on saturation ----\n")
nWinf <- nu_pre_inf(m0, a0, sig0, Delta)
cat(sprintf("nu_W^inf = %.3e, sqrt = %.3e\n", nWinf, sqrt(nWinf)))
cat(sprintf("limiting ATM gap = %.4f * p(n,m)\n", 2*pnorm(-0.5*sqrt(nWinf))))
cat(sprintf("nu_in(L=1e6) / [sigma^2 L/(a^2 Delta)] = %.4f   (should be approaching 1)\n",
            nu_in(1e6,a0,sig0,Delta) / (sig0^2*1e6/(a0^2*Delta))))
L1 <- a0^2*Delta/sig0^2
cat(sprintf("L such that nu_X ~ 1: %.2e days   (paper: order 1e6)\n", L1))
cat(sprintf("nu_X at that L = %.3f   (should be near 1)\n",
            nu_pre(m0,L1,a0,sig0,Delta) + nu_in(L1,a0,sig0,Delta)))
scl <- 1e-2/sqrt(nX0)   
vw_s <- nu_pre(m0,L0,a0,sig0*scl,Delta)
vx_s <- vw_s + nu_in(L0,a0,sig0*scl,Delta)
cat(sprintf("LO rel. error at sqrt(nu_X) = 1e-2: %.2e   (paper: around 1e-5)\n",
            abs(gap_atm_lo(vw_s,vx_s)/gap_atm_exact(vw_s,vx_s) - 1)))

cat("\n---- Base-case gap paragraph ----\n")
cat(sprintf("fr caplet = %.2f bp, br caplet = %.2f bp   (paper: 4.16, 4.94)\n",
            1e4*caplet_atm(nW0), 1e4*caplet_atm(nX0)))
cat(sprintf("ATM gap = %.2f bp   (paper: 0.78)\n", 1e4*gap_atm_exact(nW0,nX0)))
cat(sprintf("relative gap = %.1f%%   (paper: 18.8)\n", 100*gap_rel_exact(nW0,nX0)))
sHW <- seq(0.001, 0.025, length.out = 200)
vw_v <- sapply(sHW, function(s) nu_pre(m0,L0,a0,s*Delta,Delta))
vx_v <- vw_v + sapply(sHW, function(s) nu_in(L0,a0,s*Delta,Delta))
rel_v <- 100*gap_rel_exact(vw_v, vx_v)
cat(sprintf("relative gap over sigma sweep: min %.3f%%, max %.3f%%   (paper: nearly flat at 18.8)\n",
            min(rel_v), max(rel_v)))

cat("\n---- Remark on the maximum near the money ----\n")
cat(sprintf("0.5*sqrt(nu_X nu_W) = %.1e   (paper: 6.5e-7)\n", 0.5*sqrt(nX0*nW0)))
Kmx0 <- K_max(nW0, nX0, Katm0, alpha0)
cat(sprintf("K^max - K^ATM = %.3e bp   (paper: well under 0.1 bp)\n",
            1e4*(Kmx0 - Katm0)))

cat("\n---- Accrual-length paragraph ----\n")
relL <- function(L) 100*gap_rel_exact(nu_pre(m0,L,a0,sig0,Delta),
                                      nu_pre(m0,L,a0,sig0,Delta)+nu_in(L,a0,sig0,Delta))
cat(sprintf("relative gap at 1M = %.1f%%   (paper: about 6)\n", relL(30)))
cat(sprintf("relative gap at 3M = %.1f%%   (paper: 18.8)\n", relL(90)))
cat(sprintf("relative gap at 2Y = %.0f%%   (paper: roughly 133)\n", relL(720)))
L_cross <- uniroot(function(L) relL(L) - 100, c(200, 800))$root
cat(sprintf("relative gap crosses 100%% at %.2f months   (paper: near eighteen)\n",
            L_cross/30))


cat("\n---- Mean-reversion paragraph ----\n")
gap_a <- function(a) { 
  w <- nu_pre(m0,L0,a,sig0,Delta)
  1e4*gap_atm_exact(w, w + nu_in(L0,a,sig0,Delta)) 
  }
rel_a <- function(a) { 
  w <- nu_pre(m0,L0,a,sig0,Delta)
  100*gap_rel_exact(w, w + nu_in(L0,a,sig0,Delta)) 
  }

cat(sprintf("absolute gap: %.2f bp at a=0.05, %.2f bp at a=3   (paper: 0.76 to 0.82)\n",
            gap_a(0.05), gap_a(3)))
cat(sprintf("relative gap: %.0f%% at a=0.05, %.0f%% at a=3   (paper: 15 to 32)\n",
            rel_a(0.05), rel_a(3)))

