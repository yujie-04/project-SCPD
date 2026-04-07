rm(list = ls())
set.seed(2025)

M <- 5000; nu <- 50; G <- 50; D <- 50; B <- 50
S <- c(0.25, 0.5, 1, 3); sigma <- 1

ewma_params <- list(
  list(r=1.00,L=3.090), list(r=0.75,L=3.087), list(r=0.50,L=3.071),
  list(r=0.40,L=3.054), list(r=0.30,L=3.023), list(r=0.25,L=2.998),
  list(r=0.20,L=2.962), list(r=0.10,L=2.814), list(r=0.05,L=2.615),
  list(r=0.03,L=2.437))
aff_params <- list(list(alpha=0.005), list(alpha=0.008), list(alpha=0.010))
aff_eta <- 0.01; arl0_trials <- 1000; arl0_slen <- 50000

# Stream generation
generate_stream <- function() {
  xi <- rpois(M, nu)
  tau <- numeric(M); tau[1] <- G + xi[1]
  for(i in 2:M) tau[i] <- tau[i-1]+D+G+xi[i]
  n <- tau[M]+D+G+rpois(1,nu)
  mu <- numeric(M+1); mu[1] <- 0
  for(i in 2:(M+1)) mu[i] <- mu[i-1]+sample(c(-1,1),1)*sample(S,1)
  stream <- numeric(n)
  stream[1:tau[1]] <- rnorm(tau[1], mu[1], sigma)
  for(i in 1:M) {
    s <- tau[i]+1; e <- if(i<M) tau[i+1] else n
    if(s<=e) stream[s:e] <- rnorm(e-s+1, mu[i+1], sigma)
  }
  list(stream=stream, tau=tau, n=n, mu=mu)
}

# EWMA
run_ewma <- function(stream, r, L, B) {
  n <- length(stream); dets <- integer(0); pos <- 1
  while(pos+B <= n) {
    burn <- stream[pos:(pos+B-1)]
    mu_h <- mean(burn); sh <- sd(burn)
    if(is.na(sh)||sh<1e-10) sh <- 1
    Z <- mu_h; ms <- pos+B; found <- FALSE
    for(j in ms:n) {
      Z <- (1-r)*Z+r*stream[j]; t <- j-ms+1
      sZ <- sh*sqrt(r/(2-r)*(1-(1-r)^(2*t)))
      if(abs(Z-mu_h)>L*sZ) { dets <- c(dets,j); pos <- j+1; found <- TRUE; break }
    }
    if(!found) break
  }
  dets
}

ewma_arl0 <- function(r,L,B,slen) {
  s <- rnorm(slen); mu <- mean(s[1:B]); sh <- sd(s[1:B])
  if(is.na(sh)||sh<1e-10) sh <- 1
  Z <- mu
  for(j in (B+1):slen) {
    Z <- (1-r)*Z+r*s[j]; t <- j-B
    sZ <- sh*sqrt(r/(2-r)*(1-(1-r)^(2*t)))
    if(abs(Z-mu)>L*sZ) return(j)
  }
  slen
}

# AFF
run_aff <- function(stream, alpha, eta=0.01, B=50) {
  n <- length(stream); dets <- integer(0); pos <- 1; lmin <- 0.6
  zc <- qnorm(1-alpha/2)
  while(pos+B <= n) {
    burn <- stream[pos:(pos+B-1)]
    mu_h <- mean(burn); s2 <- var(burn)
    if(is.na(s2)||s2<1e-10) s2 <- 1
    m <- 0; w <- 0; lam <- 1; Dl <- 0; Om <- 0; u <- 1; xb <- 0; oc <- 0
    ms <- pos+B; found <- FALSE
    for(j in ms:n) {
      x <- stream[j]; oc <- oc+1
      if(oc>=2) {
        gl <- 2*(xb-x)*(Dl-xb*Om)/w
        lam <- max(min(lam-eta/s2*gl,1),lmin)
      }
      nD <- lam*Dl+m; nO <- lam*Om+w
      m <- lam*m+x; w <- lam*w+1; xb <- m/w
      if(oc==1) u <- 1 else { iw <- 1/w; u <- (1-iw)^2*u+iw^2 }
      Dl <- nD; Om <- nO
      if(abs(xb-mu_h) > zc*sqrt(u*s2)) { dets <- c(dets,j); pos <- j+1; found <- TRUE; break }
    }
    if(!found) break
  }
  dets
}

aff_arl0 <- function(alpha,eta=0.01,B=50,slen=50000) {
  s <- rnorm(slen); lmin <- 0.6; zc <- qnorm(1-alpha/2)
  mu <- mean(s[1:B]); s2 <- var(s[1:B])
  if(is.na(s2)||s2<1e-10) s2 <- 1
  m <- 0;w <- 0;lam <- 1;Dl <- 0;Om <- 0;u <- 1;xb <- 0;oc <- 0
  for(j in (B+1):slen) {
    x <- s[j]; oc <- oc+1
    if(oc>=2) { gl <- 2*(xb-x)*(Dl-xb*Om)/w; lam <- max(min(lam-eta/s2*gl,1),lmin) }
    nD <- lam*Dl+m; nO <- lam*Om+w; m <- lam*m+x; w <- lam*w+1; xb <- m/w
    if(oc==1) u <- 1 else { iw <- 1/w; u <- (1-iw)^2*u+iw^2 }
    Dl <- nD; Om <- nO
    if(abs(xb-mu)>zc*sqrt(u*s2)) return(j)
  }
  slen
}

# Classification
classify <- function(detections, true_cps, B) {
  nd <- length(detections); nt <- length(true_cps)
  if(nd==0) return(list(CCD=0, DNF=NA, ARL1=NA, SDRL1=NA))
  correct <- 0; arl1v <- numeric(0); prev <- 0
  for(i in 1:nd) {
    det <- detections[i]; be <- prev+B
    idx <- which(true_cps > be)
    if(length(idx)==0) { prev <- det; next }
    tau_m <- true_cps[idx[1]]
    if(det < tau_m) {
      # false
    } else {
      ac <- which(true_cps <= det & true_cps > be)
      if(length(ac)>0) {
        attr_cp <- true_cps[ac[length(ac)]]
        correct <- correct+1
        arl1v <- c(arl1v, det-attr_cp+1)
      }
    }
    prev <- det
  }
  list(CCD=correct/nt, DNF=if(nd>0) correct/nd else NA,
       ARL1=if(length(arl1v)>0) mean(arl1v) else NA,
       SDRL1=if(length(arl1v)>1) sd(arl1v) else NA)
}

# Run
data <- generate_stream()
stream <- data$stream; true_cps <- data$tau; mu_levels <- data$mu
cat(sprintf("Stream: %d obs, %d CPs\n\n", length(stream), length(true_cps)))

results <- data.frame(Algo=character(), Parameters=character(), Values=character(),
                      CCD=numeric(), DNF=numeric(), ARL1=numeric(), SDRL1=numeric(),
                      ARL0=numeric(), SDRL0=numeric(), stringsAsFactors=FALSE)

cat("Running EWMA ...\n")
for(k in seq_along(ewma_params)) {
  r <- ewma_params[[k]]$r; L <- ewma_params[[k]]$L
  cat(sprintf("  (r=%.2f, L=%.3f) ... ", r, L))
  dets <- run_ewma(stream, r, L, B)
  res <- classify(dets, true_cps, B)
  a0 <- replicate(arl0_trials, ewma_arl0(r,L,B,arl0_slen))
  cat(sprintf("done | CCD=%.2f DNF=%.2f\n", res$CCD, res$DNF))
  results <- rbind(results, data.frame(Algo="EWMA", Parameters="(r, L)",
                                       Values=sprintf("(%.2f, %.3f)",r,L),
                                       CCD=round(res$CCD,2), DNF=round(res$DNF,2),
                                       ARL1=round(res$ARL1,2), SDRL1=round(res$SDRL1,2),
                                       ARL0=round(mean(a0),2), SDRL0=round(sd(a0),2), stringsAsFactors=FALSE))
}

cat("\nRunning AFF ...\n")
for(k in seq_along(aff_params)) {
  alpha <- aff_params[[k]]$alpha
  cat(sprintf("  (alpha=%.3f) ... ", alpha))
  dets <- run_aff(stream, alpha, aff_eta, B)
  res <- classify(dets, true_cps, B)
  a0 <- replicate(arl0_trials, aff_arl0(alpha, aff_eta, B, arl0_slen))
  cat(sprintf("done | CCD=%.2f DNF=%.2f\n", res$CCD, res$DNF))
  results <- rbind(results, data.frame(Algo="AFF", Parameters="alpha",
                                       Values=sprintf("%.3f",alpha),
                                       CCD=round(res$CCD,2), DNF=round(res$DNF,2),
                                       ARL1=round(res$ARL1,2), SDRL1=round(res$SDRL1,2),
                                       ARL0=round(mean(a0),2), SDRL0=round(sd(a0),2), stringsAsFactors=FALSE))
}

# Print table
cat("\n", paste(rep("=",92),collapse=""), "\n")
cat(sprintf("%-6s %-12s %-18s %5s %5s %7s %9s %8s %10s\n",
            "Algo","Parameters","Values","CCD","DNF","ARL1","(SDRL1)","ARL0","(SDRL0)"))
cat(paste(rep("-",92),collapse=""), "\n")
for(i in 1:nrow(results)) {
  r <- results[i,]
  cat(sprintf("%-6s %-12s %-18s %5.2f %5.2f %7.2f (%7.2f) %8.2f (%8.2f)\n",
              r$Algo, r$Parameters, r$Values, r$CCD, r$DNF, r$ARL1, r$SDRL1, r$ARL0, r$SDRL0))
}

# Stream segments

tau <- true_cps
mu  <- mu_levels

ewma_seg <- function(seg, r, L, B) {
  ns <- length(seg); mu_h <- mean(seg[1:B]); sh <- sd(seg[1:B])
  if(is.na(sh) || sh < 1e-10) sh <- 1; Z <- mu_h
  for(j in (B+1):ns) {
    Z <- (1-r)*Z + r*seg[j]; t <- j - B
    sZ <- sh * sqrt(r/(2-r) * (1-(1-r)^(2*t)))
    if(abs(Z - mu_h) > L*sZ) return(j)
  }; NA
}

aff_seg <- function(seg, alpha, eta=0.01, B=50) {
  ns <- length(seg); mu_h <- mean(seg[1:B]); s2 <- var(seg[1:B])
  if(is.na(s2) || s2 < 1e-10) s2 <- 1
  zc <- qnorm(1-alpha/2); lmin <- 0.6
  m <- 0; w <- 0; lam <- 1; Dl <- 0; Om <- 0; u <- 1; xb <- 0; oc <- 0
  for(j in (B+1):ns) {
    x <- seg[j]; oc <- oc + 1
    if(oc >= 2) { gl <- 2*(xb-x)*(Dl-xb*Om)/w; lam <- max(min(lam-eta/s2*gl,1),lmin) }
    nD <- lam*Dl+m; nO <- lam*Om+w; m <- lam*m+x; w <- lam*w+1; xb <- m/w
    if(oc==1) u <- 1 else { iw <- 1/w; u <- (1-iw)^2*u+iw^2 }
    Dl <- nD; Om <- nO
    if(abs(xb-mu_h) > zc*sqrt(u*s2)) return(j)
  }; NA
}

# Select three non-overlapping 3CP windows
N_CP <- 3
win_indices <- c(1561, 3425, 4023)

wins <- list()
for(k in seq_along(win_indices)) {
  idx <- win_indices[k]
  cis <- idx:(idx + N_CP - 1)
  seg_s <- max(1, tau[idx] - B - G - 10)
  seg_e <- min(length(stream), tau[cis[N_CP]] + D + 30)
  wins[[k]] <- list(start = seg_s, end = seg_e, cps = cis)
}
w1 <- wins[[1]]; w2 <- wins[[2]]; w3 <- wins[[3]]

col_obs <- "grey70"; col_mu <- "black"; col_cp <- "#D32F2F"
col_burn <- rgb(0.26,0.65,0.96,0.15)
col_grace <- rgb(1,0.65,0.15,0.15)
col_detw <- rgb(0.4,0.73,0.42,0.15)
col_e1 <- "#1565C0"; col_e2 <- "#E65100"; col_aff <- "#6A1B9A"

plot_panel <- function(win, ttl, xlab_on=FALSE) {
  seg_s <- win$start; seg_e <- win$end; cis <- win$cps
  seg <- stream[seg_s:seg_e]; xv <- seg_s:seg_e
  yr <- range(seg) + c(-1,1)*diff(range(seg))*0.15
  
  plot(xv, seg, type='l', col=col_obs, lwd=0.4,
       xlab=if(xlab_on) expression("Observation index "*italic(t)) else "",
       ylab=expression(italic(x[t])), ylim=yr, main="",
       cex.axis=0.9, cex.lab=1.1)
  mtext(ttl, side=3, line=0.4, adj=0, font=2, cex=1.05)
  
  # Shaded regions (draw burn-in first)
  # Burn-in: [seg_s, seg_s + B]
  rect(seg_s, yr[1]-10, seg_s+B, yr[2]+10, col=col_burn, border=NA)
  
  for(k in seq_along(cis)) {
    ci <- cis[k]; cp <- tau[ci]; sv <- mu[ci+1]-mu[ci]
    
    # Grace period: [cp-G, cp]
    grace_left <- max(seg_s + B, cp - G)
    if(grace_left < cp) {
      rect(grace_left, yr[1]-10, cp, yr[2]+10, col=col_grace, border=NA)
    }
    
    # Detection window: [cp, cp+D]
    rect(cp, yr[1]-10, min(seg_e, cp+D), yr[2]+10, col=col_detw, border=NA)
    
    # CP line
    abline(v=cp, col=col_cp, lwd=2)
    
    # True mean before and after
    pre_s <- if(k==1) seg_s else tau[cis[k-1]]
    post_e <- if(k==length(cis)) seg_e else tau[cis[k+1]]
    segments(pre_s, mu[ci], cp, mu[ci], col=col_mu, lwd=2.5)
    segments(cp, mu[ci+1], post_e, mu[ci+1], col=col_mu, lwd=2.5)
    
    # Delta label
    lab_y <- if(sv > 0) mu[ci+1] + 0.5 else mu[ci+1] - 0.5
    text(cp+5, lab_y, bquote(delta==.(sprintf("%+.2f",sv))),
         cex=0.75, pos=4, col=col_cp, font=2)
  }
  # Final mean segment
  last_ci <- cis[length(cis)]
  segments(tau[last_ci], mu[last_ci+1], seg_e, mu[last_ci+1], col=col_mu, lwd=2.5)
  
  # Run detectors
  d1 <- ewma_seg(seg, 0.75, 3.087, B)
  d2 <- ewma_seg(seg, 0.10, 2.814, B)
  d3 <- aff_seg(seg, 0.005, 0.01, B)
  
  # Detection markers
  if(!is.na(d1)){dg<-seg_s+d1-1; abline(v=dg,col=col_e1,lwd=1.5,lty=2)
  points(dg,seg[d1],pch=25,col=col_e1,bg=col_e1,cex=1.5)}
  if(!is.na(d2)){dg<-seg_s+d2-1; abline(v=dg,col=col_e2,lwd=1.5,lty=4)
  points(dg,seg[d2],pch=22,col=col_e2,bg=col_e2,cex=1.5)}
  if(!is.na(d3)){dg<-seg_s+d3-1; abline(v=dg,col=col_aff,lwd=1.5,lty=3)
  points(dg,seg[d3],pch=24,col=col_aff,bg=col_aff,cex=1.5)}
  
  # Delay attribution
  mk <- function(nm, det) {
    if(is.na(det)) return(sprintf("%s: no detection", nm))
    dg <- seg_s + det - 1
    first_cp_pos <- tau[cis[1]]
    
    if(dg < first_cp_pos) {
      return(sprintf("%s: false alarm", nm))
    }
    
    # Attributed to last CP
    attr_ci <- NA
    for(ci in cis) {
      if(tau[ci] <= dg) attr_ci <- ci
    }
    delay <- dg - tau[attr_ci]
    in_w <- delay >= 0 && delay <= D
    sprintf("%s: delay = %d %s", nm, delay,
            if(in_w) "(in window)" else "(missed)")
  }
  
  legend("bottomright",
         legend=c(mk("EWMA r=0.75, L=3.087", d1),
                  mk("EWMA r=0.10, L=2.814", d2),
                  mk("AFF a=0.005, eta=0.01", d3)),
         text.col=c(col_e1, col_e2, col_aff), text.font=2,
         cex=0.75, bty="o", box.col="grey80", bg="white", inset=c(0.01,0.02))
}

png("~/Desktop/stream_segments.png", width=2200, height=2800, res=200)

par(mfrow=c(3,1), mar=c(4.5,5,3,1.5), oma=c(6,0,4,0), family="serif")

plot_panel(w1, "(a)")
plot_panel(w2, "(b)")
plot_panel(w3, "(c)", xlab_on=TRUE)

mtext("Stream segments with EWMA and AFF detection",
      outer=TRUE, side=3, line=2, cex=1.3, font=2)

par(fig=c(0,1,0,1), oma=c(0,0,0,0), mar=c(0,0,0,0), new=TRUE)
plot(0,0,type='n',bty='n',xaxt='n',yaxt='n')
legend("bottom", ncol=3, cex=1.0, bty="o", box.col="grey80",
       legend=c(expression(Observations~italic(x[t])),
                expression(True~mean~mu[t]), expression(True~changepoint~tau),
                sprintf("Burn-in (B=%d)",B), sprintf("Grace period (G=%d)",G),
                sprintf("Detection window (D=%d)",D),
                "EWMA (r=0.75, L=3.087)", "EWMA (r=0.10, L=2.814)",
                "AFF (a=0.005, eta=0.01)"),
       col=c(col_obs,col_mu,col_cp,col_burn,col_grace,col_detw,col_e1,col_e2,col_aff),
       lty=c(1,1,1,NA,NA,NA,2,4,3), lwd=c(0.6,2.5,2,NA,NA,NA,1.5,1.5,1.5),
       pch=c(NA,NA,NA,15,15,15,25,22,24), pt.cex=c(NA,NA,NA,2.2,2.2,2.2,1.5,1.5,1.5),
       pt.bg=c(NA,NA,NA,rgb(0.26,0.65,0.96,0.4),rgb(1,0.65,0.15,0.4),
               rgb(0.4,0.73,0.42,0.4),col_e1,col_e2,col_aff), xpd=TRUE)

dev.off()
cat("Saved to ~/Desktop/stream_segments.png\n")
