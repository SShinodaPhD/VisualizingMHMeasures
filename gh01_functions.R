library(ggplot2)
library(tidyr)

#--------------------------------------------------#
#-functions
##--Common
###---G_kl^(s,t)
fn.mat.G <- function(mat.p, s, t){
	R <- nrow(mat.p)
	mat.G <- matrix(0,3,3)

	for(i in 1:s){
		for(j in 1:s){
			mat.G[1,1] <- mat.G[1,1] + mat.p[i,j]
		}
		for(j in (s+1):t){
			mat.G[1,2] <- mat.G[1,2] + mat.p[i,j]
		}
		for(j in (t+1):R){
			mat.G[1,3] <- mat.G[1,3] + mat.p[i,j]
		}
	}
	for(i in (s+1):t){
		for(j in 1:s){
			mat.G[2,1] <- mat.G[2,1] + mat.p[i,j]
		}
		for(j in (s+1):t){
			mat.G[2,2] <- mat.G[2,2] + mat.p[i,j]
		}
		for(j in (t+1):R){
			mat.G[2,3] <- mat.G[2,3] + mat.p[i,j]
		}
	}
	for(i in (t+1):R){
		for(j in 1:s){
			mat.G[3,1] <- mat.G[3,1] + mat.p[i,j]
		}
		for(j in (s+1):t){
			mat.G[3,2] <- mat.G[3,2] + mat.p[i,j]
		}
		for(j in (t+1):R){
			mat.G[3,3] <- mat.G[3,3] + mat.p[i,j]
		}
	}
	return(mat.G)
}

###---H_1(k)^(s,t)
fn.H1 <- function(mat.G){
	H1 <- numeric(2)
	for(k in 1:2){
		for(i in 1:k){
			for(j in (k+1):3){
				H1[k] <- H1[k] + mat.G[i,j]
			}
		}
	}
	return(H1)
}

###---H_2(k)^(s,t)
fn.H2 <- function(mat.G){
	H2 <- numeric(2)
	for(k in 1:2){
		for(i in (k+1):3){
			for(j in 1:k){
				H2[k] <- H2[k] + mat.G[i,j]
			}
		}
	}
	return(H2)
}

###---Δ^(s,t)
fn.Delta <- function(mat.G){
	H1 <- fn.H1(mat.G); H2 <- fn.H2(mat.G)
	Delta <- (H1[1] + H2[1]) + (H1[2] + H2[2])
	return(Delta)
}

###---H_1(k)^(s,t)*
fn.H1a <- function(mat.G){
	H1 <- fn.H1(mat.G); Delta <- fn.Delta(mat.G)
	H1a <- H1/Delta
	return(H1a)
}

###---H_2(k)^(s,t)*
fn.H2a <- function(mat.G){
	H2 <- fn.H2(mat.G); Delta <- fn.Delta(mat.G)
	H2a <- H2/Delta
	return(H2a)
}

###---Q_(k)^(s,t)*
fn.Qa <- function(mat.G){
	H1a <- fn.H1a(mat.G); H2a <- fn.H2a(mat.G)
	Qa <- (H1a + H2a)/2
	return(Qa)
}

###---H_1(i)c
fn.H1c <- function(mat.G){
	H1 <- fn.H1(mat.G); H2 <- fn.H2(mat.G)
	H1c <- H1/(H1 + H2)
	return(H1c)
}

###---H_2(i)c
fn.H2c <- function(mat.G){
	H1 <- fn.H1(mat.G); H2 <- fn.H2(mat.G)
	H2c <- H2/(H1 + H2)
	return(H2c)
}

###---Δ_1^(s,t)*
fn.Delta1a <- function(mat.G){
	H1a <- fn.H1a(mat.G)
	Delta1a <- H1a[1] + H1a[2]
	return(Delta1a)
}

###---Δ_2^(s,t)*
fn.Delta2a <- function(mat.G){
	H2a <- fn.H2a(mat.G)
	Delta2a <- H2a[1] + H2a[2]
	return(Delta2a)
}

###---U_1(k)^(s,t)*
fn.U1a <- function(mat.G){
	H1a <- fn.H1a(mat.G); H2a <- fn.H2a(mat.G); Delta1a <- fn.Delta1a(mat.G)
	U1a <- Delta1a*(H1a + H2a)
	return(U1a)
}

###---U_2(k)^(s,t)*
fn.U2a <- function(mat.G){
	H1a <- fn.H1a(mat.G); H2a <- fn.H2a(mat.G); Delta2a <- fn.Delta2a(mat.G)
	U2a <- Delta2a*(H1a + H2a)
	return(U2a)
}

mylog <- function(x){
	if(x==0){
		y <- 0	
	}else{
		y <- log(x)
	}
	return(y)
}

##--Measures
###---Measure of departure from MH using collapsed table
fn.mCoMH.st <- function(mat.p, s, t){
	mat.G <- fn.mat.G(mat.p, s, t)
	H1a <- fn.H1a(mat.G); H2a <- fn.H2a(mat.G); Qa <- fn.Qa(mat.G)
	H1c <- fn.H1c(mat.G); H2c <- fn.H2c(mat.G)

	Phi.CoMH.st <- 0
	for(k in 1:2){
		Phi.CoMH.st <- Phi.CoMH.st + H1a[k]*mylog(H1c[k]) + H2a[k]*mylog(H2c[k])
	}
	Phi.CoMH.st <- Phi.CoMH.st / log(2) + 1
	return(Phi.CoMH.st)
}
fn.mCoMH <- function(mat){
	mat.p <- mat/sum(mat); R <- nrow(mat)

	Phi.CoMH <- 0
	for(s in 1:(R-2)){
		for(t in (s+1):(R-1)){
			Phi.CoMH <- Phi.CoMH + fn.mCoMH.st(mat.p, s, t)
		}
	}
	Phi.CoMH <- Phi.CoMH/choose((R-1),2)
	return(Phi.CoMH)
}
fn.ind1 <- function(k, l, s, t, R){
	res <- 0
	if( 1 <= k&k <= s & s+1 <= l&l <= t )	res <- 1
	return(res)
}
fn.ind2 <- function(k, l, s, t, R){
	res <- 0
	if( 1 <= k&k <= s & t+1 <= l&l <= R )	res <- 1
	return(res)
}
fn.ind3 <- function(k, l, s, t, R){
	res <- 0
	if( s+1 <= k&k <= t & t+1 <= l&l <= R )	res <- 1
	return(res)
}
fn.ind4 <- function(k, l, s, t, R){
	res <- 0
	if( s+1 <= k&k <= t & 1 <= l&l <= s )	res <- 1
	return(res)
}
fn.ind5 <- function(k, l, s, t, R){
	res <- 0
	if( t+1 <= k&k <= R & 1 <= l&l <= s )	res <- 1
	return(res)
}
fn.ind6 <- function(k, l, s, t, R){
	res <- 0
	if( t+1 <= k&k <= R & s+1 <= l&l <= t )	res <- 1
	return(res)
}
fn.Dkl.CoMH <- function(mat.p, k, l, R){
	Dkl <- 0
	for(s in 1:(R-2)){
		for(t in (s+1):(R-1)){
			mat.G <- fn.mat.G(mat.p, s, t)
			Delta <- fn.Delta(mat.G); H1c <- fn.H1c(mat.G); H2c <- fn.H2c(mat.G)
			mCoMH.st <- fn.mCoMH.st(mat.p, s, t)

			Dkl <- Dkl + fn.ind1(k,l,s,t,R)/log(2)/Delta*(
					mylog(H1c[1]) - (mCoMH.st - 1)*log(2)
				) + fn.ind2(k,l,s,t,R)/log(2)/Delta*(
					mylog(H1c[1]) + mylog(H1c[2]) - 2*(mCoMH.st - 1)*log(2)
				) + fn.ind3(k,l,s,t,R)/log(2)/Delta*(
					mylog(H1c[2]) - (mCoMH.st - 1)*log(2)
				)
		}
	}
	Dkl <- Dkl/choose((R-1),2)
	return(Dkl)
}
fn.Dlk.CoMH <- function(mat.p, k, l, R){
	Dlk <- 0
	for(s in 1:(R-2)){
		for(t in (s+1):(R-1)){
			mat.G <- fn.mat.G(mat.p, s, t)
			Delta <- fn.Delta(mat.G); H1c <- fn.H1c(mat.G); H2c <- fn.H2c(mat.G)
			mCoMH.st <- fn.mCoMH.st(mat.p, s, t)

			Dlk <- Dlk + fn.ind4(k,l,s,t,R)/log(2)/Delta*(
					mylog(H2c[1]) - (mCoMH.st - 1)*log(2)
				) + fn.ind5(k,l,s,t,R)/log(2)/Delta*(
					mylog(H2c[1]) + mylog(H2c[2]) - 2*(mCoMH.st - 1)*log(2)
				) + fn.ind6(k,l,s,t,R)/log(2)/Delta*(
					mylog(H2c[2]) - (mCoMH.st - 1)*log(2)
				)
		}
	}
	Dlk <- Dlk/choose((R-1),2)
	return(Dlk)
}
fn.mCoMH.sammary <- function(mat){
	Phi <- fn.mCoMH(mat)
	mat.p <- mat/sum(mat); R <- nrow(mat); N <- sum(mat)

	Phi.var1 <- Phi.var2 <- 0
	for(k in 1:(R-1)){
		for(l in (k+1):R){
			Phi.var1 <- Phi.var1 + (fn.Dkl.CoMH(mat.p, k, l, R))^2*mat.p[k,l] + (fn.Dlk.CoMH(mat.p, l, k, R))^2*mat.p[l,k]
			Phi.var2 <- Phi.var2 + fn.Dkl.CoMH(mat.p, k, l, R)*mat.p[k,l] + fn.Dlk.CoMH(mat.p, l, k, R)*mat.p[l,k]
		}
	}
	Phi.var <- Phi.var1 - (Phi.var2)^2

	Phi.se <- sqrt(Phi.var/N)
	Phi.CIl <- Phi - qnorm(0.975)*Phi.se
	Phi.CIu <- Phi + qnorm(0.975)*Phi.se
	return(round(c(Phi, Phi.se, Phi.CIl, Phi.CIu),5))
}

###---Measure of departure from CoME using collapsed table
fn.mCoME.st <- function(mat.p, s, t){
	mat.G <- fn.mat.G(mat.p, s, t)
	Delta1a <- fn.Delta1a(mat.G); Delta2a <- fn.Delta2a(mat.G)

	Phi.CoME.st <- ( Delta1a*mylog(Delta1a*2) + Delta2a*mylog(Delta2a*2) ) / log(2)
	return(Phi.CoME.st)
}
fn.mCoME <- function(mat){
	R <- nrow(mat); mat.p <- mat/sum(mat)
	Phi.CoME <- 0
	for(s in 1:(R-2)){
		for(t in (s+1):(R-1)){
			Phi.CoME <- Phi.CoME + fn.mCoME.st(mat.p, s, t)
		}
	}
	Phi.CoME <- Phi.CoME/choose((R-1),2)
	return(Phi.CoME)
}
fn.Dkl.CoME <- function(mat.p, k, l, R){
	Dkl <- 0
	for(s in 1:(R-2)){
		for(t in (s+1):(R-1)){
			mat.G <- fn.mat.G(mat.p, s, t)
			Delta <- fn.Delta(mat.G); Delta1a <- fn.Delta1a(mat.G); Delta2a <- fn.Delta2a(mat.G)
			mCoME.st <- fn.mCoME.st(mat.p, s, t)

			Dkl <- Dkl + fn.ind1(k,l,s,t,R)/log(2)/Delta*(
					mylog(2*Delta1a) - mCoME.st*log(2)
				) + fn.ind2(k,l,s,t,R)/log(2)/Delta*(
					2*mylog(2*Delta1a) - 2*mCoME.st*log(2)
				) + fn.ind3(k,l,s,t,R)/log(2)/Delta*(
					mylog(2*Delta1a) - mCoME.st*log(2)
				)
		}
	}
	Dkl <- Dkl/choose((R-1),2)
	return(Dkl)
}
fn.Dlk.CoME <- function(mat.p, k, l, R){
	Dlk <- 0
	for(s in 1:(R-2)){
		for(t in (s+1):(R-1)){
			mat.G <- fn.mat.G(mat.p, s, t)
			Delta <- fn.Delta(mat.G); Delta1a <- fn.Delta1a(mat.G); Delta2a <- fn.Delta2a(mat.G)
			mCoME.st <- fn.mCoME.st(mat.p, s, t)

			Dlk <- Dlk + fn.ind4(k,l,s,t,R)/log(2)/Delta*(
					mylog(2*Delta2a) - mCoME.st*log(2)
				) + fn.ind5(k,l,s,t,R)/log(2)/Delta*(
					2*mylog(2*Delta2a) - 2*mCoME.st*log(2)
				) + fn.ind6(k,l,s,t,R)/log(2)/Delta*(
					mylog(2*Delta2a) - mCoME.st*log(2)
				)
		}
	}
	Dlk <- Dlk/choose((R-1),2)
	return(Dlk)
}
fn.mCoME.sammary <- function(mat){
	Phi <- fn.mCoME(mat)
	mat.p <- mat/sum(mat); R <- nrow(mat); N <- sum(mat)

	Phi.var1 <- Phi.var2 <- 0
	for(k in 1:(R-1)){
		for(l in (k+1):R){
			Phi.var1 <- Phi.var1 + (fn.Dkl.CoME(mat.p, k, l, R))^2*mat.p[k,l] + (fn.Dlk.CoME(mat.p, l, k, R))^2*mat.p[l,k]
			Phi.var2 <- Phi.var2 + fn.Dkl.CoME(mat.p, k, l, R)*mat.p[k,l] + fn.Dlk.CoME(mat.p, l, k, R)*mat.p[l,k]
		}
	}
	Phi.var <- Phi.var1 - (Phi.var2)^2

	Phi.se <- sqrt(Phi.var/N)
	Phi.CIl <- Phi - qnorm(0.975)*Phi.se
	Phi.CIu <- Phi + qnorm(0.975)*Phi.se
	return(round(c(Phi, Phi.se, Phi.CIl, Phi.CIu),5))
}

###---Measure of departure from EMH using collapsed table
fn.mCoEMH.st <- function(mat.p, s, t){
	mat.G <- fn.mat.G(mat.p, s, t)
	H1a <- fn.H1a(mat.G); H2a <- fn.H2a(mat.G); U1a <- fn.U1a(mat.G); U2a <- fn.U2a(mat.G)

	Phi.CoEMH.st <- 0
	for(k in 1:2){
		Phi.CoEMH.st <- Phi.CoEMH.st + H1a[k]*mylog(H1a[k]/U1a[k]) + H2a[k]*mylog(H2a[k]/U2a[k])
	}
	Phi.CoEMH.st <- Phi.CoEMH.st / log(2)
	return(Phi.CoEMH.st)
}
fn.mCoEMH <- function(mat){
	R <- nrow(mat); mat.p <- mat/sum(mat)
	Phi.CoEMH <- 0
	for(s in 1:(R-2)){
		for(t in (s+1):(R-1)){
			Phi.CoEMH <- Phi.CoEMH + fn.mCoEMH.st(mat.p, s, t)
		}
	}
	Phi.CoEMH <- Phi.CoEMH/choose((R-1),2)
	return(Phi.CoEMH)
}
fn.Dkl.CoEMH <- function(mat.p, k, l, R){
	Dkl <- 0
	for(s in 1:(R-2)){
		for(t in (s+1):(R-1)){
			mat.G <- fn.mat.G(mat.p, s, t)
			Delta <- fn.Delta(mat.G); H1a <- fn.H1a(mat.G); H2a <- fn.H2a(mat.G)
			U1a <- fn.U1a(mat.G); U2a <- fn.U2a(mat.G)
			mCoEMH.st <- fn.mCoEMH.st(mat.p, s, t)

			Dkl <- Dkl + fn.ind1(k,l,s,t,R)/log(2)/Delta*(
					mylog(H1a[1]/U1a[1]) - mCoEMH.st*log(2)
				) + fn.ind2(k,l,s,t,R)/log(2)/Delta*(
					mylog(H1a[1]/U1a[1]) + mylog(H1a[2]/U1a[2]) - 2*mCoEMH.st*log(2)
				) + fn.ind3(k,l,s,t,R)/log(2)/Delta*(
					mylog(H1a[2]/U1a[2]) - mCoEMH.st*log(2)
				)
		}
	}
	Dkl <- Dkl/choose((R-1),2)
	return(Dkl)
}
fn.Dlk.CoEMH <- function(mat.p, k, l, R){
	Dlk <- 0
	for(s in 1:(R-2)){
		for(t in (s+1):(R-1)){
			mat.G <- fn.mat.G(mat.p, s, t)
			Delta <- fn.Delta(mat.G); H1a <- fn.H1a(mat.G); H2a <- fn.H2a(mat.G)
			U1a <- fn.U1a(mat.G); U2a <- fn.U2a(mat.G)
			mCoEMH.st <- fn.mCoEMH.st(mat.p, s, t)

			Dlk <- Dlk + fn.ind4(k,l,s,t,R)/log(2)/Delta*(
					mylog(H2a[1]/U2a[1]) - mCoEMH.st*log(2)
				) + fn.ind5(k,l,s,t,R)/log(2)/Delta*(
					mylog(H2a[1]/U2a[1]) + mylog(H2a[2]/U2a[2]) - 2*mCoEMH.st*log(2)
				) + fn.ind6(k,l,s,t,R)/log(2)/Delta*(
					mylog(H2a[2]/U2a[2]) - mCoEMH.st*log(2)
				)
		}
	}
	Dlk <- Dlk/choose((R-1),2)
	return(Dlk)
}
fn.mCoEMH.sammary <- function(mat){
	Phi <- fn.mCoEMH(mat)
	mat.p <- mat/sum(mat); R <- nrow(mat); N <- sum(mat)

	Phi.var1 <- Phi.var2 <- 0
	for(k in 1:(R-1)){
		for(l in (k+1):R){
			Phi.var1 <- Phi.var1 + (fn.Dkl.CoEMH(mat.p, k, l, R))^2*mat.p[k,l] + (fn.Dlk.CoEMH(mat.p, l, k, R))^2*mat.p[l,k]
			Phi.var2 <- Phi.var2 + fn.Dkl.CoEMH(mat.p, k, l, R)*mat.p[k,l] + fn.Dlk.CoEMH(mat.p, l, k, R)*mat.p[l,k]
		}
	}
	Phi.var <- Phi.var1 - (Phi.var2)^2

	Phi.se <- sqrt(Phi.var/N)
	Phi.CIl <- Phi - qnorm(0.975)*Phi.se
	Phi.CIu <- Phi + qnorm(0.975)*Phi.se
	return(round(c(Phi, Phi.se, Phi.CIl, Phi.CIu),5))
}


##--visualization
fn.visualization.Fig2 <- function(mat){
	mat.p <- mat/sum(mat); R <- nrow(mat)

	pt <- choose((R-1),2)
	tmp01 <- matrix(,3,pt+1); cnames <- character(pt)
	Phi.CoEMH <- Phi.CoME <- Phi.CoMH <- 0; i <- 1
	for(s in 1:(R-2)){
		for(t in (s+1):(R-1)){
			Phi.CoEMHst <- fn.mCoEMH.st(mat.p, s, t)
			Phi.CoMEst <- fn.mCoME.st(mat.p, s, t)
			Phi.CoMHst <- fn.mCoMH.st(mat.p, s, t)

			tmp01[,i] <- c(Phi.CoEMHst, Phi.CoMEst, Phi.CoMHst)
			cnames[i] <- paste0("(", s, ",", t, ")")
			i<-i+1

			Phi.CoEMH <- Phi.CoEMH + Phi.CoEMHst
			Phi.CoME <- Phi.CoME + Phi.CoMEst
			Phi.CoMH <- Phi.CoMH + Phi.CoMHst
		}
	}
	Phi.CoEMH <- Phi.CoEMH/pt; Phi.CoME <- Phi.CoME/pt; Phi.CoMH <- Phi.CoMH/pt
	tmp01[,pt+1] <- c(Phi.CoEMH, Phi.CoME, Phi.CoMH); tmp01 <- round(tmp01,3)
	colnames(tmp01) <- c(cnames, "Total"); rownames(tmp01) <- c("EMH","CoME","MH")

	tmp02 <- as.data.frame(t(tmp01[c(2,1),]))
	tmp02$Pattern <- c(cnames, "Total")
	tmp02 <- gather(tmp02, key="Measure", value="Estimate", -Pattern)
	tmp02$Measure <- factor(tmp02$Measure, levels = c("EMH","CoME"))

	tmp03 <- ggplot(tmp02, aes(x=Pattern, y=Estimate, fill=Measure)) + geom_bar(stat="identity")+
	labs(title = "Visualization of Collapsed Measures", x="Collapsed Pattern", y="Values of Measures") +
	scale_y_continuous(limits = c(0, 1)) +
 	geom_vline(xintercept=length(cnames)+0.5, linetype="dashed", color="black", size=1) + 
	theme(legend.title=element_blank())

	print(tmp03); return(tmp01)
}

fn.visualization.Fig3 <- function(mat){
	mat.p <- mat/sum(mat); R <- nrow(mat)

	pt <- choose((R-1),2)
	tmp01 <- matrix(,3,pt+1); cnames <- character(pt)
	Phi.CoEMH <- Phi.CoME <- Phi.CoMH <- 0; i <- 1
	for(s in 1:(R-2)){
		for(t in (s+1):(R-1)){
			Phi.CoEMHst <- fn.mCoEMH.st(mat.p, s, t)
			Phi.CoMEst <- fn.mCoME.st(mat.p, s, t)
			Phi.CoMHst <- fn.mCoMH.st(mat.p, s, t)

			tmp01[,i] <- c(Phi.CoEMHst, Phi.CoMEst, Phi.CoMHst)
			cnames[i] <- paste0("(", s, ",", t, ")")
			i<-i+1

			Phi.CoEMH <- Phi.CoEMH + Phi.CoEMHst
			Phi.CoME <- Phi.CoME + Phi.CoMEst
			Phi.CoMH <- Phi.CoMH + Phi.CoMHst
		}
	}
	Phi.CoEMH <- Phi.CoEMH/pt; Phi.CoME <- Phi.CoME/pt; Phi.CoMH <- Phi.CoMH/pt
	tmp01[,pt+1] <- c(Phi.CoEMH, Phi.CoME, Phi.CoMH); tmp01 <- round(tmp01,3)
	colnames(tmp01) <- c(cnames, "Total"); rownames(tmp01) <- c("EMH","CoME","MH")

	tmp02 <- as.data.frame(t(tmp01[c(2,1),]))
	tmp02$Pattern <- c(cnames, "Total")
	tmp02 <- gather(tmp02, key="Measure", value="Estimate", -Pattern)
	tmp02$Measure <- factor(tmp02$Measure, levels = c("EMH","CoME"))

	tmp03 <- ggplot(tmp02, aes(x=Pattern, y=Estimate, fill=Measure)) + geom_bar(stat="identity")+
	labs(title = "Visualization of Collapsed Measures", x="Collapsed Pattern", y="Estimates of Measures") +
	scale_y_continuous(limits = c(0, 1)) +
 	geom_vline(xintercept=length(cnames)+0.5, linetype="dashed", color="black", size=1) + 
	theme(legend.title=element_blank())

	print(tmp03); return(tmp01)
}
