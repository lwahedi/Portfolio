#Vector Autoregression Practice
#Laila Wahedi
library(MSBVAR)
library(foreign)
library(fUnitRoots)
library(plotrix)
#Read in west german financial data (e1)
wgf<-read.dta(file.choose())
wgf<-ts(wgf,start=c(1960,1),frequency=4)
plot(wgf)

#Unit root tests
#pick lag lengths for Dicky Fuller (df) test
write("Pick Lags for df",file="unitRootTests.txt",append=FALSE)
for(i in 5:7){
write(paste0("var:  ",toString(i)),file="unitRootTests.txt",append=TRUE)
for(j in 0:17){
	v<-wgf[,i]
	#For differenced data:
	#v<-diff(wgf[,i])
	bmp(file=paste0("q2-",toString(i),".bmp"),height=1275, width=1650,res=150)
	capture.output(
		urdfTest(v,lags=j,type=c("nc")),
		file="unitRootTests.txt",append=TRUE)
	urdfTest(v,lags=10,type=c("nc"))
	dev.off()
}
write("______________________",file="unitRootTests.txt",append=TRUE)
}

l=c(4,3,4)

write("Run Tests",file="unitRootTests.txt",append=TRUE)
for(i in 5:7){
write(paste0("var:  ",toString(i)),file="unitRootTests.txt",append=TRUE)
	v<-wgf[,i]
	#v<-wgf[,i]
	#dicky fuller
	bmp(file=paste0("q2-",toString(i),".bmp"),height=1275, width=1650,res=150)
	capture.output(
		urdfTest(v,lags=l[i-4],type=c("ct")),
		file="unitRootTests.txt",append=TRUE)
	dev.off()
	capture.output(
		urdfTest(v,lags=l[i-4],type=c("c")),
		file="unitRootTests.txt",append=TRUE)
	capture.output(
		urdfTest(v,lags=l[i-4],type=c("nc")),
		file="unitRootTests.txt",append=TRUE)
	#kpss long
	capture.output(
		urkpssTest(v,lags=c("long"),type=c("mu")),
		file="unitRootTests.txt",append=TRUE)
	capture.output(
		urkpssTest(v,lags=c("long"),type=c("tau")),
		file="unitRootTests.txt",append=TRUE)

	#kpss Short
	capture.output(
		urkpssTest(v,lags=c("short"),type=c("mu")),
		file="unitRootTests.txt",append=TRUE)
	capture.output(
		urkpssTest(v,lags=c("short"),type=c("tau")),
		file="unitRootTests.txt",append=TRUE)

		#phillips perron
	capture.output(
		PP.test(v,lshort=FALSE),
		file="unitRootTests.txt",append=TRUE)
	capture.output(
		PP.test(v,lshort=TRUE),
		file="unitRootTests.txt",append=TRUE)
write("______________________",file="unitRootTests.txt",append=TRUE)
}

#Var lag lengths
#start with a quarter of df. ((T/4)-1)/3= (23-1)/3=7.3. So start with around 7 lags 
#This is about right because we have quarterly data. I'll start with 8, for 2 years.
a<-var.lag.specification(wgf[,c(5,6,7)],lagmax=8)
capture.output(
	cbind(apply(a$ldets,2,rev),a$results[,-1])
	file="LagLengths.txt",append=TRUE)
mVAR<-reduced.form.var(wgf[,c(5,6,7)],3)
#Check for stability
eigen(mVAR$ar.coefs[1,,])$values
#[1]  0.64107897+0.00000000i -0.03198491+0.03140425i -0.03198491-0.03140425i
eigen(mVAR$ar.coefs[2,,])$values
#[1]  0.2386306+0.5948892i  0.2386306-0.5948892i -0.2469792+0.0000000i
eigen(mVAR$ar.coefs[3,,])$values
#[1]  0.5957801+0.3639544i  0.5957801-0.3639544i -0.3201906+0.0000000i
Box.test(mVAR$residuals[,1],40,"Ljung")
	#Fail to reject null, conclude white noise residuals
Box.test(mVAR$residuals[,2],40,"Ljung")
	#Fail to reject null, conclude white noise residuals
Box.test(mVAR$residuals[,3],40,"Ljung")
	#Fail to reject null, conclude white noise residuals
evals<-matrix(c( 
		0.64107897,0.00000000,
		-0.03198491,0.03140425,
		-0.03198491,-0.03140425,
		0.2386306,0.5948892, 
		0.2386306,-0.5948892, 
		-0.2469792,0.0000000,
		0.5957801,0.3639544,  
		0.5957801,-0.3639544, 
		-0.3201906,0.0000000
		),byrow=TRUE,nrow=9)
plot(evals,main="eigenvalues",ylab="imaginary",xlab="real",
	xlim=c(-1.1,1.1),ylim=c(-1.1,1.1))
abline(h=0)
abline(v=0)
draw.circle(0,0,1)
#everything lies within the unit circle

##with different ordering
mVAR<-reduced.form.var(wgf[,c(6,7,5)],3)
> eigen(mVAR$ar.coefs[1,,])$values
#Output:[1]  0.9489440 -0.5835523  0.1832358
> eigen(mVAR$ar.coefs[2,,])$values
#Output:[1] -0.6952129+0.0000000i  0.1595768+0.4461191i  0.1595768-0.4461191i
> eigen(mVAR$ar.coefs[3,,])$values
#Output:[1] 0.0910072+0.2160103i 0.0910072-0.2160103i 0.0234443+0.0000000i

evals<-matrix(c( 
	0.9489440,0, 
	-0.5835523,0,
	0.1832358,0,
	-0.6952129,0.0000000,
	0.1595768,0.4461191,
	0.1595768,-0.4461191,
	0.0910072,0.2160103,
	0.0910072,-0.2160103,
	0.0234443,0.0000000),byrow=TRUE,nrow=9)
plot(evals,main="eigenvalues",ylab="imaginary",xlab="real",
	xlim=c(-1.1,1.1),ylim=c(-1.1,1.1))
abline(h=0)
abline(v=0)
draw.circle(0,0,1)

#Granger Test
mGC<-granger.test(wgf[,c(5,6,7)],3)
#IRF
#Computed
mIRF<-irf(mVAR, 12)
mIRF
plot.irf(mIRF)
#Monte Carlo Estimation
mMCIRF<-mc.irf(mVAR, nsteps=6,draws=10000)
plot.mc.irf(mMCIRF,method=c("Normal Approximation"))
plot.mc.irf(mMCIRF,method=c("Sims-Zha1"))
#Boot Strapped Errors
library(VARS)
a<-VAR(wgf[,c(5,6,7)],3)
plot(irf(a))
#forecast error decomposition, 3 years
dfev(mVAR,A0=NULL,12)
print(dfev(mVAR,A0=NULL,12))

#bivariate regression between ln_consump (dv) and ln_inc (iv).
mLM<- lm(wgf[,"ln_consump"]~wgf[,"ln_inc"])
mLMR<-residuals(mLM)
PP.test(mLMR)
#reject null, no unit root. They are cointegrated which means 
#that an error correction model is probably more appropriate

#Testing the MSBVAR package
lagpad <- function(x, k) {
    c(rep(NA, k), x)[1 : length(x)] 
}
coefficients(lm(wgf[,7]~lagpad(wgf[,5],1)+lagpad(wgf[,5],2)+lagpad(wgf[,5],3)+
						lagpad(wgf[,6],1)+lagpad(wgf[,6],2)+lagpad(wgf[,6],3)+
						lagpad(wgf[,7],1)+lagpad(wgf[,7],2)+lagpad(wgf[,7],3)))