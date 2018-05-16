*linear regression substantive model with quadratic effects
clear
set obs 10000
gen x=rnormal()
gen y=x+x^2+rnormal()

summ y
gen missxb=(y-r(mean))/r(sd)
gen misspr=exp(missxb)/(1+exp(missxb))
replace x=. if runiform()<misspr

gen xsq=x^2
smcfcs regress y x xsq, reg(x) passive(xsq=x^2)


*cox substantive model with quadratic effects
clear
set obs 10000
gen x=rnormal()
gen t=-log(runiform())/exp(x)

summ t
gen missxb=(t-r(mean))/r(sd)
gen misspr=exp(missxb)/(1+exp(missxb))
replace x=. if runiform()<misspr

stset t

gen xsq=x^2
smcfcs stcox x xsq, reg(x) passive(xsq=x^2)

*competing risks substantive model with linear effect
clear
set obs 10000
gen x=rnormal()
gen t1=-log(runiform())/exp(x)
gen t2=-log(runiform())/exp(-0.5*x)

gen t=t1
replace t=t2 if t2<t1
gen d=1
replace d=2 if t2<t1

summ t
gen missxb=(t-r(mean))/r(sd)
gen misspr=exp(missxb)/(1+exp(missxb))
replace x=. if runiform()<misspr

smcfcs compet x, reg(x) time(t) failure(d)

*competing risks substantive model with linear effect
*testing savetrace
clear
set obs 10000
gen x=rnormal()
gen t1=-log(runiform())/exp(x)
gen t2=-log(runiform())/exp(-0.5*x)

gen t=t1
replace t=t2 if t2<t1
gen d=1
replace d=2 if t2<t1

summ t
gen missxb=(t-r(mean))/r(sd)
gen misspr=exp(missxb)/(1+exp(missxb))
replace x=. if runiform()<misspr

smcfcs compet x, reg(x) time(t) failure(d) savetrace(temp)
use temp, clear

*competing risks substantive model with binary covariate
clear
set obs 10000
gen x=(runiform()<0.5)
gen t1=-log(runiform())/exp(x)
gen t2=-log(runiform())/exp(-0.5*x)

gen t=t1
replace t=t2 if t2<t1
gen d=1
replace d=2 if t2<t1

summ t
gen missxb=(t-r(mean))/r(sd)
gen misspr=exp(missxb)/(1+exp(missxb))
replace x=. if runiform()<misspr

smcfcs compet x, logit(x) time(t) failure(d)

*competing risks substantive model with normal covariate, delayed entry
clear
set obs 10000
gen x=rnormal()
gen t1=-log(runiform())/exp(x)
gen t2=-log(runiform())/exp(-0.5*x)

gen t=t1
replace t=t2 if t2<t1
gen d=1
replace d=2 if t2<t1

*generate entry time
gen entry=-log(runiform())/10
drop if t<entry
stset t, failure(d==1) enter(entry)
stcox x, nohr
stset t, failure(d==2) enter(entry)
stcox x, nohr

summ t
gen missxb=(t-r(mean))/r(sd)
gen misspr=exp(missxb)/(1+exp(missxb))
replace x=. if runiform()<misspr

smcfcs compet x, reg(x) time(t) failure(d) enter(entry)

*Poisson substantive model
*first with binary missing covariate
clear
set obs 10000
gen z=rnormal()
gen pr=exp(z)/(1+exp(z))
gen x=(runiform()<pr)
gen mu=exp(x+z)
gen y=rpoisson(mu)

poisson y x z

gen misspr =exp(y-3)/(1+exp(y-3))
replace x=. if runiform()<misspr

poisson y x z

smcfcs poisson y x z, logit(x)

*now add in an exposure time
clear
set obs 10000
gen exptime = runiform()
gen z=rnormal()
gen pr=exp(z)/(1+exp(z))
gen x=(runiform()<pr)
gen mu=exp(x+z)
gen y=rpoisson(mu*exptime)

poisson y x z, exposure(exptime)

gen misspr =exp(y-3)/(1+exp(y-3))
replace x=. if runiform()<misspr

poisson y x z, exposure(exptime)

smcfcs poisson y x z, exposure(exptime) logit(x)

*now with continuous covariate missing
*note we reduce the magnitude of the x and z effects
*here in order to make the acceptance probability larger
*in rejection sampling
clear
set obs 10000
gen z=(runiform()<0.5)
gen x=z+rnormal()
gen mu=exp(0.1*x+0.1*z)
gen y=rpoisson(mu)

poisson y x z

gen misspr =exp(log(y+1)-2)/(1+exp(log(y+1)-2))
replace x=. if runiform()<misspr

poisson y x z

smcfcs poisson y x z, reg(x) rjlimit(10000)

