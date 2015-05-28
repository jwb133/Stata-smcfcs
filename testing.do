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
set obs 100
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
