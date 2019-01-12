model
{
for (i in 1:N) {
	mu[i]<-v[school[i]]+beta[1]*zage[i]
        +beta[2]*zgevocab[i]+beta[3]*gender[i]+beta[4]*zage[i]*zgevocab[i]
        +beta[5]*zage[i]*gender[i]+beta[6]*zgevocab[i]*gender[i] 	
zgeread[i]~dnorm(mu[i],prec.sigma2)
		}

	for (j in 1:N) {
	v[j]~dnorm(beta0,prec.tau2)
		}

beta[1] ~ dnorm(0,0.0001)
beta[2] ~ dnorm(0,0.0001)
beta[3] ~ dnorm(0,0.0001)
beta[4] ~ dnorm(0,0.0001)
beta[5] ~ dnorm(0,0.0001)
beta[6] ~ dnorm(0,0.0001)
beta0 ~ dnorm(0,0.0001)

prec.sigma2~dgamma(0.001,0.001)
prec.tau2~dgamma(0.001,0.001)

#Bunlari random effectleri kontrol ederken kullan
#tau <- sqrt(1/prec.tau2)
#sigma <- sqrt(1/prec.sigma2) 
}















