model {
    for( i in 1:N) {    
      num_awards[i] ~ dpois( lambda[i] )

      log(lambda[i]) <- beta0 + beta[1]*academic[i]+beta[2]*vocational[i]+beta[3]*zmath[i]+beta[4]*academic[i]*zmath[i]+beta[5]*vocational[i]*zmath[i]
    }

    beta0 ~ dnorm(0,0.0001)
    beta[1] ~ dnorm(0,0.0001)
    beta[2] ~ dnorm(0,0.0001)
    beta[3] ~ dnorm(0,0.0001)
    beta[4] ~ dnorm(0,0.0001)
    beta[5] ~ dnorm(0,0.0001)
   
  }







