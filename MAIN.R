##Required R packages
     library(MASS)
     library(pcaPP)
     library(Hmisc)
     library(TauStar)
     library(energy) 
     library( pracma)
##some parameter settings
     n          <-  100
     p          <-  50
     m          <-  1000
     a          <-  rep(0,m)

##Two extreme value distributions
     G         <-  function(s) {exp(-1/sqrt(8*pi)*exp(-s/2))}
     F         <-  function(s) {exp(-2.467/sqrt(8*pi)*exp(-s/2))}
     alpha      <-  1
     beta       <-  2
     kesai      <-  0.7
     theta      <-  0.5 

##data generation process
        y<-matrix(0,p,n)
        mu<-rnorm(p,0,0.5)
        eta<-rnorm(p*(n+51),0,1)
        eta<-matrix(eta,p,n+51)
        x<-matrix(0,p,51+n)
          for (i in 1:p){
            for (t in 2:(n+51)){
               x[i,1]<-mu[i]+eta[i,1]
               x[i,t]<-kesai*x[i,t-1]+mu[i]+eta[i,t]
         }
        }
       x<-x[,52:(n+51)]
       xbar<- apply(x, 1, sum)/n
       sigmai2<-rep(0,p)
       sigma2<-0.5/(1/p*sum((1+theta *xbar)^2))
       sigmai2<-sigma2*(1+theta *xbar)^2 
        #change the distribution
        set.seed(k)
        varepsilon<-rnorm(p*n,0,1)
        varepsilon<-matrix(varepsilon,p,n)
       #
       epsilon<-diag(sqrt(sigmai2))%*%varepsilon

##centralize and get within estimator and residuals
       yhat<-matrix(0,p,n)
       xhat<-matrix(0,p,n)
       for (i in 1:p){
            for (t in 1:n){
              y[i,t]<-alpha+x[i,t]*beta+mu[i]+epsilon[i,t]
          }
              yhat[i,]<-y[i,]-rep(1/n*sum(y[i,]),n)
              xhat[i,]<-x[i,]-rep(1/n*sum(x[i,]),n)
       }
       betahat<-1/sum(xhat^2)*sum(xhat*yhat)
       ehat<-matrix(0,p,n)
       for (i in 1:p){
              ehat[i,]   <-  yhat[i,]-xhat[i,]*betahat
       }

##Calculate statistics based on Pearson correlation
        rho     <-   matrix(0,p,p)
        for (j in 2:p){
           for (i in 1:(j-1)){
              rho[i,j]<-t(ehat[i,])%*%(ehat[j,])/sqrt(sum(ehat[i,]^2)*sum(ehat[j,]^2))
            }
        }
        r     <-   rho+t(rho)+diag(p)
        R     <-   r[which(upper.tri(r))]
        Sr   <-  (sum(R^2)-p*(p-1)/(n-1)/2)/sqrt(p^2/n^2)
        Lr   <-   n*(max(abs(R)))^2-4*log(p)+log(log(p))
        Cr   <-  min(1-pnorm(Sr),1-G(Lr))
        CDP    <-   sqrt(2*n/(p*(p-1)))*sum(R)
        LMBFK <-   sqrt(1/(p*(p-1)))*(sum(n*R^2)-p*(p-1)/2)-p/(2*(n-1))

##Calculate statistics based on Spearman's rho correlation
        S     <-   cor(t(ehat),method="spearman")
        Sp    <-   S[which(upper.tri(S))]
        Srho <-   (sum(Sp^2)-p*(p-1)/2/(n-1))*n/p
        Lrho <-  (n-1)*max(Sp^2)-4*log(p)+log(log(p))
        Crho <-  min(1-pnorm(Srho),1-G(Lrho))

##Calculate statistics based on Kendall's tau correlation
        tau   <-   cor.fk(t(ehat)) 
        Tau   <-   tau[which(upper.tri(tau))]
        Stau  <-  (sum(Tau^2)-choose(p,2)*2*(2*n+5)/9/n/(n-1))*9*n/4/p
        Ltau  <-  9*n*(n-1)*max(Tau^2)/2/(2*n+5)-4*log(p)+log(log(p))
        Ctau  <-   min(1-pnorm(Stau),1-G(Ltau))

##Calculate statistics based on Hoeffding's D correlation
        D     <-   hoeffd(t(ehat))$D 
        d     <-   D[which(upper.tri(D))] 
        SD   <-  (sum(d^2)-choose(p,2)*2*(n^2+5*n-32)/9/n/(n-1)/(n-3)/(n-4))/900*n^2/choose(5,2)^2/2/p/sqrt(1/810000^2+6/945000^2)/(1+2/n^{1/2})
        LD   <-  (pi^4*(n-1)*max(d)/30-4*log(p)+log(log(p))+pi^4/36)/(1 + 2/log(p*n))
        CD   <-   min(1-pnorm(SD),1-F(LD))

##Calculate statistics based on Bergsma-Dassios-Yanagimoto's tau* correlation
        TM      <-   matrix(0,p,p)
        for (j in 2:p){
              for (i in 1:(j-1)){
                         TM[i,j]  <-   tStar(ehat[i,],ehat[j,])  
                                 }
                       }
        tstar    <-   TM[which(upper.tri(TM))] 
        Staustar <-  n^2*(sum(tstar^2)-choose(p,2)*8/75*(3*n^2+5*n-18)/n/(n-1)/(n-2)/(n-3))/72/p/sqrt(1/225^2+24/525^2)/(1+n^{-1/2})
        Ltaustar <-  pi^4*(n-1)*max(tstar)*6/4/54-4*log(p)+log(log(p))+pi^4/36
        Ctaustar <-  min(1-pnorm(Staustar),1-F(Ltaustar))

##Calculate statistics based on Blum-Kiefer-Rosenblatt's R correlation
        BKR     <-   (5*TM*6/4-3*D)/2
        bkr      <-   BKR [which(upper.tri(BKR))] 
        SR <-    (sum(bkr^2)-choose(p,2)*(2*(n^3-3*n^2-6*n+10))/((n-4)*(n-3)*(n-2)*(n-1)*n))*n^2/choose(6,2)^2/2/p/sqrt(1/225^2+24/525^2)/(1+n^{-1/2})
        LR <-    pi^4*(n-1)*max(bkr)/90-4*log(p)+log(log(p))+pi^4/36 
        CR <-    min(1-pnorm(SR),1-F(LR))



