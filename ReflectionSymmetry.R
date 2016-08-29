library("MASS")
library("mpoly")

l0x<-mp('1 + 0.0000000000001 x')
l1x<-mp('x')
l2x<-mp('1.5 x^2 - 0.5')
l3x<-mp('2.5 x^3 - 1.5 x')
l4x<-mp('4.375 x^4 - 3.75 x^2 + 0.375 ')
l5x<-mp('7.875 x^5 - 8.75 x^3 + 1.875 x ')
polx<-mpolyList(l0x,l1x,l2x,l3x,l4x,l5x)

l0y<-mp('1 + 0.0000000000001 y')
l1y<-mp('y')
l2y<-mp('1.5 y^2 - 0.5')
l3y<-mp('2.5 y^3 - 1.5 y')
l4y<-mp('4.375 y^4 -3.75 y^2 + 0.375 ')
l5y<-mp('7.875 y^5 - 8.75 y^3 + 1.875 y ')
poly<-mpolyList(l0y,l1y,l2y,l3y,l4y,l5y)

l0z<-mp('1 + 0.0000000000001 z')
l1z<-mp('z')
l2z<-mp('1.5 z^2 - 0.5')
l3z<-mp('2.5 z^3 - 1.5 z')
l4z<-mp('4.375 z^4 - 3.75 z^2 + 0.375 ')
l5z<-mp('7.875 z^5 - 8.75 z^3 + 1.875 z ')
polz<-mpolyList(l0z,l1z,l2z,l3z,l4z,l5z)

l0v<-mp('1 + 0.0000000000001 v')
l1v<-mp('v')
l2v<-mp('1.5 v^2 - 0.5')
l3v<-mp('2.5 v^3 - 1.5 v')
l4v<-mp('4.375 v^4 - 3.75 v^2 + 0.375 ')
l5v<-mp('7.875 v^5 - 8.75 v^3 + 1.875 v ')
polv<-mpolyList(l0v,l1v,l2v,l3v,l4v,l5v)

l0w<-mp('1 + 0.0000000000001 w')
l1w<-mp('w')
l2w<-mp('1.5 w^2 - 0.5')
l3w<-mp('2.5 w^3 - 1.5 w')
l4w<-mp('4.375 w^4 - 3.75 w^2 + 0.375 ')
l5w<-mp('7.875 w^5 - 8.75 w^3 + 1.875 w ')
polw<-mpolyList(l0w,l1w,l2w,l3w,l4w,l5w)




#Two dimensional Legendre polynomials
polxy<-NULL
k=1
for(i in 1:length(polx))
{
  for(j in 1:length(poly))
  {
    if((i+j<8) && (i+j-2)%%2!=0 ) 
    {
      polxy[[k]]<-polx[[i]]*poly[[j]]
      k=k+1
    }
    
  }
}

#Three dimensional Legendre polynomials
polxy<-NULL
l=1
for(i in 1:length(polx))
{
  for(j in 1:length(poly))
  {
    for(k in 1:length(polz))
    {
      if((i+j+k<9) && (i+j+k-2)%%2==0 ) 
      {
        polxy[[l]]<-polx[[i]]*poly[[j]]*polz[[k]]
        l=l+1
      }
    }
    
    
  }
}

#4D

polxy<-NULL
l=1
for(i in 1:length(polx))
{
  for(j in 1:length(poly))
  {
    for(k in 1:length(polz))
    {
      for(m in 1:length(polz))
      {
        if((i+j+k+m<9) && (i+j+k+m-2)%%2!=0 ) 
        {
          polxy[[l]]<-polx[[i]]*poly[[j]]*polz[[k]]*polv[[m]]
          l=l+1
        }
      }
    }
    
    
  }
}

#5D

polxy<-NULL
l=1
for(i in 1:length(polx))
{
  for(j in 1:length(poly))
  {
    for(k in 1:length(polz))
    {
      for(m in 1:length(polz))
      {
        for(n in 1:length(polz))
        {
          if((i+j+k+m+n<11) && (i+j+k+m+n-2)%%2==0 ) 
          {
            polxy[[l]]<-polx[[i]]*poly[[j]]*polz[[k]]*polv[[m]]*polw[[n]]
            l=l+1
          }
        }
      }
      
    }
    
    
  }
}


#Calculates the statistic for reflection symmetry

statRef<-function(data,polxy){
  n=nrow(data)
 
  
  #plot(data)
  dataValues=NULL
  
  #integrate
  
  for(i in 1:length(polxy))
  {
    dataValues[[i]]=(1/sqrt(n))*sum(apply(data,1,as.function(polxy[[i]],varorder=c('x','y')))-apply(-data,1,as.function(polxy[[i]],varorder=c('x','y'))))
  }
  
  
  return(dataValues)
}

#Calculates the statistic for angular symmetry 
statAng<-function(data,polxy){
  n=nrow(data)
  
  
  dataValues=0
  
  
  #integrate
  dataValues=NULL
  for(i in 1:length(polxy))
  {
    dataValues[[i]]=1/n*sum(apply((data),1,as.function(polxy[[i]],varorder=c('x','y')))-apply((-data),1,as.function(polxy[[i]],varorder=c('x','y'))))
  }
  
  return(dataValues)
}

Sigma <- matrix(c(10,3,3,2),2,2)
Sigma
#number of data
n=1000
#data
data<-NULL
#Symmetric Normal
data<-mvrnorm(n, rep(0, 2), Sigma)
#Symmetric Chi Square
data<-cbind(rchisq(n, 10),rchisq(n,10))
#Symmetric Exponential
data<-cbind(rexp(n,rate=1),rexp(n,rate=1))
#Asymmetric Exponential
data<-cbind(rexp(n,rate=1),rexp(n,rate=2))
#Symmetric Uniform
data<-cbind(runif(n,-1,1),runif(n,-1,1))
#Asymmetric uniform
data<-cbind(runif(n,0,1),runif(n,0,1))
#Symmetric Weibull
data<-cbind(rweibull(n, 2, scale = 1),rweibull(n, 2, scale = 1))
#Asymmetric Weibull
data<-cbind(rweibull(n, 0.5, scale = 1),rweibull(n, 0.5, scale = 1))
#Beta
data<-cbind(rbeta(n,2,2),rbeta(n,2,2))

#3-d Symmetric Weibull
data<-cbind(rweibull(n, 2, scale = 1),rweibull(n, 2, scale = 1),rweibull(n, 2, scale = 1),rweibull(n, 2, scale = 1),rweibull(n, 2, scale = 1))

data<-cbind(runif(n,-1,1),runif(n,-1,1),runif(n,-1,1),runif(n,-1,1),runif(n,-1,1))

#Angular Symmetric but not Reflection Symmetric data

data<-cbind(runif(n,-1,1),runif(n,-1,1))
norm<-apply(data,1,norma)
data<-data/norm
data<-t(apply(data,1,proy))

#number of bootstraps
b=1000
rRef=0
rAng=0
cuantilesRef<-NULL
referenciaRef<-NULL
cuantilesAng<-NULL
referenciaAng<-NULL
for(i in 1:100)
{
#data
  
  
  data<-cbind(rweibull(n, 2, scale = 1),rweibull(n, 2, scale = 1))
  
#Data trimming
  ind<-1:n
  q<-t(apply(data,2,quantile,probs=c(0.05,0.95)))
  ind<-ind[data[ind,1]>q[1,1] & data[ind,1]<q[1,2] & data[ind,2]>q[2,1] & data[ind,2]<q[2,2]]
  data<-data[ind,]

  n1<-nrow(data)



  #####################
  #Reflection Symmetry#
  #####################

  #Bootstrap for reflection symmetry
  m1<-NULL
  m1<-apply(data, 2, mean, trim=.0)

 
  data1<-t(apply(data,1,resta,media=t(m1)))

  

  
  #Maps the data to [-1,1]
  
  data1<-normal(data1)
 
  MatrixRef<-NULL
  MatrixRef<-replicate(b,sample(c(-1,1),size=n1,replace=TRUE,prob=c(0.5,0.5))*(data1))
  bootMatrixRef<-NULL
  #Bootstrap with Legendre polynomials
  for(i in 1:b)
  {
    bootMatrixRef<-rbind(bootMatrixRef,statRef(MatrixRef[, ,i],polxy))
  }

  CRef<-NULL
  CRef<-var(bootMatrixRef)  

  #Calculate the quantiles

  QRef=NULL
  for(i in 1:b)
  {
    QRef[i]=t(bootMatrixRef[i,])%*%CRef%*%(bootMatrixRef[i,])
  }

  #Value of the statistic for the original data
  vpRef<-NULL
  vpRef=statRef(data1,polxy)


  QpRef<-NULL
  QpRef=t(vpRef)%*%CRef%*%vpRef

  #95% quantile
  qRef<-NULL
  qRef<-quantile(QRef,probs=c(0,0.5,0.9,0.95,0.975,1))
  cuantilesRef<-rbind(cuantilesRef,qRef)
  referenciaRef<-rbind(referenciaRef,QpRef)
  if(qRef[4]>QpRef) rRef=rRef+1
  

  ##################
  #Angular Symmetry#
  ##################

  #Bootstrap for angular symmetry
  m2<-NULL
  m2<-apply(data, 2, median)

  
  data2<-t(apply(data,1,resta,media=t(m2)))
  
  #Centralize data
  norm<-apply(data2,1,norma)
  #Normalization
  data2<-data2/norm
  
  MatrixAng<-NULL
  MatrixAng<-replicate(b,sample(c(-1,1),size=n1,replace=TRUE,prob=c(0.5,0.5))*(data2))
  bootMatrixAng<-NULL
  #Bootstrap with Legendre polynomials
  for(i in 1:b)
  {
    bootMatrixAng<-rbind(bootMatrixAng,statAng(MatrixAng[, ,i],polxy))
  }

  CAng<-var(bootMatrixAng)  

  #Calculate the quantiles

  QAng=NULL
  for(i in 1:b)
  {
    QAng[i]=t(bootMatrixAng[i,])%*%CAng%*%(bootMatrixAng[i,])
  }

  #Value of the statistic for the original data

  vpAng<-NULL
  vpAng<-statAng(data2,polxy)

  QpAng<-NULL
  QpAng<-t(vpAng)%*%CAng%*%vpAng

  #95% quantile
  qAng<-quantile(QAng,probs=c(0,0.5,0.9,0.95,0.975,1))
  cuantilesAng<-rbind(cuantilesAng,qAng)
  referenciaAng<-rbind(referenciaAng,QpAng)
  if(qAng[4]>QpAng) rAng=rAng+1
}

#######
#Visualize data
#######

m<-apply(data, 2, mean, trim=.0)
norma<-function(d)
{
  return(sqrt(sum(d^2)))
}



norm<-apply(data-m,1,norma)
plot((data-m)/norm)
plot(data-m)

#Auxiliary Functions
resta<-function(d,media)
{
  return(d-media)
}

normal<-function(d)
{
  d[,1]=d[,1]/max(abs(d[,1]))
  d[,2]=d[,2]/max(abs(d[,2]))            
  return(d)
}

norma<-function(d)
{
  return(sqrt(sum(d^2)))
}
