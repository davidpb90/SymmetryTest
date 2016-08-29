###########################
#Sign Symmetry
###########################


library("MASS")
library("mpoly")

#Auxiliary functions to create Legendre polynomials

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
    if((i+j<8)) #&& (i+j-2)%%2!=0 ) 
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

statRefSign<-function(data,polxy,sign){
  n=nrow(data)
  #Calculates trimmed mean
  
  
  #Maps the data to [-1,1]
  #data<-normal(data)
  

  dataValues=NULL
  
  #integrate
  
  for(i in 1:length(polxy))
  {
    dataValues[[i]]=(1/sqrt(n))*sum(apply(data,1,as.function(polxy[[i]],varorder=c('x','y')))-apply(t(sign*t(data)),1,as.function(polxy[[i]],varorder=c('x','y'))))
  }
  
  
  return(dataValues)
}

#Calculates the statistic for angular symmetry 
statAngSign<-function(data,polxy,sign){
  n=nrow(data)
  dataValues=0
  
  
  
  #integrate
  dataValues=NULL
  for(i in 1:length(polxy))
  {
    dataValues[[i]]=1/n*sum(apply((data),1,as.function(polxy[[i]],varorder=c('x','y')))-apply(t(sign*t(data)),1,as.function(polxy[[i]],varorder=c('x','y'))))
  }
  
  return(dataValues)
}


#number of data
n=1000

Sigma <- matrix(c(10,3,3,2),2,2)
Sigma
#data
data<-NULL
#Symmetric Normal
data<-mvrnorm(n, rep(0, 3), Sigma)
#Symmetric Chi Square
data<-cbind(rchisq(n, 1),rchisq(n,1))
#Symmetric Exponential
data<-cbind(rexp(n,rate=1),rexp(n,rate=1))
#Asymmetric Exponential
data<-cbind(rexp(n,rate=1),rexp(n,rate=2))
#Symmetric Uniform
data<-cbind(runif(n,-1,1),runif(n,-1,1))
#Asymmetric uniform
data<-cbind(runif(n,0,1),runif(n,0,1))
#Symmetric Weibull
data<-cbind(rweibull(n, 3, scale = 1),rweibull(n, 3, scale = 1))
#Asymmetric Weibull
data<-cbind(rweibull(n, 2, scale = 1),rweibull(n, 2, scale = 1))

#3-d Symmetric Weibull
data<-cbind(rweibull(n, 2, scale = 1),rweibull(n, 2, scale = 1),rweibull(n, 2, scale = 1),rweibull(n, 2, scale = 1),rweibull(n, 2, scale = 1))

data<-cbind(runif(n,-1,1),runif(n,-1,1),runif(n,-1,1),runif(n,-1,1),runif(n,-1,1))

#number of bootstraps
b=100

#Tail trimming
ind<-1:n
q<-t(apply(data,2,quantile,probs=c(0.05,0.95)))
ind<-ind[data[ind,1]>q[1,1] & data[ind,1]<q[1,2] & data[ind,2]>q[2,1] & data[ind,2]<q[2,2]]
data<-data[ind,]

n1<-nrow(data)



#####################
#Reflection Symmetry#
#####################

#Bootstrap for reflection symmetry
m<-NULL
m<-apply(data, 2, mean, trim=.0)
#Centering the data
resta<-function(d,media)
{
  return(d-media)
}
data<-t(apply(data,1,resta,media=t(m)))

#Normalization
normal<-function(d)
{
  d[,1]=d[,1]/max(d[,1])
  d[,2]=d[,2]/max(d[,2])
  return(d)
}


data<-normal(data)


#Calculating the statistics
QRef<-list()
QpRef<-list()
vpRef<-list()
qRef<-list()


for(i in 1:2)
{
  for(j in 1:2)
  {
    if(i+j!=2)
    {
      if(i%%2==0)
      {
        fSign=-1
      }
      else
      {
        fSign=1
      }
      if(j%%2==0)
      {
        sSign=-1
      }
      else
      {
        sSign=1
      }
      MatrixRef<-NULL
      MatrixRef<-replicate(b,t(matrix(unlist(sample(list(c(1,1),c(fSign,sSign)),size=n1,replace=TRUE,prob=c(0.5,0.5))),2))*data)
      bootMatrixRef<-NULL
      #Bootstrap with Legendre polynomials
      for(k in 1:b)
      {
        bootMatrixRef<-rbind(bootMatrixRef,statRefSign(MatrixRef[, ,k],polxy,c(fSign,sSign)))
      }  
      CRef<-NULL
      CRef<-var(bootMatrixRef)  
      QRef[[i*2+j-3]]<-matrix()
      for(k in 1:b)
      {
        QRef[[i*2+j-3]][[k]]=t(bootMatrixRef[k,])%*%CRef%*%(bootMatrixRef[k,])
      }
      
      
      #Value of the statistic for the original data
      
      vpRef[[i*2+j-3]]<-statRefSign(data,polxy,c(fSign,sSign))
      
      
      
     
      QpRef[[i*2+j-3]]=t(vpRef[[i*2+j-3]])%*%CRef%*%vpRef[[i*2+j-3]]
      
      #95% quantile
      
      
      qRef[[i*2+j-3]]<-quantile(QRef[[i*2+j-3]],probs=c(0,0.05,0.5,0.95,1))
    }
  }
}




##################
#Angular Symmetry#
##################

#Bootstrap for reflection symmetry
m<-NULL
m<-apply(data, 2, mean, trim=.0)
#Centering the data
resta<-function(d,media)
{
  return(d-media)
}
data<-t(apply(data,1,resta,media=t(m)))

#Normalization

norma<-function(d)
{
  return(sqrt(sum(d^2)))
}

norm<-apply(data,1,norma)

data<-data/norm



QAng<-list()
QpAng<-list()
vpAng<-list()
qAng<-list()


for(i in 1:2)
{
  for(j in 1:2)
  {
    if(i+j!=2)
    {
      if(i%%2==0)
      {
        fSign=-1
      }
      else
      {
        fSign=1
      }
      if(j%%2==0)
      {
        sSign=-1
      }
      else
      {
        sSign=1
      }
      
      MatrixAng<-NULL
      MatrixAng<-replicate(b,t(matrix(unlist(sample(list(c(1,1),c(fSign,sSign)),size=n1,replace=TRUE,prob=c(0.5,0.5))),2))*data)
      bootMatrixAng<-NULL
      #Bootstrap with Legendre polynomials
      for(k in 1:b)
      {
        bootMatrixAng<-rbind(bootMatrixAng,statAngSign(MatrixAng[, ,k],polxy,c(fSign,sSign)))
      }
      
      CAng<-var(bootMatrixAng)  
      
      #Calculate the quantiles
      
      QAng[[i*2+j-3]]<-matrix()
      for(k in 1:b)
      {
        QAng[[i*2+j-3]][[k]]=t(bootMatrixAng[k,])%*%CAng%*%(bootMatrixAng[k,])
      }
      
      
      #Value of the statistic for the original data
      
      vpAng[[i*2+j-3]]<-NULL
      vpAng[[i*2+j-3]]<-statAngSign(data,polxy,c(fSign,sSign))
      
      QpAng[[i*2+j-3]]<-NULL
      QpAng[[i*2+j-3]]<-t(vpAng[[i*2+j-3]])%*%CAng%*%vpAng[[i*2+j-3]]
      
      #95% quantile
      qAng[[i*2+j-3]]<-NULL
      qAng[[i*2+j-3]]<-quantile(QAng[[i*2+j-3]],probs=c(0,0.05,0.5,0.95,1))
      
    }
  }
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


