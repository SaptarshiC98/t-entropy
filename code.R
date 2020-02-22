library(DirichletReg)
setwd('C:\\Users\\User-PC\\Desktop\\Entropy TIT\\Code')
wt.euc.dist.sq=function(x1,x2,w){
  p=(x1-x2)^2
  p=w*p
  return(sum(p))
}

vec.wt.euc.dist.sq=function(x1,x2,w){
  p=(x1-x2)^2
  p=w*p
  return(p)
}

H=function(w,c=1){
  sum(w*atan(1/w^c))-pi/4
}

update_weight=function(D,c=1,lambda){
  f=function(x){
    w=exp(x)
    w=w/sum(w)
    s=sum(D*w)-lambda*H(w,c)
    return(s)
  }
  para=-D/lambda
  o=optim(par=para,fn=f,control = (maxit=1000))
  x=o$par
  w=exp(x)
  w=w/sum(w)
  return(w)
}
update_weight_MCMC=function(w_0,D,c=1,lambda){
  f=function(w){
    s=sum(D*w)-lambda*H(w,c)
    return(s)
  }
  alpha=exp(-D/lambda)
  W=rdirichlet(200,alpha)
  W[1,]=w_0
  val=numeric(200)
  for(i in 1:200){
    val[i]=f(W[i,])
  }
  ii=which.min(val)
  return(W[ii,])
}

HWkmeans=function(X,M,lambda,tmax,c=1){
  if(is.vector(M)==TRUE){
    M=as.matrix(M)
    M=t(M)
  }
  n=dim(X)[1]
  p=dim(X)[2]
  k=dim(M)[1]
  weight=matrix(1/p,nrow=k,ncol=p)
  weight=rand(k,p)
  for(j in 1:k){
    weight[j,]=weight[j,]/sum(weight[j,])
  }
  label=numeric(n)
  dist=numeric(k)
  t=0
  D=matrix(0,nrow=k,ncol=p)
  #update membership
  repeat{
    t=t+1
    
    for(i in 1 : n){
      for(j in 1 : k){
        dist[j]=wt.euc.dist.sq(X[i,],M[j,],weight[j,])
      }
      label[i]=which.min(dist)
    }
    
    #update centres
    for(j in 1:k){
      I=which(label==j)
      M[j,]=colMeans(X[I,])
    }
    
    #update weights
    #minimize the new objective function
    D=matrix(0,k,p)
#    for(l in 1:k){
#      for(j in 1:p){
#        D[l,j]=0
#      }
#    }
    for(j in 1:k){
      I=which(label==j)
      for(alu in I){
        D[j,]=D[j,]+vec.wt.euc.dist.sq(X[alu,],M[j,],rep(1,p))
      }
    }
    if(t%%1==0){
    for( j in 1:k){
      weight[j,]=update_weight_MCMC(weight[j,],D[j,],c,lambda)
    }
    }
#    cat(t)
#    cat('\n')
  #  points(M,col=2,pch=19)
    if(t>tmax){
      break
    }
    
  }
  return(list(label,M,weight))
  
}
data(iris)
X=iris

X=read.csv("hepatitis.csv")
X=read.csv('C:\\Users\\User-PC\\Desktop\\Datasets\\Classification keel\\newthyroid.csv')
X=read.table('C:\\Users\\User-PC\\Desktop\\Datasets\\microarray\\colon_x.txt')
X=data.matrix(X)
X=t(X)
toss=read.table('C:\\Users\\User-PC\\Desktop\\Datasets\\microarray\\colon_y.txt')
toss=c(data.matrix(toss))
table(toss)
X=readMat('C:\\Users\\User-PC\\Desktop\\Datasets\\ASU\\GLIOMA.mat')
toss=X$Y
X=X$X
k=2
p=dim(X)[2]
n=dim(X)[1]
toss=X[,p]
X=X[,-p]
p=dim(X)[2]
#X=X[,-c(1,81,80,79)]
#centering
for(i in 1:p){
  if(sd(X[,i])==0){
    X[,i]=0
  }else{
    X[,i]=(X[,i]-mean(X[,i]))/sd(X[,i])
  }
}

S=matrix(0,20,k)
for(i in 1:20){
  S[i,]=sample(n,k)
}

indexing=matrix(rep(0,20*2),ncol=2)
for(i in 1 : 20){
  sa=S[i,]
  M=X[sa,]
  l=HWkmeans(X,M,lambda=100,30,c=1)
  (indexing[i,1]=compare(toss,l[[1]],method='nmi'))
  (indexing[i,2]=compare(toss,l[[1]],method='adjusted.rand'))
  cat(i)
  cat('\n')
}
indexing
colMeans(indexing[1:11,])
for(i in 1 : 20){
  sa=S[i,]
  M=X[sa,]
  l=JHWkmeans(X,M,lambda=10,50)
  indexing[i,1]=compare(toss,l[[1]],method='nmi')
  indexing[i,2]=compare(toss,l[[1]],method='adjusted.rand')
  cat(i)
  cat('\n')
}
colMeans(indexing)
indexing
image(l[[3]])
plot(1:1000,l[[3]][1,])
toss=c(rep(1,100),rep(2,100))#,rep(3,100),rep(4,100),rep(5,100))
n=200
p=2
X=matrix(0,n,p)
X[1:100,1]=rnorm(100,0.5)
X[1:100,2]=rnorm(100,0,2)
X[101:200,1]=rnorm(100,10,2)
X[101:200,2]=rnorm(100,0,0.5)
plot(X,col=l[[1]])

for(i in 1:5){
  X=cbind(X,rnorm(200,0,6))
}
M=rand(2)*10
X=data_generate(1000,M,c(0.5,0.5),0.5)
plot(X$data,col=X$label)
k=5
p=20
M=rand(k,p)

prop=0.9
s=sample(p,floor(p*prop))
s=6:20
M[,s]=0
X=data_generate(1000,M,rep(1/k,k),0.02,1)
Y=tsne::tsne(X$data[,1:5])
Y=princomp(X$data[,1:5])
plot(Y$scores)
toss=X$label
X=X$data
lambda=seq(1,10000,length.out = 50)
NMII=matrix(0,50,2)
indexing=matrix(rep(0,10*2),ncol=2)
for(t in 1:100){
  for(i in 1 : 10){
    M=X[sample(n,2),]
    l=HWkmeans(X,M,lambda=lambda[t],50,c=1)
    indexing[i,1]=compare(toss,l[[1]],method='nmi')
    indexing[i,2]=compare(toss,l[[1]],method='adjusted.rand')
#   cat(i)
#    cat('\n')
  }
  NMII[t,]=colMeans(indexing)
  cat(t)
  cat('\n')
}
plot(lambda,NMII[,1],ty='l')
for(t in 1:50){
  for(i in 1 : 10){
    M=X[sample(n,2),]
    l=JHWkmeans(X,M,lambda=lambda[t],50)
    indexing[i,1]=compare(toss,l[[1]],method='nmi')
    indexing[i,2]=compare(toss,l[[1]],method='adjusted.rand')
    #    cat(i)
    #    cat('\n')
  }
  NMI[t,]=colMeans(indexing)
#  NMI[t,1]=max(indexing[,1])
  NMI[t,2]=max(indexing[,2])
  cat(t)
  cat('\n')
}
plot(lambda,NMII[,2],ty='l',ylim=c(0,1))
points(lambda,NMI[1:50,2],ty='l',col=2)
X=read.csv('C:\\Users\\User-PC\\Desktop\\Datasets\\Classification keel\\movement_libras.csv',head=FALSE)
X=data.matrix(X)
dim(X)
toss=X[,91]
X=X[,-91]
plot(X,col=toss,pch=19)
Y=tsne::tsne(X)

plot(Y,col=toss,pch=19)
D=dist(X)
dim(D)
Y=tsne::tsne(D)
library(ggplot2)
plot(Y,col=l[[1]],pch=19)
plot(Y,col=l[[1]],pch=19)
toss1=l[[1]]
toss1[l1==4]
for(i in 1:144){
  if(l[[1]][i]==1){
    toss1[i]=3
  }else if(l[[1]][i]==2){
    toss1[i]=2
  }else if(l[[1]][i]==3){
    toss1[i]=6
  }else if(l[[1]][i]==4){
    toss1[i]=1
  }else if(l[[1]][i]==5){
    toss1[i]=5
  }else if(l[[1]][i]==6){
    toss1[i]=4
  }
  
}
toss1=label_orientation(l[[1]])
Z=cbind(Y,toss1)
Z=cbind(Y,l[[1]])
Z=data.frame(Z)
names(Z)=c('X1','X2','X3')
g=ggplot(Z,aes(x=X1,y=X2,col=as.factor(X3),shape=as.factor(X3),size=1))+geom_point()#col=l$cluster,shape=l$cluster+14)
g=g+labs(x='tSNE dimension 1',y='tSNE dimension 2')
g=g+scale_color_brewer(palette="Dark2")
g=g+theme(legend.position="none")
g
pdf('libras_mwkmeans.pdf',height = 4,width = 4)
g
dev.off()
l=kmeans(X,6)
I=c(which(toss==3),which(toss==4),which(toss==5),which(toss==7),which(toss==11),which(toss==12))
X=X[I,]
toss=toss[I]
label_orientation1=function(label,m){
  m=length(label)
  u=unique(label)
  u=sort(u)
  n=length(u)
  u1=numeric(m)
  for(i in 1:n){
    I=which(label==u[i])
    u1[I]=i
  }
  return(u1)
}
toss=label_orientation(toss)
library(combinat)
a=permn(1:3)
Y1=Y
Y2=Y
library(MASS)
write.matrix(Y,'libras.csv',sep=',')
toss1[toss1==1]=10
toss1[toss1==4]=1
toss1[toss1==10]=4
