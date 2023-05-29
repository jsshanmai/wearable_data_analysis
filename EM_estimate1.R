#the function for estimating Omega_n given H.
#It uses function defined in the file identify.R.
#Arguments:
#TT: the observed time points for n curves
#YY: the observed functional data
#Group_NUM: the most updated number of clusters
#PI_int: the most updated mixing probabilities
#Beta: the most updated basis coefficient of mean functions
#gamma_int: the most updated variances of errors
#Lamda: the most updated eigenvalues
#Alpha: the most updated basis coefficient of eigenfunctions
#HYY: the most updated H
#lambda: the penalty parameter


#return:
#Group_number:  the estimated number of clusters
#PI: the estimated mixing probabilities
#beta: the estimated basis coefficient of mean functions
#sigma: the estimated variances of errors
#lamda: the estimated eigenvalues
#eigencoefficient: the estimated basis coefficient of eigenfunctions
#Sstep: the number of Iterations
#del: the estimated probabilities for individuals belong to each group
# 现在代码中beta就是那个dg啦
EM_estimate<-function(TT,YY,VV,Group_NUM,PI_int,Beta,gamma_int,Lamda,Alpha,HYY,lambda)
{
  q_nn =q_n+3
  delta<-matrix(0,nrow=length(YY),ncol=Group_NUM) 
  Beta_XY<-matrix(0,nrow=q_n+3,ncol=Group_NUM) 
  Beta_matrix<-array(0,dim=c(Group_NUM,q_n+3,q_n+3)) 
  Sigma_XY<-matrix(0,nrow=1,ncol=Group_NUM)
  Sigma_matrix<-matrix(0,nrow=1,ncol=Group_NUM)
  lamda_XY<-matrix(0,nrow=K_g,ncol=Group_NUM)
  lamda_matrix<-matrix(0,nrow=K_g,ncol=Group_NUM)
  Alpha_XY<-array(0,dim=c(q_n,K_g,Group_NUM))
  Alpha_matrix<-array(0,dim=c(q_n,q_n,K_g,Group_NUM))
  f_den<-matrix(0,nrow=length(YY),ncol=Group_NUM)
  
  SStep<-0
  Epsilon<-10
  G_1<-sum(PI_int!=0)
  Beta_1<-Beta
  gamma_1<-gamma_int
  Lamda_1<-Lamda
 
  PI_1=PI_int
    for (k in 1:Group_NUM)
    {
      Alpha[,,k]<-t(qr.Q(qr(t(Alpha[,,k])))) # 获取qr分解的正定矩阵Q
    }
    Alpha_1<-Alpha
    #cat('\n QR分解后的alpha',Alpha_1)
    
 ## Start a loop to find the estimates----
 while((SStep<50)&(Epsilon>10^(-6))){ #修改一下这里的精度Epsilon>10^(-5)
    SStep<-SStep+1
    cat('EM算法迭代了',SStep,'次\n')
    likelihood_em=0
  
    for (i in 1:length(YY)){ # 算一下f和likelihood
      XXX<-splinex1(TT[[i]])
      Phi<-matrix(NA,length(TT[[i]]),K_g)
      
      MU<-matrix(NA,length(TT[[i]]),1)
      # HYYY<-HYY[identify(YY[[i]],allpoint)] 原版是寻找最接近的
      HYYY<-YY[[i]]
      
      for (k in 1:Group_NUM){
          Phi<-XXX%*%t(Alpha[,,k])
        Delta<-matrix(nrow=length(TT[[i]]),ncol=length(TT[[i]]))
        Delta<-Phi%*%Lamda[,,k]%*%t(Phi)+gamma_int[k]*diag(1,ni[i])##
        
        MU<-VV[[i]]%*%Beta[,k]
        f_den[i,k]<-(2*pi)^(-length(TT[[i]])/2)/sqrt(det(Delta))*exp(-1/2*t(HYYY-MU)%*%solve(Delta)%*%(HYYY-MU))
      }
       likelihood_em=likelihood_em+log(f_den[i,]%*%PI_int) 
      }
    
    for (i in 1:length(YY)){
      XXX<-splinex1(TT[[i]])
      Phi<-matrix(NA,length(TT[[i]]),K_g)
      MU<-matrix(NA,length(TT[[i]]),1)
      HHYY<-YY[[i]]
      for (k in 1:Group_NUM){
          Phi<-XXX%*%t(Alpha[,,k])
        Delta<-matrix(nrow=length(TT[[i]]),ncol=length(TT[[i]]))
        Delta<-Phi%*%Lamda[,,k]%*%t(Phi)+gamma_int[k]*diag(1,ni[i])##
        #cat("Delta",Delta)
        delta[i,k]<-PI_int[k]*f_den[i,k]/((t(as.matrix(PI_int))%*%as.matrix(f_den[i,])))
        
        # 看到这，不过是把样条的12*6改成了vv的12*9
        MU<-VV[[i]]%*%Beta[,k]
        times<-matrix(HHYY-MU,ni[i],1)
        Beta_XY[,k]<-Beta_XY[,k]+delta[i,k]*t(VV[[i]])%*%solve(Delta)%*%HHYY##beta_k
        Beta_matrix[k,,]<-Beta_matrix[k,,]+delta[i,k]*t(VV[[i]])%*%solve(Delta)%*%VV[[i]]##beta_K
        
        # 这个算的应该是lambda
        Sigma_XY[k]<-Sigma_XY[k]+delta[i,k]*sum(diag(solve(Delta)%*%solve(Delta)%*%(Phi%*%Lamda[,,k]%*%t(Phi)-times%*%t(times))))
        Sigma_matrix[k]<-Sigma_matrix[k]+delta[i,k]*sum(diag(solve(Delta)%*%solve(Delta)))
        for (g in 1:K_g){
          # 下面对于R的处理很好，用Delta减去第k项的结果
          lamda_XY[g,k]<-lamda_XY[g,k]+delta[i,k]*t(Phi[,g])%*%solve(Delta)%*%(Delta-Lamda[g,g,k]*Phi[,g]%*%t(Phi[,g])-times%*%t(times))%*%solve(Delta)%*%Phi[,g]
          lamda_matrix[g,k]<-lamda_matrix[g,k]+delta[i,k]*t(Phi[,g])%*%solve(Delta)%*%Phi[,g]%*%t(Phi[,g])%*%solve(Delta)%*%Phi[,g]
          
          # cat('delta[i,k]',delta[i,k])应该是对的
          
          Alpha_XY[,g,k]<-Alpha_XY[,g,k]+delta[i,k]*t(XXX)%*%solve(Delta)%*%times%*%t(times)%*%solve(Delta)%*%Phi[,g]
          Alpha_matrix[,,g,k]<-Alpha_matrix[,,g,k]+delta[i,k]*t(XXX)%*%solve(Delta)%*%XXX
          #cat('Alpha_matrix[,,g,k]',Alpha_matrix[,,g,k],'\n')
        }
      }
    }
    PII=rep(NA,Group_NUM)
    for (k in 1:Group_NUM){
      PII[k]<-max(0,1/(1-lambda*Group_NUM)*(1/length(YY)*sum(delta[,k])-lambda)) # 判罚来了
    }
    PII<-PII/sum(PII)
 
    for (k in 1:Group_NUM){ # 算出最后的结果，beta没有下标的才是真实的值
      Beta[,k]=ginv(Beta_matrix[k,,])%*%Beta_XY[,k]
      gamma_int[k]=-Sigma_XY[k]/(Sigma_matrix[k])
      for (g in 1:K_g){
        Lamda[g,g,k]<-lamda_XY[g,k]/(-lamda_matrix[g,k])
        # 判断矩阵是否可逆
          Alpha[g,,k]<-ginv(Alpha_matrix[,,g,k])%*%Alpha_XY[,g,k]
      }
    }
    ############identification for groups小的方第一，大的放最后
    if(Group_NUM==2){
      # 这里没看懂啊,就差这了----
      #WZ=order(c(splinex1(0)%*%Beta[1:q_nn,1],splinex1(0)%*%Beta[1:q_nn,2]))
      WZ=order(c(cbind( splinex1(0),Z_i) %*%Beta[1:q_nn,1], cbind( splinex1(0),Z_i)%*%Beta[1:q_nn,2]))
      Beta=Beta[,WZ]
      gamma_int=gamma_int[WZ]
      Lamda=Lamda[,,WZ]
      Alpha=Alpha[,,WZ]
      PII=PII[WZ]
    }  

if((is.nan(sum(delta)))&(SStep==1))  {  Updated_Est<-list(Group_number=G_1,PI=PI_int,beta=Beta_1,sigma=gamma_1,lamda=Lamda_1,eigencoefficient=Alpha_1,Sstep=100,del=delta)
    break}
if(is.nan(sum(delta))) {                Updated_Est<-list(Group_number=G_1,PI=PI_int,beta=Beta_1,sigma=gamma_1,lamda=Lamda_1,eigencoefficient=Alpha_1,Sstep=100,del=delta_1)
    break}

 ############identification for eigenfunctions
    for (k in 1:Group_NUM)
    {
      Alpha[,,k]<-t(qr.Q(qr(t(Alpha[,,k]))))
    }
    eigen11<-function(t) splinex1(t)%*%Alpha[1,1:q_n,1]
    eigen12<-function(t) splinex1(t)%*%Alpha[2,1:q_n,1]
    eigen21<-function(t) splinex1(t)%*%Alpha[1,1:q_n,2]
    eigen22<-function(t) splinex1(t)%*%Alpha[2,1:q_n,2]

    
if(is.nan(fderiv(eigen11,0.6))){break}
    # <--就是赋值为负的，就是调整一下方向的意思
if(fderiv(eigen11,0.6)>0) {Alpha[1,,1]<--Alpha[1,,1]}
 if(fderiv(eigen12,0.6)>0) {Alpha[2,,1]<--Alpha[2,,1]} 

if(fderiv(eigen21,0.6)<0) {Alpha[1,,2]<--Alpha[1,,2]}
 if(fderiv(eigen22,0.6)>0) {Alpha[2,,2]<--Alpha[2,,2]} 
    #yucha<-c(max(abs(Beta-Beta_1)),max(abs(gamma_int-gamma_1)),max(abs(PII-PI_1)),max(abs(Lamda-Lamda_1)),max(abs(Alpha-Alpha_1)))
    Epsilon<-max(max(abs(Beta-Beta_1)),max(abs(gamma_int-gamma_1)),max(abs(Lamda-Lamda_1)),max(abs(PII-PI_1)),max(abs(Alpha-Alpha_1)))
    cat('此次迭代的最大变差为',Epsilon,'\n')
############keep the groups with nonzero probabilities
    G_1<-sum(PII!=0) 
    Group_NUM<-sum(PII!=0)
    Beta<-as.matrix(Beta[,PII!=0])
    Beta_1<-Beta
    gamma_int<-gamma_int[PII!=0]
    gamma_1<-gamma_int
    Lamda<-Lamda[,,PII!=0]
    Lamda_1<-Lamda
    Alpha<-Alpha[,,PII!=0]
    Alpha_1<-Alpha
    PI_1=PII[PII!=0]
    PI_int<-PII[PII!=0]
    delta_1<-delta
    df<-sum(PI_int!=0)
    Updated_Est<-list(Group_number=G_1,PI=PI_int,beta=Beta_1,sigma=gamma_1
                      ,lamda=Lamda_1,eigencoefficient=Alpha_1,Sstep=SStep,del=delta)

    #聚类数太少，重新初始化，下面不用看
     if(G_1<2){break}
    f_den<-matrix(0,nrow=length(YY),ncol=Group_NUM)
    delta<-matrix(nrow=length(YY),ncol=Group_NUM)
    Beta_XY<-matrix(0,nrow=q_nn,ncol=Group_NUM)
    Beta_matrix<-array(0,dim=c(Group_NUM,q_nn,q_nn))
    Sigma_XY<-matrix(0,nrow=1,ncol=Group_NUM)
    Sigma_matrix<-matrix(0,nrow=1,ncol=Group_NUM)
    lamda_XY<-matrix(0,nrow=K_g,ncol=Group_NUM)
    lamda_matrix<-matrix(0,nrow=K_g,ncol=Group_NUM)
    Alpha_XY<-array(0,dim=c(q_n,K_g,Group_NUM))
    Alpha_matrix<-array(0,dim=c(q_n,q_n,K_g,Group_NUM))
  }
  return(Updated_Est)
}
