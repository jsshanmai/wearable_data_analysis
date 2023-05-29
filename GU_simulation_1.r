rm(list=ls())
# 用来进行第一次数值模拟，有两类以及两个主成分
############This is the main code  for simulation 1. It uses functions defined in the files EM_estimate1.R, original function.R, first_der function.R, H_est.R, identify.R.
########the true clutser id: cluster_result[,1:n]
########the estimated clutser id: cluster_result[,(n+1):(2*n)]
setwd("C:/Users/matebook 14/Desktop/函数型数据聚类/我的论文内容/穿戴设备/code_old")
library(MASS) #ginv()
library(splines)
library(pracma)
library(birk)
#source("original function.R")
#source("first_der function.R")
#source("H_est.R")
#source("identify.R")

set.seed(192)
##----#定义正交基函数----
kontx=c(0.4,0.75)###
splinex<-function(t) bs(t,df=NULL,kontx,degree=3,intercept=T,Boundary.knots=c(0,1))
q_n <-6
work<-matrix(0,q_n,q_n)
pos<-matrix(0,q_n,1)
for (i in 1:10000)
{work=work+splinex(i/10000)[1,]%o%splinex(i/10000)[1,] # 这是向量内积
}
working=work/10000 # working 就是A
VV1<-solve(sqrtm(working)$B) #A^-1/2
splinex1<- function(t) splinex(t)%*%VV1

#开始执行----
source("EM_estimate1.R")
#----产生模拟数据的函数----
TTT <-sort(runif(12)) #总的时间点，样本产生的时间点可能有缺失

sim_data = function(n1,PI,gamma,LLamda) # 下面都是一个样本一个样本来的
{
  T1<-sort(sample(TTT,n1,replace = FALSE))
  #T1<-sort(runif(n1)) # 时间点是这样定义的
  X_i <-matrix (runif(n1),1,n1)
  Phi<-matrix(NA,n1,K_g)
  index<-runif(1)
  INDEX=1 
  G<-length(PI)-1
  for (l in 1:G){
    if ( index>=sum(PI[1:l]) )# 基于产生的随机数，确定类
    { INDEX <- INDEX+1 }
  }
  xi<-mvrnorm(1,rep(0,K_g),LLamda[,,INDEX],tol = 1e-15)
  Sigma_r=matrix(nrow=n1,ncol=n1)
  Sigma_r=gamma[INDEX]*diag(1,n1)
  error1<-mvrnorm(1,rep(0,n1),Sigma_r,tol = 1e-15)
  
  
  M_it <- function(t){ # tmpM是可以抽自变量的设计矩阵mit
    mit <-matrix(NA, nrow=length(t), ncol=q_n)
    for (t1 in t){
      t_inter <-which.closest (T1,t1)
      tmpM <-0 #matrix(NA, nrow = nrow(t), ncol=q_n)
      for (j in T1[1:t_inter]){ # 把t固定住了之后，就看s从小到大的变化
        j_inter <-which.closest (T1,j) # 第几个就积分几个
        # 产生一下X_i*Bn的积分
        tmpM <-tmpM+X_i[,1:j_inter] %*%splinex1(t1-T1[1:j_inter]) 
      }
      mit[t_inter,] <- tmpM
    }
    return(mit)}
  
  #----均值函数----
  G_igt <- function(g,t){ # 这个函数G_igt(1,T1)$tmpM是可以抽自变量的设计矩阵mit
    g_igt <-c()
    for (t1 in t){
      t_inter <-which.closest (T1,t1)
      # print (t_inter)
      tmpg <-0 
      for (j in T1[1:t_inter]){ # 把t固定住了之后，就看s从小到大的变化
        j_inter <-which.closest (T1,j) # 第几个就积分几个
      
        # 首先顺道产生一下X_i*Bn的积分
        #tmpM <-tmpM+X_i[,1:j_inter] %*%splinex1(t1-T1[1:j_inter]) 
        
        #正常算g_igt
        if (g==1){ 
          tmpg <-tmpg+X_i[,1:j_inter] %*%beta1t(t1-T1[1:j_inter])  
        }
        else {
          tmpg <-tmpg+X_i[ ,1:j_inter] %*%beta2t(t1-T1[1:j_inter])  
        }
      }
      g_igt =as.vector (c(g_igt,tmpg))
    }
    return(g_igt)}

  g1f<-function(t) G_igt(1,t) + as.vector (Z_i%*%alpha1)########mean function for cluster 1
  g2f<-function(t) G_igt(2,t) + as.vector (Z_i%*%alpha2)########mean function for cluster 2
  
  # 计算Vt,要的是行向量
  Vt <-cbind( M_it(T1),as.matrix(rep(1,length(T1) ))%*%Z_i)
  
  # 计算Y
  if(INDEX==1){
    Phi[,1]<-eigenfun11(T1)
    Phi[,2]<-eigenfun12(T1)
    Yt= g1f(T1)+xi[1]*Phi[,1]+xi[2]*Phi[,2]+error1
  }
  else{
      Phi[,1]<-eigenfun21(T1)
      Phi[,2]<-eigenfun22(T1)
      Yt=g2f(T1)+xi[1]*Phi[,1]+xi[2]*Phi[,2]+error1

      }
  data<-list(Y=Yt,VV =Vt,TT=T1,IND=INDEX)
} # sim_data 的结果是给data赋值

sim_data_all<-function(n,ni,PI,gamma,LLamda){ # 时间点是ni[i]
  data_all<-list() # list里面可以什么都不放
  for (i in 1:n){
    sample<-sim_data(ni[i],PI,gamma,LLamda)
    data_all[[i]]<-sample
  }
  return(data_all)
}
n=300######sample size
lambda<-0.05#######penalty parameter
PI<-c(1/2,1/2)###the mixing probabilities
q_n<-6###the number of spline basis functions
num=7########initial number of clusters
K_g<-2#######the number of eigengunctions
sstep<-rep(0,300)
ggamma<-c(0.01,0.04)######## the varainces of errors



# 部分定义均值函数的地方
beta1t <- function(t) (1-t^2) #时变系数 for cluster 1
beta2t <- function(t) sqrt(3)*(cos(pi/2*t)) #时变系数 for cluster 2
# plot(beta2t(seq(0,1,by=0.01)))
# plot(beta1t(seq(0,1,by=0.01)))

alpha1 <- matrix(c(1.4,5.1,3.4), ncol=1 )
alpha2 <- matrix(c(1.1,7.6,2.1), ncol=1 )
Z_i <- matrix(c(rpois(1,70),rbinom(1,1,0.5)+1,runif(1,min=20,max=30)),1,3)

LLamda<-array(0,dim=c(K_g,K_g,num)) # 相当于有聚类个数页的lambda存在

##----特征函数----
eigenfun11<-function(t) sqrt(2)*cos(pi*t)
eigenfun12<-function(t) sqrt(2)*sin(pi*t)
eigenfun21<-function(t) sqrt(2)*t
eigenfun22<-function(t) -3*t+1.5
# plot(eigenfun11(seq(0,1,by=0.01)))

# the eigenvalues类为1
LLamda[1,1,1]<-1.3
LLamda[2,2,1]<-0.3
# the eigenvalues类为2
LLamda[1,1,2]<-1.6
LLamda[2,2,2]<-0.6
for (i in 3:num) { # num是剪枝前的聚类数
  LLamda[1,1,i]<-1
  LLamda[2,2,i]<-0.5
}

#答案----
T11<-seq(0,1,by=0.001)
yita_1 <-ginv(splinex1(T11) )%*%beta1t(T11)  
yita_2 <-ginv(splinex1(T11) )%*%beta2t(T11)
eigencoe_key11 <-ginv(splinex1(T11) )%*%eigenfun11(T11)
eigencoe_key12 <-ginv(splinex1(T11) )%*%eigenfun12(T11)
eigencoe_key21 <-ginv(splinex1(T11) )%*%eigenfun21(T11)
eigencoe_key22 <-ginv(splinex1(T11) )%*%eigenfun22(T11)


# 转化函数
H<-function(t) t####the transformation function in Case 1
SQ<-function(t) t  #inverse function of H

num_group=rep(NA,300)
prob=matrix(NA,300,length(PI))
D_g=array(NA,dim=c(300,q_n+3,length(PI)))
stderror=matrix(NA,300,length(PI))
lambda_g=array(NA,dim=c(300,K_g,K_g,length(PI)))
eigencoe=array(NA,dim=c(300,K_g,q_n,length(PI)))
cluster_result<-matrix(NA,300,800)

ni<-ceiling(runif(n,8,12))####观测时点数，ceiling是超过最大值的那个整数
data<-sim_data_all(n,ni,PI,ggamma,LLamda)

YY<-list()
TT<-list() # 时间点超重要
VV<-list()
IND<-list()
  
#initialization BBeta, SSigma, LLamda, eigencoef, PI_int for the basis coefficient of mean functions, the varainces of errors, eigenvalues, the basis coefficient of eigenfunctions, the mixing probabilities
gamma_int <-rep(0.15 ,num)
BBeta_XY<-matrix(0,nrow=q_n+3,ncol=num) 
BBeta_matrix<-array(0,dim=c(num,q_n+3,q_n+3))  
BBeta<-matrix(0,nrow=q_n+3,ncol=num) 
  
Sigma_XY<-matrix(0,nrow=1,ncol=num)
Sigma_matrix<-matrix(0,nrow=1,ncol=num)
SSigma<-matrix(0,nrow=1,ncol=num)
  
lamda_XY<-matrix(0,nrow=K_g,ncol=num)
lamda_matrix<-matrix(0,nrow=K_g,ncol=num)
eigencoef<-array(0,dim=c(K_g,q_n,num))  

  ##----开始迭代simulate 300 datasets----
  
for(iter in 1:200){
  set.seed(2*iter)
  cat("\n正在重复第",iter,"次XD==================================================\n")
    
  ####################generating data
  ni<-ceiling(runif(n,10,12))####观测时点数，ceiling是超过最大值的那个整数
  data<-sim_data_all(n,ni,PI,ggamma,LLamda)
  n<-length(data)
  for (i in 1:n) {
    YY[[i]]<-data[[i]]$Y######the functional observation
    TT[[i]]<-data[[i]]$TT####obvserved time
    VV[[i]]<-data[[i]]$VV # 确保前面模拟的时候VV已经算出来了
    IND[[i]]<-data[[i]]$IND####true clutse  r id of the data with the individuals in the same order
  }
  true_index<-rep(NA,n)
  
  YYlen<-rep(NA,400)
  for (i in 1:length(YY)) YYlen[i] <-length(unlist(YY[[i]]))
  lenmin<-YYlen[which.min(YYlen[1:400])] # 取出最短和最长的行数
  lenmax<-YYlen[which.max(YYlen[1:400])]
  
  ymm<-matrix(0,nrow=length(YY),ncol=lenmin) # Y赋给矩阵
  for (i in 1:length(YY))
  { ymm[i,]=YY[[i]][1:lenmin] 
  true_index[i]=IND[[i]]
  }
  
  kgroup<-kmeans(ymm,num)
  PI_int<-kgroup$size/sum(kgroup$size)
# }  

  # 先算一个初始值
  for (k in 1:num){ # 对于每一个类考虑
    IND<-kgroup$cluster==k # 是哪一个用==不就行了吗
    Ygroup<-YY[IND] 
    Tgroup<-TT[IND]
    Vgroup<-VV[IND]
    for (i in 1:length(Ygroup)){
      Phi<-matrix(0,length(Tgroup[[i]]),K_g)
      Phi[,1]<-eigenfun11(Tgroup[[i]])
      Phi[,2]<-eigenfun12(Tgroup[[i]])
      Delta<-matrix(nrow=length(Tgroup[[i]]),ncol=length(Tgroup[[i]]))
      Delta<-Phi%*%LLamda[,,k]%*%t(Phi)+gamma_int[k]*diag(1,length(Tgroup[[i]]))
      
      BX<-splinex1(Tgroup[[i]])
      
      #cat('YY和VV的时间点是否是一致的',length(Ygroup[[i]])==nrow(Vgroup[[i]]))
      
      BBeta_XY[,k]<-BBeta_XY[,k]+t(Vgroup[[i]])%*%solve(Delta)%*%(H(Ygroup[[i]]))
      BBeta_matrix[k,,]<-BBeta_matrix[k,,]+t(Vgroup[[i]])%*%solve(Delta)%*%Vgroup[[i]]
    }
    BBeta[,k]=ginv(BBeta_matrix[k,,])%*%BBeta_XY[,k]
  } # 这是估计α的吧
  
  samp=seq(0.001,0.999,length=100)
  desig<-splinex1(samp)
  project<-solve(t(desig)%*%desig)%*%t(desig)
  
  for (k in 1:num){
    eigencoef[1,,k]<-project%*%eigenfun11(samp)
    eigencoef[2,,k]<-project%*%eigenfun12(samp)
  }
  for (k in 1:num){
    SSigma[k]<-ggamma[1]
  }
  
  # 初始化η
  #yita_g1 <-solve(t(splinex1(samp)) %*%splinex1(samp))%*%( t(splinex1(samp))%*% beta1t(samp))
  #for (k in 1:num){
  #  yita_g[1,,k]<-yita_g1 # 新的初始化yita_g
  #}
  
  # Omega里面的所有参数的初始化
  Group_NUM<-rep(NA,300)
  PI_MATRIX<-list()
  ALPHABETA<-list()
  ERRORSIGMA<-list()
  EIGENLAMDA<-list()
  EIGENCOEFFI<-list()
  HY<-list()
  Group_NUM[[1]]=num
  PI_MATRIX[[1]]=PI_int
  ALPHABETA[[1]]=BBeta
  ERRORSIGMA[[1]]=SSigma
  EIGENLAMDA[[1]]=LLamda
  EIGENCOEFFI[[1]]=eigencoef
  
  dif=10000
  likelibefore=100
  delta<-matrix(NA,n,3)
  
  #开算----
  #到这的话要把这个EM_estimate函数重新定义了
  for(i in 1:4) #10
  {
    ##估计Omega_n given H----
    EM_est<-EM_estimate(TT,YY,VV,Group_NUM[i],PI_MATRIX[[i]],ALPHABETA[[i]],ERRORSIGMA[[i]],EIGENLAMDA[[i]],EIGENCOEFFI[[i]],HY[[i]],lambda)
    Group_NUM[i+1]<-EM_est$Group_number
    PI_MATRIX[[i+1]]<-EM_est$PI
    ALPHABETA[[i+1]]<-EM_est$beta
    ERRORSIGMA[[i+1]]<-EM_est$sigma
    EIGENLAMDA[[i+1]]<-EM_est$lamda
    EIGENCOEFFI[[i+1]]<-EM_est$eigencoefficient
    delta<-EM_est$del
    sstep[iter]=i
    #print(c(EM_est$Group_number,EM_est$PI,EM_est$beta,EM_est$sigma,EM_est$lamda,EM_est$eigencoefficient,EM_est$del))
    #cat('聚类数',EM_est$Group_number,'\n')
    cat('隶属概率',EM_est$PI,'\n')
    cat('参数D_g的最大估计误差为',max(EM_est$beta[1:6,1]-yita_1,EM_est$beta[1:6,2]-yita_2),'\n')
    cat('lamda为',EM_est$lamda,'\n')
    
    #############计算各个样本属于每一类的概率
    delta<-matrix(0,nrow=length(YY),ncol=EM_est$Group_number)
    f_den<-matrix(0,nrow=length(YY),ncol=EM_est$Group_number)
    lAlpha<-EM_est$eigencoefficient
    lLamda<-EM_est$lamda
    lgamma_int<-EM_est$sigma
    lBeta<-EM_est$beta
    lPI_int<-EM_est$PI
    likeliafter=0
    
    for (j in 1:length(YY)){
      XXX<-splinex1(TT[[j]])
      Phi<-matrix(NA,length(TT[[j]]),K_g)
      MU<-matrix(NA,length(TT[[j]]),1)
      # lHYYY<-lHYY[identify(YY[[j]],allpoint)]
      lHYYY<- YY[[j]]
      for (k in 1:EM_est$Group_number){
        lPhi<-XXX%*%t(lAlpha[,,k])
        Delta<-matrix(nrow=length(TT[[j]]),ncol=length(TT[[j]]))
        Delta<-lPhi%*%lLamda[,,k]%*%t(lPhi)+lgamma_int[k]*diag(1,ni[j])##
        
        lMU<-VV[[j]]%*%lBeta[,k] # 这会出错吗----
        f_den[j,k]<-(2*pi)^(-length(TT[[j]])/2)/sqrt(det(Delta))*exp(-1/2*t(lHYYY-lMU)%*%solve(Delta)%*%(lHYYY-lMU))
      }
    }
    
    for (j in 1:length(YY)){
      for (k in 1:EM_est$Group_number){
        delta[j,k]<-lPI_int[k]*f_den[j,k]/((t(lPI_int)%*%f_den[j,]))
        likeliafter=likeliafter+delta[j,k]*log(f_den[j,k]%*%lPI_int[k])
      }
    }
    
    # cat("iter=",iter,"\n")
    #----EM参数估计完了----
    #H_diff=mean(abs(HY[[i+1]][51:120]-HY[[i]][51:120]))/mean(abs(HY[[i]][51:120]))
    if(Group_NUM[i+1]-Group_NUM[i]==0){dif=mean(c(mean(abs(sort(PI_MATRIX[[i+1]])-sort(PI_MATRIX[[i]]))),
                                                  mean(abs(ALPHABETA[[i+1]]-ALPHABETA[[i]])),
                                                  mean(abs(ERRORSIGMA[[i+1]]-ERRORSIGMA[[i]])),
                                                  mean(abs(EIGENLAMDA[[i+1]]-EIGENLAMDA[[i]])),
                                                  mean(abs(EIGENCOEFFI[[i+1]]-EIGENCOEFFI[[i]]))))} #i最多只有4？？
    # sort函数是将其排序，dif是显示两次迭代之间各个参数的平均差异
    if(dif<1e-3){ break }
  }
  ##----the end of iteration----
  ##----the estimated number of clusters----
  num_group[iter]=Group_NUM[i+1]#######the estimated number of clusters
  if(Group_NUM[i+1]!=2){
    prob[iter,]=NA
    D_g[iter,,]=NA
    stderror[iter,]=NA
    lambda_g[iter,,,]=NA
    eigencoe[iter,,,]=NA
  }
  if(Group_NUM[i+1]!=2){next}
  ##----the Bayes' optimal allocation rule----
  categorial<-rep(NA,n)
  for (j in 1:n)
  {
    categorial[j]<-which.max(delta[j,])
  }
  cluster_result[iter,1:300]<-true_index#######the true clutser id
  cluster_result[iter,301:600]<-categorial#######the estimated clutser id
  ##############the estimaed results
  prob[iter,]=PI_MATRIX[[i+1]]#######the probabilities
  D_g[iter,,]=ALPHABETA[[i+1]]#######the basis coefficient of mean functions每一个样本都一样吗？？
  stderror[iter,]=ERRORSIGMA[[i+1]]#######the varainces of errors
  lambda_g[iter,,,]=EIGENLAMDA[[i+1]]#######the diagonal matrix of eigenvalue这个估计出来怎么是错的
  eigencoe[iter,,,]=EIGENCOEFFI[[i+1]]#######the basis coefficient of eigenfunctions
}


  
  
  
  
  
#结束执行----

  
  
  
  
  
  
  
  
  
# 答案检查----
RMSE <-function(x) sqrt(mean((x)^2))

D_g[1,,1][1:6]-yita_1
D_g[2,,2][1:6]-yita_2
lambda_g[1,1,1,1]-LLamda[1,1,1]
lambda_g[1,1,1,2]-LLamda[1,1,2]
lambda_g[1,2,2,1]-LLamda[2,2,1]
lambda_g[1,2,2,2]-LLamda[2,2,2]

mean(lambda_g[1,1,1,2] - LLamda[1,1,2])
RMSE(eigencoe[1,1,,1] -t(eigencoe_key11))#相当于是查看第一个样本的，注意g和k是反过来的这

stderror[1,]-ggamma
prob[1,]
all.equal(cluster_result[iter,1:300], cluster_result[iter,301:600])

save.image(file = "SIMULATION_result_1.RData")
#5月14日，设置了估计误差，终于看懂这些是在干啥了，每一个估计十次，逐渐减小误差
#5月15日，重点盯一下lambda的误差，我严重怀疑是中间算错了.现在Alpha_matrix[,,g,k]一直是缺失值

bias <- mean(stderror[1,]-ggamma)
se <- sd(m2 - m1)
rmse <- sqrt(mean((m2 - m1)^2))

