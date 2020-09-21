library(np)
library(nlsrk)
library(splines2)
library(grpreg)
library(mvtnorm)


#### setups
n <- 100 # 100,200,300
sigma_square<-0.1 # 0.1, 0.5

lambda.min<-0.05
#lambda.min<-0.01 #for n=100
par(mfrow = c(3, 3))

#################################################################
#generate data
#################################################################
# number of covariates
pn <- 21
# number of time points
iter<-200
alpha<-1
rho<-0 
beta_index<-c(-0.5,-1,1.5,-1)
beta_cov<-c(0)
abs_bias_j<-rep(0,pn)
mse_j<-rep(0,pn)
mse_1_j<-rep(0,pn)
mse_0_j<-rep(0,pn)
mse_diff_j<-rep(0,pn)
mse<-rep()
mse_1<-rep()
mse_0<-rep()
mse_diff<-rep()
abs_bias<-rep()

count_edge_j<-rep(0,pn)
count_edge_total_j<-rep(0,pn)
count_edge_j_interaction<-rep(0,pn)
count_edge_total_j_interaction<-rep(0,pn)
count_edge<-rep()
count_edge_total<-rep()
count_edge_interaction<-rep()
count_edge_total_interaction<-rep()
main_z<-rep()
int_cept<-rep()


#construct correlation matrix
X.sigma <- matrix(0,(pn),(pn))
{
  for (i in 1:(pn))
    X.sigma[i,] <- sigma_square*rho^(abs((1:(pn))-i))
}

sigma_square_index<-1
X.sigma_index <- matrix(0,3,3)
{
  for (i in 1:(3))
    X.sigma_index[i,] <- sigma_square_index*rho^(abs((1:(3))-i))
}

#t=seq(0,20,length=n)
#t=round(t,6)
set.seed(4351)
timestart<-Sys.time()
for (ss in 1:iter)
{
  #generate the true index, disease status, and fitted index
  x_for_index<-rmvnorm(n, mean = rep(0,3), X.sigma_index)
  pi_true<-exp(cbind(1,x_for_index)%*%beta_index)/(1+exp(cbind(1,x_for_index)%*%beta_index)) #true index
  hist(pi_true)
  status<-rbinom(n,size=1,prob=pi_true) # disease status
  fit_index<-glm(status~x_for_index,family = binomial) 
  #predict(fit_index,type = "response")
  fit_value<-fitted(fit_index) #fitted index
  mean(20*abs(pi_true-fit_value))
  
  #generate true ode solutions
  t_true<-20*sort(pi_true)
  t_true[1]<-0
  t_true[n]<-20
  t<-20*sort(fit_value)
  t[1]<-0
  t[n]<-20
  length(unique(t_true))
  length(unique(t))
  min(diff(t,lag=1))
  plot(t_true)
  lines(t,col="red")
  #generate true ode solution
  fn<-function()
  { 
    h<-function(t,y,param) param[1]+param[2]*y[1]+param[3]*y[1]^2+param[4]*y[1]^3+param[5]*y[2]+param[6]*y[2]^2+param[7]*y[2]^3
    d<-function(t,y,param) param[8]+param[9]*y[1]+param[10]*y[1]^2+param[11]*y[1]^3+param[12]*y[2]+param[13]*y[2]^2+param[14]*y[2]^3
    c(h,d)
  }
  #parameter<-c(1,-0.01,1,-0.5)
  #parameter<-c(2,-0.02,1,-1)  #good for g2
  #parameter<-c(0.5,-0.01,0.5,1)  #still good for g2 not g1
  #parameter<-c(0,0,0,0,0.1,0.2,0.2,0.4,-2,0,0.4,0,0,0)
  #parameter<-c(-1,0,0,0,0.1,0.2,-0.2,0.4,-0.5,0,0.4,0,0,0)
  parameter<-c(0, 1.2,0.3,-0.6, 0.1,0.2,0.2,   0.4, -2,0,0.4, 0.5,0.2,-0.3)
  
  ode<-evrunge(t=t_true,param=parameter,y0=c(-2,2),sys=fn,dt=min(diff(t_true,lag=1))/2)
  mu1<-ode[1:n] 
  mu2<-ode[(n+1):(2*n)]
  plot(mu1)
  plot(mu2)
  
  #parameter<-c(1,-0.01,1,-0.5)
  #parameter<-c(2,-0.02,1,-1)  #good for g2
  #parameter<-c(0.5,-0.01,0.5,1)  #still good for g2 not g1
  parameter<-c(-0.2, 0,0,0, -0.3,0.4,0.1,   -0.2, 0.2,-0.1,-0.2, 0,0,0)
  ode<-evrunge(t=t_true,param=parameter,y0=c(2,-2),sys=fn,dt=min(diff(t_true,lag=1))/2)
  mu3<-ode[1:n] 
  mu4<-ode[(n+1):(2*n)]
  plot(mu3)
  plot(mu4)
  
  #parameter<-c(1,-0.01,1,-0.5)
  #parameter<-c(2,-0.02,1,-1)  #good for g2
  #parameter<-c(0.5,-0.01,0.5,1)  #still good for g2 not g1
  parameter<-c(0.05, 0,0,0, 0.1,0,-0.8,  -0.05, 0,0,0.5, 0,0,0)
  ode<-evrunge(t=t_true,param=parameter,y0=c(-1.4,1.4),sys=fn,dt=min(diff(t_true,lag=1))/2)
  mu5<-ode[1:n] 
  mu6<-ode[(n+1):(2*n)]
  plot(mu5)
  plot(mu6)
  
  #not very power
  fn<-function()
  { 
    h<-function(t,y,param) param[1]
    c(h)
  }
  #parameter<-c(1,-0.01,1,-0.5)
  #parameter<-c(2,-0.02,1,-1)  #good for g2
  #parameter<-c(0.5,-0.01,0.5,1)  #still good for g2 not g1
  parameter<-c(0)
  ode<-evrunge(t=t_true,param=parameter,y0=c(0),sys=fn,dt=min(diff(t_true,lag=1))/2)
  mu7<-ode[1:n] 
  #plot(mu7)
  
  
  
  #not very power
  fn<-function()
  { 
    h<-function(t,y,param) param[1]+param[2]*sin(0.5*y[1])*y[1]
    d<-function(t,y,param) param[3]+param[4]*cos(0.8*y[1])*abs(y[1])
    m<-function(t,y,param) param[5]+param[6]*y[1]+param[7]*y[1]^2+param[8]*y[1]^3+param[9]*y[2]+param[10]*y[2]^2+param[11]*y[2]^3
    c(h,d,m)
  }
  #parameter<-c(1,-0.01,1,-0.5)
  #parameter<-c(2,-0.02,1,-1)  #good for g2
  #parameter<-c(0.5,-0.01,0.5,1)  #still good for g2 not g1
  #parameter<-c(-0.5,1,0.5,1)
  #parameter<-c(-0.1,0.2,0.1,0.5,-1.5,0.2,0,-0.1,-0.3,0,-0.15)
  #parameter<-c(-0.1,0.2,-0.1,0.5, -0.5, 0.05,0,-0.005, -0.015,0,-0.01)
  #parameter<-c(-0.1,0.2,-0.1,0.5, -0.5, 0.05,0,-0.005, -0.5,0.6,-0.05)
  #parameter<-c(-0.1,0.2,-0.1,0.5, -0.5, 0,0,0, -0.5,0.6,-0.05)
  parameter<-c(-0.1,0.2, -0.1,0.5,  -0.5, 0.05,0,-0.005, -0.5,0.6,-0.05)
  parameter<-c(0,1.5, -0.1,0.5,  -0.5, 0.05,0,-0.005, -0.5,0.6,-0.05)
  
  
  ode<-evrunge(t=t_true,param=parameter,y0=c(0.1,-0.1,-1),sys=fn,dt=min(diff(t_true,lag=1))/2)
  mu9<-ode[1:n] 
  mu10<-ode[(n+1):(2*n)]
  mu8<-ode[(2*n+1):(3*n)]
  plot(mu8)
  plot(mu9)
  plot(mu10)
  
  #not very power
  fn<-function()
  { 
    h<-function(t,y,param) param[1]+param[2]*y[1]+param[3]*y[1]^2+param[4]*y[1]^3
    c(h)
  }
  #parameter<-c(1,-0.01,1,-0.5)
  #parameter<-c(2,-0.02,1,-1)  #good for g2
  #parameter<-c(0.5,-0.01,0.5,1)  #still good for g2 not g1
  #parameter<-c(-0.5,1,0.5,1)
  #parameter<-c(-0.1,0.2,0.1,0.5,-1.5,0.2,0,-0.1,-0.3,0,-0.15)
  #parameter<-c(-0.1,0.2,-0.1,0.5, -0.5, 0.05,0,-0.005, -0.015,0,-0.01)
  #parameter<-c(-0.1,0.2,-0.1,0.5, -0.5, 0.05,0,-0.005, -0.5,0.6,-0.05)
  #parameter<-c(-0.1,0.2,-0.1,0.5, -0.5, 0,0,0, -0.5,0.6,-0.05)
  #parameter<-c(-0.1,0.2, -0.1,0.5,  -0.5, 0.05,0,-0.005, -0.5,0.6,-0.05)
  parameter<-c(0.98,-1.56,0.91,-0.142)
  ode<-evrunge(t=t_true,param=parameter,y0=c(-2),sys=fn,dt=min(diff(t_true,lag=1))/2)
  mu11<-ode[1:n] 
  plot(mu11)
  
  #generate data used for constructing the network
  all_nodes<-rmvnorm(n, mean = rep(0,pn), X.sigma)
  x_cov<-rbinom(n,1,0.5)
  mu1_total<-mu1+as.matrix(x_cov,ncol=1)%*%beta_cov
  mu2_total<-mu2+as.matrix(x_cov,ncol=1)%*%beta_cov
  mu3_total<-mu3+as.matrix(x_cov,ncol=1)%*%beta_cov
  mu4_total<-mu4+as.matrix(x_cov,ncol=1)%*%beta_cov
  mu5_total<-mu5+as.matrix(x_cov,ncol=1)%*%beta_cov+x_cov*(mu9)
  mu6_total<-mu6+as.matrix(x_cov,ncol=1)%*%beta_cov
  mu7_total<-mu7+as.matrix(x_cov,ncol=1)%*%beta_cov+x_cov*(mu11)
  mu8_total<-mu8+as.matrix(x_cov,ncol=1)%*%beta_cov
  mu9_total<-mu9+as.matrix(x_cov,ncol=1)%*%beta_cov
  mu10_total<-mu10+as.matrix(x_cov,ncol=1)%*%beta_cov
  mu11_total<-mu11+as.matrix(x_cov,ncol=1)%*%beta_cov
  x1<-mu1_total+all_nodes[,1]
  x2<-mu2_total+all_nodes[,2]
  x3<-mu3_total+all_nodes[,3]
  x4<-mu4_total+all_nodes[,4]
  x5<-mu5_total+all_nodes[,5]
  x6<-mu6_total+all_nodes[,6]
  x7<-mu7_total+all_nodes[,7]
  x8<-mu8_total+all_nodes[,8]
  x9<-mu9_total+all_nodes[,9]
  x10<-mu10_total+all_nodes[,10]
  x11<-mu11_total+all_nodes[,11]
  
  
  
  
  data_observe<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,(all_nodes[,12:pn]+0.01*t_true))
  data_observe<-data_observe[c(1,diff(t,lag=1))>=0.0001,]
  dim(data_observe)
  x_cov<-x_cov[c(1,diff(t,lag=1))>=0.0001]
  
  mean_node<-cbind(mu1_total,mu2_total,mu3_total,mu4_total,mu5_total,mu6_total
                   ,mu7_total,mu8_total,mu9_total,mu10_total,mu11_total,matrix(0,nrow=length(mu1),ncol=(pn-11))+0.01*t_true)[c(1,diff(t,lag=1))>=0.0001,]
  mean_node_0<-cbind(mu1,mu2,mu3,mu4,mu5,mu6
                     ,mu7,mu8,mu9,mu10,mu11,matrix(0,nrow=length(mu1),ncol=(pn-11))+0.01*t_true)[c(1,diff(t,lag=1))>=0.0001,]
  mean_node_1<-cbind(mu1,mu2,mu3,mu4,mu5+mu9,mu6
                     ,mu7+mu11,mu8,mu9,mu10,mu11,matrix(0,nrow=length(mu1),ncol=(pn-11))+0.01*t_true)[c(1,diff(t,lag=1))>=0.0001,]
  mean_node_diff<-mean_node_1-mean_node_0
  
  t<-t[c(1,diff(t,lag=1))>=0.0001]
  #n<-nrow(data_observe)
  #plot(data_observe[,9])
  
  #get fitted values of observations using varying coefficient kernel regression
  data_fitted<-rep()
  data_fitted_cov<-rep()
  for(i in 1:ncol(data_observe))
  {
    bw <- npscoefbw(formula=data_observe[,i]~x_cov|t,betas=T,regtype="ll",bwtype="generalized_nn")
    bw <- npscoef(bw,betas=T,exdat=data.frame(x_cov=x_cov),ezdat=data.frame(t=t),regtype="ll")
    fitted<-coef(bw)[,1]
    data_fitted<-cbind(data_fitted,fitted)
    fitted_cov<-coef(bw)[,2]
    data_fitted_cov<-cbind(data_fitted_cov,fitted_cov)
  }
  plot(y=mean_node_0[,9],x=t)
  lines(data_fitted_cov[,5],x=t)
  plot(y=mean_node[,7],x=t)
  lines(data_fitted[,7]+x_cov*data_fitted_cov[,7],x=t)
  
  plot(y=mean_node_0[,12],x=t)
  lines(data_fitted[,12]+x_cov*data_fitted_cov[,12],x=t)
  
  #spline basis and its integration for baseline
  library(splines2)
  library(Matrix)
  library(grpreg)
  #some useful fct
  step_fct<-function(x,step)
  {
    t<-length(x)
    x_new<-x[1]
    for (i in 1:(t-1))
    {
      x_int<-seq(x[i],x[i+1],by=step)[-1]
      x_new<-c(x_new,x_int)
    }
    x_new
  }
  
  #generate B-spline and its integration for baseline
  #generate B-spline and its integration for cov
  X_big<-rep()
  X_big_cov<-rep()
  X_big_int<-rep()
  X_big_int_cov<-rep()
  degree<-3
  round_num<-attr(regexpr("(?<=\\.)0+", format(min(diff(t,lag=1)),scientific = FALSE), perl = TRUE), "match.length")+1
  t<-round(t,round_num)
  if(min(diff(t,lag=1))>=2*10^(-round_num))
  {
    step<-10^(-round_num)
  }else
  {
    step<-5*10^(-round_num-1)
    t<-round(t,round_num+1)
  }
  cluster<-ncol(data_observe)
  lambda<-seq()
  for (i in 1:cluster)
  {
    knots <- sort(data_fitted[,i])[c(round(nrow(data_fitted)/4),round(nrow(data_fitted)*3/4))]
    knots_cov <-sort(data_fitted_cov[,i])[c(round(nrow(data_fitted_cov)/4),round(nrow(data_fitted_cov)*3/4))]
    bsMat <- bSpline(data_fitted[,i], knots = knots, degree = degree, intercept=F)
    bsMat_cov <- bSpline(data_fitted_cov[,i], knots = knots_cov, degree = degree, intercept=F)
    x_int<-step_fct(as.vector(t),step)
    if(min(diff(t,lag=1))<2*10^(-round_num))
    {x_int<-round(x_int,round_num+1)}else
    {x_int<-round(x_int,round_num)}
    x_cov_int<-rep(0,length(x_int))
    bw <- npscoefbw(formula=data_observe[,i]~x_cov|t,betas=T,regtype="ll",bwtype="generalized_nn")
    bw <- npscoef(bw,betas=T,exdat=data.frame(x_cov=c(x_cov_int)),ezdat=data.frame(t=c(x_int)),regtype="ll")
    x_hat<-coef(bw)[,1]
    x_hat_cov<-coef(bw)[,2]
    
    #x_hat<-exp(cbind(1,log(x_int)) %*% beta[,i])
    basis_int<-bSpline(x_hat, knots = knots, degree = degree, intercept=F)
    basis_int_cov<-bSpline(x_hat_cov, knots = knots_cov, degree = degree, intercept=F)
    base_int<-0
    base_int_cov<-0
    for (j in 1:(length(t)-1))
    {
      int_row<-apply(basis_int[x_int>=t[j] & x_int<=t[j+1],][-1,],2,function (x) sum(x)*step)
      base_int<-rbind(base_int,int_row)
      
      int_row_cov<-apply(basis_int_cov[x_int>=t[j] & x_int<=t[j+1],][-1,],2,function (x) sum(x)*step)
      base_int_cov<-rbind(base_int_cov,int_row_cov)
    }
    base_int<-apply(base_int,2,cumsum)
    base_int_cov<-apply(base_int_cov,2,cumsum)
    #base_int<-ibs(data_fitted[,i], knots = knots, degree = degree, intercept=F)*(max(t)-min(t)) #not correct
    X_big <- cbind(X_big,bsMat)
    X_big_cov <- cbind(X_big_cov,bsMat_cov)
    
    X_big_int<-cbind(X_big_int,base_int)
    X_big_int_cov<-cbind(X_big_int_cov,base_int_cov)
  }
  dim(X_big)
  dim(X_big_int)
  
  dim(X_big_cov)
  dim(X_big_int_cov)
  
  #X_big<-apply(X_big,2,function(x) x/sum(x))
  #X_big_int<-apply(X_big_int,2,function(x) x/sum(x))
  num_cov<-degree+length(knots)
  #X_big_int_exp_intcep<-cbind(1,t,X_big_int,x_cov)
  
  #some setup for final analysis
  X_big_int_exp<-cbind(t,X_big_int,x_cov*X_big_int_cov,t*x_cov,x_cov)
  X_big_int_exp_intcep<-cbind(1,t,X_big_int,x_cov*X_big_int_cov,t*x_cov,x_cov)
  X_big_int_exp_intcep_0<-cbind(1,t,X_big_int,0*X_big_int_cov,t*0,0)
  X_big_int_exp_intcep_1<-cbind(1,t,X_big_int,1*X_big_int_cov,t*1,1)
  
  x_cov<-matrix(x_cov,ncol=1)
  group<-c(0,rep(1:(2*ncol(data_observe)),each=num_cov),rep(0,ncol(x_cov)),rep(0,ncol(x_cov)))
  tt_1<-list(c(8,3),c(3,8),18,13,c(28),23,0,c(47,52),47,47,
             53,0,0,0,0,0,0,0,0,0,0
  )
  tt_2<-list(c(0),c(0),0,0,c(128),0,138,c(0,0),0,0,
             0,0,0,0,0,0,0,0,0,0,0
  )
  for (j in 1:ncol(data_observe))  #
  {
    #cvfit <- cv.grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",nfolds=5,alpha=alpha1)
    #beta_ini<-grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",lambda = cvfit$lambda.min,alpha=alpha1)$beta
    #group_id<-c(1,2,rep(3:(2+ncol(data_observe)),each=num_cov))
    #beta_each_g<-unlist(tapply(beta_ini,group_id,function(x) x))
    #weight_ini<-rep()
    #for(i in 1:ncol(data_observe))
    #{
    #  m<-sum((beta_each_g[(3+(i-1)*num_cov):(3+(i-1)*num_cov+num_cov-1)])^2)
    #  weight_ini<-c(weight_ini,m)
    #}
    #weight<-(1/(weight_ini+0.0000001))^0.3
    cvfit <- cv.grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",nfolds=10,alpha=alpha,max.iter=100000)
    fit<-grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",lambda =cvfit$lambda.min,alpha=alpha,max.iter=100000)
    beta_group<-fit$beta
    #fit<-grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",alpha=alpha,lambda.min = lambda.min)
    #fit<-grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",alpha=alpha,lambda = cvfit$lambda.min)
    #beta_group<-select(fit,"EBIC")$beta
    #cvfit <- cv.grpreg(X_big_int_exp, data_observe[,4], group=group, penalty="grLasso",nfolds=50)
    #fit<-grpreg(X_big_int_exp, data_observe[,4], group=group, penalty="grLasso",lambda = cvfit$lambda.min)
    #fit<-grpreg(X_big_int_exp, data_observe[,2], group=group, penalty="grLasso")
    #beta_group<-select(fit,"EBIC")$beta
    #weight
    
    ###capture the beta_cov and overall 
    if(j %in% c(5,7))
    {
      main_z<-c(main_z,beta_group[length(beta_group)])
      int_cept<-c(int_cept,beta_group[1])
    }
    
    select<-which(beta_group!=0)[-c(1,2,length(which(beta_group!=0))-1,length(which(beta_group!=0)))]
    count_edge_j[j]<-sum(ifelse(select %in% tt_1[[j]],1,0))
    count_edge_total_j[j]<-sum(select %in% c(3:107))/num_cov
    count_edge_j_interaction[j]<-sum(ifelse(select %in% tt_2[[j]],1,0))
    count_edge_total_j_interaction[j]<-sum(select %in% c(108:212))/num_cov
    #beta_group[beta_group!=0]
    fitted_mean<-X_big_int_exp_intcep%*%as.matrix(beta_group,ncol=1)
    fitted_mean_1<-X_big_int_exp_intcep_1%*%as.matrix(beta_group,ncol=1)
    fitted_mean_0<-X_big_int_exp_intcep_0%*%as.matrix(beta_group,ncol=1)
    fitted_mean_diff<-fitted_mean_1-fitted_mean_0
    par(mfrow = c(1, 1))
    plot(y=mean_node_1[,j],x=t,col="red")
    lines(y=fitted_mean_1,x=t,col="red")
    points(y=mean_node_0[,j],x=t,col="blue")
    lines(y=fitted_mean_0,x=t,col="blue")
    mse_j[j]<-mean((fitted_mean-mean_node[,j])^2)
    abs_bias_j[j]<-mean(abs(fitted_mean-mean_node[,j]))
    mse_1_j[j]<-mean((fitted_mean_1-mean_node_1[,j])^2)
    mse_0_j[j]<-mean((fitted_mean_0-mean_node_0[,j])^2)
    mse_diff_j[j]<-mean((fitted_mean_diff-mean_node_diff[,j])^2)
  }
  count_edge<-cbind(count_edge,count_edge_j)
  count_edge_interaction<-cbind(count_edge_interaction,count_edge_j_interaction)
  count_edge_total<-cbind(count_edge_total,count_edge_total_j)
  count_edge_total_interaction<-cbind(count_edge_total_interaction,count_edge_total_j_interaction)
  
  mse<-cbind(mse,mse_j)
  mse_1<-cbind(mse_1,mse_1_j)
  mse_0<-cbind(mse_0,mse_0_j)
  mse_diff<-cbind(mse_diff,mse_diff_j)
  abs_bias<-cbind(abs_bias,abs_bias_j)
  print(ss)
}

####remove outliers which is larger than 5
outlier<-unique(c(which(mse>5,arr.ind = T)[2],which(mse_1>5,arr.ind = T)[2],
                  which(mse_0>5,arr.ind = T)[2],which(mse_diff>5,arr.ind = T)[2]))
if(sum(!is.na(outlier))>=1)
{
  mse<-mse[,-na.omit(outlier)]
  mse_1<-mse_1[,-na.omit(outlier)]
  mse_0<-mse_0[,-na.omit(outlier)]
  mse_diff<-mse_diff[,-na.omit(outlier)]
  count_edge<-count_edge[,-na.omit(outlier)]
  count_edge_interaction<-count_edge_interaction[,-na.omit(outlier)]
  count_edge_total<-count_edge_total[,-na.omit(outlier)]
  count_edge_total_interaction<-count_edge_total_interaction[,-na.omit(outlier)]
}

outlier

#total edges selected
apply(count_edge_total,1,mean)
apply(count_edge_total,1,sd)
sum(apply(count_edge_total,1,mean))
#total for interaction
apply(count_edge_total_interaction,1,mean)
apply(count_edge_total_interaction,1,sd)
sum(apply(count_edge_total_interaction,1,mean))

#fp
apply(count_edge_total-count_edge,1,mean)
apply(count_edge_total-count_edge,1,sd)
#overall fp
mean(apply(count_edge_total-count_edge,2,sum))
sd(na.omit(count_edge_total-count_edge))
#fp for interaction
apply(count_edge_total_interaction-count_edge_interaction,1,mean)
apply(count_edge_total_interaction-count_edge_interaction,1,sd)
#overall fp interaction
mean(apply(count_edge_total_interaction-count_edge_interaction,2,sum))
sd(na.omit(count_edge_total_interaction-count_edge_interaction))

#tp
apply(count_edge,1,mean)
apply(count_edge,1,sd)
#overall tp
mean(apply(count_edge,2,sum))
sd(na.omit(count_edge))
#tp for interaction
apply(count_edge_interaction,1,mean)
apply(count_edge_interaction,1,sd)
#overall tp interaction
mean(apply(count_edge_interaction,2,sum))
sd(na.omit(count_edge_interaction))

#fp rate
fpplustn<-pn-c(2,2,1,1,1,1,0,2,1,1,1,rep(0,pn-11))
apply((count_edge_total-count_edge)/fpplustn,1,mean)
apply((count_edge_total-count_edge)/fpplustn,1,sd)
#overall fp rate
mean(apply(na.omit(count_edge_total-count_edge),2,sum)/(sum(fpplustn)))
sd(apply(na.omit(count_edge_total-count_edge),2,sum)/(sum(fpplustn)))
#fp rate for interaction
fpplustn_interaction<-pn-c(0,0,0,0,1,0,1,0,0,0,0,rep(0,pn-11))
apply((count_edge_total_interaction-count_edge_interaction)/fpplustn_interaction,1,mean)
apply((count_edge_total_interaction-count_edge_interaction)/fpplustn_interaction,1,sd)
#overall fp rate for interaction
mean(apply(na.omit(count_edge_total_interaction-count_edge_interaction),2,sum)/(sum(fpplustn_interaction)))
sd(apply(na.omit(count_edge_total_interaction-count_edge_interaction),2,sum)/(sum(fpplustn_interaction)))


#tp rate
tpplusfn<-c(2,2,1,1,1,1,0,2,1,1,1,rep(0,pn-11))
apply(count_edge/tpplusfn,1,mean)
apply(count_edge/tpplusfn,1,sd)
#overall tp rate
mean(apply(na.omit(count_edge),2,sum)/sum(tpplusfn))
sd(apply(na.omit(count_edge),2,sum)/sum(tpplusfn))
#tp rate for interaction
tpplusfn_interaction<-c(0,0,0,0,1,0,1,0,0,0,0,rep(0,pn-11))
apply(count_edge_interaction/tpplusfn_interaction,1,mean)
apply(count_edge_interaction/tpplusfn_interaction,1,sd)
#overall tp rate for interaction
mean(apply(na.omit(count_edge_interaction),2,sum)/sum(tpplusfn_interaction))
sd(apply(na.omit(count_edge_interaction),2,sum)/sum(tpplusfn_interaction))


#mse for observed data
apply(mse,1,mean)
apply(mse,1,sd)
mean(apply(mse,1,mean))
#overallmse for observed data
mean(mse)
sd(mse)

#absolute bias
apply(abs_bias,1,mean)
apply(abs_bias,1,sd)
mean(apply(abs_bias,1,mean))
#overallbias
mean(abs_bias)
sd(abs_bias)

#for z
#mean(main_z)
#sd(main_z)

#main_z
#int_cept

#mse for individual
#mse_1
apply(mse_1,1,sd)
#mse_0
apply(mse_0,1,sd)
#mse_diff
apply(mse_diff,1,sd)

#mean mse for all
mean(mse_1)
sd(mse_1)
mean(mse_0)
sd(mse_0)
mean(mse_diff)
sd(mse_diff)

#mean mse for each individual 
apply(mse_1,1,mean)
apply(mse_1,1,sd)
apply(mse_0,1,mean)
apply(mse_0,1,sd)
apply(mse_diff,1,mean)
apply(mse_diff,1,sd)


#plot(x=count_edge_total,y=count_edge)
timeend<-Sys.time()
timeend-timestart

