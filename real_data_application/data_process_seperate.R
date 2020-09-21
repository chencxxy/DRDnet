setwd("/Users/chixiangchen/Desktop/phd_in_psu/research/Prof._Wu/Network/tissue/blood_vessel/")

#install.packages("mvtnorm")
library(mvtnorm)
library(MASS)
set.seed(1324)

###############################################
#analysis 2: do log transformation for the data
###############################################

#######################################
#useful functions in functional mapping
#######################################
#calculate P_ij in the paper
P_ij<-function (y_i,x_i,beta,sigma,omega)
{
  beta_sigma<-rbind(beta,sigma)
  whole<-omega*apply(beta_sigma,2,function (x) exp(-t(y_i-x_i%*%x[1:2])%*%(y_i-x_i%*%x[1:2])/2/x[3])/(sqrt(2*pi*x[3]))^(length(y_i))  )     
  #as.matrix(whole/sum(whole),ncol=1)
  whole
}

#calcuate P matrix
P_matrix<-function (x) {
  Y<-x
  x_i<-cbind(1,X[!is.na(Y)])
  y_i<-Y[!is.na(Y)]
  p_ij_1<-P_ij(y_i,x_i,beta,sigma,omega)
  p_ij_1
}

#calculate one part in mle
fct_A<-function (x)
{
  Y<-x[-1]
  P_j_i<-x[1]
  x_i<-cbind(1,as.numeric(X[!is.na(Y)]))
  y_i<-Y[!is.na(Y)]
  t(P_j_i*t(y_i)%*%x_i)
}

#calculate another part in mle
fct_B<-function (x)
{
  Y<-x[-1]
  P_j_i<-x[1]
  x_i<-cbind(1,as.numeric(X[!is.na(Y)]))
  y_i<-Y[!is.na(Y)]
  P_j_i*t(x_i)%*%x_i
}

#calculate sigma_square
fct_sigma<-function (x)
{
  Y<-x[-1]
  P_j_i<-x[1]
  x_i<-cbind(1,as.numeric(X[!is.na(Y)]))
  y_i<-Y[!is.na(Y)]
  numerator<-P_j_i*t(y_i-x_i%*%beta_new[,j])%*%(y_i-x_i%*%beta_new[,j])
  numerator
}

#recalculate the row mean
mean_nonzero<-function (x)
{
  xx<-x[!is.na(x)]
  mean(xx)
}


#read the data
data_vessel<-read.csv("data_vessel.csv")
head(data_vessel)
rownames(data_vessel)<-data_vessel[,1]
data_vessel<-data_vessel[,-1]
dim(data_vessel)

#delete gene with all zeros
row_mean<-apply(data_vessel,1, mean)
data_vessel<-data_vessel[row_mean>0,]
dim(data_vessel)

#let the gene expression equal to zero if it's less than 1 level:
#row_mean<-apply(data_111YS,1, function (x) sum(ifelse(x<1 & x>0,1,0)))
#summary(row_mean)
#data_111YS<-data_111YS[which(row_mean==0),]
data_vessel[data_vessel<=1]<-0
dim(data_vessel)
which(data_vessel<=1 & data_vessel>0)

#count zeros in each row:
zero_row<-apply(data_vessel,1,function (x) length(which(x==0)))
summary(zero_row)

#do log transformation
data_vessel[data_vessel!=0]<-log(data_vessel[data_vessel!=0])
dim(data_vessel)

#set threshold about zeros:
data_vessel<-data_vessel[zero_row<=0,]
dim(data_vessel)
zero_row_del<-apply(data_vessel,1,function (x) length(which(x==0)))
sum(zero_row_del)
#data_111YS[zero_row_del>0 & row_mean>5,]

#summary statistics for row means:
row_mean<-apply(data_vessel,1, mean)
summary(row_mean)
par(mfrow=c(1,1))
sum(row_mean>=8.2)/length(row_mean)
qqnorm(y=row_mean)
#data_vessel[row_mean>=8.2,1]

#read phenotype data
phenotype<-read.csv("phenotype.csv")
#phenotype[1:4,1:10]
name_pheno<-c("SUBJID","AGE","HGHT","BMI","MHASTHMA","MHCOPD","MHHRTATT","MHHRTDIS",
              "MHT1D","MHT2D","MHSMKSTS","MHDRNKSTS","SEX","MHHTN")
pheno_used<-phenotype[,colnames(phenotype) %in% name_pheno]
pheno_used<-pheno_used[-which(pheno_used$SUBJID=="K-562"),]
dim(pheno_used)
summary(pheno_used)
pheno_used$MHDRNKSTS<-ifelse(pheno_used$MHDRNKSTS %in% "Yes",1,ifelse(pheno_used$MHDRNKSTS %in% "No",0,NA))
pheno_used$MHSMKSTS<-ifelse(pheno_used$MHSMKSTS %in% "Yes",1,ifelse(pheno_used$MHSMKSTS %in% "No",0,NA))
pheno_used$SEX<-pheno_used$SEX-1
for (i in 6:ncol(pheno_used))
{
  pheno_used[pheno_used[,i]==99 & !is.na(pheno_used[,i]),i]<-NA
}

for (i in 2:ncol(pheno_used))
{
  if (sum(is.na(pheno_used[,i]))!=0)
  {
    sample_used<-pheno_used[!is.na(pheno_used[,i]),i]
    pheno_used[is.na(pheno_used[,i]),i]<-sample(sample_used,sum(is.na(pheno_used[,i])),replace = T)
  }
}
fit<-glm(MHHTN~AGE+HGHT+BMI+MHASTHMA+MHCOPD+MHHRTATT+MHHRTDIS+
           MHT1D+MHT2D+MHSMKSTS+MHDRNKSTS+SEX,data=pheno_used,family = binomial(link = "logit"))
length(unique(fitted(fit)))
summary(fitted(fit))
plot(sort(fitted(fit)))

pheno_subject_id<-substring(pheno_used$SUBJID,6,10)
vessel_subject_id<-substring(colnames(data_vessel),6,10)
vessel_subject_id<-ifelse(substring(vessel_subject_id,5,5)==".",substring(vessel_subject_id,1,4),substring(vessel_subject_id,1,5))
pheno_together<-data.frame(cbind(pheno_subject_id,fitted(fit),pheno_used$MHSMKSTS,pheno_used$SEX))
colnames(pheno_together)<-c("id","rate","smoke","sex")
vessel_together<-data.frame(cbind(vessel_subject_id,1))
colnames(vessel_together)<-c("id","haha")
index_data<-merge(x=vessel_together,y=pheno_together,by.x="id",all.x=T,sort=F)

index<-as.numeric(as.matrix(index_data$rate,ncol=1))
sort(diff(sort(index),lag=1))

#order the data by the index:
exp_index<-index[order(index)]
X<-log(exp_index)
data_vessel_order<-data_vessel[,order(index)]
index_data<-index_data[order(index),]
dim(data_vessel_order)
Y<-data_vessel_order
#Y<-data_15ER7_order[order_bind[1:nrow(data_15ER7_order),][,2],]
plot(x=exp_index,y=Y[600,])

################################################################
#select some significant genes (by sirs)
p<-nrow(data_vessel_order)
n<-ncol(data_vessel_order)
corr1<-rep()
oumiga<-matrix(0,n,p)
y<-exp_index
x<-t(Y)
for (jj in 1:n)
{
  ind=rep(0,n)
  ind[y<y[jj]]<-1
  inside<-t(rep(1/n,n))%*%(x*ind)
  oumiga[jj,]<-inside^2
}
corr1<-as.vector(t(rep(1/n,n))%*%oumiga)


#corr1<-abs(cor(x,y))
#d1<-floor(0.4*size)
size<-5000
sis<-cbind(corr1,1:p)
sis<-sis[order(sis[,1],decreasing = TRUE),2][1:size]
sis<-sis[order(sis,decreasing = FALSE)]
#lengthsis<-length(sis)
xsis<-x[,sis]
id_sub<-60
plot(x=exp_index,y=xsis[,id_sub],col="red")
plot(x=exp_index[index_data$sex==0],y=xsis[index_data$sex==0,id_sub],col="red")
points(x=exp_index[index_data$sex==0],y=xsis[index_data$sex==0,id_sub],col="blue")
################################################################



##screening by spearman correlation (seperate screening based on smoking groups)
#do smoking group 
x_original<-t(Y)
data_vessel_order_smoking<-data_vessel_order[,index_data$smoke==1]
exp_index_smoking<-exp_index[index_data$smoke==1]
Y_smoking<-Y[,index_data$smoke==1]
p<-nrow(data_vessel_order_smoking)
n<-ncol(data_vessel_order_smoking)
y<-exp_index[index_data$smoke==1]
x<-t(Y_smoking)
spearman<-rep()
for (i in 1:ncol(x))
{
  spearman_i<-abs(cor(x=exp_index_smoking, y=x[,i], method = 'spearman'))
  spearman<-c(spearman,spearman_i)
}
size<-40
sis_smoking<-cbind(spearman,1:p)
sis_smoking<-sis_smoking[order(sis_smoking[,1],decreasing = TRUE),2][1:size]
#sis<-sis[order(sis,decreasing = FALSE)]
#xsis_smoking<-x[,sis_smoking]


#do non-smoking group
data_vessel_order_nonsmoking<-data_vessel_order[,index_data$smoke==0]
exp_index_nonsmoking<-exp_index[index_data$smoke==0]
Y_nonsmoking<-Y[,index_data$smoke==0]
p<-nrow(data_vessel_order_nonsmoking)
n<-ncol(data_vessel_order_nonsmoking)
y<-exp_index[index_data$smoke==0]
x<-t(Y_nonsmoking)
spearman<-rep()
for (i in 1:ncol(x))
{
  spearman_i<-abs(cor(x=exp_index_nonsmoking, y=x[,i], method = 'spearman'))
  spearman<-c(spearman,spearman_i)
}
size<-40
sis_nonsmoking<-cbind(spearman,1:p)
sis_nonsmoking<-sis_nonsmoking[order(sis_nonsmoking[,1],decreasing = TRUE),2][1:size]
#sis<-sis[order(sis,decreasing = FALSE)]
#xsis_nonsmoking<-x[,sis]
#dim(xsis_nonsmoking)

#combine two data together
index_utilized<-unique(c(sis_smoking,sis_nonsmoking))
xsis<-x_original[,index_utilized]
dim(xsis)

#plot(x=exp_index_nonsmoking,y=xsis_nonsmoking[,id_sub],col="blue")
id_sub<-61
par(mfrow = c(1, 2))
plot(x=exp_index[index_data$smoke==1],y=xsis[index_data$smoke==1,id_sub],col="red")
plot(x=exp_index[index_data$smoke==0],y=xsis[index_data$smoke==0,id_sub],col="blue")


dim(data_observe)
write.csv(colnames(data_observe),"72_gene_information.csv")

