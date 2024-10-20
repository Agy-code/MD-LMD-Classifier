# Library needed: 
   # 1. robustbase: for MCD,
   # 2. RobStatTM: for MVE
   # 3. gam,VGAM : For gam fitting

#-----------(Global) Mahalanobis Distance (MD) Classifier-------------------


MD_classifier<-function(train.X,train.Y,test.X,test.Y,sigma_type)
{

# Input:
   # train.X<- features of the training data.
   # train.Y<- class information of the training data.
   # test.X<-  features of the test data.
   # test.Y<- class information of the training data.
   # sigma_type <- This is the type of scatter matrix, that to be used.
   

# Transformation of the classes as {0,1,...,J}, particularly required for Binary classification:
    
train.Y<-as.numeric(as.factor(train.Y))-1
test.Y<-as.numeric(as.factor(test.Y))-1

# full data:

train_data<-data.frame(class=train.Y,train.X)
test_data<-data.frame(class=test.Y,test.X)

# Dimension:
d<-ncol(train.X)

# Class informations:
  g<-sort(unique(train.Y)) # What are the classes?
  J<-length(g)             # No. of classes

# To make groups of the training data based on the classes:

  X<-lapply(g, function(j) (train.X[train.Y == j, ]))


#List of location parameters:
  means<-lapply(g, function(j) colMeans(X[[j+1]]))

#List of inverted scatter matrices:

 if (sigma_type == "varcov") {
    #Usual variance covariance matrix
    inv_sigma <- lapply(g, function(j) solve(cov(X[[j+1]])))
    
  } else if (sigma_type == "MCD") {
    # Robust covariance matrix using robustbase.
    library(robustbase)
    inv_sigma <- lapply(g, function(j) solve(covMcd(X[[j+1]],alpha=0.75)$cov)) 
    
  }  else if (sigma_type == "MVE") {
    # Robust covariance matrix using RobStatTM.
    library(RobStatTM)
    inv_sigma <- lapply(g, function(j) solve(fastmve(as.matrix(X[[j+1]]))$cov)) 
    
  } else if (sigma_type  == "diagonal") {
    # Diagonal of the covariance matrix (keeping only diagonal elements)
    inv_sigma <- lapply(g, function(j) diag(1/diag(cov(X[[j+1]]))))
    
  } else if (sigma_type  == "identity") {
   # Identity matrix of size d x d.
    inv_sigma <-  replicate(J, diag(d), simplify=FALSE)
  } else {
    stop("Please choose 'varcov', 'MCD', 'MVE', 'diagonal', or 'identity'.")
  }
  
# Mahalanobis distances based on the estimated location and scatter matrices:

  MD.train<-lapply(g, function(j) sqrt(mahalanobis(train.X,center=means[[j+1]],cov=inv_sigma[[j+1]],inverted=TRUE)))
  MD.test<- lapply(g, function(j)  sqrt(mahalanobis(test.X,center=means[[j+1]],cov=inv_sigma[[j+1]],inverted=TRUE)))

# New train and test data, considering MD's as features:

train<-data.frame(train.Y,MD.train)
colnames(train)<-c("class",paste("MD",1:J,sep=""))
test<-data.frame(class=test.Y,MD.test)
colnames(test)<-c("class",paste("MD",1:J,sep=""))

# Spline fitting for GAM:
features<-paste0("s(MD", 1:J, ")", collapse = " + ")

# Model fitting and to find the predicted classes for the test data:
if(J==2){
  library(gam); mod<-gam(as.formula(paste("as.factor(class)~", features)) , 
                          family = "binomial", data = train)  
  new<-data.frame(test[,2:(J+1)])
  prob<-predict(mod,newdata=new,type="response")
  predict<-ifelse(prob>0.50,1,0)}else{

  library(VGAM); mod<-vgam(as.formula(paste("as.factor(class) ~", features)) ,
                           family = multinomial, data = train)
  new<-data.frame(test[,2:(J+1)])
  prob<-predict(mod,newdata=new,type="response")
  predict<-apply(prob,1,which.max)-1
  }


# Confusion matrix:

  T<-table(predict,test.Y)
  accuracy<-(sum(diag(T))/sum(T))
  error<-(1-accuracy)

return(error)
}

