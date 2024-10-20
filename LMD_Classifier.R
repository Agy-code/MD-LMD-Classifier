
# This function calculates LMD of observation(s):

LMD<-function(h,data,train,inv_cov,kern=dnorm)
{
# Input:
#h<-tuning parameter(s).
#data<- a matrix or data frame such that the rows contain query points, for which LMD has to be found.
#train<-a matrix or data frame containing the x1,x2,...xn.
#cov<-the INVERSE of the scatter matrix.
#kern<-This is the kernel to be used, if not sepcified the Gaussian kernel will be considered.

# Making sure the data type is a data frame, where each row is a query point:
 data=data.frame(data)

# Dimension: 
  d<-ncol(train) 

# Finding beta for a single observation, say x:

         B<-function(x)
         {
             MD <- mahalanobis(train,x,inv_cov,inverted=TRUE)  #for all i:(x-x_i)^T \Sigma^{-1}(x-x_i)
             if(length(h)==1)                                  #for a single tuning parameter. 
               {
              k <- kern(MD/(h^2))                              #for all i:K(((x-x_i)^T \Sigma^{-1}(x-x_i))/h^2)
               beta<-mean(MD*k)            
               }else{                                          #for multiple tuning parameters.
              k <- sapply(h,function(j){kern(MD/(j^2))})
              beta <- apply(MD*k,2,mean)
               }
               return(beta)
          }
# Beta for all query points:
  beta_vec <- apply(data,1,B)

# LMD/Gamma:
         if(length(h)==1){
                  gamma<-beta_vec*ifelse(h>1,1,1/h^(d+2))}else{
                  gamma<-beta_vec*sapply(h,function(j){ifelse(j>1,1,1/j^(d+2))})
            }
  return(gamma)                                                  # returns the LMD value(s).
}



# This function provides the error rate for the LMD_classifier:

LMD_classifier<-function(train.X,train.Y,test.X,test.Y,sigma_type,boot=TRUE,n.folds=100,h_given=NULL)
{

options(warn=-1)
# Input:
   # train.X<- features of the training data.
   # train.Y<- class information of the training data.
   # test.X<-  features of the test data.
   # test.Y<- class information of the training data.
   # sigma_type <- This is the type of variance matrix, that to be used.
   # If boot is TRUE, then it performs bootstrap samling to find optimal value of h, based on n.folds bootstrap samples.
   # If boot is FALSE, then it considers the given value of h (as specified by h_given) as the optimal.
 


# The function needs a value of h_given, if boot= FALSE. Otherwise, this will stop.

if (boot==FALSE && is.null(h_given)) {
    stop("Error: No choice for the tuning parameter is given")
  }


# Tranformation of the classes as {0,1,...,J}, particularly required for Binary classification:
    
train.Y<-as.numeric(as.factor(train.Y))-1
test.Y<-as.numeric(as.factor(test.Y))-1

# full data:

train_data<-data.frame(class=train.Y,train.X)
test_data<-data.frame(class=test.Y,test.X)

# Dimension:
d<-ncol(train.X)

# Class informations:
  g<-sort(unique(train.Y)) # What are the classes
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

  
if(boot){

#To calculate the squared generalized Mahalanobis distance between
# all pairs of rows in a data frame with respect to a covariance matrix.

   M<- lapply(g, function(j) sqrt(biotools::D2.dist(train.X,cov=inv_sigma[[j+1]],inverted=TRUE)))

# List of the quantiles of the MD's:
  Q<-lapply(g, function(j) quantile(as.vector(M[[j+1]]),0.05))

# starting the gird for h:
  h.start<-min(unlist(Q))/3

# To find the range of h:

  # Mahalanobis distances:
   MD<- lapply(g, function(j) mahalanobis(train.X,center=means[[j+1]],cov=inv_sigma[[j+1]],inverted=TRUE))
   

  i <- 1
  r <- 0
  hvec <- NULL
  h.loop <- h.start  # Initial value of h

  while (r < 0.99) {
    if (length(hvec) == 15) break
    hvec[i] <- h.loop
    correlations <- numeric(J)

    for (j in g) {
     
      LMD_j <- LMD(h.loop, train.X, X[[j+1]], inv_sigma[[j+1]])
      
      # Computes the correlation between MD and LMD
      correlations [j+1] <- cor(MD[[j+1]], LMD_j)
    }
    
    # Calculate the minimum correlation across all J cases
    r <- min(correlations)
    r[is.na(r)]<-0
    h.loop <- h.loop * (1.25)^(i + 1)
    i <- i + 1
  } 

# Vector of tuning parameters:

  h<-extremefit::bandwidth.grid(hvec[1],3*hvec[length(hvec)], length=50, type = "geometric")

# A function that gives data based on a particular class:

group<-function(g,data) # The first column of data should contain the class information.
{
  class=as.vector(data[,1])
  return(as.matrix(data[class==g,])) # returns a matrix, along with class
}

# "bootstrap" is a function that performs bootstrap sampling to select the tuning parameter, from given hvec vector.

bootstrap<-function(data,hvec) # data should contain the class information, in first column.
 {
  options(warn=-1)
# making sure that given data is a data.frame

  data<-data.frame(data)

  my_list <- list()
   for(k in 1:n.folds)
        {

          # preparation of the bootstrap training and test data:


           tab.prop <- sapply(g, function(j) nrow(group(j, data)) / nrow(data))  # Proportions for each class


           boot_data_list <- lapply(g, function(j) {
            data_j <- group(j, data)
            boot_data_j <- data_j[sample(1:nrow(data_j), replace=TRUE), ]
            return(boot_data_j)
               })

             idx_list <- lapply(1:J, function(j) {
              n_j <- 0.50 * nrow(data) * tab.prop[j]  
             sample(1:nrow(boot_data_list[[j]]), size=n_j, replace=FALSE)
              })


            train_data <- do.call(rbind, lapply(1:J, function(j) {
                                    boot_data_list[[j]][idx_list[[j]], ]
                                    }))

            test_data <- do.call(rbind, lapply(1:J, function(j) {
                                    boot_data_list[[j]][-idx_list[[j]], ]
                                    }))

             #Features: .X, class information: .Y

             train.X<-train_data[,-1]
             train.Y<-train_data[,1]
             test.X<-test_data[,-1]
             test.Y<-test_data[,1]


             # To make groups of the bootstrap training data based on the classes:

               X<-lapply(g, function(j) (train.X[train.Y == j, ]))


             # List of location parameters:
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


             # Calculation of LMD for the given vector hvec, corresponding to the bootstrap training and test data.

              LMD.train<-lapply(g, function(j) LMD(h,data=train.X,train=X[[j+1]],inv_cov=inv_sigma[[j+1]]))
              LMD.test<- lapply(g, function(j) LMD(h,data=test.X,train=X[[j+1]],inv_cov=inv_sigma[[j+1]]))


             # Accuracy of the LMD classifier, across different h values.

              acc<-NULL
              for(i in 1:length(hvec))
                  {
                     LMD.train.h<-sapply(LMD.train, function(x)x[i,],simplify=TRUE)
                     LMD.test.h<-sapply(LMD.test, function(x)x[i,],simplify=TRUE)

                     train.h<-data.frame(train.Y,LMD.train.h)
                     colnames(train.h)<-c("class",paste("LMD",1:J,sep=""))
                     test.h<-data.frame(class=test.Y,LMD.test.h)
                     colnames(test.h)<-c("class",paste("LMD",1:J,sep=""))
                      
                     features<-paste0("s(LMD", 1:J, ")", collapse = " + ")

                               if(J==2){
                                library(gam); res=try(mod<-gam(as.formula(paste("as.factor(class)~", features)) , 
                                                                  family = "binomial", data = train.h),silent=TRUE)  
                                              if(inherits(res, "try-error")!=TRUE)
                                                     {
                               new1<-data.frame(test.h[,2:(J+1)])
                               r1<-predict(mod,newdata=new1,type="response")
                               predict<-ifelse(r1>0.50,1,0)
                               T<-table(predict,test.Y)
                               acc[i]=sum(diag(T))/sum(T)  
                                                    }else(next)
                               }else{
                                library(VGAM);res=try(mod<-vgam(as.formula(paste("as.factor(class) ~", 
                                                                           features)) ,family = multinomial, data = train.h),silent=TRUE)
                                               if(inherits(res, "try-error")!=TRUE)
                                                     {
                                 new1<-data.frame(test.h[,2:(J+1)])
                                r1<-predict(mod,newdata=new1,type="response")
                                predict<-apply(r1,1,which.max)-1
                                T<-table(predict,test.Y)
                                acc[i]=sum(diag(T))/sum(T)  
                                                     }else(next)
                                      } 
                       } 
        acc[is.na(acc)] <- (-100000) #penalty for bad h-values.
        my_list[[k]]=acc
      }
# Returns the best h value, which has maximum accuracy on an average.

  return(hvec[which.max(apply(matrix(unlist(my_list), ncol=n.folds, byrow=FALSE),1,mean))])  
 }

  h_opt<-bootstrap(data=train_data,hvec=h)
      }else{
  h_opt<-h_given
      }  

# New train and test data, considering LMD's as features:

  LMD.train<-lapply(g, function(j) LMD(h=h_opt,data=train.X,train=X[[j+1]],inv_cov=inv_sigma[[j+1]]))
  LMD.test<- lapply(g, function(j) LMD(h=h_opt,data=test.X,train=X[[j+1]],inv_cov=inv_sigma[[j+1]]))

  train<-data.frame(train.Y,LMD.train)
  colnames(train)<-c("class",paste("LMD",1:J,sep=""))
  test<-data.frame(class=test.Y,LMD.test)
  colnames(test)<-c("class",paste("LMD",1:J,sep=""))


# Spline fitting for GAM:

features<-paste0("s(LMD", 1:J, ")", collapse = " + ")

# Model fitting and to find the predicted classes for the test data:

if(J==2){
  library(gam); mod<-gam(as.formula(paste("class ~", features)) , 
                          family = "binomial", data = train)  
  new<-data.frame(test[,2:(J+1)])
  prob<-predict(mod,newdata=new,type="response")
  predict<-ifelse(prob>0.50,1,0)}else{

  library(VGAM);mod<-vgam(as.formula(paste("as.factor(class) ~", features)) ,
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
