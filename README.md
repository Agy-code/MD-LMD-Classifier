# MD-LMD-Classifier
The following libraries are needed:    
 robustbase: for MCD,
    
 RobStatTM: for MVE

 gam,VGAM : For gam fitting

1. The MD Classifier gives error rate for a given test dataset. It needs the following inputs:


    train.X<- features of the training data.
   
    train.Y<- class information of the training data.

    test.X<-  features of the test data.

    test.Y<- class information of the training data.

    sigma_type <- This is the type of scatter matrix, that to be used. The possible choices are:
            'varcov' = usual variance matrix,
            'MCD',
             'MVE',
            'diagonal'= diagonal matrix containing the variances only.
         or 'identity'= Identity matrix.

   
2. The LMD Classifier gives error gives the error rate for a given test dataset. For this we need to calculate the LMD first.

(i) The LMD function calculates the Local Mahalanobis distances.It needs the following inputs:

h<-tuning parameter (vector or real valued).

data<- a matrix or data frame such that the rows contain query points, for which LMD has to be found.

train<-a matrix or data frame containing {x1,x2,...xn}.

inv_cov<-the INVERSE of the sigma matrix.

kern<-This is the kernel to be used, if not specified the Gaussian kernel will be considered.



(ii) The LMD classifier needs the following inputs:

  train.X<- features of the training data.

  train.Y<- class information of the training data.
  
  test.X<-  features of the test data.
  
  test.Y<- class information of the training data.
  
  sigma_type <- This is the type of variance matrix, that to be used.The possible choices are:
            'varcov' = usual variance matrix,
            'MCD',
             'MVE',
            'diagonal'= diagonal matrix containing the variances only.
         or 'identity'= Identity matrix.
 
  If boot is TRUE, then it performs bootstrap samling to find optimal value of h, based on n.folds bootstrap samples.
  
  If boot is FALSE, then it considers the given value of h (as specified by h_given) as the optimal.
 


 
