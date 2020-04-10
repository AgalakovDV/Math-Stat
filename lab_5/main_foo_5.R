mu <- c(0,0)  

R <- c()
RS <- c()
RQ <- c()

#смесь распределений
mixDistr <- function(sizeArr = 20)
{
  Sigma1 <- matrix(c(1, 0.9, 0.9, 1), 2);
  Sigma2 <-  matrix(c(100, -90, 90, 100), 2);
  
  for (j in 1:1000)
  {
    myMatrix <- (0.9 * mvrnorm(sizeArr, mu = mu, Sigma = Sigma1)
            + 0.1 * mvrnorm(sizeArr, mu = mu, Sigma = Sigma2));
    r = cor(myMatrix[, 1], myMatrix[, 2], method = "pearson");
    rs = cor(myMatrix[, 1], myMatrix[, 2], method = "spearman");
    R = append(R, r);
    RS = append(RS, rs);
    
    moreMed_1 = myMatrix[myMatrix[,1] > median(myMatrix[,1]),];
    moreMed = moreMed_1[moreMed_1[,2] > median(myMatrix[,2]),];
    lessMed_1 = myMatrix[myMatrix[,1] < median(myMatrix[,1]),];
    lessMed = lessMed_1[lessMed_1[,2] < median(myMatrix[,2]),];
    
    numEl = length(moreMed[,1]) + length(lessMed[,1]);
    rq = 2 * numEl / sizeArr - 1;
    RQ = append(RQ, rq);
  }  
  
  vec <- c(9);
  
  vec[1] = sum(R) / 1000;
  vec[2] = sum(R^2) / 1000;
  vec[3] = vec[2] - vec[1]^2;
  vec[4] = sum(RS) / 1000;
  vec[5] = sum(RS^2) / 1000;
  vec[6] = vec[5] - vec[4]^2;
  vec[7] = sum(RQ) / 1000;
  vec[8] = sum(RQ^2) / 1000;
  vec[9] = vec[8] - vec[7]^2;
  
  return(vec)
}

#нормальное распределение
norm_distr <- function(sizeArr = 20, rho = 0)
{
  Sigma <- matrix(c(1, rho, rho, 1), 2);
  for (j in 1:1000)
  {
    myMatrix <- mvrnorm(sizeArr, mu = mu, Sigma = Sigma);
    r = cor(myMatrix[, 1], myMatrix[, 2], method = "pearson");
    rs = cor(myMatrix[, 1], myMatrix[, 2], method = "spearman");
    R = append(R, r);
    RS = append(RS, rs);
    
    moreMed_1 = myMatrix[myMatrix[,1] > median(myMatrix[,1]),];
    moreMed = moreMed_1[moreMed_1[,2] > median(myMatrix[,2]),];
    lessMed_1 = myMatrix[myMatrix[,1] < median(myMatrix[,1]),];
    lessMed = lessMed_1[lessMed_1[,2] < median(myMatrix[,2]),];
    
    numEl = length(moreMed[,1]) + length(lessMed[,1]);
    rq = 2 * numEl / sizeArr - 1;
    RQ = append(RQ, rq);
  }
  
  vec <- c(9);
  
  vec[1] = sum(R) / 1000;
  vec[2] = sum(R^2) / 1000;
  vec[3] = vec[2] - vec[1]^2;
  vec[4] = sum(RS) / 1000;
  vec[5] = sum(RS^2) / 1000;
  vec[6] = vec[5] - vec[4]^2;
  vec[7] = sum(RQ) / 1000;
  vec[8] = sum(RQ^2) / 1000;
  vec[9] = vec[8] - vec[7]^2;
  
  return(vec)
}


draw_ellipse <- function(sizeArr = 20, rho = 0){
mu <- c(0,0)  

Sigma <- matrix(c(1, rho, rho, 1), 2)

#set.seed(100)
selection <- mvrnorm(sizeArr, 
                     mu = mu, 
                     Sigma = Sigma ) 
require(ellipse)
confidence.ellipse <- ellipse(Sigma, 
                              centre = mu, 
                              level = 0.99, 
                              npoints = 100)
plot(confidence.ellipse, 
     type = "l", 
     xlim = c(-3, 3), 
     ylim = c(-3, 3))
par(new = TRUE)
plot(selection,
     axes = FALSE, 
     ann = FALSE, 
     col = "#F37F00", #orange
     pch = 1)
}