Probability <- function(){
  n = 100;
  arrX = rnorm(n, 0, 1);
  row = ceiling(1 + 3.3 * log10(n))

  
  myData <- matrix(data = 0, nrow = row + 1, ncol = 7);
  colnames(myData) = c("Левая граница",
                       "Правая граница",
                       "Вероятность попадания p_i",
                       "Число попаданий n_i",
                       "произведение вероятности
                        и числа всех точек 'n*p_i'",
                       "'n_i - n*p_i'",
                       "'Z'");
  
  rownames(myData) = c("1", "2", "3", "4", "5", "6", "7", "8", "SUM");


  left = min(arrX)
  right = max(arrX)
  arrPoints = seq(from = left, to = right, by = (right - left)/row)
  arrPoints[1] = -Inf
  arrPoints[row + 1] = Inf
  
  for (i in 1:row)
  {
    myData[i, 1] = arrPoints[i];
    myData[i, 2] = arrPoints[i + 1];
    myData[i, 3] = pnorm(arrPoints[i + 1]) - pnorm(arrPoints[i]);
    myData[i, 4] = length(arrX[arrX <= arrPoints[i + 1] & arrX >= arrPoints[i]]);
    myData[i, 5] = n * myData[i, 3];
    myData[i, 6] = myData[i, 4] - myData[i, 5];
    myData[i, 7] = (myData[i, 6])^2 / myData[i, 5];
  }
  
  for (j in 3:7)
  {
    myData[row + 1, j] = sum(myData[, j]);
  }
  View(myData);
  return(myData);
}