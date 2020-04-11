createError <- function(Arr)
{
  size_ = length(Arr);
  Arr[1] = Arr[1] + 10;
  Arr[size_] = Arr[size_] - 10;
  return(Arr);
}

lineRegr <- function(wannaEr = 0, left = -1.8, right = 2, countOfPoints = 20)
{
  step = (right - left) / (countOfPoints - 1);

  e = rnorm(n = countOfPoints, mean = 0, sd = 1);

  betaMat = matrix(data = 0, nrow = 2, ncol = 2);
  colnames(betaMat) <- c("0", "1");
  rownames(betaMat) <- c("MNK", "MNM");

  myPoints = matrix(data = 0, nrow =  2, ncol = countOfPoints);
  rownames(myPoints) <- c("X", "Y");
  myPoints["X",] = seq(from = left, to = right, by = step)
  myPoints["Y",] = 2 + 2 * myPoints["X",] + e

  if (wannaEr == 1)
  {
    myPoints["Y", ] = createError(myPoints["Y",]);
  }
  
  #Метод наименьших квадратов
  matr_stat = matrix(data = 0, nrow = 2, ncol = 3);
  colnames(matr_stat) <- c("Mean", "Square Mean", "Median");
  rownames(matr_stat) <- c("X", "Y");
  matr_stat["X", "Mean"] = mean(myPoints["X",]);
  matr_stat["Y", "Mean"] = mean(myPoints["Y",]);
  matr_stat["X", "Median"] = median(myPoints["X",]);
  matr_stat["Y", "Median"] = median(myPoints["Y",]);
  matr_stat["X", "Square Mean"] = mean(myPoints["X",]^2);
  matr_stat["Y", "Square Mean"] = mean(myPoints["Y",]^2);

  betaMat["MNK", "1"] = (mean(myPoints["X",] * myPoints["Y",]) -
                       matr_stat["X", "Mean"] * matr_stat["Y", "Mean"]) /
               (matr_stat["X", "Square Mean"] - matr_stat["X", "Mean"]^2);

  betaMat["MNK", "0"] = matr_stat["Y", "Mean"] - 
            matr_stat["X", "Mean"] * betaMat["MNK", "1"];
  

  #Метод наименьших модулей
  mas_rQ = sign(myPoints["X",] - matr_stat["X", "Median"]) * 
              sign(myPoints["Y",] - matr_stat["Y", "Median"]);
  rQ = mean(mas_rQ);
  l = countOfPoints / 4;
  j = countOfPoints - l + 1;
  xl = myPoints["X", l];
  xj = myPoints["X", j];
  yl = myPoints["Y", l];
  yj = myPoints["Y", j];
  q_x = (xj - xl) / 2;
  q_y = (yj - yl) / 2;
  betaMat["MNM", "1"] = rQ * q_y / q_x;
  betaMat["MNM", "0"] = matr_stat["Y", "Median"] - 
              matr_stat["X", "Median"] * betaMat["MNM", "1"];

  if (wannaEr == 0)
  {
    plot(myPoints["X",], myPoints["Y",], pch = 17, main = "Выборка без возмущений");
  } else {
    plot(myPoints["X",], myPoints["Y",], pch = 17, main = "Выборка с возмущениями");
  }
  abline(a = 2, b = 2, col = "#650397"); #violete
  abline(a = betaMat["MNK", "0"], b = betaMat["MNK", "1"], col = "#F37F00"); #orange
  abline(a = betaMat["MNM", "0"], b = betaMat["MNM", "1"], col = "#00E0FF"); #blue
  legend("bottomright", c("Выборка", "Модель", "МНК","МНМ"), 
       col = c("black", "#650397", "#F37F00", "#00E0FF"),
       lty = c(-1, 1, 1, 1), pch = c(17, NA, NA, NA));
  
  return(betaMat);
}