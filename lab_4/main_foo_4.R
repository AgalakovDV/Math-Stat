# распределение Гаусса
gauss_distr <- function(size) {
  randX <- rnorm(n = size,
                 mean = 0, 
                 sd = 1);
  
  return(sort(randX))
}

# распределение Коши
coshi_distr <- function(size){
  randX <- rcauchy(n = size,
                   location = 0,
                   scale = 1);
  
  return(sort(randX))
}

# распределение Пуассона
pois_distr <- function(size){
  randX <- rpois(n = size,
                 lambda = 10);
  
  return(sort(randX))
}

# распределение Лапласа
laplas_distr <- function(size){
  randX <- rlaplace(n = size, 0, 1/sqrt(2));
  
  return(sort(randX))
}

# равномерное распределение
uniform_distr <- function(size){
  randX <- runif(n = size, -sqrt(3), sqrt(3)) 
  
  return(sort(randX))
}


# 1 - Гаусс
# 2 - Коши
# 3 - Лаплас
# 5 - Пуассон
# 4 - Равномерное
make_distr <- function(size = 20, param = 1){
  column = 32*round(1 + 3.3 * log10(size));
  arrX <- array(dim = column + 1);
  custom <- array(data = 0, dim = column);
  randX <- array(dim = size);
  
  if (param > 0 & param < 5) {
    left = -4;
    right = 4;
    step = (right - left) / column;
    
    if (param == 1) {
      randX = gauss_distr(size);
    } else if (param == 2) {
      randX = coshi_distr(size);   
    } else if (param == 3) {
      randX = laplas_distr(size);
    } else if (param == 4) {
      randX = uniform_distr(size)
    }
  } else {
    left = 6;
    right = 14;
    step = (right - left) / column;
    
    randX = pois_distr(size);
  }
  
  for (i in 1:(column + 1)) {
    arrX[i] = left + step * (i - 1);
  }
  
  i = 1;
  j = 1;
  while (j <= size & i <= column) {  
    if (arrX[i] <= randX[j] & randX[j] <= arrX[i + 1]) {
      custom[i] = custom[i] + 1;
      j = j + 1;
    } else {
      if (arrX[i] > randX[j]) 
      {
        j = j + 1;
      } else {
        i = i + 1;
      }
    }
  }
  
  p <- array(data = 0, dim = column + 1);
  for (i in 2:column + 1) {
    p[i] = p[i - 1] + custom[i - 1];
  }
  
  P <- array(dim = column + 1);
  for (i in 1:column + 1) {
    P[i] = p[i] / size;
  }

  unifX = sort(runif(n = 1024, left, right))
  
  if (param == 1) {
    myLine = pnorm(unifX,
                   mean = 0,
                   sd = 1);
    }else if (param == 2) {
    myLine = pcauchy(unifX,
                     location = 0,
                     scale = 1);
  } else if (param == 3) {
    myLine = plaplace(unifX, 0, 1/sqrt(2));
    
  } else if (param == 4) {
    myLine = punif(unifX, -sqrt(3), sqrt(3));
  } else {
    myLine = ppois(unifX, lambda = 10);
  }
  
  plot(arrX, P, type = "s", 
       col = "#FF0095", #pink
       main = "Э.Ф.Р.",
       ylim = c(0,1));
  lines(unifX, myLine, 
        col = "#00E0FF"); #blue
}

nuclear_assessment <- function(size = 20, param = 1) {
  if (param > 0 & param < 5) {
    left = -4;
    right = 4;
    unifX = sort(runif(n = 1024, left, right));
    if (param == 1) {
      randX = rnorm(n = size,
                    mean = 0, 
                    sd = 1);
      myLine = dnorm(unifX);
    } else if (param == 2) {
      randX = rcauchy(n = size,
                       location = 0,
                       scale = 1);
      myLine = dcauchy(unifX);
    } else if (param == 3) {
      randX <- rlaplace(n = size, 0, 1/sqrt(2));
      myLine = dlaplace(unifX, 0, 1/sqrt(2));
    } else if (param == 4) {
      randX <- runif(n = size, -sqrt(3), sqrt(3));
      myLine = dunif(unifX, -sqrt(3), sqrt(3));
    }
  } else {
    left = 6;
    right = 14;
    unifX = sort(round(runif(n = 1024, left, right)));
    randX <- rpois(n = size,
                   lambda = 10);
    myLine = dpois(unifX, 10);
  }
  
  
  plot(density(randX, adjust = 1/2), 
       xlim = c(left, right), 
       ylim = c(0, 0.5),
       col = "#650397", #violete
       main = "Ядерная оценка, 0.5h");
  lines(unifX, myLine, 
        col = "#F37F00") #orange
  
  plot(density(randX, adjust = 1),
       xlim = c(left, right),
       ylim = c(0, 0.5),
       main = "Ядерная оценка, h",
       col = "#650397"); #violete
  lines(unifX, myLine, 
        col = "#F37F00") #orange
  
  plot(density(randX, adjust = 2),
       xlim = c(left, right),
       ylim = c(0, 0.5),
       main = "Ядерная оценка, 2h",
       col = "#650397"); #violete
  lines(unifX, myLine, 
        col = "#F37F00") #orange
}