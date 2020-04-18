#MO
MO <- function(n = 20, mo = 0, sigma = 1)
{
  if (n == 20)
  {
    t = 2.093;#хи квадрат распределение
  }else if (n == 100) {
    t = 1.984;#хи квадрат распределение
  }else {
    return(NA);
  }
  rad = sigma*t/sqrt(n - 1); #радиус, центр - mo.
  vec = c(2); #вектор, хранящий границы доверительного интервала
  vec[1] = mo - rad; #левая граница
  vec[2] = mo + rad; #правая граница
  return(vec);
}

#D 
Dis <- function(n = 20, sigma = 1)
{
  t_vec = c(2);
  if (n == 20)
  {
    t_vec[1] = 32.852;#хи квадрат распределение
    t_vec[2] = 8.907;#хи квадрат распределение
  }else if (n == 100) {
    t_vec[1] = 128.422;#хи квадрат распределение
    t_vec[2] = 73.361; #хи квадрат распределение
  }else {
    return(NA);
  }
  k = sigma * sqrt(n); #числитель дроби
  vec = c(2); #вектор, хранящий границы доверительного интервала
  vec[1] = k / sqrt(t_vec[1]); #левая граница
  vec[2] = k / sqrt(t_vec[2]); #правая граница
  return(vec);
}

#MO asimpt
MO_Asimptotic <- function(n = 20, mo = 0, sigma = 1)
{
  u = 1.96  #Квантиль нормального распределения
  rad = sigma * u / sqrt(n); # радиус, центр - mo
  vec = c(2); #вектор, хранящий границы доверительного интервала
  vec[1] = mo - rad; #левая граница
  vec[2] = mo + rad; #правая граница
  return(vec);
}

#D asimpt
Dis_Asimptotic <- function(n = 20, mo = 0, sigma = 1, selection)
{
  u = 1.96;  #Квантиль нормального распределения
  mu_4 = sum((selection - mo)^4) / n; #центральный момент 4-го порядка
  e = mu_4 / sigma^4 - 3; #четвёртый выборочный центральный момент
  vec = c(2); #вектор, хранящий границы доверительного интервала
  k = 0.5 * u * sqrt((e + 2)/n);
  vec[1] = sigma * (1 - k); #левая граница
  vec[2] = sigma * (1 + k); #правая граница
  return(vec);
}

find_interval <- function(n = 20)
{
  selection = rnorm(n, 0, 1);
  mo = mean(selection); #Математическое ожидание
  sigma = sqrt(sum((selection - mo)^2) / n);
  myData = matrix(data = 0, nrow = 2, ncol = 4);
  colnames(myData) = c("MO",
                       "D",
                       "MO asimpt",
                       "D  asimpt");
  rownames(myData) = c("left", "right");
  myData[,1] = MO(n = n, mo = mo, sigma = sigma);
  myData[,2] = Dis(n = n, sigma = sigma);
  myData[,3] = MO_Asimptotic(n = n, mo = mo, sigma = sigma);
  myData[,4] = Dis_Asimptotic(n = n, mo = mo, sigma = sigma, selection = selection);
  View(myData);
  return(myData);
}