#MO
MO <- function(n = 20, mo = 0, sigma = 1)
{
  if (n == 20)
  {
    t = 2.093;
  }else if (n == 100) {
    t = 1.984;
  }else {
    return(NA);
  }
  rad = sigma*t/sqrt(n - 1);
  vec = c(2);
  vec[1] = mo - rad;
  vec[2] = mo + rad;
  return(vec);
}

#D 
Dis <- function(n = 20, sigma = 1)
{
  t_vec = c(2);
  if (n == 20)
  {
    t_vec[1] = 32.852;
    t_vec[2] = 8.907;
  }else if (n == 100) {
    t_vec[1] = 128.422;
    t_vec[2] = 73.361;
  }else {
    return(NA);
  }
  k = sigma * sqrt(n);
  vec = c(2);
  vec[1] = k / sqrt(t_vec[1]);
  vec[2] = k / sqrt(t_vec[2]);
  return(vec);
}

#MO asimpt
MO_Asimptotic <- function(n = 20, mo = 0, sigma = 1)
{
  u = 1.96
  rad = sigma * u / sqrt(n);
  vec = c(2);
  vec[1] = mo - rad;
  vec[2] = mo + rad;
  return(vec);
}

#D asimpt
Dis_Asimptotic <- function(n = 20, mo = 0, sigma = 1, selection)
{
  u = 1.96;
  m_4 = sum((selection - mo)^4) / n;
  e = m_4 / sigma^4 - 3;
  k = 0.5 * u * sqrt((e + 2)/n);
  vec = c(2);
  vec[1] = sigma * (1 - k);
  vec[2] = sigma * (1 + k);
  return(vec);
}

find_interval <- function(n = 20)
{
  selection = rnorm(n, 0, 1);
  mo = mean(selection);
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