## noise term: epsilon

## Error by Gaussian
EGS <- function(x, v) rnorm(x, 0, sqrt(v))
EG1 <- function(x) EGV(x, 1)
EG2 <- function(x) EGV(x, 2)
EG3 <- function(x) EGV(x, 3)
EG4 <- function(x) EGV(x, 4)
EG5 <- function(x) EGV(x, 5)

## Error by T-student
EST <- function(x, v) rt(x, 2 + 2 / (v - 1))
ET2 <- function(x) ETV(x, 2)
ET3 <- function(x) ETV(x, 3)
ET4 <- function(x) ETV(x, 4)
ET5 <- function(x) ETV(x, 5)

