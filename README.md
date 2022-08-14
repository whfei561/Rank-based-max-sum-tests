# Rank-based max-sum tests

This code introduces the mainly simulation code involved in the paper “Rank-based max-sum tests for mutual independence of high-dimensional random vectors”.

The code mainly depends on the functions in the R package.
The file "MAIN.R" contains comments and program code for one simulation.

We will give examples of how to calculate the five correlation coefficients Spearman's rho, Kendall's tau, Hoeffding's D, Bergsma-Dassios-Yanagimoto's tau*, and Blum-Kiefer-Rosenblatt's R, respectively.
We also give notes on the use of functions.

## Spearman's rho

`library(pcaPP)`

#Spearman's rho statistic

#Two simple vectors

`x<-rnorm(10,0,1)`

`y<-rnorm(10,0,0.2)`

#result is a value

`cor(x,y,method="spearman")`

[1] 0.05454545

#A sample matrix composed of multiple sample vectors

`x<-rnorm(120,0,1)`

`x<-matrix(x,40,3)`

#result is a 3*3 matrix

`cor(x,method="spearman")`


           [,1]        [,2]        [,3]
[1,]  1.0000000 -0.13827392 -0.22851782
[2,] -0.1382739  1.00000000  0.01838649
[3,] -0.2285178  0.01838649  1.00000000

## Kendall's tau

`library(pcaPP)`

#Kendall's tau statistic

#Two simple vectors

`x<-rnorm(10,0,1)`

`y<-rnorm(10,0,0.2)`

#result is a value

`cor.fk(x,y)`

[1] 0.02222222
#A sample matrix composed of multiple sample vectors

`x<-rnorm(120,0,1)`

`x<-matrix(x,40,3)`

#result is a 3*3 matrix

`cor.fk(x)`


          [,1]       [,2]       [,3]
[1,] 1.0000000  0.1410256  0.1820513
[2,] 0.1410256  1.0000000 -0.1948718
[3,] 0.1820513 -0.1948718  1.0000000

## Hoeffding's D

`library(Hmisc)`

#Hoeffding's D statistic

#Two simple vectors

`x<-rnorm(10,0,1)`

`y<-rnorm(10,0,0.2)`

#result is a 2*2 matrix

`hoeffd(x,y)$D`

            x           y
x  1.00000000 -0.01984127
y -0.01984127  1.00000000

#A sample matrix composed of multiple sample vectors

`x<-rnorm(120,0,1)`

`x<-matrix(x,40,3)`

#result is a 3*3 matrix

`hoeffd(x)$D`


             [,1]         [,2]         [,3]
[1,]  1.000000000 -0.012394986 -0.001896633
[2,] -0.012394986  1.000000000 -0.001408797
[3,] -0.001896633 -0.001408797  1.000000000

##Note that  the statistic obtained by "hoeffd" function is 30 times Hoeffding's D statistic in Hoeffding W. (1948)

## Bergsma-Dassios-Yanagimoto's tau*

`library(TauStar)`

#Bergsma-Dassios-Yanagimoto's tau* statistic

#Two simple vectors

`x<-rnorm(10,0,1)`

`y<-rnorm(10,0,0.2)`

#result is a value

`tStar(x,y)`

[1] -0.04761905


##"tStar" Function can only be evaluated in pairs



##Note that  the statistic obtained by "tStar" function is 4/6 times Bergsma–Dassios–Yanagimoto‘s tau* statistic in Drton M. (2020)



## Blum-Kiefer-Rosenblatt's R

#Blum-Kiefer-Rosenblatt's R can be obtained by Hoeffding's D and Bergsma-Dassios-Yanagimoto's tau* statistics

`library(TauStar)`

`library(Hmisc)`

#Two simple vectors

`x<-rnorm(10,0,1)`

`y<-rnorm(10,0,0.2)`

#result is a value

`(30*tStar(x,y)/4-3*hoeffd(x,y)$D[1,2])/2`

[1] -0.1547619

##Note that  this statistic obtained by "tStar" and "hoeffd" function is the Blum–Kiefer–Rosenblatt's R statistic in Drton M. (2020)



## References

Hoeffding W. "A Non-Parametric Test of Independence." The Annals of Mathematical Statistics, 19(4) 546-557 December, 1948.

Drton M, Han F, Shi H. "High-dimensional consistent independence testing with maxima of rank correlations." The Annals of Statistics, 48(6) 3206-3227 December, 2020.
