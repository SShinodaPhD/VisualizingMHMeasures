source("gh01_functions.R")

#--------------------------------------------------#
#-artificial data
pmat <- matrix(c(
0.022, 0.022, 0.022, 0.022, 
0.111, 0.022, 0.022, 0.022, 
0.111, 0.111, 0.022, 0.225, 
0.111, 0.111, 0.022, 0.022
),4,4,byrow=T)
fn.visualization.Fig2(pmat)


#--------------------------------------------------#
#-real data
mat3a <- matrix(c(
29, 43, 25, 31, 4,
23, 159, 89, 38, 14,
11, 69, 184, 34, 10,
42, 147, 148, 184, 17,
42, 176, 377, 114, 298
),5,5,byrow=T)
fn.mCoMH.sammary(mat3a)
fn.mCoME.sammary(mat3a)
fn.mCoEMH.sammary(mat3a)
fn.visualization.Fig3(mat3a)


mat3b <- matrix(c(
50, 45,8,18, 8,
28, 174, 84, 154, 55,
11, 78, 110, 223, 96,
14, 150, 185, 714, 447,
3, 42, 72, 320, 411
),5,5,byrow=T)
fn.mCoMH.sammary(mat3b)
fn.mCoME.sammary(mat3b)
fn.mCoEMH.sammary(mat3b)
fn.visualization.Fig3(mat3b)

