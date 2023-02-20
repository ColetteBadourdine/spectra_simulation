##########################################
#covariance stability in hyprspectral data
#########################################
#load a test data
library(terra)
fl=rast("y:/users/ColinPrieur/Paracou/VNIR/correctedT1/PAR_FL21_20160919_162638_VNIR_1600_SN0014_PS01_IMG003_atm_slice_001_x2_post_COREG_LOCAL.tif")
plot(fl[[55]])
#crop a sample area
new_ext<-ext(c(-5890900, -5890850,587550,587600))
sel=crop(fl, new_ext)
#exclude added layers (>160)
sel <- sel[[1:160]]
plot(sel[[78]])
##
library(corrplot)
corrplot(cor(sel[]))
mean(sel[], na.rm=T)
rmean <- app(sel, mean)
plot(rmean)
dim(sel)
spec=c()
for (i in 1:dim(sel)[3])
{
  spec<-rbind(spec,mean(sel[[i]][]))
}
plot(spec)

# Cholesky decomposition
L = chol(cov(sel[]))
nvars = dim(L)[1]
nobs=1000

r = t(L) %*% matrix(rnorm(nvars*nobs), nrow=nvars, ncol=nobs)
rowMeans(r)
corrplot(cor(t(r)))



################################
L = chol(cov(test_qualea[, -c(1:3)]))
nvars = dim(L)[1]
nobs=100

r = t(L) %*% matrix(rnorm(nvars*nobs), nrow=nvars, ncol=nobs)
#test modifier var intra
r_modif = t(L) %*% matrix(rnorm(nvars*nobs, sd = 0.5), nrow=nvars, ncol=nobs)

r2 = r+matrix(colMeans(test_qualea[, -c(1:3)]),nrow = 123,ncol=100)

matplot(r2,type = 'l')

r3 = r_modif+matrix(colMeans(test_qualea[, -c(1:3)]),nrow = 125,ncol=100)

matplot(r3,type = 'l')

matplot(matrix(t(test_qualea[, -c(1:3)]),ncol = 5102),type = 'l')

dim(test_qualea[, -c(1:3)])

dim(matrix(t(test_qualea[, -c(1:3)])))


test_qualea[, -c(1:3)]

?matplot
r = matrix(data=nobs)
rowMeans(r)
corrplot(cor(t(r)))

aa=unique(ReflectanceData[, c(1,3)])
table(aa$SPID)
#Sterculia_speciosa

L = chol(cov(test_sterculia[, -c(1:3)]))
nvars = dim(L)[1]
nobs=10

r = t(L) %*% matrix(rnorm(nvars*nobs), nrow=nvars, ncol=nobs)

r2 = r+matrix(colMeans(test_sterculia[, -c(1:3)]),nrow = 125,ncol=10)

matplot(r2,type = 'l')
matplot(matrix(t(test_sterculia[, -c(1:3)]),ncol = 5102),type = 'l')

dim(test_sterculia[, -c(1:3)])

dim(matrix(t(test_sterculia[, -c(1:3)])))
