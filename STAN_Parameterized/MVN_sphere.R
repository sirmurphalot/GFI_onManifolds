library(rstan)
library(mvtnorm)

set.seed(13)
Iters = 50000
n = 4
d = 3
m1 = -sqrt(5/8)
m2 = 1/sqrt(8)
z_shift=1
m3 = 1/sqrt(16)+z_shift
sd1 = 1
sd2 = 1
sd3 = 1
# For the purposes of our experiment, our data must be the same as what we have on Matlab.
DataX = c(0.1713,-2.0252,0.8154,0.0977)
DataY = c(2.2951,0.1641,0.7312,1.5356)
DataZ = c(1.6317,1.8123,-0.4593,0.2301)

sigma=diag(3)
y=matrix(c(DataX,DataY,DataZ),ncol=3)
theta_start = 0.7854
phi_start = 0.7854

data_stan=list(d=d,n=n,y=y,sig=sigma,dataX=DataX,dataY=DataY,dataZ=DataZ)
init_fun <- function(...) list(theta = 0.7854, phi = 0.7854)
fiducial_stan = stan(file='MVN_sphere.stan',data = data_stan, init = init_fun, 
                     iter = Iters, cores = 4L) 
# FIDfit=sampling(fiducial_stan, data = data_stan, #init = init_stan, 
#                 iter = Iters, cores = 4L)

FIDfit=fiducial_stan
plot(FIDfit,pars="phi",show_density = TRUE)
plot(FIDfit,pars="theta",show_density = TRUE)

FIDdraws <- extract(FIDfit)#,permute=FALSE, inc_warmup=TRUE) 
temp_data=data.frame(Phi=FIDdraws$phi,Theta=FIDdraws$theta)
write.csv(temp_data,file="../sphereDraws.csv")


