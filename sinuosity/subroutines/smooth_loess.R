#------------------------------------------------------------------------------------------
# Script for applying temporal loess-smoothing to a NetCDF file
# [August 2017; copied from smooth.splines.R]
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# Loading R routines + packages
#------------------------------------------------------------------------------------------

# Script with base routines (netcdf, etc.)
source("/cnrm/amacs/USERS/cattiaux/APPLICATE/TOOLS/common/R/rbase.R")

#------------------------------------------------------------------------------------------
# Getting input files:
#   - in.sl.var.nc  contains the nlon x nlat x ntime variable (netcdf)
#   - in.sl.par.txt contains the span parameter (scalar)
#------------------------------------------------------------------------------------------
print("Getting input files..")
dat=myno("in.sl.var.nc")
ndim=length(dim(dat$var))
nlon=length(dat$lon)
nlat=length(dat$lat)
ntime=length(dat$time)
sp=scan(file="in.sl.par.txt")
sp=sp/ntime

#----------------------------------------------------------------------
# Computing loess with the subroutine myloessy() [defined below] to an array nlon x nlat x time
# WARNING: mind the dimension of input data and R dimension-permuting jokes with apply
#----------------------------------------------------------------------
print(paste("Computing loess with span =",sp))

myloessy=function(x,span) predict(loess(x~c(1:length(x)),span=span))

if (ndim==1) smooth.dat=array(myloessy(dat$var,span=sp),dim=c(nlon,nlat,ntime))
if (ndim==2) smooth.dat=array(t(apply(dat$var,1,myloessy,span=sp)),dim=c(nlon,nlat,ntime))
if (ndim==3) smooth.dat=aperm(apply(dat$var,1:2,myloessy,span=sp),c(2,3,1))

#------------------------------------------------------------------------------------------
# Writing output files:
#   - out.sl.nc, same as in.sl.var.nc but with smoothed values
#------------------------------------------------------------------------------------------
print("Writing ouput file..")
mync("out.sl.nc",smooth.dat,dat$varname,dat$lon,dat$lonu,dat$lat,dat$latu,dat$time,dat$timeu)

#------------------------------------------------------------------------------------------
# End of script
#------------------------------------------------------------------------------------------
print("Done.")
