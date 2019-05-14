#------------------------------------------------------------------------------------------
# Script for applying temporal spline-smoothing to a NetCDF file
# [August 2017]
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# Loading R routines + packages
#------------------------------------------------------------------------------------------

# Script with base routines (netcdf, splines, etc.)
source("/cnrm/amacs/USERS/cattiaux/APPLICATE/TOOLS/common/R/rbase.R")

#------------------------------------------------------------------------------------------
# Getting input files:
#   - in.ss.var.nc  contains the nlon x nlat x ntime variable (netcdf)
#   - in.ss.par.txt contains the number of degrees of freedom (scalar)
#------------------------------------------------------------------------------------------
print("Getting input files..")
dat=myno("in.ss.var.nc")
ndim=length(dim(dat$var))
nlon=length(dat$lon)
nlat=length(dat$lat)
ntime=length(dat$time)
dof=scan(file="in.ss.par.txt")

#----------------------------------------------------------------------
# Computing splines with the subroutine mysmoothy() to an array nlon x nlat x time
# WARNING: mind the dimension of input data and R dimension-permuting jokes with apply
#----------------------------------------------------------------------
print(paste("Computing splines with",dof,"df"))
if (ndim==1) smooth.dat=array(mysmoothy(dat$var,df=dof),dim=c(nlon,nlat,ntime))
if (ndim==2) smooth.dat=array(t(apply(dat$var,1,mysmoothy,df=dof)),dim=c(nlon,nlat,ntime))
if (ndim==3) smooth.dat=aperm(apply(dat$var,1:2,mysmoothy,df=dof),c(2,3,1))

#------------------------------------------------------------------------------------------
# Writing output files:
#   - out.ss.nc, same as in.ss.var.nc but with smoothed values
#------------------------------------------------------------------------------------------
print("Writing ouput file..")
mync("out.ss.nc",smooth.dat,dat$varname,dat$lon,dat$lonu,dat$lat,dat$latu,dat$time,dat$timeu)

#------------------------------------------------------------------------------------------
# End of script
#------------------------------------------------------------------------------------------
print("Done.")
