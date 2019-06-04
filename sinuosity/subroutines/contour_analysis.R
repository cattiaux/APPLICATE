#------------------------------------------------------------------------------------------
# Script for contour analysis, incl. sinuosity index used in Cattiaux et al. (GRL 2016)
# [last edit oct 2018]
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# Loading R routines + packages
#------------------------------------------------------------------------------------------

# Script with base routines (netcdf, time/date, etc.)
source("../../rbase.R")

# Additional package with geographic tools (spherical distance/perimeter, etc.)
library(geosphere)

#------------------------------------------------------------------------------------------
# Getting input files:
#   - in.ca.var.nc contains the nlon x nlat x ntime variable (Z500, PV, etc.)
#   - in.ca.iso.nc contains the    1 x    1 x ntime value of the iso-contour (may be constant)
# WARNING: latitudes must be in ascending order!
#------------------------------------------------------------------------------------------
print("Getting input files..")
dat=myno("in.ca.var.nc")
iso=myno("in.ca.iso.nc")$var
dts=mdy2yyyymmdd(time2mdy(dat$time,dat$timeu,caltype=dat$caltype))
ntime=length(dts)

#------------------------------------------------------------------------------------------
# List of statistics to be computed, with short names and default values
# WARNING: modify the below subroutine if you modify this list!
#------------------------------------------------------------------------------------------
LIST.OF.STATS=matrix(c(
  "total contour length"   ,"cltot" ,  0,
  "first contour length"   ,"cl1"   ,  0,
  "max contour length"     ,"clmax" ,  0,
  "min contour length"     ,"clmin" ,  0,
  "number of contours"     ,"nc"    ,  0,
  "maximal latitude"       ,"latmax", NA,
  "minimal latitude"       ,"latmin", NA,
  "mean latitude"          ,"latavg", NA,
  "median latitude"        ,"latmed", NA,
  "amplitude of latitudes" ,"latamp", NA,
  "std dev of latitudes"   ,"latstd", NA
),ncol=3,byrow=TRUE)

nstats=nrow(LIST.OF.STATS)

#------------------------------------------------------------------------------------------
# Subroutine for computing stats + flag on a lon-lat field (no time dimension)
# Returns a list:
#   - $stats: vector of contour statistics (see above list)
#   - $flag:  matrix nlon x nlat of contour flag (0 or 1)
#------------------------------------------------------------------------------------------
my.contour.analysis=function(x,lon,lat,iso){
  # Isolate contours with contourLines
  c=contourLines(lon,lat,x,levels=iso)
  # Number of separated contours
  nc=length(c)
  # Compute stats if nc > 0 ; else default values
  stats=LIST.OF.STATS[,3]
  if (nc>0){
    # Loop on nc: increment xy coordinates + vector of contour lengths
    # Note: as perimeter() assumes contour is periodic, need to compute
    #       on [ x,rev(x) | y,rev(y) ] and then divide by 2.
    cxy=c(); clen=c(); for (i in 1:nc){
      cxyi=cbind(c[[i]]$x,c[[i]]$y)
      cxy=rbind(cxy,cxyi)
      cxyi=rbind(cxyi,apply(cxyi,2,rev))
      clen=c(clen,perimeter(cxyi)/2)
    }
    # Construct stats vector (must correspond to LIST.OF.STATS above!)
    stats=c(
      sum(clen),                 # cltot
      clen[1],                   # cl1
      max(clen),                 # clmax
      min(clen),                 # clmin
      nc,                        # nc
      max(cxy[,2]),              # latmax
      min(cxy[,2]),              # latmin
      mean(cxy[,2]),             # latavg
      median(cxy[,2]),           # latmed
      max(cxy[,2])-min(cxy[,2]), # latamp
      sd(cxy[,2])                # latstd    
    )
  }
  # Construct flag matrix
  flag=matrix(0,length(lon),length(lat))
  if (nc>0){ for (i in 1:nrow(cxy)) flag[which.min(abs(lon-cxy[i,1])),which.min(abs(lat-cxy[i,2]))]=1 }
  # Output
  return(list(stats=stats,flag=flag))
}

#------------------------------------------------------------------------------------------
# Computing contour analysis and writing output files:
#   - out.ca.stats.txt is a ASCII file containing contour stats for each time step
#   - out.ca.flag.nc   is a nlon x nlat x ntime NetCDF file containing contour flag
#------------------------------------------------------------------------------------------

print("Starting contour analysis..")

# Initialization
STATS=matrix(NA,ntime,nstats)
FLAG=array(NA,dim=dim(dat$var))

# Loop on time steps
for (i in 1:ntime){
  if (i/50==trunc(i/50)) print(paste(".. time step",i,"/",ntime,".."))
  dum=my.contour.analysis(dat$var[,,i],dat$lon,dat$lat,iso[i])
  STATS[i,]=dum$stats
  FLAG[,,i]=dum$flag
}

# Writing ouput..
print("Writing output files..")

# .. stats file
STATS=cbind(dts,iso,STATS)
STATS=as.data.frame(STATS)
names(STATS)=c("date","iso",LIST.OF.STATS[,2])
write.table(STATS,"out.ca.stats.txt",quote=F,row.names=F,col.names=T)

# .. flag file
mync("out.ca.flag.nc",FLAG,"flag",dat$lon,dat$lonu,dat$lat,dat$latu,dat$time,dat$timeu,dat$caltype)

#------------------------------------------------------------------------------------------
# End of script
#------------------------------------------------------------------------------------------
print("Done.")


