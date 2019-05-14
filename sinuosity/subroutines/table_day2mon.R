#------------------------------------------------------------------------------------------
# Script for computing monthly mean on a ASCII table with dates in yyyymmdd format
# [August 2017]
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# Loading R routines + packages
#------------------------------------------------------------------------------------------

# Script with base routines (time/date, etc.)
source("/cnrm/amacs/USERS/cattiaux/APPLICATE/TOOLS/common/R/rbase.R")

#------------------------------------------------------------------------------------------
# Getting input file
#   - in.d2m.txt contains the daily table
# WARNING: columns must be named (header=T), dates must be under $date with yyyymmdd format
#------------------------------------------------------------------------------------------
print("Getting input files..")
dat=read.table("in.d2m.txt",header=T)
dts=yyyymmdd2mdy(dat$date)

#------------------------------------------------------------------------------------------
# Computing monthly means and writing output file:
#   - out.d2m.txt contains the monthly-mean table
#------------------------------------------------------------------------------------------
print("Computing monthly averages..")

out=c(); for (y in unique(dts$y)){ for (m in 1:12){
  iym=which(dts$y==y & dts$m==m)
  outym=rep(NA,ncol(dat))
  if (length(iym)>0) outym=apply(dat[iym,],2,mean,na.rm=T)
  outym[1]=10000*y+100*m+1
  out=rbind(out,outym)
}}

print("Writing output files..")
out=as.data.frame(out)
names(out)=names(dat)
write.table(out,"out.d2m.txt",quote=F,row.names=F,col.names=T)

#------------------------------------------------------------------------------------------
# End of script
#------------------------------------------------------------------------------------------
print("Done.")


