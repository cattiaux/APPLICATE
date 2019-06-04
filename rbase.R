###########################################
# Basic functions used for other scripts  #
# + loading of usual packages             #
###########################################

#-------------------------------------------------------------------------------------
# Packages
#-------------------------------------------------------------------------------------
#library(ncdf,quietly=TRUE,verbose=FALSE)
library(ncdf4,quietly=TRUE,verbose=FALSE)
#library(clim.pact,quietly=TRUE,verbose=FALSE)
library(fields,quietly=TRUE,verbose=FALSE)
library(abind,quietly=TRUE,verbose=FALSE)
library(maps,quietly=TRUE,verbose=FALSE)
library(mapproj,quietly=TRUE,verbose=FALSE)
library(sn,quietly=TRUE,verbose=FALSE)
library(class,quietly=TRUE,verbose=FALSE)
#library(mclust,quietly=TRUE,verbose=FALSE)
#library(grImport,quietly=TRUE,verbose=FALSE)
library(chron,quietly=TRUE,verbose=FALSE)


#-------------------------------------------------------------------------------------
# Aliases on existing functions
#-------------------------------------------------------------------------------------
big=toupper#.case
small=tolower#.case


#-------------------------------------------------------------------------------------
# Vectors/time-series related functions (x)
#-------------------------------------------------------------------------------------
# Basic functions to get/remove first/last/penul(n-1) element of a vector
first=function(x)   x[1]
nofirst=function(x) x[-1]
last=function(x)    x[length(x)]
nolast=function(x)  x[-length(x)]
penul=function(x)   x[length(x)-1]

# Difference between consecutive elements of a vector (eg t+1/t for time series)
autodiff=function(x) {
  out=numeric(0); if (length(x)>1) out=x[2:length(x)]-x[1:(length(x)-1)]
  out
}

# Function replacing vector values persisting less than n times by NAs (default)
mypersist=function(x,n,replace=NA){
  chg=c(0,which(is.na(autodiff(x)) | autodiff(x)!=0),length(x))
  na=which(autodiff(chg) < n)
  for (i in 1:length(na)) x[(chg[na[i]]+1):chg[na[i]+1]]=replace
  x
}

# Normalization of a time series
mynorm=function(x,mean=NULL,sd=NULL){
  if (is.null(mean)) mean=mean(x,na.rm=TRUE); if (is.null(sd)) sd=sd(x,na.rm=TRUE)
  (x-mean)/sd
}

# Biased formula for variance
mybiasedvar=function(x,na.rm=TRUE){
  if (na.rm) out=mybiasedvar(na.omit(x),na.rm=FALSE) else out=mean((x-mean(x))^2)
  out
}

# Computing percentages
mypercent=function(x) x*100/sum(x)

# Slope of the linear regression of a vector
mytrend=function(x){
  out=NA; if (any(!is.na(x))) out=lm(x~c(1:length(x)))$coef[2]
  out
}

# Linear fit of a vector
mylmfit=function(x,ref=NULL,predict.na=FALSE){
  if (is.null(ref)) ref=1:length(x)
  if (any(!is.na(x))){ fit=lm(x~ref,na.action=na.exclude)
    if (predict.na)  out=as.vector(predict(fit,as.data.frame(ref)))
    if (!predict.na) out=as.vector(predict(fit))
  }
  out
}

# Detrending a time series (wrt ref, only if significant wrt pval)
mydetrend=function(x,ref=NULL,pval=0.05){
  if (is.null(ref)) ref=1:length(x)
  if (mypvalue(x,ref)<pval) out=x-mylmfit(x~ref)
  else out=x-mean(x)
  out
}

# Splines accounting for NAs + possibility of deducing df from length(x) (sdf=length(x)/df)
# (Use smooth.na=F for not predicting at NA values of x)
mysmoothy=function(x,df=NULL,sdf=30,smooth.na=TRUE,periodic=FALSE){
  if (is.null(df) & is.null(sdf)) sdf=length(x)
  if (is.null(df)) df=length(x)/sdf
  if (length(unique(na.omit(x))) < max(4,df)) out=x#rep(NA,length(x))
  else {
    if (df>1){
      if (periodic) { x=rep(x,3) ; df=df*3 }
      sms=smooth.spline(which(!is.na(x)),na.omit(x),df=df)
      out=predict(sms$fit,1:length(x))$y
      if (!smooth.na) out[which(is.na(x))]=NA
      if (periodic) { out=matrix(out,ncol=3)[,2] ; df=df/3}
    }
    if (df<=1) out=mylmfit(x,predict.na=smooth.na)
  }
  out
}

# Running averages
myrunmean=function(x,n,na.rm=T,fill.with.NAs=F,align="center"){
  out=c(); for (i in 1:(length(x)-n+1)) out=c(out,mean(x[i:(i+n-1)],na.rm=na.rm))
  if (fill.with.NAs){
    if (align=="center") out=c(rep(NA,floor((n-1)/2)),out,rep(NA,ceiling((n-1)/2)))
    if (align=="left")   out=c(out,rep(NA,n-1))
    if (align=="right")  out=c(rep(NA,n-1),out)
  }
  out
}

# Average by blocks (useful for seasonal averages, for instance)
myblockmean=function(x,n,na.rm=TRUE){
  if (length(x)/n != trunc(length(x)/n)) print("Error: length(x) must be a multiple of n.")
  else return(apply(matrix(x,nrow=n),2,mean,na.rm=na.rm))
}

myblockmax=function(x,n,na.rm=TRUE){
  if (length(x)/n != trunc(length(x)/n)) print("Error: length(x) must be a multiple of n.")
  else return(apply(matrix(x,nrow=n),2,max,na.rm=na.rm))
}

myblockmin=function(x,n,na.rm=TRUE){
  if (length(x)/n != trunc(length(x)/n)) print("Error: length(x) must be a multiple of n.")
  else return(apply(matrix(x,nrow=n),2,min,na.rm=na.rm))
}

#-------------------------------------------------------------------------------------
# Functions for comparing two vectors/time-series (x,y)
#-------------------------------------------------------------------------------------
# Euclidean distance between two vectors
myeuclid=function(x,y,na.rm=TRUE) sqrt(sum((x-y)^2,na.rm=na.rm))

# RMS difference between two vectors
myrmsd=function(x,y,na.rm=TRUE) sqrt(mean((x-y)^2,na.rm=na.rm))

# Mahalanobis distance between two populations (samples in lines)
mymaha=function(x,y,sym=TRUE){
  mux=apply(x,2,mean,na.rm=TRUE)
  muy=apply(y,2,mean,na.rm=TRUE)
  if (!sym)  out=(t(mux-muy)) %*% (solve(cov(x))) %*% (mux-muy)
  if (sym) out=0.5*(mymaha(x,y,sym=FALSE)+mymaha(y,x,sym=FALSE))
  out
}

# Correlation function accounting for NAs
mycor=function(x,y) {
  ix=which(is.na(x)); iy=which(is.na(y)); ixy=unique(c(ix,iy))
  if (length(ixy)>0) {x=x[-ixy];y=y[-ixy]}
  cor(x,y)
}

# P.value of a linear trend (if y=NULL) or a correlation test with y
mypvalue=function(x,y=NULL){
  if (is.null(y)) y=1:length(x)
  ix=which(is.na(x)); iy=which(is.na(y)); ixy=unique(c(ix,iy))
  if (length(ixy)>0) {x=x[-ixy];y=y[-ixy]}
  if (length(x)>2) return(cor.test(x,y)$p.value) else return(NA)
}


#-------------------------------------------------------------------------------------
# Weighted functions (eg for cos(latitude) weighting)
# Edit: IC 95% with bootstrap
#-------------------------------------------------------------------------------------
# Weighted mean sd var and skewness
mywmean=function(x,w,na.rm=TRUE,bootstrap=FALSE){
  nas=which(is.na(x)); if (length(nas)>0 & na.rm) {x=x[-nas]; w=w[-nas]}
  if (length(x)==0) return(NA)
  if (length(w)==1) w=rep(w,length(x))
  w=w/sum(w); out=sum(x*w)
  if (bootstrap) {tmp=c(); for (i in 1:1000){
    ii=sample(1:length(x),length(x),replace=TRUE); ww=w[ii]; xx=x[ii]
    ww=ww/sum(ww); tmp=c(tmp,sum(xx*ww))
  }; out=c(out,quantile(tmp,0.025),quantile(tmp,0.975))}
  out
}
mywsd=function(x,w,na.rm=TRUE,unbiased=TRUE,bootstrap=FALSE){
  nas=which(is.na(x)); if (length(nas)>0 & na.rm) {x=x[-nas]; w=w[-nas]}
  if (length(x)==0) return(NA)
  n=length(x); unbias=1; if (unbiased) unbias=n/(n-1)
  mx=mywmean(x,w); out=sqrt(unbias*mywmean((x-mx)^2,w))
  if (bootstrap) {tmp=c(); for (i in 1:1000){
    ii=sample(1:length(x),length(x),replace=TRUE); ww=w[ii]; xx=x[ii]
    mxx=mywmean(xx,ww); tmp=c(tmp,sqrt(unbias*mywmean((xx-mxx)^2,ww)))
  }; out=c(out,quantile(tmp,0.025),quantile(tmp,0.975))}
  out
}
mywvar=function(x,w,na.rm=TRUE,unbiased=TRUE,bootstrap=FALSE){
  nas=which(is.na(x)); if (length(nas)>0 & na.rm) {x=x[-nas]; w=w[-nas]}
  if (length(x)==0) return(NA)
  n=length(x); unbias=1; if (unbiased) unbias=n/(n-1)
  mx=mywmean(x,w); out=unbias*mywmean((x-mx)^2,w)
  if (bootstrap) {tmp=c(); for (i in 1:1000){
    ii=sample(1:length(x),length(x),replace=TRUE); ww=w[ii]; xx=x[ii]
    mxx=mywmean(xx,ww); tmp=c(tmp,unbias*mywmean((xx-mxx)^2,ww))
  }; out=c(out,quantile(tmp,0.025),quantile(tmp,0.975))}
  out
}
mywskew=function(x,w,na.rm=TRUE,unbiased=TRUE,bootstrap=FALSE){
  nas=which(is.na(x)); if (length(nas)>0 & na.rm) {x=x[-nas]; w=w[-nas]}
  if (length(x)==0) return(NA)
  n=length(x); unbias=1; if (unbiased) unbias=n^2/((n-1)*(n-2))
  mx=mywmean(x,w); sx=mywsd(x,w,unbiased=unbiased); out=unbias*mywmean(((x-mx)/sx)^3,w)
  if (bootstrap) {tmp=c(); for (i in 1:1000){
    ii=sample(1:length(x),length(x),replace=TRUE); ww=w[ii]; xx=x[ii]
    mxx=mywmean(xx,ww); sxx=mywsd(xx,ww,unbiased=unbiased);
    tmp=c(tmp,unbias*mywmean(((xx-mxx)/sxx)^3,ww))
  }; out=c(out,quantile(tmp,0.025),quantile(tmp,0.975))}
  out
}

# Weighted quantiles
mywqs=function(x,w,p=myquantiles,na.rm=TRUE){
  nas=which(is.na(x)); if (length(nas)>0 & na.rm) {x=x[-nas]; w=w[-nas]}
  if (length(x)==0) return(rep(NA,length(p)))
  x=x[order(x)]; w=w[order(x)]; out=c()
  for (pp in p) out=c(out,x[which.max(cumsum(w)/sum(w)>=pp)])
  out
}

# Weighted Euclidean distance
myweuclid=function(x,y,w,na.rm=TRUE) sqrt(sum(w*(x-y)^2,na.rm=na.rm))

# Weighted RMS difference
mywrmsd=function(x,y,w,na.rm=TRUE) sqrt(mywmean((x-y)^2,w,na.rm=na.rm))

# Weighted correlation
mywcor=function(x,y,w,na.rm=T){
  nas=which(is.na(x) | is.na(y))
  if (length(nas)>0 & na.rm) {x=x[-nas]; y=y[-nas]; w=w[-nas]}
  w=w/sum(w) ; x=x-sum(w*x) ; y=y-sum(w*y)
  sum(w*x*y)/sqrt(sum(w*x*x)*sum(w*y*y))
}

# Weighted coordinates for Taylor diagrams (x0 = reference)
mywtaylorcoords=function(x0,x,w) {
  sig0=mywsd(x0,w); sig=mywsd(x,w); r=mywcor(x0,x,w)
  sig/sig0*c(r,sqrt(1-r^2))
}

# Weighted nearest neighbours of dat among cls - (eg for weather regimes attribution step)
# Input: dat matrix ntime*nspace - cls matrix ncluster*nspace - weigths vector nspace
mywnn=function(dat,cls,weights,method=c("rmsd","cor","euclid")){
  meth=method[1]; nt=nrow(dat); ns=ncol(dat); nc=nrow(cls)
  if (ncol(cls)!=ns | length(weights)!=ns) return("Error in input dimensions")
  if (meth=="rmsd")   { fun=function(x,y) apply(x,1,mywrmsd,y,weights)   ; best=min }
  if (meth=="cor")    { fun=function(x,y) apply(x,1,mywcor,y,weights)    ; best=max }
  if (meth=="euclid") { fun=function(x,y) apply(x,1,myweuclid,y,weights) ; best=min }
  out=c(); for (i in 1:nt){
    vals=fun(cls,dat[i,]); icls=first(which(vals==best(vals)))
    out=c(out,icls,vals)
  }
  out=as.data.frame(matrix(out,nrow=nt,byrow=TRUE)); names(out)=c("class",paste(meth,1:nc,sep=""))
  out
}

#-------------------------------------------------------------------------------------
# Bootstrap functions
#-------------------------------------------------------------------------------------
# Mean of a time series with associated IC
mybootstrap=function(x,n=1000,p=0.05,nmin=length(x),na.rm=TRUE,replace=TRUE){
  if (nmin==length(x)) size=rep(length(x),n) else size=sample(nmin:length(x),n,replace=TRUE)
  boot=c(); for (i in 1:n){
    xi=x[sample(1:length(x),size[i],replace=replace)]
    boot=c(boot,mean(xi,na.rm=na.rm))
  }
  c(mean(boot),quantile(boot,p/2),quantile(boot,1-p/2))
}

# Mean of a difference between 2 time series with associated IC
mybootstrapdiff=function(x,y,n=1000,p=0.05,nmin=length(x),na.rm=TRUE,replace=TRUE){
  if (nmin==length(x)) size=rep(length(x),n) else size=sample(nmin:length(x),n,replace=TRUE)
  boot=c();  for (i in 1:n){
    xi=x[sample(1:length(x),size[i],replace=replace)]
    yi=y[sample(1:length(y),size[i],replace=replace)]
    boot=c(boot,mean(xi-yi,na.rm=na.rm))
  }
  c(mean(boot),quantile(boot,p/2),quantile(boot,1-p/2))
}


#-------------------------------------------------------------------------------------
# Array-related functions
#-------------------------------------------------------------------------------------
# Repetition of an array (x,y,t) along t, until length.out (if 2D-array, force 3D)
rep.abind=function(tab,length.out){
  if (length(dim(tab))==2){
    tmp=array(dim=c(dim(tab),1)); tmp[,,1]=tab; tab=tmp; rm(tmp)}
  d=dim(tab)[3]
  while (d<length.out){
    tab=abind(tab,tab,along=3)
    d=dim(tab)[3]
  }
  tab[,,1:length.out]
}


#-------------------------------------------------------------------------------------
# Time/Date-related functions - Dates are data.frames $m $d $y or yyyymmdd numbers
# Warning to CDO or NCO conventions for dates of monthly means (CDO: last day - NCO: median day)
#-------------------------------------------------------------------------------------
# Is year y leap (bissextile)?
bissex=function(y,caltype="leap"){
  if (caltype %in% c("noleap","360_day")) biss=FALSE
  else  biss=(trunc(y/400)*400==y | (trunc(y/4)*4==y & trunc(y/100)*100!=y))
  biss
}

# Nb of days of month (m,y)
ndays=function(m,y,caltype="leap"){
  if (caltype=="360_day") n=rep(30,12)[m]
  else n=c(31,28+bissex(y,caltype),31,30,31,30,31,31,30,31,30,31)[m]
  n
}

# Nb of days of season (s,y)
ndays.season=function(s,y,caltype="leap"){
  n=0; for (m in season2months(s)) n=n+ndays(m,y,caltype)
  n
}

# My modulo (mymod(n,n)=n) - default: modulo 12 (for months)
mymod=function(x,n=12){
  mm=c(); for (i in 1:length(x)){
    if (x[i]>=1) mm=c(mm,x[i]-n*trunc((x[i]-1)/n))
    if (x[i]<1)  mm=c(mm,x[i]+n*(trunc(-x[i]/n)+1)) }
  mm
}

# Season (eg "JJA") to index of first month ($i=6) and nb of months ($n=3), to months (6:8)
# or to cdo-months ("6,7,8")
season2i=function(season){
  ii=c(); nn=c();
  mm=c("JF","FM","MA","AM","MJ","JJ","JA","AS","SO","ON","ND","DJ")
  for (s in season){
    split=strsplit(s,"")[[1]]  
    ii=c(ii,which(mm==paste(split[1:2],collapse="")))
    nn=c(nn,length(split))
  }
  data.frame(i=ii,n=nn)
}
season2months=function(season)    { ii=season2i(season); mymod(seq(ii$i,by=1,le=ii$n)) }
season2cdomonths=function(season) { mm=season2months(season); paste(mm,collapse=",")}

# Function similar to julday() but dealing with all caltypes (leap/standard,noleap,360_day)
# (Odd reference to 1/1/-4713 in julday...)
# EDIT 2014: switch from julday to julian BUT reference left to -4713 for no-leap cases
myjulday=function(m,d,y,caltype="leap"){
  refy=-4713
  if (caltype %in% c("leap","standard")) out=julian(m,d,y)#julday(m,d,y)
  if (caltype %in% c("noleap","360_day")){
    Ndyr=sum(ndays(1:12,1,caltype))
    dts=matrix(c(m,d,y),ncol=3)
    out=apply(dts,1,function(x) Ndyr*(x[3]-refy-(x[3]>0))+(x[1]!=1)*sum(ndays(1:(x[1]-1),x[3],caltype))+x[2])
  }
  out
}

# Function similar to caldat() but dealing with all caltypes (inverse function of myjulday)
# EDIT 2014: switch from caldat to month.day.year BUT reference left to -4713 for no-leap cases
mycaldat=function(jday,caltype="leap"){
  refy=-4713
  if (caltype %in% c("leap","standard")) {out=month.day.year(jday) #out=caldat(jday)
                        names(out)=c("m","d","y")}
  if (caltype %in% c("noleap","360_day")){
    Ndyr=sum(ndays(1:12,1,caltype)); CNdyr=cumsum(ndays(1:12,1,caltype))
    y=floor(refy+jday/Ndyr)+(floor(refy+jday/Ndyr)>=0)
    diny=jday-Ndyr*(y-refy-(y>0))
    y=y-(diny==0); diny=diny+(diny==0)*Ndyr
    dts=matrix(c(diny,y),ncol=2)
    m=apply(dts,1,function(x) min(which(sort(c(CNdyr,x[1]))==x[1])))
    dts=cbind(m,dts)
    d=apply(dts,1,function(x) x[2]-(x[1]!=1)*sum(ndays(1:(x[1]-1),x[3],caltype)))
    out=list(m=m,d=d,y=y)
  }
  out
}

# Function creating dates (data.frame) from m,d,y
makedates=function(m=NULL,d=NULL,y=NULL,avg="day",conv="cdo",caltype="leap"){
  out=c()
  if (length(y)!=0){ for (yy in y){ M=m; if (length(m)==0) M=1:12; for (mm in M){
      D=d; if (length(d)==0) D=1:ndays(mm,yy,caltype)
      if (avg=="mon" & conv=="cdo") D=last(D)
      if (avg=="mon" & conv=="nco") D=ceiling(median(D))
      if (avg=="mon" & conv=="1")   D=1
      for (dd in D) out=c(out,mm,dd,yy)
  }}}
  if (length(out)!=0) {
    out=data.frame(matrix(out,ncol=3,byrow=TRUE)); names(out)=c("m","d","y") }
  out
}

# Function creating dates corresponding to season over period
# (if inside=TRUE and season=DJF, JF of first year and D of last year will be included)
makedates.season=function(season,period,avg="day",inside=FALSE,caltype="leap",conv="cdo"){
  # Getting indices corresp. to season
  is=season2i(season)
  msta=min(is$i); mend=max(is$i)+is$n[which(is$i==max(is$i))]-1
  # Computing of dates
  dts=c(); for (y in period){ for (s in 1:length(season)){
    issta=is$i[s]; isend=is$i[s]+is$n[s]-1
    if (isend<=12)
      dts=c(dts,c(t(as.matrix(makedates(m=issta:isend,y=y,avg=avg,caltype=caltype,conv=conv)))))
    else {
      if (inside) dts=c(dts,c(t(as.matrix(makedates(m=1:mymod(isend),y=y,avg=avg,caltype=caltype,conv=conv)))),c(t(as.matrix(makedates(m=issta:12,y=y,avg=avg,caltype=caltype,conv=conv)))))
      else dts=c(dts,c(t(as.matrix(makedates(m=issta:12,y=y,avg=avg,caltype=caltype,conv=conv)))),c(t(as.matrix(makedates(m=1:mymod(isend),y=y+1,avg=avg,caltype=caltype,conv=conv)))))
    }
  }}
  dts=data.frame(matrix(dts,ncol=3,byrow=TRUE)); names(dts)=c("m","d","y")
  # Removing duplicated dates
  if (any(duplicated(dts))) dts=dts[-which(duplicated(dts)),]
  dts
}

# Convert date $m$d$y to time index wrt date0
mdy2i=function(date,date0,avg="day",caltype="leap"){
  if (avg=="day") {
    fun=function(x,y) myjulday(x[1],x[2],x[3],caltype)-myjulday(y[1],y[2],y[3],caltype)+1
    out=apply(as.matrix(date),1,fun,as.matrix(date0))
  }
  if (avg=="mon") out=apply(as.matrix(date),1,function(x,y) 12*(x[3]-y[3])+x[1],as.matrix(date0))
  out
}

# Convert date $m$d$y to time values wrt units (and inverse function)
mdy2time=function(date,units,caltype="leap"){
  what=strsplit(units," ")[[1]][1]
  date0=as.numeric(strsplit(strsplit(units," ")[[1]][3],"-")[[1]])[c(2,3,1)]
  fun=function(x,y)
    myjulday(x[1],x[2],x[3],caltype)-myjulday(y[1],y[2],y[3],caltype)
  diff=apply(as.matrix(date),1,fun,date0)
  if (what=="days")    out=diff
  if (what=="hours")   out=24*diff
  if (what=="minutes") out=1440*diff
  if (what=="seconds") out=86400*diff
  out
}
time2mdy=function(time,units,caltype="leap"){
  what=strsplit(units," ")[[1]][1]
  date0=as.numeric(strsplit(strsplit(units," ")[[1]][3],"-")[[1]])[c(2,3,1)]
  m0=date0[1]; d0=date0[2]; y0=date0[3]
  fun=function(x)
    as.data.frame(mycaldat(x+myjulday(m0,d0,y0,caltype),caltype))
  if (what=="days")    out=fun(time)
  if (what=="hours")   out=fun(time/24)
  if (what=="minutes") out=fun(time/1440)
  if (what=="seconds") out=fun(time/86400)
  names(out)=c("m","d","y")
  out
}

# Convert date $m$d$y to yyyymmdd (and inverse function)
mdy2yyyymmdd=function(date) 10000*date$y+100*date$m+date$d
yyyymmdd2mdy=function(x){
  yy=trunc(x/10000); mm=trunc((x-10000*yy)/100)
  data.frame(m=mm,d=x-100*trunc(x/100),y=yy)
}

# Convert date $m$d$y to cdo-date "yyyy-mm-dd"
mdy2cdo=function(date){
  out=c(); for (i in 1:nrow(date)){
    yyyy=paste(c(rep(0,4-nchar(date$y[i])),date$y[i]),collapse="")
    mm=paste(c(rep(0,2-nchar(date$m[i])),date$m[i]),collapse="")
    dd=paste(c(rep(0,2-nchar(date$d[i])),date$d[i]),collapse="")
    out=c(out,paste(yyyy,mm,dd,sep="-"))
  }
  out
}


#----------------------------------------------------------------------
# NETCDF functions (new version with package ncdf4)
#----------------------------------------------------------------------
# Create netcdf
mync=function(file,dat,var,lon,lon.units,lat,lat.units,time,time.units,caltype=NA){
  if (length(lon)==0)  lon=0  ; if (length(lon.units)==0)  lon.units=""
  if (length(lat)==0)  lat=0  ; if (length(lat.units)==0)  lat.units=""
  if (length(time)==0) time=0 ; if (length(time.units)==0) time.units=""
  LON=ncdim_def("lon",lon.units,lon,create_dimvar=TRUE)
  LAT=ncdim_def("lat",lat.units,lat,create_dimvar=TRUE)
  TIME=ncdim_def("time",time.units,time,unlim=TRUE,create_dimvar=TRUE,calendar=caltype)
  varnc=ncvar_def(var,"",list(LON,LAT,TIME),missval=1.e+30)
  ncnew=nc_create(file,varnc)
  if (length(dim(dat))==2) dat=array(dat,dim=c(dim(dat),1))
  ncvar_put(ncnew,var,dat,start=NA,count=dim(dat))
  nc_close(ncnew); rm(ncnew); gc()
  system(paste("ncatted -a _FillValue,",var,",c,f,1.e+30 ",file,sep=""))
}

# Open netcdf > list (dealing with fill value when missing value not specified)
myno=function(file,varid=NA){
  nc=nc_open(file)
  if (system(paste("ncdump -h",file,"| grep 'missing_value ='"))==1
      & system(paste("ncdump -h",file,"| grep '_FillValue ='"))!=1 ){
    system(paste("ncdump -h",file,"| grep '_FillValue ='  > tmp"))
    tmp=scan("tmp",what="character"); system("rm -f tmp")
    varname=strsplit(tmp[1],split=":")[[1]][1]
    fv=strsplit(tmp[3],split="")[[1]]
    fillval=paste(fv[-length(fv)],collapse="")
    nc$var[[paste(varname)]]$missval=as.numeric(fillval)
  }
  out=list(varname=last(names(nc$var)))
  out$lon=nc$dim$lon$vals; out$lonu=nc$dim$lon$units
  out$lat=nc$dim$lat$vals; out$latu=nc$dim$lat$units
  if (any(names(nc$dim)=="plev")) {out$lev=nc$dim$plev$vals; out$levu=nc$dim$plev$units}
  out$time=floor(as.numeric(nc$dim$time$vals)); out$timeu=nc$dim$time$units
  if (any(names(nc$dim$time)=="calendar")) out$caltype=nc$dim$time$calendar
  out$var=ncvar_get(nc,varid=varid); nc_close(nc); rm(nc); gc()
  out
}

# Only getting dimensions of a netcdf (not the var)
mynd=function(file){
  nc=nc_open(file)
  out=list(varname=names(nc$var))
  out$lon=nc$dim$lon$vals; out$lonu=nc$dim$lon$units
  out$lat=nc$dim$lat$vals; out$latu=nc$dim$lat$units
  if (any(names(nc$dim)=="plev")) {out$lev=nc$dim$plev$vals; out$levu=nc$dim$plev$units}
  out$time=floor(as.numeric(nc$dim$time$vals)); out$timeu=nc$dim$time$units
  if (any(names(nc$dim$time)=="calendar")) out$caltype=nc$dim$time$calendar
  nc_close(nc); rm(nc); gc()
  out
}

# View netcdf with ncview
mynv=function(file) system(paste("ncview",file,"&"))

# Applying a cdo command to one (or more) ifile(s) -> ofile
mycdo=function(cdocmd,ifiles,ofile)
  system(paste("cdo",cdocmd,ifiles,ofile))

