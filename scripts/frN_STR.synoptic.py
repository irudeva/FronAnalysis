from netCDF4 import Dataset
import datetime as datetime
import calendar
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import scipy

import os.path

# import xlwt,xlrd
# from xlutils.copy import copy
# import os

# ************************************************
#  Excel
# ************************************************

def output(ireg,reg,filename, sheet, title,tperiod, cc):
    pad = nt+5

    if os.path.exists(filename):
        rb = xlrd.open_workbook(filename, formatting_info=True)
        sheets = rb.sheet_names()
        wb = copy(rb)
        if sheet in sheets:
            # r_sheet = rb.sheet_by_name(sheet)
            # print r_sheet
            sh = wb.get_sheet(sheets.index(sheet))
        else:
            sh = wb.add_sheet(sheet)
            for it in np.arange(nt):
                sh.write(it+3,0,tperiod[it+1])
                sh.write(pad+it,0,tperiod[it+1])
            sh.write(0,0,title)
            sh.write(2,0,"corr")
            sh.write(pad,0,"corr_dt")

    else:
        wb = xlwt.Workbook()
        sh = wb.add_sheet(sheet)
        for it in np.arange(nt):
            sh.write(it+3,0,tperiod[it+1])
            sh.write(pad+it+1,0,tperiod[it+1])
        sh.write(0,0,title)
        sh.write(2,0,"corr")
        sh.write(pad,0,"corr_dt")


    sh.write(2,ireg+1,reg)
    sh.write(pad,ireg+1,reg)

    for it in np.arange(nt):
        # sh.write(it,ireg+0,tperiod[it+1])
        sh.write(it+3,ireg+1,cc[0,it])
        sh.write(pad+it+1,ireg+1,cc[1,it])

    wb.save(filename)

def output_trend(ireg,reg,filename, sheet, title,tperiod, trend):
    pad = nt+5

    if os.path.exists(filename):
        rb = xlrd.open_workbook(filename, formatting_info=True)
        sheets = rb.sheet_names()
        wb = copy(rb)
        if sheet in sheets:
            # r_sheet = rb.sheet_by_name(sheet)
            # print r_sheet
            sh = wb.get_sheet(sheets.index(sheet))
        else:
            sh = wb.add_sheet(sheet)
            for it in np.arange(nt):
                sh.write(it+3,0,tperiod[it+1])
                sh.write(pad+it+1,0,tperiod[it+1])
            sh.write(0,0,title)
            sh.write(2,0,"slope")
            sh.write(pad,0,"p value")

    else:
        wb = xlwt.Workbook()
        sh = wb.add_sheet(sheet)
        for it in np.arange(nt):
            sh.write(it+3,0,tperiod[it+1])
            sh.write(pad+it+1,0,tperiod[it+1])
        sh.write(0,0,title)
        sh.write(2,0,"slope")
        sh.write(pad,0,"p value")


    sh.write(2,ireg+1,reg)
    sh.write(pad,ireg+1,reg)

    for it in np.arange(nt):
        # sh.write(it,ireg+0,tperiod[it+1])
        sh.write(it+3,ireg+1,trend[0,it])
        # if trend[1,it-1] < .05:
        sh.write(pad+it+1,ireg+1,trend[1,it])

    wb.save(filename)


def trend(ireg,reg,nt,yrs,var,var_tr,title,tperiod,filetrend):

    for it in np.arange(nt):
        mask = ~np.isnan(var[0,it,:])
        masknan = np.isnan(var[0,it,:])
        slope, intercept, r_value, p_value, std_err = stats.linregress(yrs[mask],var[0,it,mask])
        var[1,it,mask] = intercept + (slope * yrs[mask])
        var[1,it,masknan] = None
        var_tr[0,it] = slope
        var_tr[1,it] = p_value

    output_trend(ireg,reg,filetrend, "trend", title,tperiod,var_tr)


# for multiple regression
def fn(x, a, b, c):
    return a + b*x[0] + c*x[1]



# ************************************************
#  Selection
# ************************************************

# time
year = [1979,2016]
yrs = np.arange(year[0],year[1],1)
nyrs = np.size(yrs)
x = np.arange(0,nyrs-1,1)

mm = [1,12]  #[ 9,9]
dd = [1,31]
hr = [0,18]

lag = 5 #lag in days

#choose time scale
tscale = "ssn"
# tscale = "mon"

# ssn   / use "YYY" for whole year
ssn=["","YYY","DJF","MAM","JJA","SON"]
month_abbr = ["","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"]
if tscale == "ssn":
    nt = 5+1      #for ssn
    tperiod = ssn
if tscale == "mon":    #for monthly analysis
    nt = 12+1
    tperiod = month_abbr

# nf
nf1 = 0# 5518 ;0
nf2 = 20000# 5518 ;20000

nfig = 1  # possible combinations of vars
fig1 = 0
fig2 = 1


hs = "SH"
# reg
nreg = 5

cc = np.zeros((nreg+1,nfig+1,2,nt),dtype=np.float) # correlation


for ireg in np.arange(3,nreg+1):
    if ireg == 0:
        lon = [ -90,361]
        lat = [ -40, -20]
        reg = "%d_%dS"%(np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 1:
        lon = [ 30,90]
        lat = [ -40, -20]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 2:
        lon = [ 90,150]
        lat = [ -40, -20]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 3:
        lon = [ 150,210]
        lat = [ -40, -20]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 4:
        lon = [ 210,285]
        lat = [ -40, -20]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 5:
        lon = [ 300,340]
        lat = [ -40, -20]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))



    maxnf = 200 # max # of fronts per timestep
    maxnp = 100 # max # of frontal points


    # ************************************************
    # read in netCDF files
    # ************************************************
    # 1. slp
    # ************************************************

    STRlat = np.zeros((nt,nyrs,1464+lag*4*2+31*4),dtype=np.float); STRlat.fill(np.nan)
    STRslp = np.zeros_like(STRlat); STRslp.fill(np.nan)
    T_hr   = np.zeros((nt,nyrs,1464+lag*4*2+31*4),dtype=np.float); T_hr.fill(np.nan)

    STRslp_clim = np.zeros((12,31,4,nyrs),dtype=np.float);STRslp_clim.fill(np.nan)
    STRlat_clim = np.zeros_like(STRslp_clim);STRlat_clim.fill(np.nan)

    for yr in yrs:


        for slpyr in [yr,yr-1,yr+1]:
            print "STR %d for %d "%(slpyr,yr)
            fin = "/Users/irudeva/work/DATA/ERAint/Mslp_highres/erain_hr.mslp.%d.10_70S.nc"%slpyr
            print fin, os.path.isfile(fin)

            if os.path.isfile(fin):
                print "Slp from", fin
                slpnc = Dataset(fin, 'r')

                print slpnc.variables.keys()
                if slpnc.variables.keys()[0]=="lon":
                  lonslp = slpnc.variables['lon'][:]
                  latslp = slpnc.variables['lat'][:]
                elif slpnc.variables.keys()[0]=="longitude":
                  lonslp = slpnc.variables['longitude'][:]
                  latslp = slpnc.variables['latitude'][:]
                else:
                 print"Check lon/lat names in the netcdf file, %d "%slpyr
                time = slpnc.variables['time'][:]
                mslp = slpnc.variables['msl'][:]/100.

                if np.any(lonslp<0):
                    print "!!!!WARNING!!!! Check lons"
                    print lonslp
                    print "Converting from [-180,180) to [0,360)"

                    lonslp_old = lonslp
                    lonslp = lonslp_old[lonslp_old>=0]
                    lonslp = np.append(lonslp,lonslp_old[lonslp_old<0]+360.)


                    mslp_old = mslp
                    mslp = mslp_old[:,:,lonslp_old>=0]
                    mslp = np.append(mslp,mslp_old[:,:,lonslp_old<0],axis=2)



                dt_slp = [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
                   for t in time]
                if slpyr == yr :
                    ndt = len(dt_slp)

                STRlat_yr = np.zeros(len(dt_slp),dtype=np.float); STRlat_yr.fill(np.nan)
                STRslp_yr = np.zeros_like(STRlat_yr); STRslp_yr.fill(np.nan)

                latslpSH  = latslp[np.logical_and(latslp<-20,latslp>-65)]
                mslpSH    = mslp[:,np.logical_and(latslp<-20,latslp>-65),:]
                print mslpSH.shape
                print lonslp
                mslpSHlon = mslpSH[:,:,np.logical_and(lonslp<lon[1],lonslp>lon[0])]
                print mslpSHlon
                # mslpSHlonDec = mslpSHlon[-1,:,:]

                # zonal average
                slp_z = np.mean(mslpSHlon,axis=2)

                for it in np.arange(len(dt_slp)):
                    STRslp_yr[it] = np.amax(slp_z[it,:])
                    STRlat_yr[it] = latslpSH[np.argmax(slp_z[it,:])]

                    if slpyr <= year[1]-1:
                        STRslp_clim[dt_slp[it].month-1,dt_slp[it].day-1,dt_slp[it].hour/6-1,slpyr-year[0]]=STRslp_yr[it]
                        STRlat_clim[dt_slp[it].month-1,dt_slp[it].day-1,dt_slp[it].hour/6-1,slpyr-year[0]]=STRlat_yr[it]

                if slpyr == yr:
                    tarr = np.arange((lag+31)*4,(lag+31)*4+ndt,1)
                    STRlat[0,yr-yrs[0],tarr] = STRlat_yr[:ndt]
                    STRslp[0,yr-yrs[0],tarr] = STRslp_yr[:ndt]
                    T_hr[0,yr-yrs[0],tarr] = time[:len(dt_slp)]


                if slpyr == yr-1:
                    tarr = np.arange(0,(lag+31)*4,1)
                    STRlat[0,yr-yrs[0],tarr] = STRlat_yr[len(dt_slp)-(lag+31)*4:]
                    STRslp[0,yr-yrs[0],tarr] = STRslp_yr[len(dt_slp)-(lag+31)*4:]
                    T_hr[0,yr-yrs[0],tarr] = time[len(dt_slp)-(lag+31)*4:]


                if slpyr == yr+1:
                    tarr = np.arange(ndt+(lag+31)*4,ndt+lag*4*2+31*4,1)
                    STRlat[0,yr-yrs[0],tarr] = STRlat_yr[:lag*4]
                    STRslp[0,yr-yrs[0],tarr] = STRslp_yr[:lag*4]
                    T_hr[0,yr-yrs[0],tarr] = time[:lag*4]

        # print [t for t in T_hr]
        # print [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
        #    for t in T_hr[~np.isnan(T_hr)]]
        # quit()

        for iss,issn in enumerate(ssn[1:]):
            if issn == 'YYY':
              mons = np.arange(1,13)
              dt1 = (datetime.datetime(yr, mons[0], 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
              dt2 = (datetime.datetime(yr+1, 1, 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
            if issn == 'MAM':
              mons = np.array([3,4,5])
              dt1 = (datetime.datetime(yr, mons[0], 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
              dt2 = (datetime.datetime(yr, mons[2]+1, 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
            if issn == 'JJA':
              mons = np.array([6,7,8])
              dt1 = (datetime.datetime(yr, mons[0], 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
              dt2 = (datetime.datetime(yr, mons[2]+1, 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
            if issn == 'SON':
              mons = np.array([9,10,11])
              dt1 = (datetime.datetime(yr, mons[0], 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
              dt2 = (datetime.datetime(yr, mons[2]+1, 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
            if issn == 'DJF':
              mons = np.array([12,1,2])
              dt1 = (datetime.datetime(yr-1, mons[0], 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
              dt2 = (datetime.datetime(yr, mons[2]+1, 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()

            print ssn[iss+1]

            t1 = int(dt1 // (3600))
            t2 = int(dt2 // (3600))
            # print t1,t2

            t1_lag = t1-lag*24
            t2_lag = t2+lag*24
            # print "lagged t1,ts:",t1_lag,t2_lag

            mask = np.logical_and(T_hr[0,yr-yrs[0],:]<t2_lag,T_hr[0,yr-yrs[0],:]>=t1_lag)
            # print T_hr[mask]

            T_ssn  = T_hr[0,yr-yrs[0],mask]
            # print T_ssn.size
            # print T_ssn
            # print [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
            #    for t in T_ssn]

            m = T_ssn.size
            # STRlat[iss+1,yr-yrs[0],:m] = STRlat[0,yr-yrs[0],mask]-np.mean(STRlat[0,yr-yrs[0],mask])
            # STRslp[iss+1,yr-yrs[0],:m] = STRslp[0,yr-yrs[0],mask]-np.mean(STRslp[0,yr-yrs[0],mask])
            STRlat[iss+1,yr-yrs[0],:m] = STRlat[0,yr-yrs[0],mask]
            STRslp[iss+1,yr-yrs[0],:m] = STRslp[0,yr-yrs[0],mask]
            T_hr[iss+1,yr-yrs[0],:m]   = T_hr[0,yr-yrs[0],mask]



        slpnc.close()



    # to subtract the daily average

    for yr in yrs:
        for iss,issn in enumerate(ssn[1:]):
            for ct,it in enumerate(T_hr[iss+1,yr-yrs[0],:]):
                if ~np.isnan(it):
                    date = datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(it))
                    # if date.month==2 and date.day==29:
                    #     print date
                    #     print STRslp[iss+1,yr-yrs[0],ct]
                    #     print STRslp_clim[date.month-1,date.day-1,date.hour/6-1,:]
                    #     print STRslp[iss+1,yr-yrs[0],ct] - np.mean(STRslp_clim[date.month-1,date.day-1,date.hour/6-1,~np.isnan(STRslp_clim[date.month-1,date.day-1,date.hour/6-1,:])])
                    #     quit()
                    STRslp[iss+1,yr-yrs[0],ct] = STRslp[iss+1,yr-yrs[0],ct] - np.mean(STRslp_clim[date.month-1,date.day-1,date.hour/6-1,~np.isnan(STRslp_clim[date.month-1,date.day-1,date.hour/6-1,:])])
                    STRlat[iss+1,yr-yrs[0],ct] = STRlat[iss+1,yr-yrs[0],ct] - np.mean(STRlat_clim[date.month-1,date.day-1,date.hour/6-1,~np.isnan(STRlat_clim[date.month-1,date.day-1,date.hour/6-1,:])])



    # write max/min STR years to a file
    # for it in np.arange(nt):
    #     print tperiod[it+1],STR1slp[0,it,:]
    #     fSTRint = "../output/STRint.%s.%s%d_%d.txt"%(reg,tperiod[it+1],year[0],year[1]-1)
    #     with open(fSTRint, "w") as text_file:
    #         text_file.write("{:>7}{:>7}{:>7}{:>7}\n".format("STRmax",  "Yrmax",  "STRmin",   "Yrmin"))
    #
    #         nsel = 10
    #         top = (-STR1slp[0,it,:]).argsort()[:nsel]
    #         bottom = (STR1slp[0,it,:]).argsort()[:nsel]
    #
    #         for i in range(nsel):
    #             text_file.write("{0[0]:.2f} {0[1]:6d} {0[2]:.2f} {0[3]:6d}\n".format([STR1slp[0,it,top][i],yrs[top][i],STR1slp[0,it,bottom][i],yrs[bottom][i]]))
    #
    #     fSTRlat = "../output/STRloc.%s.%s%d_%d.txt"%(reg,tperiod[it+1],year[0],year[1]-1)
    #     with open(fSTRlat, "w") as text_file:
    #         text_file.write("{:>7}{:>7}{:>7}{:>7}\n".format("STRmax",  "Yrmax",  "STRmin",   "Yrmin"))
    #
    #         nsel = 10
    #         top = (STR1lat[0,it,:]).argsort()[:nsel]
    #         bottom = (-STR1lat[0,it,:]).argsort()[:nsel]
    #
    #         for i in range(nsel):
    #             text_file.write("{0[0]:.2f} {0[1]:6d} {0[2]:.2f} {0[3]:6d}\n".format([STR1lat[0,it,top][i],yrs[top][i],STR1lat[0,it,bottom][i],yrs[bottom][i]]))
    #
        # quit()

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #  Fronts
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    for yr in yrs:
      fin = "../frontgrd/frline.%d.nc"%yr
      print "fronts from %s"%fin
      ncf = Dataset(fin, 'r')

      print ncf.variables.keys()
      flon = ncf.variables['flon'][:]
      flat = ncf.variables['flat'][:]
      time = ncf.variables['time'][:]
      dv = ncf.variables['dv'][:]
      nf = ncf.variables['nfr'][:]
      npts = ncf.variables['npts'][:]

      fdt = [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
           for t in time]


      if np.any(nf>maxnf) :
        print np.amax(nf)
        print "ERROR: nf > maxnf"
        quit()

      if np.any(npts > maxnp) :
        print np.amax(npts)
        print "ERROR:  np > maxnp"
        quit()

      # ************************************************
      # masking
      # ************************************************
      # create mask

      #   from 1 Dec Yr-1 to 31 Dec Yr

      for ind,t in enumerate(fdt):
          if yr==year[0]:
              if t.year == yr and t.month==1 and t.day == 1 and t.hour == 0:
                  tfr1 = ind
                  tmp = (datetime.datetime(yr, 1, 1, 0) -  \
                             datetime.datetime(yr-1, 12, 1, 0)).total_seconds()
                  tfr0 = int(tmp//(3600))/6
          else:
              if t.year == yr-1 and t.month==12 and t.day == 1 and t.hour == 0:
                  tfr1 = ind
                  tfr0 = 0
          if t.year == yr and t.month==12 and t.day == 31 and t.hour == 18:
               tfr2 = ind
               break

      if yr == year[0]:
        # find the size of the array:
        yrtmp = 1980
        if calendar.isleap(yrtmp):
            dt_delta = (datetime.datetime(yrtmp+1, 1, 31, 18) -  \
                       datetime.datetime(yrtmp-1, 12, 1, 0)).total_seconds()
            dt_delta_hr = int(dt_delta // (3600))
            # dt_tmp = [ datetime.datetime(yrtmp-1, 12, 1, 0)+i*datetime.timedelta(hours=6) for i in range(dt_delta_hr/6+1)]
            stfr = dt_delta_hr/6+1
        else:
            print "%d is not a leap year"%yrtmp
            quit()


        #   if calendar.isleap(yr):
        #       frmask = np.zeros((nyrs,tfr2-tfr1+1,maxnf),dtype=np.int)
        #       frdv  = np.zeros((nyrs,tfr2-tfr1+1,maxnf),dtype=np.float)
        #       frmask_it = np.zeros((nt,nyrs,tfr2-tfr1+1,maxnf),dtype=np.int)
        #       frdv_it  = np.zeros((nt,nyrs,tfr2-tfr1+1,maxnf),dtype=np.float)
        #   else:
        frmask = np.zeros((nyrs,stfr,maxnf),dtype=np.int)
        frdv   = np.zeros((nyrs,stfr,maxnf),dtype=np.float32)
        frmask_it = np.zeros((nt,nyrs,stfr,maxnf),dtype=np.int)
        frdv_it   = np.zeros((nt,nyrs,stfr,maxnf),dtype=np.float32)

        frNPlat = np.zeros_like(frdv)
        frNPlat_it = np.zeros_like(frdv_it)

        frN = np.zeros([nt,nyrs,stfr],dtype=np.float); frN.fill(np.nan)
        frT = np.zeros([nt,nyrs,stfr],dtype=np.float); frT.fill(np.nan)

        frN_clim = np.zeros((12,31,4,nyrs),dtype=np.float);frN_clim.fill(np.nan)



      for ct in np.arange(tfr1,tfr2+1):
          if fdt[ct].year>=year[0] and fdt[ct].year<=year[1]:
                for ifr in range(nf[ct]):
                  if ifr+1>=nf1 and ifr+1<=nf2:
                      for ip in range(npts[ct,ifr]):
                          if flat[ct,ifr,ip]>=lat[0] and flat[ct,ifr,ip]<=lat[1]:
                              if flon[ct,ifr,ip]>=lon[0] and flon[ct,ifr,ip]<=lon[1]:
                                  frmask[yr-year[0],tfr0+(ct-tfr1),ifr] = time[ct]
                                  frdv[yr-year[0],tfr0+(ct-tfr1),ifr] = np.sum(dv[ct,ifr,:npts[ct,ifr]])
                                  frNPlat[yr-year[0],tfr0+(ct-tfr1),ifr] = np.amax(flat[ct,ifr,:npts[ct,ifr]])
                                  break


      for ct in np.arange(tfr1,tfr2+1):
          frN[0,yr-year[0],tfr0+(ct-tfr1)] = np.sum(np.where(frmask[yr-year[0],:,:]==time[ct],1,0))
          print yr, ct, frN[0,yr-year[0],tfr0+(ct-tfr1)]
          frT[0,yr-year[0],tfr0+(ct-tfr1)] = time[ct]

          frN_clim[fdt[ct].month-1,fdt[ct].day-1,fdt[ct].hour/6-1,yr-year[0]]=frN[0,yr-year[0],tfr0+(ct-tfr1)]

        #   print ct,tfr0+(ct-tfr1),fdt[ct],time[ct],frT[0,yr-year[0],tfr0+(ct-tfr1)],frN[0,yr-year[0],tfr0+(ct-tfr1)], tfr0

    #   if yr-year[0]+1 < nyrs:
    #       for ct in np.arange(len(time)-31*4,len(time)):
    #           frN[0,yr-year[0]+1,ct-len(time)+31*4] = np.sum(np.where(frmask[yr-year[0],:,:]==time[ct],1,0))
    #           frT[0,yr-year[0]+1,ct-len(time)+31*4] = time[ct]

      for iss,issn in enumerate(ssn[1:]):
        if issn == 'YYY':
          mons = np.arange(1,13)
          dt1 = (datetime.datetime(yr, mons[0], 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
          dt2 = (datetime.datetime(yr+1, 1, 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
        if issn == 'MAM':
          mons = np.array([3,4,5])
          dt1 = (datetime.datetime(yr, mons[0], 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
          dt2 = (datetime.datetime(yr, mons[2]+1, 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
        if issn == 'JJA':
          mons = np.array([6,7,8])
          dt1 = (datetime.datetime(yr, mons[0], 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
          dt2 = (datetime.datetime(yr, mons[2]+1, 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
        if issn == 'SON':
          mons = np.array([9,10,11])
          dt1 = (datetime.datetime(yr, mons[0], 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
          dt2 = (datetime.datetime(yr, mons[2]+1, 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
        if issn == 'DJF':
          mons = np.array([12,1,2])
          dt1 = (datetime.datetime(yr-1, mons[0], 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()
          dt2 = (datetime.datetime(yr, mons[2]+1, 1, 0)-datetime.datetime(1900, 1, 1, 0)).total_seconds()

        print ssn[iss+1]

        t1 = int(dt1 // (3600))
        t2 = int(dt2 // (3600))
        # print t1,t2

        # print t1, t2
        # print frT[0,yr-year[0],tfr0:tfr2+tfr0]

        mask = np.logical_and(frT[0,yr-year[0],:]<t2,frT[0,yr-year[0],:]>=t1)
        # print T_hr[mask]

        T_ssn  = frT[0,yr-year[0],mask]
        # print T_ssn
        # print [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
        #    for t in T_ssn]

        m = T_ssn.size
        # frN[iss+1,yr-yrs[0],:m] = frN[0,yr-year[0],mask]-np.mean(frN[0,yr-year[0],mask])
        frN[iss+1,yr-yrs[0],:m] = frN[0,yr-year[0],mask]
        frT[iss+1,yr-yrs[0],:m ]= frT[0,yr-year[0],mask]

    # to subtract daily means
    for yr in yrs:
        for iss,issn in enumerate(ssn[1:]):
            for ct,it in enumerate(frT[iss+1,yr-yrs[0],:]):
                if ~np.isnan(it):
                    date = datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(it))
                    frN[iss+1,yr-yrs[0],ct] = frN[iss+1,yr-yrs[0],ct] - np.mean(frN_clim[date.month-1,date.day-1,date.hour/6-1,:])

    # ************************************************
    # flattening and correlation
    # ************************************************
    corr= np.zeros([2,nt,lag*8+1],dtype = np.float)

    for iv,STRsel in enumerate(["STRlat","STRslp"]):
        if STRsel == "STRlat":
            STR = STRlat
        elif STRsel == "STRslp":
            STR = STRslp


        for iss,issn in enumerate(ssn[1:]) :
            print ssn[iss+1]
            for il,ilag in enumerate(np.arange(-lag*6*4,lag*6*4+1,6)):
                if issn == "YYY":
                    STR_yy = np.zeros_like(STR[iss+1,:,:],dtype=float); STR_yy.fill(np.nan)
                    T_hr_yy = np.zeros_like(T_hr[iss+1,:,:],dtype=float); T_hr_yy.fill(np.nan)
                    for yy in np.arange(nyrs):
                        mask = np.in1d(T_hr[iss+1,yy,:],frT[iss+1,yy,:]+ilag)
                        ms = STR[iss+1,yy,mask].size
                        STR_yy[yy,:ms] = STR[iss+1,yy,mask]
                        T_hr_yy[yy,:ms]   = T_hr[iss+1,yy,mask]
                else:
                    STR_yy = STR[iss+1,:,:]
                    T_hr_yy   = T_hr[iss+1,:,:]

                frNflat = np.ndarray.flatten(frN[iss+1,:,:])
                frTflat = np.ndarray.flatten(frT[iss+1,:,:])
                mask = ~np.isnan(frNflat)

                # STRflat = np.ndarray.flatten(STR[iss+1,:,:])
                STRflat = np.ndarray.flatten(STR_yy)
                T_hrflat = np.ndarray.flatten(T_hr_yy)
                # mask2 = np.in1d(T_hr[iss+1,:,:],frT[iss+1,:,: ]+ilag)
                mask2 = np.in1d(T_hrflat,frTflat[mask]+ilag)


                if frNflat[mask].size != STRflat[mask2].size:
                    print "!!!! WARNING: mask sizes are not the same %d vs %d"%(frNflat[mask].size,STRflat[mask2].size)
                    tmp = frTflat[mask]+ilag
                    tmpfrN = frNflat[mask]
                    mask3 = np.in1d(tmp, T_hrflat[mask2])
                    # print frTflat[mask3].size
                    # print "mask2:",T_hrflat[mask2]
                    # print "mask+ilag:",frTflat[mask]+ilag
                    frT_masked = tmp[mask3]
                    frN_masked = tmpfrN[mask3]
                else:
                    frT_masked = frTflat[mask]+ilag
                    frN_masked = frNflat[mask]

                corr[iv,iss+1,il] = np.corrcoef(frN_masked,STRflat[mask2])[1,0]
                # slope, intercept, r_value, p_value, std_err = stats.linregress(frN_masked,STRflat[mask2])

                print ssn[iss+1],il," ",ilag,"h ",corr[iv,iss+1,il]


    # quit()
    # #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # # Plotting
    # # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    print "start plotting..."

    # majorLocator = MultipleLocator(5)
    # majorFormatter = FormatStrFormatter('%d')
    # minorLocator = MultipleLocator(1)

    if tscale == "mon":
        nc = 2
        nr = 6
        wdth = 3.5
    if tscale == "ssn":
        nc = 1
        nr = nt-1
        # wdth = 3.
        wdth = 2.



    # plt.ion()
    plt.close('all')
    for ifig in np.arange(fig1,fig2+1):
        print ifig, nr,nc
        # fig = plt.figure(ifig)
        # f, ax = plt.subplots(nr, nc,  sharex='col', sharey='row',figsize=(wdth*3.13,4.2*3.13))
        f, ax = plt.subplots(nr, nc,  sharex='col', figsize=(wdth*3.13,4.2*3.13))

        if ifig ==0:
            print ifig
            title = 'lag correlation: STR location vs number of fronts, %s'%reg
            fout = "corrlag.frN_STRloc"
        if ifig ==1:
            print ifig
            title = 'lag correlation: STR intensity vs number of fronts, %s'%reg
            fout = "corrlag.frN_STRint"
        plt.suptitle(title,fontsize=14)
        print title


        for ic in np.arange(nc):
        # for ic in [0,1]:
            for ir in np.arange(nr):
            # for ir in [2,3,4,5]:
                it = ir+ic*6
                print tperiod[it+1]," start"

                if tscale == 'mon':
                    a = ax[ir, ic]
                if tscale == 'ssn':
                    a = ax[it]
                # a2 = a.twinx()

                var1 = corr[ifig,1+it,:]
                a.plot(np.arange(-lag*6*4,lag*6*4+1,6),var1,color='b',lw=1)
                # if ifig ==11:
                #     a2.plot(yrs,var2[it,:],color='g',lw=1,label=var2legend)
                # else:
                #     a2.plot(yrs,var2[0,it,:],color='g',lw=1,label=var2legend)

                # ax[ir, ic].set_title('STR intensity, %s, SH'%tperiod[it+1])

                a.set_ylim(-.7, .7)
                a.set_xlim(-lag*24, lag*24)
                # a2.axis((year[0]-1, year[1], yax2[0], yax2[1]))
                # ax2.set_ylim(0, 35)
                # print tperiod[it+1],np.amin(var1[0,it,:]),np.amax(var1[0,it,:])
                # a.set_ylim(0.9*np.amin(var1[0,it,:]),1.1*np.amax(var1[0,it,:]))
                # a2.set_ylim(0.9*np.amin(var1[0,it,:]),1.1*np.amax(var1[0,it,:]))

                # a.xaxis.set_ticks(yrs)
                # xa = a.get_xaxis()
                # a.xaxis.set_tick_params(axis='x',which='major',length=6,width =2)
                # a.xaxis.set_tick_params(axis='x',which='minor',length=3,width =2)

                # a.xaxis.set_major_locator(majorLocator)
                # a.xaxis.set_major_formatter(majorFormatter)

                # for the minor ticks, use no labels; default NullFormatter
                # a.xaxis.set_minor_locator(minorLocator)

                a.grid()

                if tscale == "mon":
                    # if it == 3:
                    #     a.set_ylabel(yax1label,color='b')
                    # if it == 9:
                    #     a2.set_ylabel(yax2label,color='g')
                    if ir == 5 :
                        a.set_xlabel('lag')
                if tscale == "ssn":
                    # if it == 3 :
                    #     a.set_ylabel(yax1label,color='b')
                    #     # a2.set_ylabel(yax2label,color='g')
                    if it == nt-1 :
                        a.set_xlabel('lag')
                a.set_title(' %s, %s'%(tperiod[it+1],reg))
                print tperiod[it+1]," done"

                #  trend & correlation
                # if ifig != 11:
                #     a.plot(yrs,var1[1,it,:],color='b',lw=2)
                #     a2.plot(yrs,var2[1,it,:],color='g',lw=2)
                #
                #     mask  = ~np.isnan(var1[0,it,:]) & ~np.isnan(var2[0,it,:])
                #     cctmp = np.corrcoef(var1[0,it,mask],var2[0,it,mask])
                # else:
                #     mask  = ~np.isnan(var1[0,it,:]) & ~np.isnan(var2[it,:])
                #     cctmp = np.corrcoef(var1[0,it,mask],var2[it,mask])
                # cc[ireg,ifig,0,it] = cctmp[1,0]
                # if ifig != 11:
                #     cctmp = np.corrcoef(var1[0,it,mask]-var1[1,it,mask],var2[0,it,mask]-var2[1,it,mask])
                #     cc[ireg,ifig,1,it] = cctmp[1,0]
                #     a.set_title(' %s, %s, r = %.2f, r_dt = %.2f'%(tperiod[it+1],reg,cc[ireg,ifig,0,it],cc[ireg,ifig,1,it]))
                # else:
                #     a.set_title(' %s, %s, r = %.2f'%(tperiod[it+1],reg,cc[ireg,ifig,0,it]))
                #
                # if it == nt-1 :
                #     # ask matplotlib for the plotted objects and their labels
                #     lines, labels = a.get_legend_handles_labels()
                #     lines2, labels2 = a2.get_legend_handles_labels()
                #     a2.legend(lines + lines2, labels + labels2, loc='lower right')

                    # ax[ir, ic].legend(loc='lower right')


        # plt.setp([a.set_xlabel('Year') for a in ax[5, :]])


        f.subplots_adjust(hspace=0.3,wspace = 0.3)
        plt.draw()
        f.savefig("../output/%s.%s.mon%d_%d.png"%(fout,reg,year[0],year[1]-1))
        plt.close(f)

    #     # ************************************************
    #     #  correlations to Excel
    #     # ************************************************
    #     fxls = "../output/%s.%d_%d.xls"%(fout,year[0],year[1]-1)
    #     if ireg == 0 and os.path.exists(fxls):
    #             os.remove(fxls)
    #     output(ireg,reg,fxls, "corr", title, tperiod, cc[ireg,ifig,:,:])
