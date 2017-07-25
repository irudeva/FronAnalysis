from netCDF4 import Dataset
import datetime as datetime
import calendar
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import scipy

import xlwt,xlrd
from xlutils.copy import copy
import os

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

#choose time scale
tscale = "ssn"
# tscale = "mon"

#choose type of SAM index
# SAMscale = "SH"
SAMscale = "reg"


# ssn   / use "YYY" for whole year
ssn=["","YYY","DJF","MAM","JJA","SON"]
month_abbr = ["   ","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"]
if tscale == "ssn":
    nt = 5      #for ssn
    tperiod = ssn
if tscale == "mon":    #for monthly analysis
    nt = 12
    tperiod = month_abbr

# nf
nf1 = 0# 5518 ;0
nf2 = 20000# 5518 ;20000

nfig = 3  # possible combinations of vars
fig1 = 0
fig2 = 0


hs = "SH"
# reg
nreg = 5

cc = np.zeros((nreg+1,nfig+1,2,nt),dtype=np.float) # correlation

# ************************************************
# read ONI
# ************************************************
ONI = np.zeros((nyrs,12)); ONI.fill(np.nan)

fONI = "../output/SAM_STR_ENSO/oni.txt" # ONI is a 3-months average ENSO3.4
with open(fONI, "r") as text_file:
    header1 = text_file.readline()
    for line in text_file:
        # print (line)
        columns = line.split()
        # print (columns)
        yr = int(columns[0])

        if np.any(yrs==yr) :
            ONI[yr-year[0],:] = columns[1:13]
        if(yr==year[1]):
            break

    text_file.close

ONIssn = np.zeros((2,nt,nyrs),dtype=np.float);ONIssn.fill(np.nan)
ONIssn_tr = np.zeros((2,nt),dtype=np.float);ONIssn.fill(np.nan)

if tscale=='ssn':
  for it,issn in enumerate(ssn[1:]) :
      print issn
      cy = yr-year[0]
      if issn == 'YYY':
          mons = np.array([2,5,8,11])
      if issn == 'MAM':
          mons = np.array([4])
      if issn == 'JJA':
          mons = np.array([7])
      if issn == 'SON':
          mons = np.array([10])
      if issn == 'DJF':
          mons = np.array([1])

      mons = mons-1

      if issn != "YYY":
          for yr in yrs:
            #   print yr, SAM[yr-year[0],mons]
              if issn != "DJF":
                  ONIssn[0,it,yr-year[0]] = ONI[yr-year[0],mons]
              if yr != year[0] and issn == "DJF":
                #   ONIssn[0,it,yr-year[0]] = ONI[yr-year[0]-1,-1]
                  ONIssn[0,it,yr-year[0]] = ONI[yr-year[0]-1,mons]
      else:
          ONIssn[0,it,:] = np.mean(ONI[:,mons],axis=1)

#       print ONIssn[0,it,:]
#       print ONI[:,mons]
#
# print(ONIssn[0,:,:])
# quit()
# ************************************************
# set region
# ************************************************

for ireg in np.arange(nreg+1):
# for ireg in np.arange(1):

    if ireg == 0:
        lon = [ -90,361]
        lat = [ -90, 0]
        reg = "%d_%dS"%(np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 1:
        lon = [ 30,90]
        lat = [ -90, 0]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 2:
        lon = [ 90,150]
        lat = [ -90, 0]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 3:
        lon = [ 150,210]
        lat = [ -90, 0]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 4:
        lon = [ 210,285]
        lat = [ -90, 0]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 5:
        lon = [ 300,340]
        lat = [ -90, 0]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))



    maxnf = 200 # max # of fronts per timestep
    maxnp = 100 # max # of frontal points


    # ************************************************
    # read in netCDF files
    # ************************************************
    # 1. slp
    # ************************************************

    mslp40 = np.zeros((nyrs,12)); mslp40.fill(np.nan)
    mslp65 = np.zeros((nyrs,12)); mslp65.fill(np.nan)

    SAM = np.zeros((nyrs,12)); SAM.fill(np.nan)
    SAM40 = np.zeros((nyrs,12)); SAM40.fill(np.nan)
    SAM65 = np.zeros((nyrs,12)); SAM65.fill(np.nan)

    for yr in yrs:
      slpyr = yr
      fin = "/Users/irudeva/work/DATA/ERAint/Mslp_highres/erain.mslp.avm.%d.nc"%slpyr
      print "Slp from", fin
      nc = Dataset(fin, 'r')

      # varnam=['longitude','latitude','time','msl']
      # v=0
      # for var in varnam:
      #   if nc.variables[varnam[v]].name != var:
      #       print "Variables don't agree", var, nc.variables[varnam[v]].name, v
      #       exit()
      #   v += 1

      # nc_attrs, nc_dims, nc_vars = ncdump(nc)

      # print nc_attrs, nc_dims, nc_vars
      print nc.variables.keys()
      if nc.variables.keys()[0]=="lon":
          lonslp = nc.variables['lon'][:]
          latslp = nc.variables['lat'][:]
      elif nc.variables.keys()[0]=="longitude":
          lonslp = nc.variables['longitude'][:]
          latslp = nc.variables['lattitude'][:]
      else:
         print"Check lon/lat names in the netcdf file, %d "%yr
      time = nc.variables['time'][:]
      mslp = nc.variables['msl'][:]/100.

      dt_nc = [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
           for t in time]
      print dt_nc[1].year


      latslpSH  = latslp[np.logical_and(latslp<-20,latslp>-65)]
      mslpSH    = mslp[:,np.logical_and(latslp<-20,latslp>-65),:]
      mslpSHlon = mslpSH[:,:,np.logical_and(lonslp<lon[1],lonslp>lon[0])]
      mslpSHlonDec = mslpSHlon[-1,:,:]
      if SAMscale == "SH":
          mslp40[yr-year[0],:]    = np.squeeze(np.mean( mslp[:,latslp==-40,:],axis=2))
          mslp65[yr-year[0],:]    = np.squeeze(np.mean( mslp[:,latslp==-65,:],axis=2))
      elif SAMscale == "reg":
          tmp40 = np.squeeze(mslp[:,latslp==-40,:])
          tmp65 = np.squeeze(mslp[:,latslp==-65,:])
          mslp40[yr-year[0],:]    = np.mean( tmp40[:,np.logical_and(lonslp<lon[1],lonslp>lon[0])],axis=1)
          mslp65[yr-year[0],:]    = np.mean( tmp65[:,np.logical_and(lonslp<lon[1],lonslp>lon[0])],axis=1)


      # zonal average
      slp_z = np.mean(mslpSHlon,axis=2)
      slp_zDec = np.mean(mslpSHlon[-1,:,:],axis=1)

      if yr==year[0]:
          STR1slp = np.zeros((2,nt,nyrs),dtype=np.float)
          STR1lat = np.zeros_like(STR1slp)
          STR1slp_tr = np.zeros((2,nt),dtype=np.float)
          STR1lat_tr = np.zeros_like(STR1slp_tr)

          STR2slp = np.zeros_like(STR1slp)
          STR2lat = np.zeros_like(STR2slp)
          STR2slp_tr = np.zeros_like(STR1slp_tr)
          STR2lat_tr = np.zeros_like(STR1slp_tr)


      if tscale == 'mon' :
          for im in np.arange(slp_z[:,0].size):
            #   print dt_nc[im].month
              STR1slp[0,im,yr-year[0]] = np.amax(slp_z[im,:])
              STR1lat[0,im,yr-year[0]] = latslpSH[np.argmax(slp_z[im,:])]

              STR2slp[0,im,yr-year[0]] = np.mean(np.amax(mslpSH[im,:,np.logical_and(lonslp<lon[1],lonslp>lon[0])],axis=1))
              STR2lat[0,im,yr-year[0]] = np.mean(latslpSH[np.argmax(mslpSH[im,:,np.logical_and(lonslp<lon[1],lonslp>lon[0])],axis=1)])

      if tscale == 'ssn' :
          for it,issn in enumerate(ssn[1:]) :
              cy = yr-year[0]
              if issn == 'YYY':
                  mons = np.arange(1,13)
              if issn == 'MAM':
                  mons = np.array([3,4,5])
              if issn == 'JJA':
                  mons = np.array([6,7,8])
              if issn == 'SON':
                  mons = np.array([9,10,11])
              if issn == 'DJF':
                  mons = np.array([1,2])

              mons = mons-1

              if not (issn == 'DJF' and cy==0) :

                  if not issn == 'DJF':
                      slp_z_mean = np.mean(slp_z[mons,:],axis=0)
                      slpSHlon_mean = np.mean(mslpSHlon[mons,:,:],axis=0)
                  else:
                      slp_z_tmp = np.array([slp_zDec_pyr,slp_z[1,:],slp_z[2,:]])
                      slp_z_mean = np.mean(slp_z_tmp,axis=0)
                      slpSHlon_tmp = np.array([mslpSHlonDec_pyr,mslpSHlon[1,:,:],mslpSHlon[2,:,:]])
                      slpSHlon_mean = np.mean(slpSHlon_tmp,axis=0)

                  STR1slp[0,it,cy] = np.amax(slp_z_mean)
                  STR1lat[0,it,cy] = latslpSH[np.argmax(slp_z_mean)]

                  STR2slp[0,it,cy] = np.mean(np.amax(slpSHlon_mean,axis=0))
                  STR2lat[0,it,cy] = np.mean(latslpSH[np.argmax(np.mean(mslpSHlon[mons,:,:],axis=0),axis=0)])

              else : # for DJF in first year
                  STR1slp[0,it,cy] = None
                  STR1lat[0,it,cy] = None

                  STR2slp[0,it,cy] = None
                  STR2lat[0,it,cy] = None

      slp_zDec_pyr = slp_zDec
      mslpSHlonDec_pyr = mslpSHlonDec


      nc.close()

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

    # SAM index
    nyr1 = 2000-1979
    for im in np.arange(12):
        for yr in yrs:
            # print yr, mslp40[yr-year[0],im],mslp65[yr-year[0],im]
            # print mslp40[yr-year[0],im]-mslp65[yr-year[0],im]
            # print (mslp40[yr-year[0],im]-np.mean(mslp40[:nyr1,im]))/np.std(mslp40[:nyr1,im])-(mslp65[yr-year[0],im]-np.mean(mslp65[:nyr1,im]))/np.std(mslp65[:nyr1,im])
            # print mslp40[:,im]-np.mean(mslp40[:nyr1,im])
            # print np.std(mslp40[:nyr1,im])
            # print mslp65[:,im]-np.mean(mslp65[:nyr1,im])
            # print np.std(mslp65[:nyr1,im])

            SAM[yr-year[0],im] = (mslp40[yr-year[0],im]-np.mean(mslp40[:nyr1,im]))/ \
            np.std(mslp40[:nyr1,im])-(mslp65[yr-year[0],im]-np.mean(mslp65[:nyr1,im]))/np.std(mslp65[:nyr1,im])
            SAM40[yr-year[0],im] = (mslp40[yr-year[0],im]-np.mean(mslp40[:nyr1,im]))/ \
            np.std(mslp40[:nyr1,im])
            SAM65[yr-year[0],im] = (mslp65[yr-year[0],im]-np.mean(mslp65[:nyr1,im]))/np.std(mslp65[:nyr1,im])
            print im, yr, SAM[yr-year[0],im]


    fSAM = "../output/SAM.%s.%d_%d.txt"%(reg,year[0],year[1]-1)
    with open(fSAM, "w") as text_file:
        text_file.write("    {:>7}\n".format('     '.join(month_abbr)))
        # print  ' '.join(month_abbr[1:])

        yrs2d = np.zeros ((nyrs,1))
        yrs2d[:,0] = yrs
        np.savetxt(text_file, np.append(yrs2d,SAM,axis=1), fmt='%7.2f')

        text_file.close()


    SAMssn = np.zeros_like(STR1slp);SAMssn.fill(np.nan)
    SAMssn_tr = np.zeros_like(STR1slp_tr);SAMssn.fill(np.nan)
    SAM40ssn = np.zeros_like(STR1slp);SAM40ssn.fill(np.nan)
    SAM40ssn_tr = np.zeros_like(STR1slp_tr);SAM40ssn.fill(np.nan)
    SAM65ssn = np.zeros_like(STR1slp);SAM65ssn.fill(np.nan)
    SAM65ssn_tr = np.zeros_like(STR1slp_tr);SAM65ssn.fill(np.nan)

    if tscale=='ssn':
      for it,issn in enumerate(ssn[1:]) :
          print issn
          cy = yr-year[0]
          if issn == 'YYY':
              mons = np.arange(1,13)
          if issn == 'MAM':
              mons = np.array([3,4,5])
          if issn == 'JJA':
              mons = np.array([6,7,8])
          if issn == 'SON':
              mons = np.array([9,10,11])
          if issn == 'DJF':
              mons = np.array([1,2])

          mons = mons-1

          SAMssn[0,it,:] = np.mean(SAM[:,mons],axis=1)
          if issn != "YYY":
              for yr in yrs:
                #   print yr, SAM[yr-year[0],mons]
                  if yr != year[0] and issn != "DJF":
                      SAMssn[0,it,yr-year[0]] = np.mean(SAM[yr-year[0],mons])
                      SAM40ssn[0,it,yr-year[0]] = np.mean(SAM40[yr-year[0],mons])
                      SAM65ssn[0,it,yr-year[0]] = np.mean(SAM65[yr-year[0],mons])
                  if yr != year[0] and issn == "DJF":
                      SAMssn[0,it,yr-year[0]] = np.mean([SAM[yr-year[0]-1,-1],SAM[yr-year[0],0],SAM[yr-year[0],1]])
                      SAM40ssn[0,it,yr-year[0]] = np.mean([SAM40[yr-year[0]-1,-1],SAM40[yr-year[0],0],SAM40[yr-year[0],1]])
                      SAM65ssn[0,it,yr-year[0]] = np.mean([SAM65[yr-year[0]-1,-1],SAM65[yr-year[0],0],SAM65[yr-year[0],1]])
            #   print SAMssn[0,it,:]
    #
    # print SAMssn[0,:,:]
    # print SAMssn.shape


    # quit()

    # # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # #  Fronts
    # # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #
    # for yr in yrs:
    #   fin = "../frontgrd/frline.%d.nc"%yr
    #   print "fronts from %s"%fin
    #   ncf = Dataset(fin, 'r')
    #
    #   print ncf.variables.keys()
    #   flon = ncf.variables['flon'][:]
    #   flat = ncf.variables['flat'][:]
    #   time = ncf.variables['time'][:]
    #   dv = ncf.variables['dv'][:]
    #   nf = ncf.variables['nfr'][:]
    #   npts = ncf.variables['npts'][:]
    #
    #   fdt = [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
    #        for t in time]
    #
    #
    #   if np.any(nf>maxnf) :
    #     print np.amax(nf)
    #     print "ERROR: nf > maxnf"
    #     quit()
    #
    #   if np.any(npts > maxnp) :
    #     print np.amax(npts)
    #     print "ERROR:  np > maxnp"
    #     quit()
    #
    #   # ************************************************
    #   # masking
    #   # ************************************************
    #   # create mask
    #
    #   for ind,t in enumerate(fdt):
    #     if t.year == yr and t.month==1 and t.day == 1 and t.hour == 0:
    #        tfr1 = ind
    #     if t.year == yr and t.month==12 and t.day == 31 and t.hour == 18:
    #        tfr2 = ind
    #        break
    #
    #   if yr == year[0]:
    #       if calendar.isleap(yr):
    #           frmask = np.zeros((nyrs,tfr2-tfr1+1,maxnf),dtype=np.int)
    #           frdv  = np.zeros((nyrs,tfr2-tfr1+1,maxnf),dtype=np.float)
    #           frmask_it = np.zeros((nt,nyrs,tfr2-tfr1+1,maxnf),dtype=np.int)
    #           frdv_it  = np.zeros((nt,nyrs,tfr2-tfr1+1,maxnf),dtype=np.float)
    #       else:
    #           frmask = np.zeros((nyrs,tfr2-tfr1+5,maxnf),dtype=np.int)
    #           frdv   = np.zeros((nyrs,tfr2-tfr1+5,maxnf),dtype=np.float32)
    #           frmask_it = np.zeros((nt,nyrs,tfr2-tfr1+5,maxnf),dtype=np.int)
    #           frdv_it   = np.zeros((nt,nyrs,tfr2-tfr1+5,maxnf),dtype=np.float32)
    #       frNPlat = np.zeros_like(frdv)
    #       frNPlat_it = np.zeros_like(frdv_it)
    #
    #
    #   for ct in np.arange(tfr1,tfr2+1):
    #       if fdt[ct].year>=year[0] and fdt[ct].year<=year[1]:
    #             for ifr in range(nf[ct]):
    #               if ifr+1>=nf1 and ifr+1<=nf2:
    #                   for ip in range(npts[ct,ifr]):
    #                       if flat[ct,ifr,ip]>=lat[0] and flat[ct,ifr,ip]<=lat[1]:
    #                           if flon[ct,ifr,ip]>=lon[0] and flon[ct,ifr,ip]<=lon[1]:
    #                               frmask[yr-year[0],ct-tfr1,ifr] = fdt[ct].month
    #                               frdv[yr-year[0],ct-tfr1,ifr] = np.sum(dv[ct,ifr,:npts[ct,ifr]])
    #                               frNPlat[yr-year[0],ct-tfr1,ifr] = np.amax(flat[ct,ifr,:npts[ct,ifr]])
    #                               break
    #
    #
    #
    # # ************************************************
    # # front stats
    # # ************************************************
    # if tscale == 'mon':
    #     for it in np.arange(nt):
    #         frmask_it[it,:,:,:] = np.where(frmask==it+1,1,0)
    #         frdv_it[it,:,:,:]   = np.where(frmask==it+1,frdv,0)
    #         frNPlat_it[it,:,:,:]= np.where(frmask==it+1,frNPlat,0)
    #
    # if tscale == 'ssn':
    #     for it,issn in enumerate(ssn[1:]):
    #       print issn
    #       if issn == 'YYY':
    #           mons = np.arange(1,13)
    #       if issn == 'MAM':
    #           mons = [3,4,5]
    #       if issn == 'JJA':
    #           mons = [6,7,8]
    #       if issn == 'SON':
    #           mons = [9,10,11]
    #       if issn == 'DJF':
    #           mons = [12,1,2]
    #
    #       if issn is not 'DJF' and issn is not 'YYY':
    #           frmask_it[it,:,:,:] = np.where(frmask==mons[0],1,frmask_it[it,:,:,:] )
    #           frmask_it[it,:,:,:] = np.where(frmask==mons[1],1,frmask_it[it,:,:,:] )
    #           frmask_it[it,:,:,:] = np.where(frmask==mons[2],1,frmask_it[it,:,:,:] )
    #           frdv_it[it,:,:,:]   = np.where(frmask_it[it,:,:,:]==1,frdv,0)
    #           frNPlat_it[it,:,:,:]= np.where(frmask_it[it,:,:,:]==1,frNPlat,0)
    #
    #       elif issn is 'DJF' :
    #           for yr in np.arange(1,nyrs):
    #               frmask_it[it,yr,:,:] = np.where(frmask[yr-1,:,:]==mons[0],1,frmask_it[it,yr,:,:] )
    #               frdv_it[it,yr,:,:]   = np.where(frmask[yr-1,:,:]==mons[0],frdv[yr-1,:,:],frdv_it[it,yr,:,:] )
    #               frNPlat_it[it,yr,:,:]= np.where(frmask[yr-1,:,:]==mons[0],frNPlat[yr-1,:,:],frNPlat_it[it,yr,:,:] )
    #
    #               frmask_it[it,yr,:,:] = np.where(frmask[yr,:,:]==mons[1],1,frmask_it[it,yr,:,:] )
    #               frmask_it[it,yr,:,:] = np.where(frmask[yr,:,:]==mons[2],1,frmask_it[it,yr,:,:] )
    #               frdv_it[it,yr,:,:]   = np.where(frmask[yr,:,:]==mons[1],frdv[yr,:,:],frdv_it[it,yr,:,:] )
    #               frNPlat_it[it,yr,:,:]= np.where(frmask[yr,:,:]==mons[1],frNPlat[yr,:,:],frNPlat_it[it,yr,:,:] )
    #               frdv_it[it,yr,:,:]   = np.where(frmask[yr,:,:]==mons[2],frdv[yr,:,:],frdv_it[it,yr,:,:] )
    #               frNPlat_it[it,yr,:,:]= np.where(frmask[yr,:,:]==mons[2],frNPlat[yr,:,:],frNPlat_it[it,yr,:,:] )
    #
    #       elif issn is  'YYY' :
    #               frmask_it[it,:,:,:] = np.where(frmask!=0,1,frmask_it[it,:,:,:])
    #               frdv_it[it,:,:,:]   = np.where(frmask!=0,frdv,0 )
    #               frNPlat_it[it,:,:,:]= np.where(frmask!=0,frNPlat,0 )
    #
    #
    #
    # dvcrit           = np.zeros(nt,dtype = np.float)
    # nfr_my           = np.zeros([ 3,nt,nyrs],dtype=np.float)
    # nfr_my3tr        = np.zeros_like(nfr_my)
    # frNPlat_my       = np.zeros_like(nfr_my)
    # frNPlat_my3tr    = np.zeros_like(nfr_my)
    #
    # nfr_my_tr        = np.zeros([ 2,nt],dtype=np.float)
    # nfr_my3tr_tr     = np.zeros_like(nfr_my_tr)
    # frNPlat_my_tr    = np.zeros_like(nfr_my_tr)
    # frNPlat_my3tr_tr = np.zeros_like(nfr_my_tr)
    #
    # for it in range(nt):
    #
    #     frdvtmp = frdv_it[it,:,:,:][np.where(frmask_it[it,:,:,:] == 1)]
    #     print frdvtmp.shape
    #
    #     dvcrit[it] =  np.percentile(frdvtmp,67)
    #     print  tperiod[it+1], " 67th %% threshold=",dvcrit[it]
    #
    #     # hist, bin_edges = np.histogram(frdvtmp,range=(0.1,25))
    #     # print "pdf: ",hist
    #     # print bin_edges
    #
    #
    #     for yr in yrs:
    #         for it in range(nt):
    #             if yr==year[0] and tperiod[it+1]=="DJF":
    #                 nfr_my[0,it,yr-year[0]]        = None
    #                 frNPlat_my[0,it,yr-year[0]]    = None
    #                 nfr_my3tr[0,it,yr-year[0]]     = None
    #                 frNPlat_my3tr[0,it,yr-year[0]] = None
    #             else:
    #                 nfr_my[0,it,yr-year[0]]=np.sum(frmask_it[it,yr-year[0],:,:]==1)
    #
    #                 tmp0 = frNPlat_it[it,yr-year[0],:,:][np.where(frmask_it[it,yr-year[0],:,:] == 1)]
    #                 frNPlat_my[0,it,yr-year[0]]=np.mean(tmp0)
    #
    #                 tmp  = frdv_it[it,yr-year[0],:,:][np.where(frmask_it[it,yr-year[0],:,:] == 1)]
    #                 nfr_my3tr[0,it,yr-year[0]]=np.count_nonzero(np.where(tmp>=dvcrit[it], 1,0))
    #                 frNPlat_my3tr[0,it,yr-year[0]]=np.mean(tmp0[np.where(tmp>=dvcrit[it])])
    #                 # nfr_my3tr[it,yr-year[0]]=np.sum(frmask[yr-year[0],:,:]==im and frdv[yr-year[0],:,:]>=dvcrit[im-1])
    #                 # print yr, im, month_abbr[im], frNPlat_my[0,im-1,yr-year[0]],frNPlat_my3tr[0,im-1,yr-year[0]]
    #                 # print np.count_nonzero(np.where(frmask[yr-year[0],:,:]== im, 1,0))
    #
    # # ************************************************
    # #  Record Numbers for multiple regression
    # # ************************************************
    #
    # for it in np.arange(nt):
    #     # print tperiod[it+1]
    #     fN = "../output/frN_STR.%s.%s%d_%d.txt"%(reg,tperiod[it+1],year[0],year[1]-1)
    #     with open(fN, "w") as text_file:
    #         text_file.write("{:>4}{:>7}{:>7}{:>7}\n".format("year","frN",  "STRloc",  "STRint"))
    #
    #         for yr in yrs:
    #             iyr = yr - year[0]
    #             # print "{:4d} {:.2f} {:.2f} {:.2f}\n".format(yr, nfr_my[0,it,iyr],STR1lat[0,it,yr-year[0]],STR1slp[0,it,yr-year[0]])
    #             text_file.write("{:4d} {:.2f} {:.2f} {:.2f}\n".format(yr, nfr_my[0,it,iyr],STR1lat[0,it,yr-year[0]],STR1slp[0,it,yr-year[0]]))
    #
    #         mask = ~np.isnan(nfr_my[0,it,:])
    #         x = np.array((np.ones(len(STR1lat[0,it,mask])),STR1lat[0,it,mask],STR1slp[0,it,mask]))
    #         X = x.T
    #         # print X
    #         # print nfr_my[0,it,mask]
    #         rmult = np.linalg.lstsq(X,nfr_my[0,it,mask])[0]
    #         # nfr_my[2,it,:] =  a[0]+a[1]*STR1lat[0,it,:]+a[2]*STR1slp[0,it,:]
    #         nfr_my[2,it,mask] = np.dot(X,rmult)
    #         # print "true: ",nfr_my[0,it,mask]
    #         # print mask
    #         nfr_my[2,it,~mask] = np.nan
    #         # print "fitted: ",nfr_my[2,it,:]
    #
    #
    # ************************************************
    # trends
    # ************************************************

    # fronts

    # title = "Trends in N of fronts"
    # filetrend = "../output/trend.Nfr.%d_%d.xls"%(year[0],year[1]-1)
    # trend(ireg,reg,nt,yrs,nfr_my,nfr_my_tr,title,tperiod,filetrend)
    #
    # title = "Trends in N of STRONG fronts"
    # filetrend = "../output/trend.Nfrstr.%d_%d.xls"%(year[0],year[1]-1)
    # trend(ireg,reg,nt,yrs,nfr_my3tr,nfr_my3tr_tr,title,tperiod,filetrend)
    #
    # title = "Trends in frontal NP lat"
    # filetrend = "../output/trend.frNPlat.%d_%d.xls"%(year[0],year[1]-1)
    # trend(ireg,reg,nt,yrs,frNPlat_my,frNPlat_my_tr,title,tperiod,filetrend)
    #
    # title = "Trends in STRONG fronts NP lat"
    # filetrend = "../output/trend.frNPlat_str.%d_%d.xls"%(year[0],year[1]-1)
    # trend(ireg,reg,nt,yrs,frNPlat_my3tr,frNPlat_my3tr_tr,title,tperiod,filetrend)
    #

    # STR

    title = "Trends in STR1 location"
    filetrend = "../output/trend.STR1loc.%d_%d.xls"%(year[0],year[1]-1)
    trend(ireg,reg,nt,yrs,STR1lat,STR1lat_tr,title,tperiod,filetrend)

    title = "Trends in STR1 intensity"
    filetrend = "../output/trend.STR1int.%d_%d.xls"%(year[0],year[1]-1)
    trend(ireg,reg,nt,yrs,STR1slp,STR1slp_tr,title,tperiod,filetrend)

    title = "Trends in STR2 location"
    filetrend = "../output/trend.STR2loc.%d_%d.xls"%(year[0],year[1]-1)
    trend(ireg,reg,nt,yrs,STR2lat,STR2lat_tr,title,tperiod,filetrend)

    title = "Trends in STR2 intensity"
    filetrend = "../output/trend.STR2int.%d_%d.xls"%(year[0],year[1]-1)
    trend(ireg,reg,nt,yrs,STR2slp,STR2slp_tr,title,tperiod,filetrend)

    # SAM

    title = "Trends in SAM"
    filetrend = "../output/trend.SAM.%d_%d.xls"%(year[0],year[1]-1)
    trend(ireg,reg,nt,yrs,SAMssn,SAMssn_tr,title,tperiod,filetrend)

    print SAMssn

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Plotting
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    majorLocator = MultipleLocator(5)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(1)

    if tscale == "mon":
        nc = 2
        nr = 6
        wdth = 3.5
    if tscale == "ssn":
        nc = 2
        nr = nt
        wdth = 3.


    # plt.ion()
    plt.close('all')
    for ifig in np.arange(fig1,fig2+1):
        # fig = plt.figure(ifig)
        # f, ax = plt.subplots(nr, nc,  sharex='col', sharey='row',figsize=(wdth*3.13,4.2*3.13))
        f, ax = plt.subplots(nr, nc,  sharex='col', figsize=(wdth*3.13,4.2*3.13))

        if ifig ==0:

            title = 'SAM%s vs (STR loc, STR int), %s'%(SAMscale,reg)
            fout = "SAM%s_STRloc_STRlat"%(SAMscale)

            plt.suptitle(title,fontsize=14)

        # if any ([0,2,11] == ifig ):
            print 'set var1 == SAMssn'
            var1 = SAMssn
            yax1 = [-5, 5]
            yax1label = 'SAM'
            var1legend = 'SAM'

            var2 = STR1lat
            yax2 = [-40, -25]
            yax2label = 'latitude'
            var2legend = 'STR location'

            var3 = STR1slp
            yax3 = [1010, 1025]
            yax3label = 'MSL Pressure (Pa)'
            var3legend = 'STR intensity'


        if ifig ==1:
            title = 'SAM%s vs (mslp40, mslp65), %s'%(SAMscale,reg)
            fout = "SAM%s_mslp"%(SAMscale)

            plt.suptitle(title,fontsize=14)

        # if any ([0,2,11] == ifig ):
            print 'set var1 == SAMssn'
            var1 = SAMssn
            yax1 = [-5, 5]
            yax1label = 'SAM'
            var1legend = 'SAM'

            var2 = SAM40ssn
            yax2 = [-5, 5]
            yax2label = 'mslp40'
            var2legend = 'mslp40'

            var3 = -SAM65ssn
            yax3 = [-5, 5]
            yax3label = '-mslp65'
            var3legend = '(-1)mslp65'

        if ifig ==2:
            title = 'mslp40 vs mslp65, %s'%reg
            fout = "mslp40_mslp65"

            plt.suptitle(title,fontsize=14)

        # if any ([0,2,11] == ifig ):

            var1 = SAM40ssn
            yax1 = [-5, 5]
            yax1label = 'mslp40'
            var1legend = 'mslp40'

            var2 = -SAM65ssn
            yax2 = [-5, 5]
            yax2label = '-mslp65'
            var2legend = '(-1)mslp65'

        if ifig ==3:
            title = 'SAM%s vs ENSO3.4, %s'%(SAMscale,reg)
            fout = "SAM%s_ONI"%(SAMscale)

            plt.suptitle(title,fontsize=14)

        # if any ([0,2,11] == ifig ):

            var1 = SAMssn
            yax1 = [-5, 5]
            yax1label = 'SAMreg'
            var1legend = 'SAMreg'

            var2 = ONIssn
            yax2 = [-5, 5]
            yax2label = 'ENSO3.4'
            var2legend = 'ENSO3.4'


        # for ic in np.arange(nc):
        for ic in np.arange(2):
            for ir in np.arange(nr):
            # for ir in [2,3,4,5]:
                it = ir+ic*6

                # if tscale == 'mon':
                a = ax[ir, ic]
                # if tscale == 'ssn':
                #     a = ax[it]
                a2 = a.twinx()

                a.plot(yrs,var1[0,ir,:],color='b',lw=1,label=var1legend)
                if ic == 0:
                    var = var2
                    yax = yax2
                    yaxlabel = yax2label
                    varlegend = var2legend
                else:
                    var = var3
                    yax = yax3
                    yaxlabel = yax3label
                    varlegend = var3legend

                a2.plot(yrs,var[0,ir,:],color='g',lw=1,label=varlegend)


                # ax[ir, ic].set_title('STR intensity, %s, SH'%tperiod[it+1])

                # a.axis((year[0]-1, year[1], yax1[0], yax1[1]))
                # a2.axis((year[0]-1, year[1], yax2[0], yax2[1]))
                # ax2.set_ylim(0, 35)
                print tperiod[ir+1],np.amin(var1[0,ir,:]),np.amax(var1[0,ir,:])
                # a.set_ylim(0.9*np.amin(var1[0,it,:]),1.1*np.amax(var1[0,it,:]))
                # a2.set_ylim(0.9*np.amin(var1[0,it,:]),1.1*np.amax(var1[0,it,:]))

                # a.xaxis.set_ticks(yrs)
                # xa = a.get_xaxis()
                # a.xaxis.set_tick_params(axis='x',which='major',length=6,width =2)
                # a.xaxis.set_tick_params(axis='x',which='minor',length=3,width =2)

                a.xaxis.set_major_locator(majorLocator)
                a.xaxis.set_major_formatter(majorFormatter)
                # a.yaxis.set_major_formatter(majorFormatter)
                # a2.yaxis.set_major_formatter(majorFormatter)
                # a.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                # a2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

                # for the minor ticks, use no labels; default NullFormatter
                a.xaxis.set_minor_locator(minorLocator)

                a.grid()

                if tscale == "mon":
                    # if it == 3:
                    #     a.set_ylabel(yax1label,color='b')
                    # if it == 9:
                    #     a2.set_ylabel(yax2label,color='g')
                    if ir == 5 :
                        a.set_xlabel('Year')
                if tscale == "ssn":
                    if ir == 3 :
                        a.set_ylabel(yax1label,color='b')
                        if nc == 0:
                            a2.set_ylabel(yax2label,color='g')
                        if nc == 1:
                            a2.set_ylabel(yax3abel,color='g')

                    if ir == nt-1 :
                        a.set_xlabel('Year')

                #  trend & correlation

                a.plot(yrs,var1[1,ir,:],color='b',lw=2)
                a2.plot(yrs,var[1,ir,:],color='g',lw=2)

                mask  = ~np.isnan(var1[0,ir,:]) & ~np.isnan(var[0,ir,:])
                cctmp = np.corrcoef(var1[0,ir,mask],var[0,ir,mask])

                cc[ireg,ifig,0,ir] = cctmp[1,0]

                cctmp = np.corrcoef(var1[0,ir,mask]-var1[1,ir,mask],var[0,ir,mask]-var[1,ir,mask])
                cc[ireg,ifig,1,ir] = cctmp[1,0]
                a.set_title(' %s, %s, r = %.2f, r_dt = %.2f'%(tperiod[ir+1],reg,cc[ireg,ifig,0,ir],cc[ireg,ifig,1,ir]))

                if ir == nt-1 :
                    # ask matplotlib for the plotted objects and their labels
                    lines, labels = a.get_legend_handles_labels()
                    lines2, labels2 = a2.get_legend_handles_labels()
                    a2.legend(lines + lines2, labels + labels2, loc='lower right')

                    # ax[ir, ic].legend(loc='lower right')


                # ************************************************
                #  correlations to Excel
                # ************************************************
                if ifig == 0:
                    if ic == 0:
                        fxls = "../output/SAM%s_STRlat.%d_%d.xls"%(SAMscale,year[0],year[1]-1)
                    elif  ic ==1:
                        fxls = "../output/SAM%s_STRslp.%d_%d.xls"%(SAMscale,year[0],year[1]-1)
                    if ireg == 0 and os.path.exists(fxls):
                            os.remove(fxls)
                    output(ireg,reg,fxls, "corr", title, tperiod, cc[ireg,ifig,:,:])


        f.subplots_adjust(hspace=0.3,wspace = 0.3)
        plt.draw()
        f.savefig("../output/%s.%s.mon%d_%d.png"%(fout,reg,year[0],year[1]-1))
        plt.close(f)
