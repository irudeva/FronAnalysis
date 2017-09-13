from netCDF4 import Dataset
import datetime as datetime
import calendar
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np
from scipy import stats
from scipy import signal
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

nfig = 12  # possible combinations of vars
fig1 = 12
fig2 = 12


hs = "SH"
# reg
nreg = 0

cc = np.zeros((nreg+1,nfig+1,2,nt),dtype=np.float) # correlation



for ireg in np.arange(0,nreg+1):
    lat = [-40, -20]
    if ireg == 0:
        lon = [ -90,361]
        # lat = [ -40,-20]
        reg = "%d_%dS"%(np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 1:
        lon = [ 30,90]
        # lat = [ -40,-20]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 2:
        lon = [ 90,150]
        # lat = [ -40,-20]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 3:
        lon = [ 150,210]
        # lat = [ -40,-20]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 4:
        lon = [ 210,285]
        # lat = [ -40,-20]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))

    if ireg == 5:
        lon = [ 300,20]
        # lat = [ -40,-20]
        reg = "%d_%dE.%d_%dS"%(np.abs(lon[0]),np.abs(lon[1]),np.abs(lat[1]),np.abs(lat[0]))



    maxnf = 200 # max # of fronts per timestep
    maxnp = 100 # max # of frontal points


    # ************************************************
    # read in netCDF files
    # ************************************************
    # 1. slp
    # ************************************************

    # mslp40 = np.zeros((nyrs,12)); mslp40.fill(np.nan)
    # mslp65 = np.zeros((nyrs,12)); mslp65.fill(np.nan)
    #
    # SAM = np.zeros((nyrs,12)); SAM.fill(np.nan)
    # SAM40 = np.zeros((nyrs,12)); SAM40.fill(np.nan)
    # SAM65 = np.zeros((nyrs,12)); SAM65.fill(np.nan)

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


      latslpSH  = latslp[np.logical_and(latslp<=-20,latslp>=-65)]
      mslpSH    = mslp[:,np.logical_and(latslp<=-20,latslp>=-65),:]
      if lon[0]<lon[1]:
          mslpSHlon = mslpSH[:,:,np.logical_and(lonslp<=lon[1],lonslp>=lon[0])]
      elif lon[0]>lon[1]:
          mslpSHlon = mslpSH[:,:,lonslp>=lon[0]]
          mslpSHlon = np.append(mslpSHlon,mslpSH[:,:,lonslp<=lon[1]],axis=2)

      mslpSHlonDec = mslpSHlon[-1,:,:]


      # 5 deg longitude bands average:
      wlon = 1.
      lonslp5deg = np.arange(0,360,wlon)
      slp_z = np.zeros((time.size, latslpSH.size,lonslp5deg.size));slp_z.fill(np.nan)
      ilon = 0
      slp_zDec = np.zeros((latslpSH.size,lonslp5deg.size));slp_zDec.fill(np.nan)
      for i,ilon in enumerate(lonslp5deg):
          slp_z[:,:,i] = np.mean(mslpSHlon[:,:,np.logical_and(lonslp>=ilon,lonslp<ilon+wlon)],axis=2)
          slp_zDec[:,i] = slp_z[-1,:,i]

      print slp_z.shape

    #   quit()

      if yr==year[0]:
          STR1slp = np.zeros((nt,nyrs,lonslp5deg.size),dtype=np.float)
          STR1lat = np.zeros_like(STR1slp)


    #   if tscale == 'mon' :
    #       for im in np.arange(slp_z[:,0].size):
    #         #   print dt_nc[im].month
    #           STR1slp[0,im,yr-year[0]] = np.amax(slp_z[im,:])
    #           STR1lat[0,im,yr-year[0]] = latslpSH[np.argmax(slp_z[im,:])]
      #
    #           if lon[0]<lon[1]:
    #               STR2slp[0,im,yr-year[0]] = np.mean(np.amax(mslpSH[im,:,np.logical_and(lonslp<=lon[1],lonslp>=lon[0])],axis=1))
    #               STR2lat[0,im,yr-year[0]] = np.mean(latslpSH[np.argmax(mslpSH[im,:,np.logical_and(lonslp<=lon[1],lonslp>=lon[0])],axis=1)])
    #           elif lon[0]>lon[1]:
    #               STR2slp[0,im,yr-year[0]] = np.mean(np.amax(np.append(mslpSH[im,:,lonslp>=lon[0]],mslpSH[im,:,lonslp<=lon[1]],axis=2),axis=1))
    #               STR2lat[0,im,yr-year[0]] = np.mean(latslpSH[np.argmax(np.append(mslpSH[im,:,lonslp>=lon[0]],mslpSH[im,:,lonslp<=lon[1]],axis=2),axis=1)])


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
                      slp_z_mean = np.mean(slp_z[mons,:,:],axis=0)
                  else:
                      slp_z_tmp = np.array([slp_zDec_pyr,slp_z[1,:,:],slp_z[2,:,:]])
                      slp_z_mean = np.mean(slp_z_tmp,axis=0)

                  STR1slp[it,cy,:] = np.amax(slp_z_mean,axis = 0)
                  STR1lat[it,cy,:] = latslpSH[np.argmax(slp_z_mean,axis = 0)]


              else : # for DJF in first year
                  STR1slp[it,cy,:] = None
                  STR1lat[it,cy,:] = None

      slp_zDec_pyr = slp_zDec


      nc.close()

    # write STR years to a file
    for it in np.arange(nt):
        fSTRloc5deg = "../output/STRloc_5deg.%s.%d_%d.txt"%(reg,year[0],year[1]-1)
        STRlat_ssn= np.zeros((nt,lonslp5deg.size),dtype=np.float)
        with open(fSTRloc5deg, "w") as text_file:
            text_file.write("       {:>7}\n".format('     '.join(ssn)))
            # print  ' '.join(month_abbr[1:])

            lon5d = np.zeros ((lonslp5deg.size,1),float)

            for it,issn in enumerate(ssn[1:]) :
              y1 = 0
              if issn == 'DJF':
                  y1 = 1

              print np.mean(STR1lat[it,y1:,:],axis=0).shape
              STRlat_ssn[it,:]=(np.mean(STR1lat[it,y1:,:],axis=0))

            lon5d[:,0] = lonslp5deg+wlon/2
            np.savetxt(text_file, np.append(lon5d,np.transpose(STRlat_ssn),axis=1), fmt='%7.2f')

            text_file.close()
